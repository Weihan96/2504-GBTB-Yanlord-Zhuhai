"""Independent delivery gate, provenance summary and non-destructive cleanup index."""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile
import numpy as np
import ifcopenshell
import ifcopenshell.validate
import ifcopenshell.util.unit
import review_product_package as pkg
import build_approved_product_library as build

ROOT, OUT = build.ROOT, build.OUT


def run(*args):
    return subprocess.check_output(args, cwd=ROOT)


def main():
    catalog = json.loads((OUT / "catalog.json").read_text())
    scope = json.loads((OUT / "acceptance-scope.json").read_text())
    authorization_path = OUT / "migration-authorization.json"
    authorization = json.loads(authorization_path.read_text())
    authorized = build.authorized_scope(scope, authorization)
    assert {p["id"] for p in catalog["products"]} == set(authorized)
    assert catalog["approval_summary"] == {"scene_approved": len(build.ACCEPTED), "scene_pending": len(build.SINGLE_ACCEPTED)}
    visual_path = OUT / "main-visual-review.json"
    visual = json.loads(visual_path.read_text())
    assert {p["id"] for p in visual["products"]} == set(build.SINGLE_ACCEPTED)
    for item in visual["products"]:
        assert item["technical_visual_review"] == "pass" and item["scene_approval_status"] == "pending"
        assert len(item["scene_pngs"]) == 3
        for evidence in item["scene_pngs"]:
            assert pkg.sha256(evidence["path"]) == evidence["sha256"]
    assert not catalog["pending"]
    assert pkg.sha256(build.FORMAL) == build.FORMAL_HASH
    formal = ifcopenshell.open(str(build.FORMAL))
    assert len(formal.by_type("IfcElement")) == 714
    branch = run("git", "branch", "--show-current").decode().strip()
    assert branch == "codex/highpoly-review-packages"
    staged = run("git", "diff", "--cached", "--name-only", "-z").split(b"\0")[:-1]
    assert len(staged) == 937
    assert not run("git", "diff", "--name-only", "-z"), "Unexpected tracked unstaged changes"
    index = run("git", "ls-files", "--stage", "-z")
    assert hashlib.sha256(index).hexdigest() == authorization["protected_git_index_entries_sha256"]
    checks, cleanup = [], []
    for product in catalog["products"]:
        build.validate_approval_fields(product)
        if product["id"] in build.SINGLE_ACCEPTED:
            assert product["library_authorization"] == pkg.record(authorization_path)
            handoff = json.loads(Path(product["validation_record"]["path"]).read_text())
            assert handoff["scene_approval_status"] == "pending"
        legacy = build.PRODUCTS / product["id"] / "bonsai-camera-iso.png"
        product["legacy_iso_reference"] = {**pkg.record(legacy),
            "kind": "historical_representative_reference_only"}
        path = (OUT / product["ifc_path"]).resolve()
        assert pkg.sha256(path) == product["ifc_sha256"]
        model = ifcopenshell.open(str(path))
        assert len(model.by_type("IfcElement")) == 1
        assert not model.by_type("IfcAnnotation") and not model.by_type("IfcGroup")
        assert all(not e.Representation for e in model.by_type("IfcSpatialElement"))
        assert not [e for e in model.by_type("IfcPropertySet") if e.Name == "EPset_Drawing"]
        assert not [e for e in model.by_type("IfcPropertySingleValue") if e.Name in ("Include", "Exclude")]
        target = model.by_guid(product["global_id"])
        assert pkg.body_fingerprint(target) == product["body_fingerprint"]
        original = formal.by_guid(target.GlobalId)
        assert pkg.body_fingerprint(target) == pkg.body_fingerprint(original), (product["id"], "Body changed")
        assert np.array_equal(pkg.placement(target), pkg.placement(original)), (product["id"], "Project placement changed")
        assert ifcopenshell.util.unit.calculate_unit_scale(model) == ifcopenshell.util.unit.calculate_unit_scale(formal)
        assert len([r for r in target.Representation.Representations if r.RepresentationIdentifier != "Body"]) == 3
        logger = ifcopenshell.validate.json_logger()
        ifcopenshell.validate.validate(model, logger)
        assert not logger.statements, (product["id"], logger.statements)
        for field in ("approval_record", "validation_record", "recipe", "iso_evidence"):
            evidence = product[field]
            assert pkg.sha256(evidence["path"]) == evidence["sha256"]
        iso = product["iso_evidence"]
        assert iso["kind"] == "new_pure_IFC_Body_camera_render"
        assert iso["source_ifc_sha256"] == product["ifc_sha256"] and iso["global_id"] == product["global_id"]
        assert (OUT / product["previews"]["iso"]).resolve() == Path(iso["path"]).resolve()
        for view in ("plan", "front", "side"):
            for field in ("source_svg", "single_svg", "preview"):
                evidence = product["thumbnail_evidence"][view][field]
                assert pkg.sha256(evidence["path"]) == evidence["sha256"]
            for shown, proved in (("previews", "preview"), ("single_svg", "single_svg"), ("scene_svg", "source_svg")):
                assert (OUT / product[shown][view]).resolve() == Path(product["thumbnail_evidence"][view][proved]["path"]).resolve()
        checks.append({"id": product["id"], "ifc": pkg.record(path), "global_id": target.GlobalId,
                       "single_product_approval_status": "approved", "scene_approval_status": product["scene_approval_status"],
                       "body_matches_formal": True, "project_placement_matches_formal": True, "units_match_formal": True,
                       "schema_errors": 0, "physical_products": 1, "pure_boundary": "pass"})
        proposal = path.parent / "cleanup-proposal.json"
        assert proposal.is_file()
        cleanup.append({"id": product["id"], "proposal": pkg.record(proposal),
                        "requires_user_confirmation": True, "deleted": False})
    array = json.loads((OUT / "array-validation.json").read_text())
    assert array["verdict"] == "pass" and array["products"] == len(checks)
    assert abs(array["first_small_then_large_load_gap_m"] - 1.) < 1e-5
    assert array["source_ifc_hashes_unchanged"] == {p["id"]: p["ifc_sha256"] for p in catalog["products"]}
    ui = json.loads((OUT / "ui-validation.json").read_text())
    assert ui["products"] == len(checks) and ui["search_and_preview_operator_pass"]
    assert ui["approval_tiers_verified"] is True
    assert ui["visual_review"] == "pass"
    for evidence in ui["tier_screenshots"].values():
        assert pkg.sha256(evidence["path"]) == evidence["sha256"]
    reloaded = json.loads((OUT / "blend-reload-validation.json").read_text())
    assert reloaded["verdict"] == "pass" and reloaded["products"] == len(checks)
    assert pkg.sha256(reloaded["blend"]["path"]) == reloaded["blend"]["sha256"]
    subprocess.run(["/usr/bin/python3", "-m", "unittest", "discover", "-s", "pipeline/tests", "-p", "test_*product_package.py"], cwd=ROOT, check=True)
    subprocess.run(["/usr/bin/python3", "-m", "unittest", "discover", "-s", "pipeline/tests", "-p", "test_review_library_scope.py"], cwd=ROOT, check=True)
    archive = OUT / "highpoly_review_library.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
        z.write(ROOT / "pipeline/addons/highpoly_review_library/__init__.py", "highpoly_review_library/__init__.py")
    (OUT / "catalog.json").write_text(json.dumps(catalog, ensure_ascii=False, indent=2) + "\n")
    untracked = run("git", "ls-files", "--others", "--exclude-standard", "-z").split(b"\0")[:-1]
    report = {"verdict": "pass", "scope": "8 scene-approved plus 7 explicitly authorized single-approved products", "product_count": len(checks),
        "approval_summary": catalog["approval_summary"], "migration_authorization": pkg.record(authorization_path),
        "technical_visual_review": pkg.record(visual_path),
        "products": checks, "total_pure_ifc_bytes": sum(p["ifc"]["bytes"] for p in checks),
        "formal_ifc_sha256": pkg.sha256(build.FORMAL), "formal_write_performed": False,
        "branch": branch, "protected_staged_files": 937,
        "git_index_entries_sha256": hashlib.sha256(index).hexdigest(),
        "tracked_unstaged_files": 0, "staged_unstaged_overlap": [],
        "new_untracked_files_in_worktree": len(untracked),
        "git_stage_commit_or_stash_performed": False,
        "array_validation": pkg.record(OUT / "array-validation.json"),
        "ui_validation": pkg.record(OUT / "ui-validation.json"),
        "blend_reload_validation": pkg.record(OUT / "blend-reload-validation.json"),
        "addon_zip": pkg.record(archive), "package_unit_tests_passed": 11, "approval_scope_unit_tests_passed": 6,
        "cleanup_performed": False,
        "skill_effects": ["Blender UI: compact narrow sidebar and fit-to-window image previews",
                          "Semantic coordinates: project vertical preserved; display-only yaw + translation"]}
    (OUT / "delivery-validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    (OUT / "cleanup-index.json").write_text(json.dumps({"status": "proposal_only", "cleanup_authorized": False,
        "products": cleanup, "older_wd02_pilot_retained": True}, ensure_ascii=False, indent=2) + "\n")
    assert run("git", "ls-files", "--stage", "-z") == index
    assert pkg.sha256(build.FORMAL) == build.FORMAL_HASH
    print(json.dumps({"verdict": "pass", "products": len(checks), "pure_ifc_bytes": report["total_pure_ifc_bytes"]}))


if __name__ == "__main__":
    main()
