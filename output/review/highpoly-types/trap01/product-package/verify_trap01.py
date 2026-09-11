"""Independent persisted IFC and regenerated SVG checks; no legacy changes."""
from pathlib import Path
import json
import subprocess
import xml.etree.ElementTree as ET
import ifcopenshell
import ifcopenshell.validate
from shapely.geometry import LineString
from shapely.ops import unary_union
import migrate_trap01 as m
import review_product_package as pkg
import pure_product_package as pure_pkg


def line_geometry(path, guid):
    root = ET.parse(path).getroot()
    segments = [e for e in root.iter() if f"GlobalId-{guid}" in e.get("class", "").split()]
    assert segments and all(e.tag.endswith("line") for e in segments)
    assert not any(e.tag.endswith("image") for e in root.iter())
    return unary_union([LineString([(float(e.get("x1")), float(e.get("y1"))), (float(e.get("x2")), float(e.get("y2")))]) for e in segments]), len(segments)


def main():
    report = json.loads(m.REPORT.read_text())
    assert report["stage"] == "bonsai_roundtrip_and_scene_complete", report
    assert report["scene_reload_verified"]
    audit = json.loads((m.OUT / "source-audit.json").read_text())
    pure = ifcopenshell.open(str(m.SINGLE))
    state = m.pure_state(pure)
    assert state == report["pure_roundtrip"]
    assert state["body_fingerprint"] == audit["target_body_fingerprint"]
    assert state["placement"] == audit["target_placement"]
    logger = ifcopenshell.validate.json_logger()
    ifcopenshell.validate.validate(pure, logger)
    assert not logger.statements, logger.statements
    recipe = json.loads(m.RECIPE.read_text())
    graph = pure_pkg.graph_from_json(recipe["graph"])
    assert all(a.Representation is None for a in graph.by_type("IfcAnnotation") if a.ObjectType != "DRAWING")
    assert all(e.Representation is None for e in graph.by_type("IfcElement"))
    comparisons = []
    for old in audit["approved_scene_views"]:
        view = old["view"]
        guid = old["annotation_global_id"]
        path = m.OUT / f"TRAP01-SCENE-{view.upper()}.svg"
        old_lines, old_count = line_geometry(m.ROOT / old["svg"]["path"], guid)
        new_lines, new_count = line_geometry(path, guid)
        error_mm = old_lines.hausdorff_distance(new_lines) * 5
        # Native serializer may merge microscopic SVG segments. Persisted IFC
        # curve content is checked separately and must remain exactly equal.
        assert error_mm < .01, (view, error_mm)
        preview = path.with_suffix(".png")
        subprocess.run(["/Applications/Inkscape.app/Contents/MacOS/inkscape", str(path),
                        "--export-area-page", "--export-background=white", "--export-background-opacity=1", "--export-width=1600", f"--export-filename={preview}"], check=True, capture_output=True)
        comparisons.append({"view": view, "approved_svg": old["svg"], "new_svg": pkg.record(path), "preview": pkg.record(preview),
                            "annotation_guid": guid, "old_line_segments": old_count, "new_line_segments": new_count,
                            "maximum_world_linework_error_mm": error_mm, "tolerance_mm": .01,
                            "segment_count_note": "Native Bonsai output consolidates sub-0.01 mm segments; pure IFC geometry fingerprint remains exact. Compare physical line geometry, not serialization segment count."})
    assert all(pkg.sha256(r["path"]) == r["sha256"] for r in audit["protected_files"])
    assert pkg.sha256(m.FORMAL) == m.FORMAL_SHA
    report.update({"independent_ifc_schema_errors": 0, "independent_pure_state": state, "approved_scene_comparison": comparisons,
                   "external_recipe_contains_no_linework_or_body_geometry": True, "linework_geometry_source": "pure_product_ifc_representations",
                   "protected_source_and_approval_files_unchanged": True, "formal_sha256_after_independent_verification": pkg.sha256(m.FORMAL),
                   "verdict": "validated_pending_visual_inspection"})
    m.write(m.REPORT, report)
    legacy = [m.OLD, m.PRODUCT / "Geberit-151.116.11.1-TRAP01-official-detail-v4.blend", m.PRODUCT / "Geberit-151.116.11.1-TRAP01-official-detail-v4.blend1"]
    m.write(m.OUT / "cleanup-proposal.json", {"status": "proposal_only_no_deletion", "requires_user_confirmation": True,
        "legacy_files": [pkg.record(p) for p in legacy if p.is_file()], "temporary_directories": [report["temporary_directory"]],
        "retained_older_superseded_versions": "v3 and earlier records are not selected for cleanup by this migration",
        "keep": [str(m.SINGLE), str(m.RECIPE), "scene SVG/PNG", "verification and audit records", "official-source files", str(m.APPROVAL)]})
    m.write(m.OUT / "package-manifest.json", {"schema_version": 1, "profile_key": "trap01", "status": "validated_pending_visual_inspection",
        "target_global_id": m.GUID, "persistent_ifc": pkg.record(m.SINGLE), "scene_recipe": pkg.record(m.RECIPE),
        "pure_boundary": "one product Body + approved official v4 Plan/Front/Side + necessary coordinate/unit/property/style dependencies; zero annotations, drawing cameras, scene filters or other product geometry",
        "source_kind": "native_dwg_configured_installation", "source_label_zh": "基于官方原生DWG固定接头刚性移位及直管长度调整的二维安装细节表达",
        "approval": pkg.record(m.APPROVAL), "approval_scope": "official-detail-v4 full installed detail, not old proxy outlines", "source_audit": pkg.record(m.OUT / "source-audit.json"),
        "line_semantics": audit["line_semantics"], "dashed_reference_policy": audit["blue_dashed_reference_policy"],
        "legacy_full_project_ifc": pkg.record(m.OLD), "size_reduction_percent": (1 - m.SINGLE.stat().st_size / m.OLD.stat().st_size) * 100,
        "scene_outputs": comparisons, "validation_record": "validation.json", "cleanup_proposal": "cleanup-proposal.json",
        "formal_authoritative_ifc_write_allowed": False, "formal_ifc_sha256": m.FORMAL_SHA, "legacy_copy_retained": True})
    print(json.dumps({"schema_errors": 0, "comparisons": comparisons, "formal_unchanged": True}, ensure_ascii=False))


if __name__ == "__main__":
    main()
