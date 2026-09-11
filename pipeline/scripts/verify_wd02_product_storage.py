"""Independent disk/SVG validation and a non-destructive cleanup proposal."""
import json
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
import ifcopenshell
import ifcopenshell.validate
from shapely.geometry import LineString
from shapely.ops import unary_union
import review_product_package as pkg
import pilot_wd02_product_storage as pilot


def line_geometry(svg, guid):
    root = ET.parse(svg).getroot()
    lines = [e for e in root.iter() if guid in e.get("class", "")]
    assert lines and all(e.tag.endswith("line") for e in lines)
    for e in root.iter():
        assert not e.tag.endswith("image"), "Unexpected external raster dependency"
    return unary_union([LineString([(float(e.get("x1")), float(e.get("y1"))),
                                    (float(e.get("x2")), float(e.get("y2")))]) for e in lines])


def main():
    report = json.loads(pilot.REPORT.read_text())
    assert report["stage"] == "bonsai_roundtrip_and_scene_complete"
    assert report["scene_reload_verified"]
    source = ifcopenshell.open(str(pilot.OLD))
    package = ifcopenshell.open(str(pilot.SINGLE))
    result = pkg.validate(source, package, pilot.GUID, pilot.annotation_guids())
    logger = ifcopenshell.validate.json_logger()
    ifcopenshell.validate.validate(package, logger)
    assert not logger.statements, logger.statements
    result["ifc_schema_validation_errors"] = 0
    comparisons = []
    for v in pilot.view_records():
        path = pilot.OUT / f"WD02-SCENE-{v['view'].upper()}.svg"
        old = line_geometry(Path(v["svg"]["path"]), v["annotation_global_id"])
        new = line_geometry(path, v["annotation_global_id"])
        error_mm = old.segmentize(.1).hausdorff_distance(new.segmentize(.1)) * 25
        assert error_mm < .01, f"{v['view']} differs from approved scene: {error_mm} mm"
        png = path.with_suffix(".png")
        subprocess.run(["/Applications/Inkscape.app/Contents/MacOS/inkscape", str(path),
                        "--export-area-page", "--export-background=white", "--export-background-opacity=1",
                        "--export-width=1400", f"--export-filename={png}"], check=True, capture_output=True)
        drawing = package.by_guid(v["drawing_global_id"])
        references = [r.RelatingDocument for r in drawing.HasAssociations if r.is_a("IfcRelAssociatesDocument")]
        assert len(references) == 1 and (pilot.OUT / references[0].Location).resolve() == path.resolve()
        comparisons.append({"view": v["view"], "max_linework_hausdorff_error_world_mm": error_mm,
                            "tolerance_mm": .01, "svg": pkg.record(path), "preview": pkg.record(png),
                            "document_reference_resolves": True})
    protected = report["preState"]["protected_files"]
    assert all(pkg.sha256(r["path"]) == r["sha256"] for r in protected)
    for resource in report.get("temporary_external_style_resources", []):
        if resource.get("source"):
            assert pkg.sha256(resource["source"]["path"]) == resource["source"]["sha256"]
    assert pkg.sha256(pilot.FORMAL) == pilot.FORMAL_HASH
    report["independent_validation"] = result
    report["approved_scene_comparison"] = comparisons
    report["formal_sha256_after_independent_validation"] = pkg.sha256(pilot.FORMAL)
    report["verdict"] = "validated_pending_final_visual_check"
    report["known_limitations"] = [
        "Pilot only: one physical product without nested parts/openings/ports; unsupported relationship scopes stop extraction.",
        "OD_Textures/Materials.blend is absent at the formal source, not lost by packaging; current 2D Create Drawing does not use it.",
        "Second bridge call exceeded its response deadline; do not infer failure/success from transport. Disk report, reload queries and independent verification confirmed completion.",
        "Old full-project IFC and Blend remain until explicit cleanup approval. No batch migration yet."
    ]
    pilot.write(pilot.REPORT, report)
    cleanup = []
    for path, reason in [(pilot.OLD, "Replaced by validated product IFC for review storage"),
                         (pilot.PRODUCT / "Poliform-Senzafine-WD02-project-drawings.blend", "Legacy complete-project work session; IFC is reopened for acceptance")]:
        if path.is_file():
            cleanup.append({**pkg.record(path), "reason": reason, "requires_user_confirmation": True})
    temporary_dirs = {str(Path(report["temporary_project_ifc"]).parent)}
    temporary_dirs.update(str(Path(r["path"]).parent) for r in report.get("previous_temporary_attempts", []))
    temporary_records = [{"directory": d, "files": [pkg.record(p) for p in sorted(Path(d).rglob("*")) if p.is_file()],
                          "requires_user_confirmation": True} for d in sorted(temporary_dirs)]
    pilot.write(pilot.OUT / "cleanup-proposal.json", {"status": "proposal_only_no_deletion", "legacy_files": cleanup,
                "temporary_scene_directories": temporary_records,
                "keep": [str(pilot.SINGLE), "scene SVG/PNG", "validation records", "official-source/**", str(pilot.APPROVAL)],
                "batch_cleanup_authorized": False})
    pilot.write(pilot.OUT / "package-manifest.json", {
        "schema_version": 1, "profile_key": "wd02", "status": "pilot_validated_pending_visual_signoff",
        "persistent_ifc": pkg.record(pilot.SINGLE), "represented_physical_elements": 1,
        "representative_global_id": pilot.GUID, "unit_scale_to_m": result["unit_scale_to_m"],
        "product_placement_project_units": result["placement_matrix"], "coordinate_policy": "retain project coordinates; never center or rotate package",
        "source_kind": "geometry_derived_simplified_proxy", "source_label_zh": "基于原始高模几何生成的简化图纸表达",
        "approval_record": pkg.record(pilot.APPROVAL),
        "source_records": [r for r in protected if "/official-source/" in r["path"]],
        "scene_recipe": {"matching": "existing GlobalId + identical Body + identical placement + identical units",
                         "body_append": False, "repeat_attachment": "idempotent",
                         "save_target": "unique temporary project copy only", "drawing": "Bonsai Create Drawing"},
        "formal_ifc_sha256": pilot.FORMAL_HASH, "formal_authoritative_ifc_write_allowed": False,
        "scene_outputs": comparisons, "validation_record": "validation.json", "cleanup_proposal": "cleanup-proposal.json",
        "legacy_copy_retained": True, "batch_conversion_started": False})
    print(json.dumps({"independent_validation": result, "comparisons": comparisons}, ensure_ascii=False))


if __name__ == "__main__":
    main()
