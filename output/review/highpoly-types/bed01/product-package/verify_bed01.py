"""Disk-only verification of BED01 package and regenerated approved blue lines."""
from pathlib import Path
import json
import subprocess
import xml.etree.ElementTree as ET
import ifcopenshell
import ifcopenshell.validate
from shapely.geometry import LineString
from shapely.ops import unary_union
import migrate_bed01 as task

pkg = task.pkg


def line_geometry(path, guid):
    root = ET.parse(path).getroot()
    lines = [e for e in root.iter() if guid in e.get("class", "") and e.tag.endswith("line")]
    assert lines
    assert not [e for e in root.iter() if e.tag.endswith("image")]
    return unary_union([LineString([(float(e.get("x1")), float(e.get("y1"))),
                                    (float(e.get("x2")), float(e.get("y2")))]) for e in lines])


def main():
    report = json.loads(task.REPORT.read_text())
    assert report["pure_saved_reloaded"] and report["scene_saved_reloaded"]
    model = ifcopenshell.open(str(task.SINGLE))
    state = task.check_pure(model)
    assert state == report["pure_pre_bonsai"]
    logger = ifcopenshell.validate.json_logger()
    ifcopenshell.validate.validate(model, logger)
    assert not logger.statements, logger.statements
    assert all(pkg.sha256(r["path"]) == r["sha256"] for r in report["protected_files"])
    for dependency in report["necessary_external_style_dependencies"]:
        assert pkg.sha256(dependency["source"]["path"]) == dependency["source"]["sha256"]
        assert pkg.sha256(dependency["packaged"]["path"]) == dependency["source"]["sha256"]
    comparisons = []
    for v in report["views"]:
        svg = task.OUT / f"BED01-SCENE-{v['view'].upper()}.svg"
        old = line_geometry(task.ROOT / v["old_svg"]["path"], v["annotation_guid"])
        new = line_geometry(svg, v["annotation_guid"])
        error = old.segmentize(.1).hausdorff_distance(new.segmentize(.1)) * 25
        assert error < .01, {"view": v["view"], "world_error_mm": error}
        png = svg.with_suffix(".png")
        subprocess.run(["/Applications/Inkscape.app/Contents/MacOS/inkscape", str(svg),
                        "--export-area-page", "--export-background=white", "--export-background-opacity=1",
                        "--export-width=1400", f"--export-filename={png}"], check=True, capture_output=True)
        comparisons.append({"view": v["view"], "max_linework_hausdorff_error_world_mm": error,
                            "tolerance_mm": .01, "svg": pkg.record(svg), "preview": pkg.record(png)})
    report["independent_validation"] = {"pure_state": state, "schema_errors": 0,
                                        "approved_scene_comparison": comparisons,
                                        "protected_files_unchanged": True,
                                        "formal_sha256": pkg.sha256(task.FORMAL)}
    assert pkg.sha256(task.FORMAL) == task.FORMAL_HASH
    report["verdict"] = "validated_pending_visual_check"
    task.write(task.REPORT, report)
    source_record = task.PRODUCT / "official-source/source-access-record.json"
    source_data = json.loads(source_record.read_text())
    task.write(task.OUT / "manifest.json", {
        "schema_version": 1, "profile_key": "bed01", "display_name": "Baxter Casablanca 180 / BED01",
        "status": "validated_pending_visual_check", "target_global_id": task.GUID,
        "persistent_ifc": pkg.record(task.SINGLE), "source_kind": "native_dwg_review_reference",
        "source_label_zh": "Baxter Casablanca 官方独立 2D DWG 原生蓝线",
        "scope": "approved 180x200 family linework; not an exact project Body fit; no geometry rescaling",
        "source_dwg": pkg.record(task.ROOT / "output/review/highpoly-types/bed01/official-source/2d-3d-download/CASABLANCA/2D/Casablanca_Letto.dwg"),
        "official_download_url": source_data["official_product_2d_3d_cad"]["download_url"],
        "source_access_record": pkg.record(source_record), "approval_record": pkg.record(task.APPROVAL),
        "geometry": {"physical_products": 1, "original_body_preserved": True, "approved_views": ["plan", "front", "side"],
                     "placement_matrix_project_units": state["placement"], "units_to_m": state["unit_scale_to_m"],
                     "annotations": 0, "drawing_cameras": 0, "spatial_geometry": 0},
        "scene_recipe": pkg.record(task.RECIPE), "scene_outputs": comparisons,
        "external_style_dependencies": report["necessary_external_style_dependencies"],
        "validation": "validation.json", "legacy_source_retained": True, "formal_write_allowed": False,
        "formal_sha256": task.FORMAL_HASH,
        "known_limitations": ["Approved official family CAD does not exactly fit the project's modified highpoly Body.",
                              "Missing OD_Textures/Materials.blend was inherited from formal baseline; 2D linework succeeds."]})
    task.write(task.OUT / "cleanup-proposal.json", {
        "status": "proposal_only_no_deletion", "requires_user_confirmation": True,
        "legacy_files": [pkg.record(task.OLD)],
        "temporary_directories": sorted(set([report["temporary_project_directory"], *report.get("previous_temporary_attempts", [])])),
        "keep": [str(task.SINGLE), str(task.RECIPE), "scene SVG/PNG", "validation and source records", "official-source/**"]})
    task.write(task.OUT / "handoff.json", {
        "product": "bed01", "status": "validated_pending_visual_check", "package": "manifest.json",
        "ifc": "BED01-product.ifc", "scene_recipe": "scene-recipe.json", "validation": "validation.json",
        "new_files_only": True, "staged_baseline_modified": False, "legacy_deleted": False,
        "formal_modified": False, "bridge_port": 9885})
    print(json.dumps({"state": state, "comparisons": comparisons}, ensure_ascii=False))


if __name__ == "__main__":
    main()
