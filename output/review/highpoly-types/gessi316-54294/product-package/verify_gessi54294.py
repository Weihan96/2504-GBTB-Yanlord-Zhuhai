"""Disk-only verification of GESSI54294 package and regenerated approved blue lines."""
from pathlib import Path
import json
import subprocess
import xml.etree.ElementTree as ET
import copy
import re
import ifcopenshell
import ifcopenshell.validate
from shapely.geometry import LineString
from shapely.ops import unary_union
import migrate_gessi54294 as task

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
        svg = task.OUT / f"GESSI54294-SCENE-{v['view'].upper()}.svg"
        old = line_geometry(task.ROOT / v["old_svg"]["path"], v["annotation_guid"])
        new = line_geometry(svg, v["annotation_guid"])
        error = old.segmentize(.1).hausdorff_distance(new.segmentize(.1)) * 25
        assert error < .01, {"view": v["view"], "world_error_mm": error}
        if v["view"] == "side":
            root = ET.parse(svg).getroot()
            groups = [e for e in root.iter() if any(k.endswith("}guid") and val == "3WekjaeUn1qfL6oR1KYD_2" for k, val in e.attrib.items())]
            assert groups
            wall_x = []
            for group in groups:
                for e in group.iter():
                    if e.tag.endswith("path"):
                        nums = [float(x) for x in re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", e.get("d", ""))]
                        wall_x.extend(nums[::2])
                    elif e.tag.rsplit("}", 1)[-1] == "line":
                        wall_x.extend([float(e.get("x1")), float(e.get("x2"))])
                    elif e.tag.endswith("polygon") or e.tag.endswith("polyline"):
                        wall_x.extend([float(x) for x in re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", e.get("points", ""))][::2])
            assert wall_x
            anchor_x = [float(e.get(k)) for e in root.iter() if e.tag.endswith("line") and v["annotation_guid"] in e.get("class", "") for k in ("x1", "x2")]
            wall_face = max(wall_x)
            residual = (min(anchor_x, key=lambda x: abs(x-wall_face)) - wall_face) * 25
            assert abs(residual) < .2, residual
            excluded = {"0nlHIdaVrV7xIKsYXraYCa", "04rs0EDjn2EvxytEQSxWRB", "26PwmC1AX6ZROC65DohUg5"}
            assert not [e for e in root.iter() if any(k.endswith("}guid") and val in excluded for k, val in e.attrib.items())]
            report["side_finished_wall_recheck"] = {"svg_residual_mm": residual, "tolerance_mm": .2,
                "occluding_substrate_absent": True, "annotation_translation_mm": [0, 0, 0], "pass": True}
        png = svg.with_suffix(".png")
        subprocess.run(["/Applications/Inkscape.app/Contents/MacOS/inkscape", str(svg),
                        "--export-area-page", "--export-background=white", "--export-background-opacity=1",
                        "--export-width=1400", f"--export-filename={png}"], check=True, capture_output=True)
        comparisons.append({"view": v["view"], "max_linework_hausdorff_error_world_mm": error,
                            "tolerance_mm": .01, "svg": pkg.record(svg), "preview": pkg.record(png)})
        # Same persisted linework, isolated and re-framed for the library only.
        thumb = ET.Element("{http://www.w3.org/2000/svg}svg", width="1000", height="760")
        minx, miny, maxx, maxy = new.bounds
        margin = max(maxx-minx, maxy-miny) * .08
        thumb.set("viewBox", f"{minx-margin} {miny-margin} {maxx-minx+2*margin} {maxy-miny+2*margin}")
        for e in ET.parse(svg).getroot().iter():
            if e.tag.endswith("line") and v["annotation_guid"] in e.get("class", ""):
                preview_line = copy.deepcopy(e)
                width = max(maxx-minx, maxy-miny) * .002
                preview_line.set("style", f"stroke:#1677c8;stroke-width:{width};stroke-linecap:round;stroke-linejoin:round;fill:none")
                thumb.append(preview_line)
        single_svg = task.OUT / f"GESSI54294-SINGLE-{v['view'].upper()}.svg"
        ET.ElementTree(thumb).write(single_svg, encoding="utf-8", xml_declaration=True)
        single_png = single_svg.with_suffix(".png")
        subprocess.run(["/Applications/Inkscape.app/Contents/MacOS/inkscape", str(single_svg),
                        "--export-area-page", "--export-background=white", "--export-background-opacity=1",
                        "--export-width=1000", f"--export-filename={single_png}"], check=True, capture_output=True)
        comparisons[-1]["single_product_thumbnail"] = pkg.record(single_png)
        comparisons[-1]["single_product_svg"] = pkg.record(single_svg)
        comparisons[-1]["single_product_display_style"] = "Thumbnail-only 2px-equivalent strokes and round caps; source coordinates and all scene SVG styles unchanged"
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
        "schema_version": 1, "profile_key": "gessi316-54294", "display_name": "Gessi316 54294",
        "status": "validated_pending_visual_check", "target_global_id": task.GUID,
        "persistent_ifc": pkg.record(task.SINGLE), "source_kind": "native_dwg_review_simplification",
        "source_label_zh": "基于官方 Gessi 54294 原生 DWG 轮廓的去纹审核简化表达",
        "scope": "Approved main bathroom de-textured Plan/Front/Side; zero annotation translation; Side finished-wall filter retained",
        "source_dwg": pkg.record(task.PRODUCT / "official-source/GPF5429400000G000_3.dwg"),
        "official_download_url": source_data["official_product_cad"]["native_dwg_zip"]["url"],
        "source_access_record": pkg.record(source_record), "approval_record": pkg.record(task.APPROVAL),
        "geometry": {"physical_products": 1, "original_body_preserved": True, "approved_views": ["plan", "front", "side"],
                     "placement_matrix_project_units": state["placement"], "units_to_m": state["unit_scale_to_m"],
                     "annotations": 0, "drawing_cameras": 0, "spatial_geometry": 0},
        "scene_recipe": pkg.record(task.RECIPE), "scene_outputs": comparisons,
        "external_style_dependencies": report["necessary_external_style_dependencies"],
        "validation": "validation.json", "legacy_source_retained": True, "formal_write_allowed": False,
        "formal_sha256": task.FORMAL_HASH,
        "known_limitations": ["Blue linework is a de-textured simplification, not untouched official CAD; original DWG and its source record remain preserved."]})
    task.write(task.OUT / "cleanup-proposal.json", {
        "status": "proposal_only_no_deletion", "requires_user_confirmation": True,
        "legacy_files": report["legacy_full_projects"],
        "temporary_directories": sorted(set([report["temporary_project_directory"], *report.get("previous_temporary_attempts", [])])),
        "keep": [str(task.SINGLE), str(task.RECIPE), "scene SVG/PNG", "validation and source records", "official-source/**"]})
    task.write(task.OUT / "handoff.json", {
        "product": "gessi54294", "status": "validated_pending_visual_check", "package": "manifest.json",
        "ifc": "GESSI54294-product.ifc", "scene_recipe": "scene-recipe.json", "validation": "validation.json",
        "new_files_only": True, "staged_baseline_modified": False, "legacy_deleted": False,
        "formal_modified": False, "bridge_port": 9889})
    print(json.dumps({"state": state, "comparisons": comparisons}, ensure_ascii=False))


if __name__ == "__main__":
    main()
