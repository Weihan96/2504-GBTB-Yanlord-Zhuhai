#!/usr/bin/env python3
"""Overlay Miami Soft F03 on its complete project plan and R20 side elevation."""

import xml.etree.ElementTree as ET

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from render_sis04_project_context import add_white_review_background, suppress_plan_elevation_markers, write_uncached_png_preview
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/miamisoft-f03"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "1y3lLYK$H9fPbuh1R1vAKX"
shared.OVERLAY_ID_PREFIX = "miamisoft-f03"
shared.GENERATOR = "pipeline/scripts/render_miamisoft_f03_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-01-02-R20-PY.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 1.3
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 1),
    ("side", "side", shared.SIDE_SOURCE, "project-context-r20-side-elevation.svg", (8.0, 6.0), 0),
)


def projection_group_points(source):
    groups = []
    for element in ET.parse(source).getroot().iter():
        if not element.tag.endswith("g") or shared.guid(element) != shared.GLOBAL_ID or "projection" not in element.attrib.get("class", "").split():
            continue
        points = [point for child in element.iter() if child.tag.endswith("path") for point in shared.path_points(child)]
        if points:
            groups.append(points)
    if len(groups) != 1:
        raise RuntimeError(f"expected one complete Miami Soft F03 projection in {source}, found {len(groups)}")
    return groups[0]


shared.group_points = projection_group_points
original_add_overlay = shared.add_overlay


def add_overlay_with_opaque_mask(source, target, view, path):
    original_add_overlay(source, target, view, path)
    content = target.read_text(encoding="utf-8").replace('class="geometry-derived-proxy-mask"', 'class="geometry-derived-proxy-mask" fill="#ffffff" fill-rule="evenodd"', 1)
    target.write_text(content, encoding="utf-8")


shared.add_overlay = add_overlay_with_opaque_mask


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    for record in manifest["views"]:
        crop = ROOT / record["review_crop"]
        preview = ROOT / record["review_preview"]
        if record["view"] == "plan":
            suppress_plan_elevation_markers(crop)
        add_white_review_background(crop)
        write_uncached_png_preview(crop, preview)
        record["review_crop_sha256"] = sha256(crop)
        record["review_preview_sha256"] = sha256(preview)
    side = next(record for record in manifest["views"] if record["view"] == "side")
    side_delta = max(side["overlay"]["fit"]["bbox_absolute_delta_svg_units"])
    manifest["semantic_view_mapping"] = {
        "plan": {"candidate_axes": [0, 1], "source": relative(shared.PLAN_SOURCE), "rotate_quarter_turns": 1, "project_direction": "+Z"},
        "side": {"candidate_axes": [1, 2], "source": relative(shared.SIDE_SOURCE), "rotate_quarter_turns": 0, "project_direction": "+Y", "complete_project_projection_group_only": True},
    }
    manifest["context_view_scope"] = {"included": ["plan", "side"], "excluded": ["front"], "reason": "The R20 +Y drawing contains a complete YZ side projection. The +X product projection is height-occluded and is not presented as a mechanically complete front context."}
    manifest["mechanical_side_fit_tolerance_svg_units"] = 0.08
    manifest["mechanical_side_fit_pass"] = side_delta <= 0.08
    manifest["review_annotation_suppression"] = {"plan": "official-elevation-anchor groups only", "side": "none", "scope": "review crops only; complete project SVGs retain their annotations", "walls_furniture_and_ifc_geometry_removed": False}
    manifest["review_preview_background"] = {"plan": "opaque white review-only paper background", "side": "opaque white review-only paper background", "full_project_svgs_remain_transparent": True, "drawing_geometry_changed": False}
    manifest["project_projection_note"] = {"official_nominal_mm": [1700.0, 1080.0, 400.0], "ifc_body_local_xyz_mm": [1747.539001, 1103.159912, 464.938965], "interpretation": "The complete geometry-derived proxy is overlaid at the drawing's native 1:50 scale and centre-aligned; no axis is stretched to force the nominal upholstered envelope.", "geometry_stretched": False}
    manifest["pass"] = manifest["pass"] and manifest["mechanical_side_fit_pass"]
    write_json(manifest_path, manifest)
