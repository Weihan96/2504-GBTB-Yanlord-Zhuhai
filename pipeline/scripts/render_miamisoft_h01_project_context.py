#!/usr/bin/env python3
"""Overlay Miami Soft H01 on its complete project plan and R20 side elevation."""

import re
import xml.etree.ElementTree as ET

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from render_sis04_project_context import add_white_review_background, suppress_plan_elevation_markers, write_uncached_png_preview
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/miamisoft-h01"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "0B8rO47U14sgJ_ydZxqvs6"
shared.OVERLAY_ID_PREFIX = "miamisoft-h01"
shared.GENERATOR = "pipeline/scripts/render_miamisoft_h01_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-01-02-R20-PY.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 3.1
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 0),
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
        raise RuntimeError(f"expected one complete Miami Soft H01 projection in {source}, found {len(groups)}")
    return groups[0]


shared.group_points = projection_group_points
original_add_overlay = shared.add_overlay


def add_overlay_with_opaque_mask(source, target, view, path):
    original_add_overlay(source, target, view, path)
    content = target.read_text(encoding="utf-8")
    projection_pattern = re.compile(
        rf'<g(?=[^>]*\bifc:guid="{re.escape(shared.GLOBAL_ID)}")(?=[^>]*\bclass="[^"]*\bprojection\b[^"]*")'
    )
    content, hidden_count = projection_pattern.subn(
        '<g style="display:none" data-review-replaced-by="geometry-derived-body-proxy"',
        content,
        count=1,
    )
    if hidden_count != 1:
        raise RuntimeError(f"expected one representative H01 projection group to replace in {source}, found {hidden_count}")
    content = content.replace(
        'class="geometry-derived-proxy-mask"',
        'class="geometry-derived-proxy-mask" fill="#ffffff" fill-rule="evenodd"',
        1,
    )
    target.write_text(content, encoding="utf-8")


shared.add_overlay = add_overlay_with_opaque_mask


def suppress_noninteger_diagnostic_highlights(path):
    """Remove generated residual-error boxes from the review crop only."""
    content = path.read_text(encoding="utf-8")
    content, count = re.subn(
        r'<g id="noninteger-highlights".*?</g>',
        "",
        content,
        count=1,
        flags=re.DOTALL,
    )
    if count != 1:
        raise RuntimeError(f"expected noninteger diagnostic highlights in {path}")
    path.write_text(content, encoding="utf-8")


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    for record in manifest["views"]:
        crop = ROOT / record["review_crop"]
        preview = ROOT / record["review_preview"]
        if record["view"] == "plan":
            suppress_plan_elevation_markers(crop)
        else:
            suppress_noninteger_diagnostic_highlights(crop)
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
    manifest["mechanical_side_fit_tolerance_svg_units"] = 0.15
    manifest["mechanical_side_fit_pass"] = side_delta <= 0.15
    plan = next(record for record in manifest["views"] if record["view"] == "plan")
    plan_delta = max(plan["overlay"]["fit"]["bbox_absolute_delta_svg_units"])
    manifest["plan_source_representation_comparison"] = {
        "source_representation": "existing Bonsai PLAN_VIEW projection group",
        "replacement_representation": "actual MODEL_VIEW Body projection",
        "maximum_bbox_delta_svg_units": plan_delta,
        "maximum_bbox_delta_mm_at_1_50": round(plan_delta / shared.PROJECT_SCALE_SVG_UNITS_PER_MM, 6),
        "representative_source_group_hidden_in_review_copy": True,
        "other_project_elements_retained": True,
        "geometry_stretched_to_match_source": False,
        "reason": "The existing PLAN_VIEW H01 group is taller than the actual Body projection. The review copy replaces only the representative group so the contextual product matches the candidate SVG without scaling distortion.",
    }
    manifest["review_annotation_suppression"] = {"plan": "official-elevation-anchor groups only", "side": "noninteger residual-error diagnostic highlights only", "scope": "review crops only; complete project SVGs retain their annotations", "walls_furniture_and_ifc_geometry_removed": False}
    manifest["review_product_projection_replacement"] = {"representative_global_id": shared.GLOBAL_ID, "source_projection_groups_hidden_per_view": 1, "replacement_proxy_top_layer": True, "second_type_instance_retained": "0b9rDAHDD2qAIQqZ8m9INg"}
    manifest["review_preview_background"] = {"plan": "opaque white review-only paper background", "side": "opaque white review-only paper background", "full_project_svgs_remain_transparent": True, "drawing_geometry_changed": False}
    manifest["project_projection_note"] = {"official_nominal_flexible_face_width_height_mm": [800.0, 800.0], "ifc_body_local_xyz_mm": [809.757111, 564.975739, 403.963976], "interpretation": "The complete geometry-derived proxy is overlaid at the drawing's native 1:50 scale and centre-aligned. The installed flexible cushion remains leaning and compressed; no axis is stretched or flattened to force the catalogue 80 x 80 cm square.", "geometry_stretched": False}
    manifest["pass"] = manifest["pass"] and manifest["mechanical_side_fit_pass"]
    write_json(manifest_path, manifest)
