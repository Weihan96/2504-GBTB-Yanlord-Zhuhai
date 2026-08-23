#!/usr/bin/env python3
"""Overlay the Geberit 154.154.00.1.F Body proxy on real project drawings."""

from falper_sorgente_linework import ROOT, load_json, sha256, write_json
from render_sis04_project_context import add_white_review_background, write_uncached_png_preview
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-154-154-00-1-f"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "2jjNIn9gHBYwNSWwlM5T_i"
shared.OVERLAY_ID_PREFIX = "geberit-154-154-00-1-f"
shared.GENERATOR = "pipeline/scripts/render_geberit_154_154_00_1_f_project_context.py"
shared.ALIGNMENT_MODE = "centre"
shared.BBOX_TOLERANCE_SVG_UNITS = 2.0
shared.PLAN_SOURCE = ROOT / "drawings/Sanitary Plan.svg"
shared.FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-06-20-R12-NY.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-06-21-R12-PX.svg"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-sanitary-plan.svg", (8.0, 8.0), 0),
    ("front", "front", shared.FRONT_SOURCE, "project-context-r12-front-elevation.svg", (5.0, 4.0), 0),
    ("side", "side", shared.SIDE_SOURCE, "project-context-r12-side-elevation.svg", (5.0, 4.0), 0),
)


def suppress_plan_elevation_markers(path):
    """Remove oversized elevation symbols from the review crop, not the project SVG."""
    content = path.read_text(encoding="utf-8")
    lines = content.splitlines()
    filtered = [line for line in lines if 'class="official-elevation-anchor"' not in line]
    if len(filtered) == len(lines):
        raise RuntimeError("expected official elevation anchors in the plan review crop")
    path.write_text("\n".join(filtered) + "\n", encoding="utf-8")


if __name__ == "__main__":
    shared.main()
    path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(path)
    for record in manifest["views"]:
        crop = ROOT / record["review_crop"]
        preview = ROOT / record["review_preview"]
        if record["view"] == "plan":
            suppress_plan_elevation_markers(crop)
        add_white_review_background(crop)
        write_uncached_png_preview(crop, preview)
        record["review_crop_sha256"] = sha256(crop)
        record["review_preview_sha256"] = sha256(preview)
    manifest["semantic_view_mapping"] = {
        "plan": {"candidate_axes": [0, 1], "source": "drawings/Sanitary Plan.svg", "project_direction": "+Z"},
        "front": {"candidate_axes": [0, 2], "source": "drawings/elevations/native/EL-06-20-R12-NY.svg", "project_direction": "-Y"},
        "side": {"candidate_axes": [1, 2], "source": "drawings/elevations/native/EL-06-21-R12-PX.svg", "project_direction": "+X"},
    }
    manifest["project_projection_note"] = {
        "ifc_body_local_xyz_mm": [88.0, 210.0, 104.999998],
        "official_standalone_component_dimensions_available": False,
        "native_plan_projection_bbox_mm": [88.0, 210.0],
        "native_front_projection_bbox_mm": [63.100047, 57.55335],
        "native_side_projection_bbox_mm": [210.0, 18.3831],
        "body_projection_bbox_mm": {"plan": [88.0, 210.0], "front": [88.0, 104.999998], "side": [210.0, 104.999998]},
        "interpretation": "The plan projection matches the actual Body. Both native R12 elevations clip concealed flange geometry; the review overlay restores the complete Body at the fixed 1:50 scale and aligns it to the original projection centre.",
        "geometry_stretched": False,
        "native_project_plan_and_elevations_found": True,
        "component_boundary": "The black overlay represents only the project .F Body; the complete official parent assembly is not substituted."
    }
    manifest["review_preview_background"] = {
        "all_views": "opaque white review-only paper background",
        "full_project_svg_remains_transparent": True,
        "drawing_geometry_changed": False,
    }
    manifest["review_annotation_suppression"] = {
        "plan": ["official-elevation-anchor"],
        "scope": "review crop only; full project SVG retains all elevation markers",
        "walls_furniture_and_ifc_geometry_removed": False,
    }
    write_json(path, manifest)
