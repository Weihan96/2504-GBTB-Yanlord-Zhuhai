#!/usr/bin/env python3
"""Overlay BST01 geometry-derived views on its real project plan and two elevations."""

from falper_sorgente_linework import ROOT, load_json, sha256, write_json
from render_sis04_project_context import add_white_review_background, suppress_plan_elevation_markers, write_uncached_png_preview
import render_hima01_project_context as shared

shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/bst01"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "3eic1dzkn5heTIn4PhF37v"
shared.OVERLAY_ID_PREFIX = "bst01"
shared.GENERATOR = "pipeline/scripts/render_bst01_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-05-17-R09-NX.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-05-14-R09-PY.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.002
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 1),
    ("front", "front", shared.FRONT_SOURCE, "project-context-front-elevation.svg", (9.0, 7.0), 0),
    ("side", "side", shared.SIDE_SOURCE, "project-context-side-elevation.svg", (9.0, 7.0), 0),
)

if __name__ == "__main__":
    shared.main()
    path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(path)
    for record in manifest["views"]:
        crop = ROOT / record["review_crop"]
        preview = ROOT / record["review_preview"]
        if record["view"] == "plan": suppress_plan_elevation_markers(crop)
        add_white_review_background(crop)
        write_uncached_png_preview(crop, preview)
        record["review_crop_sha256"] = sha256(crop)
        record["review_preview_sha256"] = sha256(preview)
    manifest["semantic_view_mapping"] = {
        "plan": {"candidate_axes": [0, 1], "source": "drawings/Furniture Plan.svg", "rotate_quarter_turns": 1, "project_direction": "+Z"},
        "front": {"candidate_axes": [0, 2], "source": "drawings/elevations/native/EL-05-17-R09-NX.svg", "rotate_quarter_turns": 0, "project_direction": "-X world / local front after +90 degree placement"},
        "side": {"candidate_axes": [1, 2], "source": "drawings/elevations/native/EL-05-14-R09-PY.svg", "rotate_quarter_turns": 0, "project_direction": "+Y world / local side after +90 degree placement"}
    }
    manifest["context_view_scope"] = {"included": ["plan", "front", "side"], "excluded": [], "reason": "The BST01 GlobalId occurs in Furniture Plan and two native R09 elevations. The object placement is a +90 degree XY rotation, so NX maps local X-Z Front and PY maps local Y-Z Side."}
    manifest["review_annotation_suppression"] = {"plan": "official-elevation-anchor groups only", "front": "none", "side": "none", "scope": "review crops only; full project SVGs retain all annotations", "walls_furniture_and_ifc_geometry_removed": False}
    manifest["project_projection_note"] = {"project_plan_bbox_mm": [399.026099, 399.026099], "project_front_bbox_mm": [399.030905, 449.99999], "project_side_bbox_mm": [399.031897, 449.99999], "ifc_body_local_xyz_mm": [399.030914, 399.031906, 450.0], "official_nominal_mm": [420.0, 420.0, 450.0], "geometry_stretched": False, "opening_side_claimed": False}
    write_json(path, manifest)
