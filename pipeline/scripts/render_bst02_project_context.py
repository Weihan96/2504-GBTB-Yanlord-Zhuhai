#!/usr/bin/env python3
"""Overlay BST02 geometry-derived views on its real project plan and side elevation."""

from falper_sorgente_linework import ROOT, load_json, sha256, write_json
from render_sis04_project_context import add_white_review_background, suppress_plan_elevation_markers, write_uncached_png_preview
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/bst02"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "3hA0vKpcn44u4Tsx4tqiUz"
shared.OVERLAY_ID_PREFIX = "bst02"
shared.GENERATOR = "pipeline/scripts/render_bst02_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-05-16-R09-NY.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.21
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 1),
    ("side", "side", shared.SIDE_SOURCE, "project-context-side-elevation.svg", (10.0, 7.0), 0),
)


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    plan = manifest["views"][0]
    plan_crop = ROOT / plan["review_crop"]
    plan_preview = ROOT / plan["review_preview"]
    suppress_plan_elevation_markers(plan_crop)
    add_white_review_background(plan_crop)
    write_uncached_png_preview(plan_crop, plan_preview)
    plan["review_crop_sha256"] = sha256(plan_crop)
    plan["review_preview_sha256"] = sha256(plan_preview)
    side = manifest["views"][1]
    side_crop = ROOT / side["review_crop"]
    side_preview = ROOT / side["review_preview"]
    add_white_review_background(side_crop)
    write_uncached_png_preview(side_crop, side_preview)
    side["review_crop_sha256"] = sha256(side_crop)
    side["review_preview_sha256"] = sha256(side_preview)
    manifest["semantic_view_mapping"] = {
        "plan": {"candidate_axes": [0, 1], "source": "drawings/Furniture Plan.svg", "rotate_quarter_turns": 1, "project_direction": "+Z"},
        "side": {"candidate_axes": [1, 2], "source": "drawings/elevations/native/EL-05-16-R09-NY.svg", "rotate_quarter_turns": 0, "project_direction": "-X"},
    }
    manifest["context_view_scope"] = {
        "included": ["plan", "side"],
        "excluded": ["front"],
        "reason": "The BST02 GlobalId occurs in Furniture Plan and the R09 NY native elevation. The R09 bbox is 577.889 x 250 mm and mechanically matches the candidate Y-Z side view; no native project front elevation exists, so one is not invented."
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "side": "none",
        "scope": "review crops only; full project SVGs retain all annotations",
        "walls_furniture_and_ifc_geometry_removed": False,
    }
    manifest["project_projection_note"] = {
        "plan_original_project_bbox_mm": [577.999987, 549.999988],
        "side_original_project_bbox_mm": [577.889239, 249.999994],
        "ifc_body_local_xyz_mm": [560.0, 577.889252, 250.0],
        "interpretation": "The plan is a nominal 550 x 578 mm project projection while the complete Body is 560 x 577.889 mm. The R09 NY elevation is the exact 577.889 x 250 mm Y-Z side projection. Both overlays retain uniform 1:50 scale and the Body is not stretched to the official 350 mm height.",
        "geometry_stretched": False,
        "official_height_discrepancy_mm": 100.0,
    }
    write_json(manifest_path, manifest)
