#!/usr/bin/env python3
"""Overlay the BST03 geometry-derived plan proxy on the project furniture plan."""

from falper_sorgente_linework import ROOT, load_json, sha256, write_json
from render_sis04_project_context import (
    add_white_review_background,
    suppress_plan_elevation_markers,
    write_uncached_png_preview,
)
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/bst03"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "1HZoxe$df4cBb4UXXH5J2S"
shared.OVERLAY_ID_PREFIX = "bst03"
shared.GENERATOR = "pipeline/scripts/render_bst03_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.3
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 1),
)


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    plan = manifest["views"][0]
    crop = ROOT / plan["review_crop"]
    preview = ROOT / plan["review_preview"]
    suppress_plan_elevation_markers(crop)
    add_white_review_background(crop)
    write_uncached_png_preview(crop, preview)
    plan["review_crop_sha256"] = sha256(crop)
    plan["review_preview_sha256"] = sha256(preview)
    manifest["semantic_view_mapping"] = {
        "plan": {
            "candidate_axes": [0, 1],
            "source": "drawings/Furniture Plan.svg",
            "rotate_quarter_turns": 1,
            "project_direction": "+Z"
        }
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "scope": "review crop only; full project SVG retains all annotations",
        "walls_furniture_and_ifc_geometry_removed": False
    }
    manifest["review_preview_background"] = {
        "plan": "opaque white review-only HTML paper background",
        "full_project_svg_remains_transparent": True,
        "drawing_geometry_changed": False
    }
    manifest["project_projection_note"] = {
        "original_project_bbox_mm": [450.0, 450.0],
        "ifc_body_plan_bbox_mm": [462.460968, 462.47998],
        "interpretation": "The project plan uses the nominal 450 mm catalogue footprint. The overlay preserves the complete high-poly Body envelope at uniform 1:50 scale and aligns it to the original projection centre.",
        "geometry_stretched": False
    }
    write_json(manifest_path, manifest)
