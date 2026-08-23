#!/usr/bin/env python3
"""Overlay the TAB02 geometry-derived plan proxy on the project furniture plan."""

from falper_sorgente_linework import ROOT, load_json, sha256, write_json
from render_sis04_project_context import (
    add_white_review_background,
    suppress_plan_elevation_markers,
    write_uncached_png_preview,
)
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/tab02"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "2eX84IyLr8_e34nuOoRPLQ"
shared.OVERLAY_ID_PREFIX = "tab02"
shared.GENERATOR = "pipeline/scripts/render_tab02_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.02
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 0),
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
            "rotate_quarter_turns": 0,
            "project_direction": "+Z",
        }
    }
    manifest["context_view_scope"] = {
        "included": ["plan"],
        "excluded": ["front", "side"],
        "reason": "The TAB02 representative has an actual projection group in Furniture Plan.svg and none in the native project elevation SVGs; project elevations are not invented.",
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "scope": "review crop only; full project SVG retains all annotations",
        "walls_furniture_and_ifc_geometry_removed": False,
    }
    manifest["review_preview_background"] = {
        "plan": "opaque white review-only HTML paper background",
        "full_project_svg_remains_transparent": True,
        "drawing_geometry_changed": False,
    }
    manifest["project_projection_note"] = {
        "original_project_bbox_mm": [500.0, 500.0],
        "ifc_body_plan_bbox_mm": [500.0, 500.0],
        "interpretation": "The project plan and complete MODEL_VIEW Body share the exact official 500 x 500 mm footprint. The overlay keeps uniform 1:50 scale and aligns by centre without stretching.",
        "geometry_stretched": False,
    }
    write_json(manifest_path, manifest)
