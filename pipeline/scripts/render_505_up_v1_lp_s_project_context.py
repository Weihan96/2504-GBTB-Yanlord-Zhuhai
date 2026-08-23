#!/usr/bin/env python3
"""Overlay the 505 UP project proxy on its real furniture-plan context."""

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from render_sis04_project_context import (
    add_white_review_background,
    suppress_plan_elevation_markers,
    write_uncached_png_preview,
)
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/505-up-v1-lp-s"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "19MpdkWqXC7uhUNhLQgrce"
shared.OVERLAY_ID_PREFIX = "molteni-505-up-v1-lp-s"
shared.GENERATOR = "pipeline/scripts/render_505_up_v1_lp_s_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.1
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
            "source": relative(shared.PLAN_SOURCE),
            "rotate_quarter_turns": 0,
            "project_direction": "+Z",
        }
    }
    manifest["context_view_scope"] = {
        "included": ["plan"],
        "excluded": ["front", "side"],
        "reason": "The representative occurs in Furniture Plan.svg, but no native project elevation SVG contains this IFC GlobalId; no elevation context was invented.",
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "scope": "review crop only; the complete project SVG retains its annotations",
        "walls_furniture_and_ifc_geometry_removed": False,
    }
    manifest["review_preview_background"] = {
        "plan": "opaque white review-only paper background",
        "full_project_svg_remains_transparent": True,
        "drawing_geometry_changed": False,
    }
    manifest["project_projection_note"] = {
        "original_project_bbox_mm": [2911.999935, 436.499441],
        "ifc_body_plan_bbox_mm": [2911.714112, 436.499451],
        "interpretation": "The complete geometry-derived proxy is overlaid at the drawing's native 1:50 scale and aligned to the original 505 UP projection centre.",
        "geometry_stretched": False,
    }
    write_json(manifest_path, manifest)
