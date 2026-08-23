#!/usr/bin/env python3
"""Overlay one complete domestic-custom STREET Body proxy on the project sanitary plan."""

from falper_sorgente_linework import ROOT, load_json, sha256, write_json
from render_sis04_project_context import add_white_review_background, suppress_plan_elevation_markers, write_uncached_png_preview
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/street"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "1FgLPMw$5B4wBH2ySMkXE1"
shared.OVERLAY_ID_PREFIX = "street"
shared.GENERATOR = "pipeline/scripts/render_street_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Sanitary Plan.svg"
shared.ALIGNMENT_MODE = "centre"
shared.BBOX_TOLERANCE_SVG_UNITS = 2.5
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-sanitary-plan.svg", (10.0, 10.0), 1),
)


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    for record in manifest["views"]:
        crop = ROOT / record["review_crop"]
        preview = ROOT / record["review_preview"]
        suppress_plan_elevation_markers(crop)
        add_white_review_background(crop)
        write_uncached_png_preview(crop, preview)
        record["review_crop_sha256"] = sha256(crop)
        record["review_preview_sha256"] = sha256(preview)
        record["overlay"]["official_cad_acquired"] = True
        record["overlay"]["official_family_cad_archived_as_rejected_reference"] = True
        record["overlay"]["official_family_paths_used_as_project_representation"] = False
    manifest["official_cad_acquired"] = True
    manifest["official_family_cad_archived_as_rejected_reference"] = True
    manifest["official_family_paths_used_as_project_representation"] = False
    manifest["semantic_view_mapping"] = {
        "plan": {
            "candidate_axes": [0, 1],
            "source": "drawings/Sanitary Plan.svg",
            "rotate_quarter_turns": 1,
            "project_direction": "+Z"
        }
    }
    manifest["context_view_scope"] = {
        "included": ["plan"],
        "excluded": ["front", "side"],
        "reason": "The representative STREET GlobalId has one actual projection group in Sanitary Plan.svg and none in the native project elevation SVGs; project elevations are not invented."
    }
    manifest["project_projection_note"] = {
        "original_visible_projection_bbox_mm": [353.188457, 960.0],
        "complete_ifc_body_plan_bbox_mm": [1000.0, 470.0],
        "interpretation": "The original PLAN_VIEW projection is smaller than the complete MODEL_VIEW Body envelope. The review overlay preserves the complete Body at uniform project scale and aligns it to the original projection centre.",
        "geometry_stretched": False
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "scope": "derived review crop only; complete original project SVG remains unchanged",
        "walls_furniture_and_ifc_geometry_removed": False
    }
    manifest["review_preview_background"] = {
        "plan": "opaque white review-only paper background",
        "full_project_svg_remains_transparent": True,
        "drawing_geometry_changed": False
    }
    write_json(manifest_path, manifest)
