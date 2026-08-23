#!/usr/bin/env python3
"""Overlay the 154.154.00.1 geometry-derived proxy on BATHM project drawings."""

import re

from falper_sorgente_linework import ROOT, load_json, sha256, write_json
from render_sis04_project_context import add_white_review_background, write_uncached_png_preview
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-154-154-00"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "14EazrLgP8whYZWY_yCuKy"
shared.OVERLAY_ID_PREFIX = "geberit-154-154-00"
shared.GENERATOR = "pipeline/scripts/render_geberit_154_154_00_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Sanitary Plan-P202-candidate.svg"
shared.FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-06-20-R12-NY.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-06-21-R12-PX.svg"
shared.ALIGNMENT_MODE = "centre"
shared.BBOX_TOLERANCE_SVG_UNITS = 1.0
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-sanitary-plan.svg", (8.0, 8.0), 0),
    ("front", "front", shared.FRONT_SOURCE, "project-context-r12-front-elevation.svg", (5.0, 4.0), 0),
    ("side", "side", shared.SIDE_SOURCE, "project-context-r12-side-elevation.svg", (5.0, 4.0), 0),
)


def suppress_plan_review_markers(path):
    """Hide review-obscuring annotations without removing project geometry."""
    content = path.read_text(encoding="utf-8")
    content, placement_count = re.subn(
        r'<g class="p202-location-marker[^"]*"[^>]*>.*?</g>',
        "",
        content,
        flags=re.DOTALL,
    )
    lines = [line for line in content.splitlines() if 'class="official-elevation-anchor"' not in line]
    if placement_count == 0 or len(lines) == len(content.splitlines()):
        raise RuntimeError("expected P-202 placement and official elevation markers in plan review crop")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    for record in manifest["views"]:
        crop = ROOT / record["review_crop"]
        preview = ROOT / record["review_preview"]
        if record["view"] == "plan":
            suppress_plan_review_markers(crop)
        add_white_review_background(crop)
        write_uncached_png_preview(crop, preview)
        record["review_crop_sha256"] = sha256(crop)
        record["review_preview_sha256"] = sha256(preview)
    manifest["semantic_view_mapping"] = {
        "plan": {"candidate_axes": [0, 1], "source": "drawings/Sanitary Plan-P202-candidate.svg", "project_direction": "+Z"},
        "front": {"candidate_axes": [0, 2], "source": "drawings/elevations/native/EL-06-20-R12-NY.svg", "project_direction": "-Y"},
        "side": {"candidate_axes": [1, 2], "source": "drawings/elevations/native/EL-06-21-R12-PX.svg", "project_direction": "+X"}
    }
    manifest["project_projection_note"] = {
        "front_original_visible_bbox_width_mm": 320.70652,
        "ifc_body_front_width_mm": 361.700024,
        "interpretation": "The original -Y project elevation clips concealed installation-set geometry; the review overlay preserves the complete IFC Body projection at uniform project scale and aligns it to the original projection centre.",
        "geometry_stretched": False
    }
    manifest["review_preview_background"] = {
        "all_views": "opaque white review-only HTML paper background",
        "full_project_svg_remains_transparent": True,
        "drawing_geometry_changed": False
    }
    manifest["review_annotation_suppression"] = {
        "plan": ["p202-location-marker", "official-elevation-anchor"],
        "scope": "review crop only; full project SVG retains all annotations",
        "walls_furniture_and_ifc_geometry_removed": False
    }
    write_json(manifest_path, manifest)
