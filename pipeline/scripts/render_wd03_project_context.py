#!/usr/bin/env python3
"""Overlay WD03 on its real furniture plan and R14 side-oriented elevation."""

import re
import xml.etree.ElementTree as ET

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from render_sis04_project_context import (
    add_white_review_background,
    suppress_plan_elevation_markers,
    write_uncached_png_preview,
)
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/wd03"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "3cmikd9MTB$egM5KQaNgUf"
shared.OVERLAY_ID_PREFIX = "poliform-senzafine-wd03"
shared.GENERATOR = "pipeline/scripts/render_wd03_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-01-01-R14-PY.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.4
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 0),
    ("side", "side", shared.SIDE_SOURCE, "project-context-side-elevation.svg", (10.0, 8.0), 0),
)


original_group_points = shared.group_points
original_add_overlay = shared.add_overlay


def group_points(source):
    if source != shared.SIDE_SOURCE:
        return original_group_points(source)
    matches = []
    for element in ET.parse(source).getroot().iter():
        if not element.tag.endswith("rect") or element.attrib.get("data-guid") != shared.GLOBAL_ID:
            continue
        x = float(element.attrib["x"])
        y = float(element.attrib["y"])
        width = float(element.attrib["width"])
        height = float(element.attrib["height"])
        matches.append([(x, y), (x + width, y), (x + width, y + height), (x, y + height)])
    if len(matches) != 1:
        raise RuntimeError(f"expected one WD03 diagnostic bbox in {source}, found {len(matches)}")
    return matches[0]


shared.group_points = group_points


def add_overlay(source, target, view, path, overlay_bbox):
    original_add_overlay(source, target, view, path, overlay_bbox)
    content = target.read_text(encoding="utf-8")
    content = re.sub(
        r'<g id="noninteger-highlights"[^>]*>.*?</g>',
        "",
        content,
        flags=re.DOTALL,
    )
    target.write_text(content, encoding="utf-8")


shared.add_overlay = add_overlay


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    for view in manifest["views"]:
        crop = ROOT / view["review_crop"]
        preview = ROOT / view["review_preview"]
        if view["view"] == "plan":
            suppress_plan_elevation_markers(crop)
        add_white_review_background(crop)
        write_uncached_png_preview(crop, preview)
        view["review_crop_sha256"] = sha256(crop)
        view["review_preview_sha256"] = sha256(preview)
    manifest["semantic_view_mapping"] = {
        "plan": {
            "candidate_axes": [0, 1],
            "source": relative(shared.PLAN_SOURCE),
            "rotate_quarter_turns": 0,
            "project_direction": "+Z",
        },
        "side": {
            "candidate_axes": [1, 2],
            "source": relative(shared.SIDE_SOURCE),
            "rotate_quarter_turns": 0,
            "project_direction": "R14 native elevation side-oriented project diagnostic bbox",
            "target_geometry_kind": "same-GlobalId diagnostic bounding rectangle; original project projection group absent",
        },
    }
    manifest["context_view_scope"] = {
        "included": ["plan", "side"],
        "excluded": ["front"],
        "reason": "Furniture Plan.svg contains the actual WD03 projection. EL-01-01-R14-PY.svg retains only a same-GlobalId diagnostic bbox whose dimensions map to the WD03 side projection. No native project front-elevation SVG contains the product projection, so no front context was invented.",
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "side": "none",
        "scope": "review crop only; complete project SVGs retain their annotations",
        "walls_furniture_and_ifc_geometry_removed": False,
    }
    manifest["review_preview_background"] = {
        "plan": "opaque white review-only paper background",
        "side": "opaque white review-only paper background",
        "full_project_svg_remains_transparent": True,
        "drawing_geometry_changed": False,
    }
    write_json(manifest_path, manifest)
