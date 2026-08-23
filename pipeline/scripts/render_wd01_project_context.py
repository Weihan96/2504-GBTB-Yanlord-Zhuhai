#!/usr/bin/env python3
"""Overlay WD01 on its real Furniture Plan and R09 +X/+Y elevations."""

import re
import xml.etree.ElementTree as ET

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from render_sis04_project_context import (
    add_white_review_background,
    suppress_plan_elevation_markers,
    write_uncached_png_preview,
)
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/wd01"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "2mPTt7$nvBSAn28I9NWt6T"
shared.OVERLAY_ID_PREFIX = "poliform-pivot-senzafine-wd01"
shared.GENERATOR = "pipeline/scripts/render_wd01_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-05-15-R09-PX.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-05-14-R09-PY.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.30
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 1),
    ("front", "front", shared.FRONT_SOURCE, "project-context-front-elevation.svg", (10.0, 8.0), 0),
    ("side", "side", shared.SIDE_SOURCE, "project-context-side-elevation.svg", (10.0, 8.0), 0),
)


def projection_group_points(source):
    """Use the native projection group when Bonsai also emitted a same-GUID cut group."""
    groups = []
    for element in ET.parse(source).getroot().iter():
        if (
            not element.tag.endswith("g")
            or shared.guid(element) != shared.GLOBAL_ID
            or "projection" not in element.attrib.get("class", "").split()
        ):
            continue
        points = [
            point
            for child in element.iter()
            if child.tag.endswith("path")
            for point in shared.path_points(child)
        ]
        if points:
            groups.append(points)
    if len(groups) != 1:
        raise RuntimeError(f"expected one WD01 projection group in {source}, found {len(groups)}")
    return groups[0]


shared.group_points = projection_group_points
original_add_overlay = shared.add_overlay


def add_overlay(source, target, view, path):
    original_add_overlay(source, target, view, path)
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
            "rotate_quarter_turns": 1,
            "project_direction": "+Z",
        },
        "front": {
            "candidate_axes": [0, 2],
            "source": relative(shared.FRONT_SOURCE),
            "rotate_quarter_turns": 0,
            "project_direction": "+X",
        },
        "side": {
            "candidate_axes": [1, 2],
            "source": relative(shared.SIDE_SOURCE),
            "rotate_quarter_turns": 0,
            "project_direction": "+Y",
        },
    }
    manifest["context_view_scope"] = {
        "included": ["plan", "front", "side"],
        "excluded": [],
        "reason": "Furniture Plan.svg and the R09 +X/+Y native elevation SVGs contain an actual projection group for the same WD01 GlobalId. The alternate R09 -Y elevation is retained unchanged but is not duplicated in this three-view review.",
    }
    side = next(item for item in manifest["views"] if item["view"] == "side")
    manifest["project_projection_note"] = {
        "side_original_visible_projection_bbox_mm": [
            round((side["overlay"]["fit"]["target_bbox"][1][0] - side["overlay"]["fit"]["target_bbox"][0][0]) / 0.02, 6),
            round((side["overlay"]["fit"]["target_bbox"][1][1] - side["overlay"]["fit"]["target_bbox"][0][1]) / 0.02, 6),
        ],
        "complete_ifc_body_side_bbox_mm": [673.52434, 2390.022879],
        "interpretation": "The +Y project elevation clips or occludes about 14 mm of the complete Body side width. The review overlay preserves the complete Body at native 1:50 scale and aligns it to the original projection centre.",
        "geometry_stretched": False,
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "front": "noninteger-highlights diagnostic group only",
        "side": "noninteger-highlights diagnostic group only",
        "scope": "derived review copies only; complete original project SVG sources remain unchanged",
        "walls_furniture_and_ifc_geometry_removed": False,
    }
    manifest["review_preview_background"] = {
        "plan": "opaque white review-only paper background",
        "front": "opaque white review-only paper background",
        "side": "opaque white review-only paper background",
        "full_project_svg_remains_transparent": True,
        "drawing_geometry_changed": False,
    }
    write_json(manifest_path, manifest)
