#!/usr/bin/env python3
"""Overlay Baxter Marilyn native-DWG views on complete project drawings."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from render_sis04_project_context import (
    add_white_review_background,
    suppress_plan_elevation_markers,
    write_uncached_png_preview,
)


PRODUCT_DIR = ROOT / "output/review/highpoly-types/marilyn-01"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-09-33-R22-PY.svg"
SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-09-36-R22-NX.svg"
GLOBAL_ID = "3l2Ji4k2H9oOTDGWYUq7uV"
SOURCE_KIND = "native_dwg"
SOURCE_LABEL_ZH = "Baxter 官方原生 DWG"
SOURCE_DWG_SHA256 = "cb9825eecab7c78ee28e7668e933c2fe413395a63bf75f3afab8cc4959864724"
BLUE = "#1677c8"
NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")
PROJECT_SCALE_SVG_UNITS_PER_MM = 0.02
BBOX_TOLERANCE_SVG_UNITS = 1.2
ALIGNMENT_MODE = "baseline_centre"
OVERLAY_ID_PREFIX = "marilyn-01"
GENERATOR = "pipeline/scripts/render_marilyn_01_project_context.py"
CONTEXT_VIEWS = (
    ("plan", "plan", PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 0),
    ("front", "front", FRONT_SOURCE, "project-context-r22-front-elevation.svg", (9.0, 7.0), 0),
    ("side", "side", SIDE_SOURCE, "project-context-r22-side-elevation.svg", (10.0, 7.0), 0),
)


def guid(element) -> str | None:
    return next((value for key, value in element.attrib.items() if key.endswith("guid")), None)


def path_points(element) -> list[tuple[float, float]]:
    values = [float(value) for value in NUMBER.findall(element.attrib.get("d", ""))]
    return list(zip(values[0::2], values[1::2]))


def group_points(source: Path) -> list[tuple[float, float]]:
    groups = []
    for element in ET.parse(source).getroot().iter():
        if not element.tag.endswith("g") or guid(element) != GLOBAL_ID or "projection" not in element.attrib.get("class", "").split():
            continue
        points = [
            point
            for child in element.iter()
            if child.tag.endswith("path")
            for point in path_points(child)
        ]
        if points:
            groups.append(points)
    if not groups:
        raise RuntimeError(f"expected a {GLOBAL_ID} projection in {source}")
    # R22 side contains a second small/clipped projection group. Select the
    # mechanically complete group with the largest projected bounding area.
    return max(groups, key=lambda points: (bbox(points)[1][0] - bbox(points)[0][0]) * (bbox(points)[1][1] - bbox(points)[0][1]))


def bbox(points):
    return (
        (min(point[0] for point in points), min(point[1] for point in points)),
        (max(point[0] for point in points), max(point[1] for point in points)),
    )


def transform_for_context(paths, target_points, rotate_quarter_turns=0, alignment_mode=ALIGNMENT_MODE):
    turns = rotate_quarter_turns % 4

    def rotate(point):
        x, y = point
        for _ in range(turns):
            x, y = -y, x
        return x, y

    original_raw = [tuple(point) for path in paths for point in path]
    raw = [rotate(point) for point in original_raw]
    source_min, source_max = bbox(raw)
    target_min, target_max = bbox(target_points)
    source_size = [source_max[axis] - source_min[axis] for axis in range(2)]
    target_size = [target_max[axis] - target_min[axis] for axis in range(2)]
    expected_size = [value * PROJECT_SCALE_SVG_UNITS_PER_MM for value in source_size]
    deltas = [abs(expected_size[axis] - target_size[axis]) for axis in range(2)]
    if max(deltas) > BBOX_TOLERANCE_SVG_UNITS:
        raise RuntimeError(f"Marilyn native DWG/project projection bbox mismatch: {deltas}")

    # Project drawings use downward-positive SVG Y. Reflect the second projected
    # coordinate and align the complete proxy bbox with the original Hima group.
    def apply(point):
        point = rotate(point)
        if alignment_mode in {"centre", "baseline_centre"}:
            source_centre = (
                (source_min[0] + source_max[0]) / 2.0,
                (source_min[1] + source_max[1]) / 2.0,
            )
            target_centre = (
                (target_min[0] + target_max[0]) / 2.0,
                (target_min[1] + target_max[1]) / 2.0,
            )
            return (
                target_centre[0] + (point[0] - source_centre[0]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
                (
                    target_max[1] - (point[1] - source_min[1]) * PROJECT_SCALE_SVG_UNITS_PER_MM
                    if alignment_mode == "baseline_centre"
                    else target_centre[1] - (point[1] - source_centre[1]) * PROJECT_SCALE_SVG_UNITS_PER_MM
                ),
            )
        return (
            target_min[0] + (point[0] - source_min[0]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
            target_max[1] - (point[1] - source_min[1]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
        )

    transformed_bbox = bbox([apply(point) for point in original_raw])
    return apply, {
        "rotate_quarter_turns": turns,
        "alignment_mode": alignment_mode,
        "flip_projected_y": True,
        "scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "target_bbox": [[round(value, 6) for value in target_min], [round(value, 6) for value in target_max]],
        "proxy_bbox": [[round(value, 6) for value in transformed_bbox[0]], [round(value, 6) for value in transformed_bbox[1]]],
        "bbox_absolute_delta_svg_units": [round(value, 6) for value in deltas],
        "bbox_tolerance_svg_units": BBOX_TOLERANCE_SVG_UNITS,
        "uniform_scale_preserved": True,
        "transformation_mode": "rigid_view_orientation_reflection_translation_and_native_1_50_scale_only",
        "pass": max(deltas) <= BBOX_TOLERANCE_SVG_UNITS,
    }


def svg_path(paths, transform):
    commands = []
    for path in paths:
        if len(path) < 2:
            continue
        first = transform(tuple(path[0]))
        commands.append(f"M {first[0]:.6f},{first[1]:.6f}")
        commands.extend(
            f"L {point[0]:.6f},{point[1]:.6f}"
            for point in (transform(tuple(item)) for item in path[1:])
        )
    return " ".join(commands)


def add_overlay(source: Path, target: Path, view: str, official_path: str, mask_path: str):
    content = source.read_text(encoding="utf-8")
    product_pattern = re.compile(rf'<g(?=[^>]*\bifc:guid="{re.escape(GLOBAL_ID)}")(?=[^>]*\bclass="[^"]*\bprojection\b[^"]*")')
    content, hidden_count = product_pattern.subn('<g style="display:none" data-review-replaced-by="official-native-dwg"', content)
    if hidden_count < 1:
        raise RuntimeError(f"Marilyn project projection missing from {source}")
    group = f'''<g id="{OVERLAY_ID_PREFIX}-{view}-{GLOBAL_ID}" class="official-native-dwg project-context-overlay" data-source-kind="{SOURCE_KIND}" data-source-label-zh="{SOURCE_LABEL_ZH}" data-source-sha256="{SOURCE_DWG_SHA256}" data-official-cad-used="true" data-third-party-cad-used="false" data-ifc-guid="{GLOBAL_ID}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <path class="official-native-dwg-mask" d="{mask_path}" fill="#ffffff" fill-rule="evenodd" stroke="#ffffff" stroke-width="0.18"/>
  <path class="official-native-dwg-blue" d="{official_path}" stroke="{BLUE}" stroke-width="0.064"/>
</g>'''
    if "</svg>" not in content:
        raise RuntimeError(f"invalid project SVG: {source}")
    target.write_text(content.rsplit("</svg>", 1)[0] + group + "\n</svg>\n", encoding="utf-8")


def write_review_crop(source: Path, target: Path, target_bbox, padding_x, padding_y):
    minimum, maximum = target_bbox
    x, y = minimum[0] - padding_x, minimum[1] - padding_y
    width = maximum[0] - minimum[0] + 2.0 * padding_x
    height = maximum[1] - minimum[1] + 2.0 * padding_y
    content = source.read_text(encoding="utf-8")
    content, count = re.subn(
        r'viewBox="[^"]+"',
        f'viewBox="{x:.6f} {y:.6f} {width:.6f} {height:.6f}"',
        content,
        count=1,
    )
    if count != 1:
        raise RuntimeError(f"could not crop {source}")
    target.write_text(content, encoding="utf-8")


def write_png_preview(source: Path, target: Path):
    qlmanage = shutil.which("qlmanage")
    if qlmanage is None:
        raise RuntimeError("qlmanage is required to render the project-context preview")
    with tempfile.TemporaryDirectory(prefix=f"{OVERLAY_ID_PREFIX}-context-preview-") as temporary:
        process = subprocess.run(
            [qlmanage, "-t", "-s", "1800", "-o", temporary, str(source)],
            capture_output=True,
            text=True,
        )
        generated = Path(temporary) / f"{source.name}.png"
        if process.returncode != 0 or not generated.is_file():
            raise RuntimeError(f"project-context preview failed: {process.stderr.strip()}")
        generated.replace(target)


def suppress_noninteger_diagnostic_highlights(path: Path):
    """Remove residual-error diagnostics from a review crop only."""
    content = path.read_text(encoding="utf-8")
    content, count = re.subn(
        r'<g id="noninteger-highlights".*?</g>',
        "",
        content,
        count=1,
        flags=re.DOTALL,
    )
    if count != 1:
        raise RuntimeError(f"expected noninteger diagnostic highlights in {path}")
    path.write_text(content, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=PRODUCT_DIR)
    args = parser.parse_args()
    output = args.output.resolve()
    candidate = load_json(CANDIDATE)
    if (
        candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("official_cad_used") is not True
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("source_dwg_sha256") != SOURCE_DWG_SHA256
    ):
        raise RuntimeError("Marilyn native-DWG source gate failed")
    records = []
    for view, candidate_view, source, filename, crop_padding, rotate_quarter_turns in CONTEXT_VIEWS:
        target_points = group_points(source)
        paths = candidate["views"][candidate_view]["official_native_dwg_paths_mm"]
        proxy_paths = candidate["views"][candidate_view]["proxy_paths_mm"]
        transform, fit = transform_for_context(
            paths,
            target_points,
            rotate_quarter_turns,
            "centre" if view == "plan" else "baseline_centre",
        )
        overlay_path = svg_path(paths, transform)
        mask_path = svg_path(proxy_paths, transform)
        target = output / filename
        add_overlay(source, target, view, overlay_path, mask_path)
        crop = output / filename.replace(".svg", "-review.svg")
        write_review_crop(target, crop, fit["target_bbox"], *crop_padding)
        preview = output / filename.replace(".svg", "-review-preview.png")
        write_png_preview(crop, preview)
        records.append({
            "view": view,
            "candidate_view": candidate_view,
            "source": relative(source),
            "output": relative(target),
            "output_sha256": sha256(target),
            "review_crop": relative(crop),
            "review_crop_sha256": sha256(crop),
            "review_preview": relative(preview),
            "review_preview_sha256": sha256(preview),
            "overlay": {
                "ifc_guid": GLOBAL_ID,
                "source_kind": SOURCE_KIND,
                "source_label_zh": SOURCE_LABEL_ZH,
                "official_cad_used": True,
                "third_party_cad_used": False,
                "source_dwg_sha256": SOURCE_DWG_SHA256,
                "path_count": len(paths),
                "fit": fit,
            },
        })
    manifest = {
        "schema_version": 1,
        "generator": GENERATOR,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "project_context_retained": True,
        "walls_and_surrounding_project_elements_retained": True,
        "overlay_top_layer_with_white_mask": True,
        "blue_line_present": True,
        "source_dwg_sha256": SOURCE_DWG_SHA256,
        "project_scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "views": records,
        "pass": all(record["overlay"]["fit"]["pass"] for record in records),
    }
    target = output / "project-context-manifest.json"
    write_json(target, manifest)
    print(json.dumps({"manifest": relative(target), "outputs": [record["output"] for record in records], "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
    manifest_path = PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    for record in manifest["views"]:
        crop = ROOT / record["review_crop"]
        preview = ROOT / record["review_preview"]
        if record["view"] == "plan":
            suppress_plan_elevation_markers(crop)
        else:
            suppress_noninteger_diagnostic_highlights(crop)
        add_white_review_background(crop)
        write_uncached_png_preview(crop, preview)
        record["review_crop_sha256"] = sha256(crop)
        record["review_preview_sha256"] = sha256(preview)
    manifest["semantic_view_mapping"] = {
        "plan": {
            "candidate_axes": [0, 1],
            "source": relative(PLAN_SOURCE),
            "rotate_quarter_turns": 0,
            "project_direction": "+Z",
        },
        "front": {
            "candidate_axes": [0, 2],
            "source": relative(FRONT_SOURCE),
            "rotate_quarter_turns": 0,
            "project_direction": "+Y",
            "complete_project_projection_group_only": True,
        },
        "side": {
            "candidate_axes": [1, 2],
            "source": relative(SIDE_SOURCE),
            "rotate_quarter_turns": 0,
            "project_direction": "-X",
            "complete_project_projection_group_only": True,
        },
    }
    manifest["context_view_scope"] = {
        "included": ["plan", "front", "side"],
        "excluded": [],
        "reason": "The Furniture Plan and selected R22 elevations contain mechanically complete Marilyn projections.",
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "front": "noninteger residual-error diagnostic highlights only",
        "side": "noninteger residual-error diagnostic highlights only",
        "scope": "review crops only; complete project SVGs retain their annotations",
        "walls_furniture_and_ifc_geometry_removed": False,
    }
    manifest["review_preview_background"] = {
        "plan": "opaque white review-only paper background",
        "front": "opaque white review-only paper background",
        "side": "opaque white review-only paper background",
        "full_project_svgs_remain_transparent": True,
        "drawing_geometry_changed": False,
    }
    manifest["project_projection_note"] = {
        "official_nominal_xyz_mm": [860.0, 1000.0, 940.0],
        "ifc_body_local_xyz_mm": [877.174957, 1013.911865, 970.476962],
        "project_ifc_description": "Bergère armchair with swivel base W86D100H78",
        "description_conflict": "The 78 cm description is stale; the current official bergere and IFC Body match the 94 cm family variant.",
        "interpretation": "Official native-DWG linework is overlaid at the drawing's native 1:50 scale and aligned by centre/floor baseline; no view is stretched to force-fit the project projection.",
        "geometry_stretched": False,
    }
    write_json(manifest_path, manifest)
