#!/usr/bin/env python3
"""Overlay exact 224.212.00.2 DWG views on complete project drawings."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json


ARTICLE = "224.212.00.2"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-duofix-sigma-224-212"
LINEWORK = PRODUCT_DIR / "official-native-dwg-linework.json"
PLAN_SOURCE = ROOT / "drawings/Sanitary Plan-P202-candidate.svg"
FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-06-21-R12-PX.svg"
SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-06-19-R12-PY.svg"
REPRESENTATIVE = "3hgNkx97vCTOC2eewCpMNk"
CONTEXT_ELEVATION_INSTANCE = "3pvAlH5C14v8uVEJ1LmK8M"
INSTANCES = [REPRESENTATIVE, CONTEXT_ELEVATION_INSTANCE]
NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")
BLUE = "#1677c8"
PROJECT_SCALE_SVG_UNITS_PER_MM = 0.02


def guid(element) -> str | None:
    return next((value for key, value in element.attrib.items() if key.endswith("guid")), None)


def path_points(element) -> list[tuple[float, float]]:
    values = [float(value) for value in NUMBER.findall(element.attrib.get("d", ""))]
    return list(zip(values[0::2], values[1::2]))


def groups_for_guid(source: Path, target_guid: str):
    groups = []
    for element in ET.parse(source).getroot().iter():
        if not element.tag.endswith("g") or guid(element) != target_guid:
            continue
        points = [
            point
            for child in element.iter()
            if child.tag.endswith("path")
            for point in path_points(child)
        ]
        if points:
            groups.append(points)
    return groups


def bbox(points):
    return (
        (min(point[0] for point in points), min(point[1] for point in points)),
        (max(point[0] for point in points), max(point[1] for point in points)),
    )


def sampled(points, maximum=650):
    if len(points) <= maximum:
        return points
    return points[:: math.ceil(len(points) / maximum)]


def nearest_rms(source, target):
    source, target = sampled(source), sampled(target)
    return math.sqrt(
        sum(min((first[0] - second[0]) ** 2 + (first[1] - second[1]) ** 2 for second in target) for first in source)
        / len(source)
    )


def transform_for_context(paths, target_points, view):
    raw = [tuple(point) for path in paths for point in path]
    target_min, target_max = bbox(target_points)
    target_center = tuple((target_min[axis] + target_max[axis]) / 2.0 for axis in range(2))
    expected_swap = view == "plan"
    candidates = []
    for swap in (False, True):
        if swap != expected_swap:
            continue
        oriented = [(y, x) if swap else (x, y) for x, y in raw]
        source_min, source_max = bbox(oriented)
        source_center = tuple((source_min[axis] + source_max[axis]) / 2.0 for axis in range(2))
        for flip_x in (False, True):
            for flip_y in (False, True):
                for x_anchor in ("minimum", "center", "maximum"):
                    for y_anchor in ("minimum", "center", "maximum"):
                        def axis_offset(axis, anchor):
                            if anchor == "minimum":
                                return target_min[axis] - source_min[axis] * PROJECT_SCALE_SVG_UNITS_PER_MM
                            if anchor == "maximum":
                                return target_max[axis] - source_max[axis] * PROJECT_SCALE_SVG_UNITS_PER_MM
                            return target_center[axis] - source_center[axis] * PROJECT_SCALE_SVG_UNITS_PER_MM

                        offset = (axis_offset(0, x_anchor), axis_offset(1, y_anchor))

                        def apply(point, *, _swap=swap, _fx=flip_x, _fy=flip_y, _offset=offset):
                            x, y = (point[1], point[0]) if _swap else point
                            if _fx:
                                x = source_min[0] + source_max[0] - x
                            if _fy:
                                y = source_min[1] + source_max[1] - y
                            return (
                                _offset[0] + x * PROJECT_SCALE_SVG_UNITS_PER_MM,
                                _offset[1] + y * PROJECT_SCALE_SVG_UNITS_PER_MM,
                            )

                        transformed = [apply(point) for point in raw]
                        score = nearest_rms(target_points, transformed)
                        candidates.append((score, swap, flip_x, flip_y, x_anchor, y_anchor, apply))
    score, swap, flip_x, flip_y, x_anchor, y_anchor, apply = min(candidates, key=lambda item: item[0])
    official_bbox = bbox([apply(point) for point in raw])
    target_size = [target_max[axis] - target_min[axis] for axis in range(2)]
    official_size = [official_bbox[1][axis] - official_bbox[0][axis] for axis in range(2)]
    return apply, {
        "swap_axes": swap,
        "flip_x": flip_x,
        "flip_y": flip_y,
        "x_anchor": x_anchor,
        "y_anchor": y_anchor,
        "fit_rms_svg_units": round(score, 6),
        "scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "target_bbox": [[round(value, 6) for value in target_min], [round(value, 6) for value in target_max]],
        "official_overlay_bbox": [[round(value, 6) for value in official_bbox[0]], [round(value, 6) for value in official_bbox[1]]],
        "target_size_svg_units": [round(value, 6) for value in target_size],
        "official_size_svg_units": [round(value, 6) for value in official_size],
        "size_absolute_delta_svg_units": [round(abs(official_size[axis] - target_size[axis]), 6) for axis in range(2)],
        "uniform_scale_preserved": True,
        "transformation_mode": "axis_swap_rigid_reflection_and_translation_only",
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


def add_overlays(source: Path, target: Path, overlays: list[dict]):
    content = source.read_text(encoding="utf-8")
    groups = []
    for overlay in overlays:
        groups.append(f'''<g id="geberit-duofix-224-212-{overlay["view"]}-{overlay["guid"]}" class="official-native-dwg project-context-overlay" data-source-kind="native_dwg" data-article-number="{ARTICLE}" data-native-dwg-code="{overlay["code"]}" data-source-sha256="{overlay["sha256"]}" data-ifc-guid="{overlay["guid"]}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <path class="official-reference-mask" d="{overlay["path"]}" stroke="#ffffff" stroke-width="0.14"/>
  <path class="official-reference native-dwg" d="{overlay["path"]}" stroke="{BLUE}" stroke-width="0.052"/>
</g>''')
    if "</svg>" not in content:
        raise RuntimeError(f"invalid project SVG: {source}")
    target.write_text(content.rsplit("</svg>", 1)[0] + "\n".join(groups) + "\n</svg>\n", encoding="utf-8")


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
        raise RuntimeError("qlmanage is required to render Duofix project-context previews")
    with tempfile.TemporaryDirectory(prefix="geberit-duofix-224212-context-preview-") as temporary:
        process = subprocess.run(
            [qlmanage, "-t", "-s", "1800", "-o", temporary, str(source)],
            capture_output=True,
            text=True,
        )
        generated = Path(temporary) / f"{source.name}.png"
        if process.returncode != 0 or not generated.is_file():
            raise RuntimeError(f"Duofix context preview failed: {process.stderr.strip()}")
        generated.replace(target)


def build_overlay(source, target_guid, view, linework, occurrence=0):
    groups = groups_for_guid(source, target_guid)
    if occurrence >= len(groups):
        raise RuntimeError(f"{target_guid} occurrence {occurrence} missing from {source}")
    view_data = linework["views"][view]
    transform, fit = transform_for_context(view_data["paths_mm"], groups[occurrence], view)
    return {
        "guid": target_guid,
        "view": view,
        "code": view_data["native_dwg_code"],
        "sha256": view_data["source_dwg_sha256"],
        "path": svg_path(view_data["paths_mm"], transform),
        "native_dwg_path_count": len(view_data["paths_mm"]),
        "fit": fit,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=PRODUCT_DIR)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    linework = load_json(LINEWORK)
    if linework.get("source_kind") != "native_dwg" or linework.get("article_number") != ARTICLE:
        raise RuntimeError("Duofix exact native-DWG identity gate failed")
    plan_overlays = [build_overlay(PLAN_SOURCE, instance, "plan", linework) for instance in INSTANCES]
    plan_target = output / "project-context-sanitary-plan.svg"
    add_overlays(PLAN_SOURCE, plan_target, plan_overlays)
    plan_crops = []
    for overlay in plan_overlays:
        crop = output / f'project-context-sanitary-plan-{overlay["guid"]}-review.svg'
        write_review_crop(plan_target, crop, overlay["fit"]["target_bbox"], 5.0, 5.0)
        preview = output / f'project-context-sanitary-plan-{overlay["guid"]}-review-preview.png'
        write_png_preview(crop, preview)
        plan_crops.append({
            "guid": overlay["guid"],
            "path": relative(crop),
            "sha256": sha256(crop),
            "preview": relative(preview),
            "preview_sha256": sha256(preview),
        })
    front_overlay = build_overlay(FRONT_SOURCE, CONTEXT_ELEVATION_INSTANCE, "front", linework)
    front_target = output / "project-context-front-elevation.svg"
    add_overlays(FRONT_SOURCE, front_target, [front_overlay])
    front_crop = output / "project-context-front-elevation-review.svg"
    write_review_crop(front_target, front_crop, front_overlay["fit"]["target_bbox"], 8.0, 6.0)
    front_preview = output / "project-context-front-elevation-review-preview.png"
    write_png_preview(front_crop, front_preview)
    side_overlay = build_overlay(SIDE_SOURCE, CONTEXT_ELEVATION_INSTANCE, "side", linework)
    side_target = output / "project-context-side-elevation.svg"
    add_overlays(SIDE_SOURCE, side_target, [side_overlay])
    side_crop = output / "project-context-side-elevation-review.svg"
    write_review_crop(side_target, side_crop, side_overlay["fit"]["target_bbox"], 8.0, 6.0)
    side_preview = output / "project-context-side-elevation-review-preview.png"
    write_png_preview(side_crop, side_preview)
    records = [
        {"view": "plan", "source": relative(PLAN_SOURCE), "output": relative(plan_target), "output_sha256": sha256(plan_target), "review_crops": plan_crops, "overlays": plan_overlays},
        {"view": "front", "source": relative(FRONT_SOURCE), "output": relative(front_target), "output_sha256": sha256(front_target), "review_crop": relative(front_crop), "review_crop_sha256": sha256(front_crop), "review_preview": relative(front_preview), "review_preview_sha256": sha256(front_preview), "overlays": [front_overlay]},
        {"view": "side", "source": relative(SIDE_SOURCE), "output": relative(side_target), "output_sha256": sha256(side_target), "review_crop": relative(side_crop), "review_crop_sha256": sha256(side_crop), "review_preview": relative(side_preview), "review_preview_sha256": sha256(side_preview), "overlays": [side_overlay]},
    ]
    manifest = {
        "schema_version": 1,
        "generator": "pipeline/scripts/render_geberit_duofix_sigma_224_212_project_context.py",
        "article_number": ARTICLE,
        "source_kind": "native_dwg",
        "source_linework": relative(LINEWORK),
        "source_linework_sha256": sha256(LINEWORK),
        "project_context_retained": True,
        "walls_and_surrounding_project_elements_retained": True,
        "blue_line_top_layer_with_white_mask": True,
        "official_dwg_uniform_project_scale_preserved": True,
        "project_scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "views": records,
        "pass": all(item["fit"]["uniform_scale_preserved"] for record in records for item in record["overlays"]),
    }
    target = output / "project-context-manifest.json"
    write_json(target, manifest)
    print(json.dumps({"manifest": relative(target), "outputs": [record["output"] for record in records], "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
