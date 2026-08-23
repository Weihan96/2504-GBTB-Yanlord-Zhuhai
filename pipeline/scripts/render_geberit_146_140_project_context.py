#!/usr/bin/env python3
"""Overlay exact Geberit 146.140 native-DWG views on project-context SVGs."""

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


PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-146-140"
LINEWORK = PRODUCT_DIR / "official-native-dwg-linework.json"
PLAN_SOURCE = ROOT / "drawings/Sanitary Plan-P202-candidate.svg"
FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-08-29-R17-PX.svg"
SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-08-31-R17-NY.svg"
REPRESENTATIVE = "1rhZG98PPCSxaLeMFLTYb9"
INSTANCES = ["0UtU7yPb10ku4gsbGoM_sp", REPRESENTATIVE]
NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")
BLUE = "#1677c8"
PROJECT_SCALE_SVG_UNITS_PER_MM = 0.02
BBOX_TOLERANCE_SVG_UNITS = 0.001


def guid(element) -> str | None:
    return next((value for key, value in element.attrib.items() if key.endswith("guid")), None)


def path_points(element) -> list[tuple[float, float]]:
    values = [float(value) for value in NUMBER.findall(element.attrib.get("d", ""))]
    return list(zip(values[0::2], values[1::2]))


def groups_for_guid(source: Path, target_guid: str):
    root = ET.parse(source).getroot()
    groups = []
    for element in root.iter():
        if not element.tag.endswith("g") or guid(element) != target_guid:
            continue
        points = [point for child in element.iter() if child.tag.endswith("path") for point in path_points(child)]
        if points:
            groups.append(points)
    return groups


def bbox(points):
    return (
        (min(point[0] for point in points), min(point[1] for point in points)),
        (max(point[0] for point in points), max(point[1] for point in points)),
    )


def sample(points, maximum=450):
    if len(points) <= maximum:
        return points
    step = math.ceil(len(points) / maximum)
    return points[::step]


def nearest_rms(source, target):
    source = sample(source)
    target = sample(target)
    return math.sqrt(
        sum(min((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for b in target) for a in source)
        / len(source)
    )


def fitted_transform(paths, target_points, swaps, flip_ys=(False, True)):
    target_min, target_max = bbox(target_points)
    target_center = tuple((target_min[index] + target_max[index]) / 2.0 for index in range(2))
    candidates = []
    raw_points = [tuple(point) for path in paths for point in path]
    for swap in swaps:
        oriented = [(y, x) if swap else (x, y) for x, y in raw_points]
        source_min, source_max = bbox(oriented)
        source_center = tuple((source_min[index] + source_max[index]) / 2.0 for index in range(2))
        source_size = [source_max[i] - source_min[i] for i in range(2)]
        target_size = [target_max[i] - target_min[i] for i in range(2)]
        official_size = [value * PROJECT_SCALE_SVG_UNITS_PER_MM for value in source_size]
        size_deltas = [abs(official_size[index] - target_size[index]) for index in range(2)]
        if max(size_deltas) > BBOX_TOLERANCE_SVG_UNITS:
            continue
        for flip_x in (False, True):
            for flip_y in flip_ys:
                def apply(
                    point,
                    *,
                    _swap=swap,
                    _flip_x=flip_x,
                    _flip_y=flip_y,
                    _source_center=source_center,
                ):
                    x, y = (point[1], point[0]) if _swap else point
                    if _flip_x:
                        x = 2.0 * _source_center[0] - x
                    if _flip_y:
                        y = 2.0 * _source_center[1] - y
                    return (
                        target_center[0] + (x - _source_center[0]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
                        target_center[1] + (y - _source_center[1]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
                    )
                transformed = [apply(point) for point in raw_points]
                score = nearest_rms(transformed, target_points)
                candidates.append((score, swap, flip_x, flip_y, size_deltas, official_size, apply))
    if not candidates:
        raise RuntimeError("official Geberit DWG/project projection bbox mismatch")
    score, swap, flip_x, flip_y, size_deltas, official_size, apply = min(candidates, key=lambda item: item[0])
    return apply, {
        "swap_axes": swap,
        "flip_x": flip_x,
        "flip_y": flip_y,
        "fit_rms_svg_units": round(score, 6),
        "scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "uniform_scale_preserved": True,
        "target_size_svg_units": [round(value, 6) for value in target_size],
        "official_overlay_size_svg_units": [round(value, 6) for value in official_size],
        "size_absolute_delta_svg_units": [round(value, 6) for value in size_deltas],
        "bbox_tolerance_svg_units": BBOX_TOLERANCE_SVG_UNITS,
        "target_bbox": [[round(value, 6) for value in target_min], [round(value, 6) for value in target_max]],
    }


def svg_path(paths, transform):
    commands = []
    for path in paths:
        if len(path) < 2:
            continue
        first = transform(tuple(path[0]))
        commands.append(f"M {first[0]:.6f},{first[1]:.6f}")
        commands.extend(f"L {point[0]:.6f},{point[1]:.6f}" for point in (transform(tuple(item)) for item in path[1:]))
    return " ".join(commands)


def add_overlays(source: Path, target: Path, overlays: list[dict]):
    content = source.read_text(encoding="utf-8")
    groups = []
    for overlay in overlays:
        path = overlay["path"]
        groups.append(f'''<g id="geberit-146-140-{overlay["view"]}-{overlay["guid"]}" class="official-native-dwg project-context-overlay" data-source-kind="native_dwg" data-article-number="146.140.11.1" data-native-dwg-code="{overlay["code"]}" data-source-sha256="{overlay["sha256"]}" data-ifc-guid="{overlay["guid"]}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <path class="official-reference-mask" d="{path}" stroke="#ffffff" stroke-width="0.13"/>
  <path class="official-reference native-dwg" d="{path}" stroke="{BLUE}" stroke-width="0.052"/>
</g>''')
    insertion = "\n".join(groups) + "\n"
    if "</svg>" not in content:
        raise RuntimeError(f"invalid project SVG: {source}")
    target.write_text(content.rsplit("</svg>", 1)[0] + insertion + "</svg>\n", encoding="utf-8")


def write_review_crop(source: Path, target: Path, target_bbox, padding_x=18.0, padding_y=14.0):
    (minimum, maximum) = target_bbox
    x = minimum[0] - padding_x
    y = minimum[1] - padding_y
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
        raise RuntimeError(f"project context review crop could not replace viewBox: {source}")
    target.write_text(content, encoding="utf-8")


def write_png_preview(source: Path, target: Path):
    qlmanage = shutil.which("qlmanage")
    if qlmanage is None:
        raise RuntimeError("qlmanage is required to render Geberit project-context previews")
    with tempfile.TemporaryDirectory(prefix="geberit-146140-context-preview-") as temporary:
        process = subprocess.run(
            [qlmanage, "-t", "-s", "1800", "-o", temporary, str(source)],
            capture_output=True,
            text=True,
        )
        generated = Path(temporary) / f"{source.name}.png"
        if process.returncode != 0 or not generated.is_file():
            raise RuntimeError(f"Geberit context preview failed: {process.stderr.strip()}")
        generated.replace(target)


def build_overlay(source, target_guid, view, linework, occurrence=0, swaps=(False,), flip_ys=(False, True)):
    groups = groups_for_guid(source, target_guid)
    if occurrence >= len(groups):
        raise RuntimeError(f"{target_guid} occurrence {occurrence} missing from {source}")
    view_data = linework["views"][view]
    transform, fit = fitted_transform(view_data["paths_mm"], groups[occurrence], swaps, flip_ys)
    if not fit["uniform_scale_preserved"] or max(fit["size_absolute_delta_svg_units"]) > BBOX_TOLERANCE_SVG_UNITS:
        raise RuntimeError(f"official {view} fixed-scale transform drifted in {source}: {fit}")
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
    if linework.get("source_kind") != "native_dwg" or linework.get("article_number") != "146.140.11.1":
        raise RuntimeError("Geberit native-DWG linework identity gate failed")
    plan_overlays = [
        build_overlay(PLAN_SOURCE, instance, "plan", linework, swaps=(False, True))
        for instance in INSTANCES
    ]
    plan_target = output / "project-context-sanitary-plan.svg"
    add_overlays(PLAN_SOURCE, plan_target, plan_overlays)
    plan_crops = []
    for overlay in plan_overlays:
        crop = output / f'project-context-sanitary-plan-{overlay["guid"]}-review.svg'
        write_review_crop(plan_target, crop, overlay["fit"]["target_bbox"])
        preview = output / f'project-context-sanitary-plan-{overlay["guid"]}-review-preview.png'
        write_png_preview(crop, preview)
        plan_crops.append({
            "guid": overlay["guid"],
            "path": relative(crop),
            "sha256": sha256(crop),
            "preview": relative(preview),
            "preview_sha256": sha256(preview),
        })
    front_overlay = build_overlay(FRONT_SOURCE, REPRESENTATIVE, "front", linework, flip_ys=(True,))
    front_target = output / "project-context-front-elevation.svg"
    add_overlays(FRONT_SOURCE, front_target, [front_overlay])
    front_crop = output / "project-context-front-elevation-review.svg"
    write_review_crop(front_target, front_crop, front_overlay["fit"]["target_bbox"], 12.0, 8.0)
    front_preview = output / "project-context-front-elevation-review-preview.png"
    write_png_preview(front_crop, front_preview)
    side_overlay = build_overlay(SIDE_SOURCE, REPRESENTATIVE, "side", linework, flip_ys=(True,))
    side_target = output / "project-context-side-elevation.svg"
    add_overlays(SIDE_SOURCE, side_target, [side_overlay])
    side_crop = output / "project-context-side-elevation-review.svg"
    write_review_crop(side_target, side_crop, side_overlay["fit"]["target_bbox"], 12.0, 8.0)
    side_preview = output / "project-context-side-elevation-review-preview.png"
    write_png_preview(side_crop, side_preview)
    records = [
        {"view": "plan", "source": relative(PLAN_SOURCE), "output": relative(plan_target), "output_sha256": sha256(plan_target), "review_crops": plan_crops, "overlays": plan_overlays},
        {"view": "front", "source": relative(FRONT_SOURCE), "output": relative(front_target), "output_sha256": sha256(front_target), "review_crop": relative(front_crop), "review_crop_sha256": sha256(front_crop), "review_preview": relative(front_preview), "review_preview_sha256": sha256(front_preview), "overlays": [front_overlay]},
        {"view": "side", "source": relative(SIDE_SOURCE), "output": relative(side_target), "output_sha256": sha256(side_target), "review_crop": relative(side_crop), "review_crop_sha256": sha256(side_crop), "review_preview": relative(side_preview), "review_preview_sha256": sha256(side_preview), "overlays": [side_overlay]},
    ]
    manifest = {
        "schema_version": 1,
        "generator": "pipeline/scripts/render_geberit_146_140_project_context.py",
        "article_number": "146.140.11.1",
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
