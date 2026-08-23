#!/usr/bin/env python3
"""Overlay both CHA02 geometry-derived plan proxies on the furniture plan."""

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


PRODUCT_DIR = ROOT / "output/review/highpoly-types/cha02"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
GLOBAL_IDS = ("1GkE2YZ316gQw10Ev5DZOp", "3YoxxZCgbF3Ap7gztAKkcs")
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
PROJECT_SCALE_SVG_UNITS_PER_MM = 0.02
BBOX_TOLERANCE_SVG_UNITS = 0.02
NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")


def guid(element) -> str | None:
    return next((value for key, value in element.attrib.items() if key.endswith("guid")), None)


def path_points(element) -> list[tuple[float, float]]:
    values = [float(value) for value in NUMBER.findall(element.attrib.get("d", ""))]
    return list(zip(values[0::2], values[1::2]))


def group_points(source: Path, target_guid: str) -> list[tuple[float, float]]:
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
    if len(groups) != 1:
        raise RuntimeError(f"expected one {target_guid} projection in {source}, found {len(groups)}")
    return groups[0]


def bbox(points):
    return (
        (min(point[0] for point in points), min(point[1] for point in points)),
        (max(point[0] for point in points), max(point[1] for point in points)),
    )


def sampled(points, maximum=600):
    if len(points) <= maximum:
        return points
    return points[:: math.ceil(len(points) / maximum)]


def nearest_rms(source, target):
    source, target = sampled(source), sampled(target)
    return math.sqrt(
        sum(
            min((first[0] - second[0]) ** 2 + (first[1] - second[1]) ** 2 for second in target)
            for first in source
        )
        / len(source)
    )


def transform_for_context(paths, target_points):
    original = [tuple(point) for path in paths for point in path]
    target_min, target_max = bbox(target_points)
    target_center = tuple((target_min[index] + target_max[index]) / 2.0 for index in range(2))
    candidates = []
    for turns in range(4):
        def rotate(point):
            x, y = point
            for _ in range(turns):
                x, y = -y, x
            return x, y

        rotated = [rotate(point) for point in original]
        source_min, source_max = bbox(rotated)
        source_center = tuple((source_min[index] + source_max[index]) / 2.0 for index in range(2))
        source_size = [source_max[index] - source_min[index] for index in range(2)]
        target_size = [target_max[index] - target_min[index] for index in range(2)]
        deltas = [
            abs(source_size[index] * PROJECT_SCALE_SVG_UNITS_PER_MM - target_size[index])
            for index in range(2)
        ]
        if max(deltas) > BBOX_TOLERANCE_SVG_UNITS:
            continue

        def apply(point, *, _rotate=rotate, _source_center=source_center):
            x, y = _rotate(point)
            return (
                target_center[0] + (x - _source_center[0]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
                target_center[1] - (y - _source_center[1]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
            )

        transformed = [apply(point) for point in original]
        candidates.append((nearest_rms(transformed, target_points), turns, deltas, apply, bbox(transformed)))
    if not candidates:
        raise RuntimeError("CHA02 plan proxy/project projection bbox mismatch")
    score, turns, deltas, apply, transformed_bbox = min(candidates, key=lambda item: item[0])
    return apply, {
        "rotate_quarter_turns": turns,
        "alignment_mode": "centre",
        "flip_projected_y": True,
        "scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "nearest_rms_svg_units": round(score, 6),
        "target_bbox": [[round(value, 6) for value in target_min], [round(value, 6) for value in target_max]],
        "proxy_bbox": [[round(value, 6) for value in transformed_bbox[0]], [round(value, 6) for value in transformed_bbox[1]]],
        "bbox_absolute_delta_svg_units": [round(value, 6) for value in deltas],
        "bbox_tolerance_svg_units": BBOX_TOLERANCE_SVG_UNITS,
        "uniform_scale_preserved": True,
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
        commands.append("Z")
    return " ".join(commands)


def add_overlays(source: Path, target: Path, overlays: list[dict]):
    content = source.read_text(encoding="utf-8")
    groups = []
    for overlay in overlays:
        groups.append(f'''<g id="cha02-plan-{overlay["ifc_guid"]}" class="geometry-derived-simplified-proxy project-context-overlay" data-source-kind="{SOURCE_KIND}" data-source-label-zh="{SOURCE_LABEL_ZH}" data-official-cad-used="false" data-ifc-guid="{overlay["ifc_guid"]}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <path class="geometry-derived-proxy-mask" d="{overlay["path"]}" stroke="#ffffff" stroke-width="0.15"/>
  <path class="geometry-derived-proxy" d="{overlay["path"]}" stroke="#111820" stroke-width="0.058"/>
</g>''')
    if "</svg>" not in content:
        raise RuntimeError(f"invalid project SVG: {source}")
    target.write_text(content.rsplit("</svg>", 1)[0] + "\n".join(groups) + "\n</svg>\n", encoding="utf-8")


def write_review_crop(source: Path, target: Path, boxes, padding_x=14.0, padding_y=14.0):
    minimum = (min(box[0][0] for box in boxes), min(box[0][1] for box in boxes))
    maximum = (max(box[1][0] for box in boxes), max(box[1][1] for box in boxes))
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
        raise RuntimeError("qlmanage is required to render the CHA02 project-context preview")
    with tempfile.TemporaryDirectory(prefix="cha02-context-preview-") as temporary:
        process = subprocess.run(
            [qlmanage, "-t", "-s", "1800", "-o", temporary, str(source)],
            capture_output=True,
            text=True,
        )
        generated = Path(temporary) / f"{source.name}.png"
        if process.returncode != 0 or not generated.is_file():
            raise RuntimeError(f"CHA02 context preview failed: {process.stderr.strip()}")
        generated.replace(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=PRODUCT_DIR)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    candidate = load_json(CANDIDATE)
    if (
        candidate.get("profile_key") != "cha02"
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
    ):
        raise RuntimeError("CHA02 geometry-derived source gate failed")
    paths = candidate["views"]["plan"]["proxy_paths_mm"]
    overlays = []
    for target_guid in GLOBAL_IDS:
        target_points = group_points(PLAN_SOURCE, target_guid)
        transform, fit = transform_for_context(paths, target_points)
        overlays.append({
            "ifc_guid": target_guid,
            "source_kind": SOURCE_KIND,
            "source_label_zh": SOURCE_LABEL_ZH,
            "official_cad_used": False,
            "third_party_cad_used": False,
            "path_count": len(paths),
            "path": svg_path(paths, transform),
            "fit": fit,
        })
    target = output / "project-context-furniture-plan.svg"
    add_overlays(PLAN_SOURCE, target, overlays)
    review = output / "project-context-furniture-plan-review.svg"
    write_review_crop(target, review, [overlay["fit"]["target_bbox"] for overlay in overlays])
    preview = output / "project-context-furniture-plan-review-preview.png"
    write_png_preview(review, preview)
    manifest = {
        "schema_version": 1,
        "generator": "pipeline/scripts/render_cha02_project_context.py",
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_used": False,
        "third_party_cad_used": False,
        "project_context_retained": True,
        "walls_and_surrounding_project_elements_retained": True,
        "overlay_top_layer_with_white_mask": True,
        "blue_line_present": False,
        "project_scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "views": [{
            "view": "plan",
            "candidate_view": "plan",
            "source": relative(PLAN_SOURCE),
            "output": relative(target),
            "output_sha256": sha256(target),
            "review_crop": relative(review),
            "review_crop_sha256": sha256(review),
            "review_preview": relative(preview),
            "review_preview_sha256": sha256(preview),
            "overlays": overlays,
        }],
        "pass": all(overlay["fit"]["pass"] for overlay in overlays),
    }
    manifest_path = output / "project-context-manifest.json"
    write_json(manifest_path, manifest)
    print(json.dumps({"manifest": relative(manifest_path), "output": relative(target), "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
