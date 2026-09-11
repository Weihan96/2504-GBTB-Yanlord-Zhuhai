#!/usr/bin/env python3
"""Overlay exact archived CleanLine50 DWG views on project-context SVGs."""

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


ARTICLE = "154.446.KS.1"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-154-446-ks-1"
LINEWORK = PRODUCT_DIR / "official-native-dwg-linework.json"
PLAN_SOURCE = ROOT / "drawings/Sanitary Plan-P202-candidate.svg"
FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-06-21-R12-PX.svg"
SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-06-20-R12-NY.svg"
GLOBAL_ID = "2S2c498tb7$gzdukjhCGVQ"
NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")
BLUE = "#1677c8"
PROJECT_SCALE_SVG_UNITS_PER_MM = 0.02


def guid(element) -> str | None:
    return next((value for key, value in element.attrib.items() if key.endswith("guid")), None)


def path_points(element) -> list[tuple[float, float]]:
    values = [float(value) for value in NUMBER.findall(element.attrib.get("d", ""))]
    return list(zip(values[0::2], values[1::2]))


def group_points(source: Path) -> list[tuple[float, float]]:
    groups = []
    for element in ET.parse(source).getroot().iter():
        if not element.tag.endswith("g") or guid(element) != GLOBAL_ID:
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
        raise RuntimeError(f"expected one {GLOBAL_ID} projection in {source}, found {len(groups)}")
    return groups[0]


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
    source = sampled(source)
    target = sampled(target)
    return math.sqrt(
        sum(min((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for b in target) for a in source)
        / len(source)
    )


def transform_for_context(paths, target_points, view):
    raw = [tuple(point) for path in paths for point in path]
    target_min, target_max = bbox(target_points)
    target_center = tuple((target_min[axis] + target_max[axis]) / 2.0 for axis in range(2))
    candidates = []
    expected_swap = view == "plan"
    for swap in (False, True):
        if swap != expected_swap:
            continue
        oriented = [(y, x) if swap else (x, y) for x, y in raw]
        oriented_min, oriented_max = bbox(oriented)
        oriented_center = tuple((oriented_min[axis] + oriented_max[axis]) / 2.0 for axis in range(2))
        flip_x_values = (True,) if view == "side" else (False, True)
        for flip_x in flip_x_values:
            for flip_y in (False, True):
                x_anchors = ("native_origin",) if view == "side" else ("minimum", "center", "maximum")
                for x_anchor in x_anchors:
                    for y_anchor in ("minimum", "center", "maximum"):
                        def axis_offset(axis, anchor):
                            if anchor == "native_origin":
                                native_origin = oriented_min[axis] + oriented_max[axis] if flip_x else 0.0
                                return target_center[axis] - native_origin * PROJECT_SCALE_SVG_UNITS_PER_MM
                            if anchor == "minimum":
                                return target_min[axis] - oriented_min[axis] * PROJECT_SCALE_SVG_UNITS_PER_MM
                            if anchor == "maximum":
                                return target_max[axis] - oriented_max[axis] * PROJECT_SCALE_SVG_UNITS_PER_MM
                            return target_center[axis] - oriented_center[axis] * PROJECT_SCALE_SVG_UNITS_PER_MM

                        offset = (axis_offset(0, x_anchor), axis_offset(1, y_anchor))

                        def apply(point, *, _swap=swap, _fx=flip_x, _fy=flip_y, _offset=offset):
                            x, y = (point[1], point[0]) if _swap else point
                            if _fx:
                                x = oriented_min[0] + oriented_max[0] - x
                            if _fy:
                                y = oriented_min[1] + oriented_max[1] - y
                            return (
                                _offset[0] + x * PROJECT_SCALE_SVG_UNITS_PER_MM,
                                _offset[1] + y * PROJECT_SCALE_SVG_UNITS_PER_MM,
                            )

                        transformed = [apply(point) for point in raw]
                        score = nearest_rms(target_points, transformed)
                        candidates.append((score, swap, flip_x, flip_y, x_anchor, y_anchor, apply))
    score, swap, flip_x, flip_y, x_anchor, y_anchor, apply = min(candidates, key=lambda item: item[0])
    transformed_bbox = bbox([apply(point) for point in raw])
    return apply, {
        "swap_axes": swap,
        "flip_x": flip_x,
        "flip_y": flip_y,
        "x_anchor": x_anchor,
        "y_anchor": y_anchor,
        "fit_rms_svg_units": round(score, 6),
        "scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "uniform_scale_preserved": True,
        "transformation_mode": "axis_swap_rigid_reflection_and_translation_only",
        "target_bbox": [[round(value, 6) for value in target_min], [round(value, 6) for value in target_max]],
        "official_overlay_bbox": [
            [round(value, 6) for value in transformed_bbox[0]],
            [round(value, 6) for value in transformed_bbox[1]],
        ],
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


def add_overlay(source: Path, target: Path, overlay: dict):
    content = source.read_text(encoding="utf-8")
    path = overlay["path"]
    group = f'''<g id="geberit-154-446-ks-1-{overlay["view"]}-{GLOBAL_ID}" class="official-native-dwg project-context-overlay" data-source-kind="native_dwg" data-article-number="{ARTICLE}" data-native-dwg-code="{overlay["code"]}" data-source-sha256="{overlay["sha256"]}" data-ifc-guid="{GLOBAL_ID}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <path class="official-reference-mask" d="{path}" stroke="#ffffff" stroke-width="0.13"/>
  <path class="official-reference native-dwg" d="{path}" stroke="{BLUE}" stroke-width="0.052"/>
</g>'''
    if "</svg>" not in content:
        raise RuntimeError(f"invalid project SVG: {source}")
    target.write_text(content.rsplit("</svg>", 1)[0] + group + "\n</svg>\n", encoding="utf-8")


def write_review_crop(source: Path, target: Path, target_bbox, padding_x, padding_y):
    minimum, maximum = target_bbox
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
        raise RuntimeError(f"could not crop {source}")
    target.write_text(content, encoding="utf-8")


def write_png_preview(source: Path, target: Path):
    qlmanage = shutil.which("qlmanage")
    if qlmanage is None:
        raise RuntimeError("qlmanage is required to render CleanLine50 project-context previews")
    with tempfile.TemporaryDirectory(prefix="geberit-154446ks1-context-preview-") as temporary:
        process = subprocess.run(
            [qlmanage, "-t", "-s", "1800", "-o", temporary, str(source)],
            capture_output=True,
            text=True,
        )
        generated = Path(temporary) / f"{source.name}.png"
        if process.returncode != 0 or not generated.is_file():
            raise RuntimeError(f"CleanLine50 context preview failed: {process.stderr.strip()}")
        generated.replace(target)


def build_overlay(source: Path, view: str, linework: dict):
    target_points = group_points(source)
    view_data = linework["views"][view]
    transform, fit = transform_for_context(view_data["paths_mm"], target_points, view)
    return {
        "guid": GLOBAL_ID,
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
        raise RuntimeError("CleanLine50 native-DWG identity gate failed")
    if linework.get("retirement_identity", {}).get("replacement_cad_used") is not False:
        raise RuntimeError("replacement article CAD must not be used")
    records = []
    for view, source, filename, crop_padding in (
        ("plan", PLAN_SOURCE, "project-context-sanitary-plan.svg", (3.0, 3.0)),
        ("front", FRONT_SOURCE, "project-context-front-elevation.svg", (6.0, 4.0)),
        ("side", SIDE_SOURCE, "project-context-side-elevation.svg", (3.0, 3.0)),
    ):
        overlay = build_overlay(source, view, linework)
        target = output / filename
        add_overlay(source, target, overlay)
        crop = output / filename.replace(".svg", "-review.svg")
        write_review_crop(target, crop, overlay["fit"]["target_bbox"], *crop_padding)
        preview = output / filename.replace(".svg", "-review-preview.png")
        write_png_preview(crop, preview)
        records.append({
            "view": view,
            "source": relative(source),
            "output": relative(target),
            "output_sha256": sha256(target),
            "review_crop": relative(crop),
            "review_crop_sha256": sha256(crop),
            "review_preview": relative(preview),
            "review_preview_sha256": sha256(preview),
            "overlays": [overlay],
        })
    manifest = {
        "schema_version": 1,
        "generator": "pipeline/scripts/render_geberit_154_446_ks_1_project_context.py",
        "article_number": ARTICLE,
        "source_kind": "native_dwg",
        "source_linework": relative(LINEWORK),
        "source_linework_sha256": sha256(LINEWORK),
        "project_context_retained": True,
        "walls_and_surrounding_project_elements_retained": True,
        "blue_line_top_layer_with_white_mask": True,
        "replacement_cad_used": False,
        "official_dwg_uniform_project_scale_preserved": True,
        "project_scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "views": records,
        "pass": all(
            record["overlays"][0]["native_dwg_path_count"] > 0
            and record["overlays"][0]["fit"]["uniform_scale_preserved"]
            for record in records
        ),
    }
    target = output / "project-context-manifest.json"
    write_json(target, manifest)
    print(json.dumps({"manifest": relative(target), "outputs": [record["output"] for record in records], "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
