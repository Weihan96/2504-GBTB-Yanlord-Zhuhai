#!/usr/bin/env python3
"""Overlay the geometry-derived Hima proxy on complete project drawings."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json


PRODUCT_DIR = ROOT / "output/review/highpoly-types/hima01"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-05-14-R09-PY.svg"
SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-05-17-R09-NX.svg"
GLOBAL_ID = "2xmcLzu1rDTeMzRuNxPDyE"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
SOURCE_DWG_SHA256 = ""
OFFICIAL_CAD_USED = False
PATH_KEY = "proxy_paths_mm"
CLOSE_PATHS = True
BLUE = "#1677c8"
NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")
PROJECT_SCALE_SVG_UNITS_PER_MM = 0.02
BBOX_TOLERANCE_SVG_UNITS = 0.08
ALIGNMENT_MODE = "minimum"
OVERLAY_ID_PREFIX = "hima01"
GENERATOR = "pipeline/scripts/render_hima01_project_context.py"
REVIEW_HIDDEN_CSS = ""
CONTEXT_VIEWS = (
    ("plan", "plan", PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 0),
    ("front", "front", FRONT_SOURCE, "project-context-front-elevation.svg", (8.0, 6.0), 0),
    ("side", "side", SIDE_SOURCE, "project-context-side-elevation.svg", (10.0, 6.0), 0),
)


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


def transform_for_context(paths, target_points, rotate_quarter_turns=0):
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
        raise RuntimeError(f"Hima proxy/project projection bbox mismatch: {deltas}")

    # Project drawings use downward-positive SVG Y. Reflect the second projected
    # coordinate and align the complete proxy bbox with the original Hima group.
    def apply(point):
        point = rotate(point)
        if ALIGNMENT_MODE == "centre":
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
                target_centre[1] - (point[1] - source_centre[1]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
            )
        return (
            target_min[0] + (point[0] - source_min[0]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
            target_max[1] - (point[1] - source_min[1]) * PROJECT_SCALE_SVG_UNITS_PER_MM,
        )

    transformed_bbox = bbox([apply(point) for point in original_raw])
    return apply, {
        "rotate_quarter_turns": turns,
        "alignment_mode": ALIGNMENT_MODE,
        "flip_projected_y": True,
        "scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "target_bbox": [[round(value, 6) for value in target_min], [round(value, 6) for value in target_max]],
        "proxy_bbox": [[round(value, 6) for value in transformed_bbox[0]], [round(value, 6) for value in transformed_bbox[1]]],
        "bbox_absolute_delta_svg_units": [round(value, 6) for value in deltas],
        "bbox_tolerance_svg_units": BBOX_TOLERANCE_SVG_UNITS,
        "uniform_scale_preserved": True,
        "transformation_mode": "axis_swap_rigid_reflection_and_translation_only",
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
        if CLOSE_PATHS:
            commands.append("Z")
    return " ".join(commands)


def rebase_external_resources(content: str, source: Path, target: Path) -> str:
    def replace(match: re.Match[str]) -> str:
        attribute, reference = match.groups()
        if reference.startswith(("#", "data:", "http://", "https://")):
            return match.group(0)
        resolved = (source.parent / reference).resolve()
        if not resolved.is_file():
            raise RuntimeError(f"missing project SVG resource: {resolved}")
        rebased = Path(os.path.relpath(resolved, target.parent)).as_posix()
        return f'{attribute}="{rebased}"'

    return re.sub(r'((?:xlink:)?href)="([^"]+)"', replace, content)


def add_overlay(source: Path, target: Path, view: str, path: str, overlay_bbox):
    content = source.read_text(encoding="utf-8")
    if OFFICIAL_CAD_USED:
        minimum, maximum = overlay_bbox
        padding = 0.18
        mask_x = minimum[0] - padding
        mask_y = minimum[1] - padding
        mask_width = maximum[0] - minimum[0] + 2 * padding
        mask_height = maximum[1] - minimum[1] + 2 * padding
        group = f'''<g id="{OVERLAY_ID_PREFIX}-{view}-{GLOBAL_ID}" class="official-reference native-dwg project-context-overlay" data-source-kind="{SOURCE_KIND}" data-source-label-zh="{SOURCE_LABEL_ZH}" data-source-dwg-sha256="{SOURCE_DWG_SHA256}" data-official-cad-used="true" data-ifc-guid="{GLOBAL_ID}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <rect class="official-reference-envelope-mask" x="{mask_x:.6f}" y="{mask_y:.6f}" width="{mask_width:.6f}" height="{mask_height:.6f}" fill="#ffffff" stroke="none"/>
  <path class="official-reference-mask" d="{path}" stroke="#ffffff" stroke-width="0.15"/>
  <path class="official-reference native-dwg" d="{path}" stroke="{BLUE}" stroke-width="0.058"/>
</g>'''
    else:
        group = f'''<g id="{OVERLAY_ID_PREFIX}-{view}-{GLOBAL_ID}" class="geometry-derived-simplified-proxy project-context-overlay" data-source-kind="{SOURCE_KIND}" data-source-label-zh="{SOURCE_LABEL_ZH}" data-official-cad-used="false" data-ifc-guid="{GLOBAL_ID}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <path class="geometry-derived-proxy-mask" d="{path}" stroke="#ffffff" stroke-width="0.15"/>
  <path class="geometry-derived-proxy" d="{path}" stroke="#111820" stroke-width="0.058"/>
</g>'''
    if "</svg>" not in content:
        raise RuntimeError(f"invalid project SVG: {source}")
    content = content.rsplit("</svg>", 1)[0] + group + "\n</svg>\n"
    target.write_text(rebase_external_resources(content, source, target), encoding="utf-8")


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
    if REVIEW_HIDDEN_CSS:
        root_end = content.find(">")
        if root_end < 0:
            raise RuntimeError(f"invalid SVG root: {source}")
        content = content[: root_end + 1] + f"<style>{REVIEW_HIDDEN_CSS}</style>" + content[root_end + 1 :]
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=PRODUCT_DIR)
    args = parser.parse_args()
    output = args.output.resolve()
    candidate = load_json(CANDIDATE)
    if (
        candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("official_cad_used") is not OFFICIAL_CAD_USED
        or candidate.get("third_party_cad_used") is not False
        or (OFFICIAL_CAD_USED and candidate.get("source_dwg_sha256") != SOURCE_DWG_SHA256)
    ):
        raise RuntimeError(f"{OVERLAY_ID_PREFIX} project-context source gate failed")
    records = []
    for view, candidate_view, source, filename, crop_padding, rotate_quarter_turns in CONTEXT_VIEWS:
        target_points = group_points(source)
        paths = candidate["views"][candidate_view][PATH_KEY]
        transform, fit = transform_for_context(paths, target_points, rotate_quarter_turns)
        overlay_path = svg_path(paths, transform)
        target = output / filename
        add_overlay(source, target, view, overlay_path, fit["proxy_bbox"])
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
                "official_cad_used": OFFICIAL_CAD_USED,
                "third_party_cad_used": False,
                "source_dwg_sha256": SOURCE_DWG_SHA256 if OFFICIAL_CAD_USED else None,
                "path_count": len(paths),
                "fit": fit,
            },
        })
    manifest = {
        "schema_version": 1,
        "generator": GENERATOR,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_used": OFFICIAL_CAD_USED,
        "third_party_cad_used": False,
        "project_context_retained": True,
        "walls_and_surrounding_project_elements_retained": True,
        "overlay_top_layer_with_white_mask": True,
        "blue_line_present": OFFICIAL_CAD_USED,
        "project_scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "views": records,
        "pass": all(record["overlay"]["fit"]["pass"] for record in records),
    }
    target = output / "project-context-manifest.json"
    write_json(target, manifest)
    print(json.dumps({"manifest": relative(target), "outputs": [record["output"] for record in records], "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
