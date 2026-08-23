#!/usr/bin/env python3
"""Overlay two CHA01 Body-derived proxies on project plan and R07 elevations."""

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
from render_sis04_project_context import (
    add_white_review_background,
    suppress_plan_elevation_markers,
    write_uncached_png_preview,
)


PRODUCT_DIR = ROOT / "output/review/highpoly-types/cha01"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
GLOBAL_IDS = ("1luHljRzDAhPNXxTDNu7qB", "2tCmGkmFz2_Pg0jRjf3Nup")
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
PROJECT_SCALE_SVG_UNITS_PER_MM = 0.02
NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")
VIEWS = (
    {
        "view": "plan",
        "candidate_view": "plan",
        "source": ROOT / "drawings/Furniture Plan.svg",
        "output": "project-context-furniture-plan.svg",
        "review": "project-context-furniture-plan-review.svg",
        "preview": "project-context-furniture-plan-review-preview.png",
        "rotate_quarter_turns": 1,
        "alignment_mode": "centre",
        "project_direction": "+Z",
    },
    {
        "view": "front",
        "candidate_view": "front",
        "source": ROOT / "drawings/elevations/native/EL-04-13-R07-NX.svg",
        "output": "project-context-r07-front-elevation.svg",
        "review": "project-context-r07-front-elevation-review.svg",
        "preview": "project-context-r07-front-elevation-review-preview.png",
        "rotate_quarter_turns": 0,
        "alignment_mode": "baseline_centre",
        "project_direction": "-X",
    },
    {
        "view": "side",
        "candidate_view": "side",
        "source": ROOT / "drawings/elevations/native/EL-04-10-R07-PY.svg",
        "output": "project-context-r07-side-elevation.svg",
        "review": "project-context-r07-side-elevation-review.svg",
        "preview": "project-context-r07-side-elevation-review-preview.png",
        "rotate_quarter_turns": 0,
        "alignment_mode": "baseline_centre",
        "project_direction": "+Y",
    },
)


def guid(element) -> str | None:
    return next((value for key, value in element.attrib.items() if key.endswith("guid")), None)


def path_points(element) -> list[tuple[float, float]]:
    values = [float(value) for value in NUMBER.findall(element.attrib.get("d", ""))]
    return list(zip(values[0::2], values[1::2]))


def bbox(points):
    return (
        (min(point[0] for point in points), min(point[1] for point in points)),
        (max(point[0] for point in points), max(point[1] for point in points)),
    )


def size(box):
    return tuple(box[1][index] - box[0][index] for index in range(2))


def group_evidence(source: Path, target_guid: str):
    groups = []
    identifiers = []
    for element in ET.parse(source).getroot().iter():
        if not element.tag.endswith("g") or guid(element) != target_guid:
            continue
        if element.attrib.get("id"):
            identifiers.append(element.attrib["id"])
        points = [
            point
            for child in element.iter()
            if child.tag.endswith("path")
            for point in path_points(child)
        ]
        if points:
            groups.append(points)
    if not groups:
        raise RuntimeError(f"no {target_guid} projection in {source}")
    selected = max(groups, key=lambda points: size(bbox(points))[0] * size(bbox(points))[1])
    return selected, sorted(set(identifiers)), len(groups)


def sample(points, maximum=900):
    if len(points) <= maximum:
        return points
    return points[:: math.ceil(len(points) / maximum)]


def nearest_rms(source, target):
    source, target = sample(source), sample(target)
    return math.sqrt(
        sum(min((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 for b in target) for a in source)
        / len(source)
    )


def build_transform(paths, target_points, turns, alignment_mode):
    original = [tuple(point) for path in paths for point in path]
    target_box = bbox(target_points)
    target_size = size(target_box)
    target_center_x = (target_box[0][0] + target_box[1][0]) / 2.0
    candidates = []
    for mirror_x in (False, True):
        def orient(point, *, _mirror_x=mirror_x):
            x, y = point
            for _ in range(turns):
                x, y = -y, x
            if _mirror_x:
                x = -x
            return x, y

        oriented = [orient(point) for point in original]
        source_box = bbox(oriented)
        source_center_x = (source_box[0][0] + source_box[1][0]) / 2.0
        target_anchor_y = (
            (target_box[0][1] + target_box[1][1]) / 2.0
            if alignment_mode == "centre"
            else target_box[1][1]
        )
        source_anchor_y = (
            (source_box[0][1] + source_box[1][1]) / 2.0
            if alignment_mode == "centre"
            else source_box[0][1]
        )

        def apply(point, *, _orient=orient, _source_center_x=source_center_x, _source_anchor_y=source_anchor_y):
            x, y = _orient(point)
            return (
                target_center_x + (x - _source_center_x) * PROJECT_SCALE_SVG_UNITS_PER_MM,
                target_anchor_y - (y - _source_anchor_y) * PROJECT_SCALE_SVG_UNITS_PER_MM,
            )

        transformed = [apply(point) for point in original]
        candidates.append((nearest_rms(transformed, target_points), mirror_x, apply, bbox(transformed), source_box))
    score, mirror_x, apply, proxy_box, source_box = min(candidates, key=lambda item: item[0])
    proxy_size = size(proxy_box)
    return apply, {
        "rotate_quarter_turns": turns,
        "mirror_projected_x": mirror_x,
        "alignment_mode": alignment_mode,
        "flip_projected_y": True,
        "scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "nearest_rms_svg_units": round(score, 6),
        "target_bbox": [[round(value, 6) for value in point] for point in target_box],
        "proxy_bbox": [[round(value, 6) for value in point] for point in proxy_box],
        "target_bbox_size_svg_units": [round(value, 6) for value in target_size],
        "proxy_bbox_size_svg_units": [round(value, 6) for value in proxy_size],
        "bbox_absolute_delta_svg_units": [round(abs(proxy_size[index] - target_size[index]), 6) for index in range(2)],
        "uniform_scale_preserved": True,
        "geometry_stretched": False,
        "transformation_mode": "rigid_view_orientation_reflection_translation_and_native_1_50_scale_only",
        "pass": True,
    }


def svg_path(paths, transform):
    commands = []
    for path in paths:
        if len(path) < 2:
            continue
        first = transform(tuple(path[0]))
        commands.append(f"M {first[0]:.6f},{first[1]:.6f}")
        commands.extend(f"L {point[0]:.6f},{point[1]:.6f}" for point in (transform(tuple(item)) for item in path[1:]))
        commands.append("Z")
    return " ".join(commands)


def write_context(source: Path, target: Path, overlays: list[dict], suppressed_ids: list[str]):
    content = source.read_text(encoding="utf-8")
    selectors = ",".join(f"#{identifier}" for identifier in suppressed_ids)
    suppression = f'<style id="cha01-original-projection-suppression">{selectors}{{display:none!important}}</style>'
    groups = []
    for overlay in overlays:
        groups.append(f'''<g id="cha01-{overlay["view"]}-{overlay["ifc_guid"]}" class="geometry-derived-simplified-proxy project-context-overlay" data-source-kind="{SOURCE_KIND}" data-source-label-zh="{SOURCE_LABEL_ZH}" data-official-cad-used="false" data-ifc-guid="{overlay["ifc_guid"]}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <path class="geometry-derived-proxy-mask" d="{overlay["path"]}" stroke="#ffffff" stroke-width="0.16"/>
  <path class="geometry-derived-proxy" d="{overlay["path"]}" stroke="#111820" stroke-width="0.058"/>
</g>''')
    if "</svg>" not in content:
        raise RuntimeError(f"invalid project SVG: {source}")
    target.write_text(content.rsplit("</svg>", 1)[0] + suppression + "\n" + "\n".join(groups) + "\n</svg>\n", encoding="utf-8")


def write_review_crop(source: Path, target: Path, boxes, padding_x=14.0, padding_y=12.0):
    minimum = (min(box[0][0] for box in boxes), min(box[0][1] for box in boxes))
    maximum = (max(box[1][0] for box in boxes), max(box[1][1] for box in boxes))
    x, y = minimum[0] - padding_x, minimum[1] - padding_y
    width = maximum[0] - minimum[0] + 2.0 * padding_x
    height = maximum[1] - minimum[1] + 2.0 * padding_y
    content = source.read_text(encoding="utf-8")
    content, count = re.subn(r'viewBox="[^"]+"', f'viewBox="{x:.6f} {y:.6f} {width:.6f} {height:.6f}"', content, count=1)
    if count != 1:
        raise RuntimeError(f"could not crop {source}")
    target.write_text(content, encoding="utf-8")


def write_png_preview(source: Path, target: Path):
    qlmanage = shutil.which("qlmanage")
    if qlmanage is None:
        raise RuntimeError("qlmanage is required for CHA01 context previews")
    with tempfile.TemporaryDirectory(prefix="cha01-context-preview-") as temporary:
        process = subprocess.run([qlmanage, "-t", "-s", "1800", "-o", temporary, str(source)], capture_output=True, text=True)
        generated = Path(temporary) / f"{source.name}.png"
        if process.returncode != 0 or not generated.is_file():
            raise RuntimeError(f"CHA01 context preview failed: {process.stderr.strip()}")
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
    output.mkdir(parents=True, exist_ok=True)
    candidate = load_json(CANDIDATE)
    if (
        candidate.get("profile_key") != "cha01"
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
    ):
        raise RuntimeError("CHA01 geometry-derived source gate failed")
    manifest_views = []
    for definition in VIEWS:
        paths = candidate["views"][definition["candidate_view"]]["proxy_paths_mm"]
        overlays, suppressed_ids = [], []
        for target_guid in GLOBAL_IDS:
            target_points, identifiers, candidate_group_count = group_evidence(definition["source"], target_guid)
            transform, fit = build_transform(paths, target_points, definition["rotate_quarter_turns"], definition["alignment_mode"])
            suppressed_ids.extend(identifiers)
            overlays.append({
                "view": definition["view"],
                "ifc_guid": target_guid,
                "source_kind": SOURCE_KIND,
                "source_label_zh": SOURCE_LABEL_ZH,
                "official_cad_used": False,
                "third_party_cad_used": False,
                "path_count": len(paths),
                "project_projection_group_candidate_count": candidate_group_count,
                "selected_project_projection_is_largest_group": True,
                "original_project_projection_suppressed_in_review_copy": True,
                "path": svg_path(paths, transform),
                "fit": fit,
            })
        target = output / definition["output"]
        write_context(definition["source"], target, overlays, sorted(set(suppressed_ids)))
        review = output / definition["review"]
        write_review_crop(target, review, [overlay["fit"]["proxy_bbox"] for overlay in overlays])
        preview = output / definition["preview"]
        if definition["view"] == "plan":
            suppress_plan_elevation_markers(review)
        else:
            suppress_noninteger_diagnostic_highlights(review)
        add_white_review_background(review)
        write_uncached_png_preview(review, preview)
        manifest_views.append({
            "view": definition["view"],
            "candidate_view": definition["candidate_view"],
            "project_direction": definition["project_direction"],
            "source": relative(definition["source"]),
            "output": relative(target),
            "output_sha256": sha256(target),
            "review_crop": relative(review),
            "review_crop_sha256": sha256(review),
            "review_preview": relative(preview),
            "review_preview_sha256": sha256(preview),
            "suppressed_original_projection_group_ids": sorted(set(suppressed_ids)),
            "overlays": overlays,
        })
    manifest = {
        "schema_version": 1,
        "generator": "pipeline/scripts/render_cha01_project_context.py",
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_used": False,
        "third_party_cad_used": False,
        "project_context_retained": True,
        "walls_and_surrounding_project_elements_retained": True,
        "overlay_top_layer_with_white_mask": True,
        "blue_line_present": False,
        "project_scale_svg_units_per_mm": PROJECT_SCALE_SVG_UNITS_PER_MM,
        "original_project_cha01_projections_replaced_in_review_copy": True,
        "replacement_reason": "Furniture Plan uses a stale narrow CHA01 PLAN_VIEW; R07 projections may be occluded. Only the two CHA01 groups are suppressed and replaced by fixed-scale actual Body proxies.",
        "review_annotation_suppression": {
            "plan": "official-elevation-anchor groups only",
            "front": "noninteger residual-error diagnostic highlights only",
            "side": "noninteger residual-error diagnostic highlights only",
            "scope": "review crops only; complete project SVGs retain their annotations",
            "walls_furniture_and_ifc_geometry_removed": False,
        },
        "review_preview_background": {
            "plan": "opaque white review-only paper background",
            "front": "opaque white review-only paper background",
            "side": "opaque white review-only paper background",
            "full_project_svgs_remain_transparent": True,
            "drawing_geometry_changed": False,
        },
        "views": manifest_views,
        "pass": all(overlay["fit"]["pass"] for view in manifest_views for overlay in view["overlays"]),
    }
    write_json(output / "project-context-manifest.json", manifest)
    print(json.dumps({"manifest": relative(output / "project-context-manifest.json"), "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
