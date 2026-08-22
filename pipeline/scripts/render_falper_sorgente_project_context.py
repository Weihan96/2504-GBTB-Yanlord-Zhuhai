#!/usr/bin/env python3
"""Overlay approved Falper WFB drawing representations on project-context SVGs."""

from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

import ifcopenshell

from falper_sorgente_linework import EXPECTED, ROOT, load_json, relative, sha256


GLOBAL_ID = "350tdaubr8QP3Cu2YMQZIN"
BLUE = "#1677c8"
IFC_NAMESPACE = "http://www.ifcopenshell.org/ns"
SOURCE_PLAN = ROOT / "drawings/Sanitary Plan.svg"
SOURCE_ELEVATION = ROOT / "drawings/elevations/native/EL-08-30-R16-NY.svg"
REPRESENTATION_IFC = (
    ROOT
    / "output/review/highpoly-types/falper-sorgente/"
    "Falper-Sorgente-WFB-bonsai-isolated.ifc"
)
REVIEW_MANIFEST = (
    ROOT
    / "output/review/highpoly-types/falper-sorgente/"
    "bonsai-review-manifest.json"
)
OUTPUT_DIR = ROOT / "output/review/highpoly-types/falper-sorgente"


def representation_paths(product, identifier: str) -> list[list[tuple[float, float]]]:
    representation = next(
        (
            item
            for item in product.Representation.Representations
            if item.RepresentationIdentifier == identifier
        ),
        None,
    )
    if representation is None or representation.RepresentationType != "GeometricCurveSet":
        raise RuntimeError(f"missing approved {identifier} GeometricCurveSet")
    axes = {
        "FalperWFBPlan": (0, 1),
        "FalperWFBFront": (0, 2),
        "FalperWFBSide": (1, 2),
    }[identifier]
    return [
        [
            (
                float(point.Coordinates[axes[0]]),
                float(point.Coordinates[axes[1]]),
            )
            for point in polyline.Points
        ]
        for curve_set in representation.Items
        for polyline in curve_set.Elements
    ]


def product_group_pattern() -> re.Pattern[str]:
    return re.compile(
        rf'<g\b(?=[^>]*\bifc:guid="{re.escape(GLOBAL_ID)}")[^>]*>.*?</g>',
        re.DOTALL,
    )


def group_bounds(group: str) -> tuple[float, float, float, float]:
    points = [
        (float(x), float(y))
        for x, y in re.findall(
            r"[ML]\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)", group
        )
    ]
    if not points:
        raise RuntimeError("could not recover Bonsai product placement from SVG group")
    return (
        min(point[0] for point in points),
        min(point[1] for point in points),
        max(point[0] for point in points),
        max(point[1] for point in points),
    )


def svg_path(points: list[tuple[float, float]]) -> str:
    return " ".join(
        f"{'M' if index == 0 else 'L'} {x:.6f},{y:.6f}"
        for index, (x, y) in enumerate(points)
    )


def transform_plan(
    paths: list[list[tuple[float, float]]], center_x: float, center_y: float
) -> list[list[tuple[float, float]]]:
    return [
        [(center_x + x / 50.0, center_y - y / 50.0) for x, y in path]
        for path in paths
    ]


def transform_elevation(
    paths: list[list[tuple[float, float]]], center_x: float, baseline_y: float
) -> list[list[tuple[float, float]]]:
    return [
        [(center_x + x / 50.0, baseline_y - y / 50.0) for x, y in path]
        for path in paths
    ]


def elevation_mask(paths: list[list[tuple[float, float]]]) -> str:
    profiles = sorted(paths, key=len, reverse=True)[:2]
    if len(profiles) != 2:
        raise RuntimeError("official elevation does not contain two side profiles")
    profiles.sort(key=lambda path: sum(point[0] for point in path) / len(path))
    left, right = profiles
    if left[0][1] < left[-1][1]:
        left = list(reversed(left))
    if right[0][1] > right[-1][1]:
        right = list(reversed(right))
    return svg_path(left + right + [left[0]]) + " Z"


def replace_and_append(source: Path, output: Path, overlay: str) -> dict:
    text = source.read_text(encoding="utf-8")
    match = product_group_pattern().search(text)
    if match is None:
        raise RuntimeError(f"Falper {GLOBAL_ID} is not present in {source}")
    bounds = group_bounds(match.group(0))
    marker = (
        f'<g id="falper-original-body-projection-removed" '
        f'data-guid="{GLOBAL_ID}" data-replaced-by="official-native-dwg"/>'
    )
    text = text[: match.start()] + marker + text[match.end() :]
    if "</svg>" not in text:
        raise RuntimeError(f"invalid SVG: {source}")
    text = text.replace("</svg>", overlay + "\n</svg>", 1)
    output.write_text(text, encoding="utf-8")
    return {
        "source": source.relative_to(ROOT).as_posix(),
        "output": relative(output),
        "source_product_bounds": [round(value, 6) for value in bounds],
        "output_sha256": sha256(output),
    }


def write_review_crop(
    source: Path, output: Path, view_box: tuple[float, float, float, float]
) -> dict:
    text = source.read_text(encoding="utf-8")
    x, y, width, height = view_box
    root_end = text.find(">")
    if root_end < 0:
        raise RuntimeError(f"invalid SVG root: {source}")
    root = text[:root_end]
    root = re.sub(
        r'viewBox="[^"]+"',
        f'viewBox="{x:.6f} {y:.6f} {width:.6f} {height:.6f}"',
        root,
        count=1,
    )
    root = re.sub(r'width="[^"]+"', f'width="{width:.6f}mm"', root, count=1)
    root = re.sub(r'height="[^"]+"', f'height="{height:.6f}mm"', root, count=1)
    output.write_text(root + text[root_end:], encoding="utf-8")
    return {
        "path": relative(output),
        "sha256": sha256(output),
        "view_box": [round(value, 6) for value in view_box],
        "purpose": "visual_review_crop_with_project_context",
    }


def render_plan(paths: list[list[tuple[float, float]]]) -> dict:
    source_text = SOURCE_PLAN.read_text(encoding="utf-8")
    match = product_group_pattern().search(source_text)
    if match is None:
        raise RuntimeError("Falper is missing from the project sanitary plan")
    min_x, min_y, max_x, max_y = group_bounds(match.group(0))
    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    transformed = transform_plan(paths, center_x, center_y)
    path_elements = "\n".join(
        f'  <path d="{svg_path(path)}"/>' for path in transformed
    )
    overlay = f"""
<g id="falper-wfb-plan-approved-overlay" class="IfcSanitaryTerminal official-native-dwg" ifc:name="Falper Sorgente WFB" ifc:guid="{GLOBAL_ID}" xmlns:ifc="{IFC_NAMESPACE}" data-representation="FalperWFBPlan" data-source-kind="native_dwg" data-source-sha256="{EXPECTED['wfb_2d']}" fill="none" stroke="{BLUE}" stroke-width="0.32" stroke-linecap="round" stroke-linejoin="round">
  <circle cx="{center_x:.6f}" cy="{center_y:.6f}" r="5.28" fill="white" stroke="none"/>
{path_elements}
</g>"""
    output = OUTPUT_DIR / "project-context-sanitary-plan.svg"
    record = replace_and_append(SOURCE_PLAN, output, overlay)
    review_crop = write_review_crop(
        output,
        OUTPUT_DIR / "project-context-sanitary-plan-review.svg",
        (center_x - 30.0, center_y - 30.0, 60.0, 60.0),
    )
    record.update(
        {
            "view": "plan",
            "representation": "FalperWFBPlan",
            "native_dwg_path_count": len(paths),
            "placement_center": [round(center_x, 6), round(center_y, 6)],
            "white_mask": True,
            "blue_top_layer": True,
            "review_crop": review_crop,
        }
    )
    return record


def render_elevation(paths: list[list[tuple[float, float]]]) -> dict:
    source_text = SOURCE_ELEVATION.read_text(encoding="utf-8")
    match = product_group_pattern().search(source_text)
    if match is None:
        raise RuntimeError("Falper is missing from project elevation EL-08-30-R16-NY")
    min_x, _min_y, max_x, max_y = group_bounds(match.group(0))
    center_x = (min_x + max_x) / 2.0
    # The Body projection includes the high-poly mesh's -8 mm local base. Moving
    # that amount upward preserves the IFC placement while using the WFB 0 mm base.
    baseline_y = max_y - 8.0 / 50.0
    transformed = transform_elevation(paths, center_x, baseline_y)
    path_elements = "\n".join(
        f'  <path d="{svg_path(path)}"/>' for path in transformed
    )
    overlay = f"""
<g id="falper-wfb-front-approved-overlay" class="IfcSanitaryTerminal official-native-dwg" ifc:name="Falper Sorgente WFB" ifc:guid="{GLOBAL_ID}" xmlns:ifc="{IFC_NAMESPACE}" data-representation="FalperWFBFront" data-source-kind="native_dwg" data-source-sha256="{EXPECTED['wfb_2d']}" fill="none" stroke="{BLUE}" stroke-width="0.32" stroke-linecap="round" stroke-linejoin="round">
  <path d="{elevation_mask(transformed)}" fill="white" stroke="none"/>
{path_elements}
</g>"""
    output = OUTPUT_DIR / "project-context-elevation-EL-08-30-R16-NY.svg"
    record = replace_and_append(SOURCE_ELEVATION, output, overlay)
    review_crop = write_review_crop(
        output,
        OUTPUT_DIR / "project-context-elevation-EL-08-30-R16-NY-review.svg",
        (center_x - 8.0, 29.0, 20.0, 25.0),
    )
    record.update(
        {
            "view": "front_elevation",
            "representation": "FalperWFBFront",
            "native_dwg_path_count": len(paths),
            "placement_center_x": round(center_x, 6),
            "placement_baseline_y": round(baseline_y, 6),
            "white_mask": True,
            "blue_top_layer": True,
            "review_crop": review_crop,
        }
    )
    return record


def main() -> None:
    global OUTPUT_DIR
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    OUTPUT_DIR = args.output_dir.resolve()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    manifest = load_json(REVIEW_MANIFEST)
    if sha256(ROOT / "2504 GBTB Yanlord Zhuhai.ifc") != EXPECTED["formal_ifc"]:
        raise RuntimeError("formal IFC bytes changed")
    if (
        manifest.get("pass") is not True
        or manifest.get("review_status") != "approved"
        or manifest.get("proxy_geometry_included_in_drawing_representations") is not False
        or manifest.get("source_kind") != "native_dwg"
        or sha256(REPRESENTATION_IFC) != manifest.get("isolated_ifc", {}).get("sha256")
    ):
        raise RuntimeError("approved isolated representation IFC gate failed")

    model = ifcopenshell.open(REPRESENTATION_IFC)
    product = model.by_guid(GLOBAL_ID)
    plan = representation_paths(product, "FalperWFBPlan")
    front = representation_paths(product, "FalperWFBFront")
    if len(plan) != 5 or len(front) != 4:
        raise RuntimeError("approved native DWG representation path counts changed")
    records = [render_plan(plan), render_elevation(front)]

    manifest_output = OUTPUT_DIR / "project-context-manifest.json"
    payload = {
        "schema_version": 1,
        "status": "approved",
        "formal_ifc_write": False,
        "formal_ifc_sha256": EXPECTED["formal_ifc"],
        "representation_ifc": REPRESENTATION_IFC.relative_to(ROOT).as_posix(),
        "representation_ifc_sha256": manifest["isolated_ifc"]["sha256"],
        "representative_global_id": GLOBAL_ID,
        "source_kind": "native_dwg",
        "source_dwg_sha256": EXPECTED["wfb_2d"],
        "scope": "family_reference_not_project_shop_drawing",
        "rendering_rule": "project context retained; original Body projection replaced by approved native-DWG representation; white mask below official blue top layer",
        "approval": {
            "reviewer": "project_owner",
            "review_date": "2026-08-23",
            "status": "approved",
        },
        "views": records,
    }
    from falper_sorgente_linework import write_json

    write_json(manifest_output, payload)
    print(html.escape(relative(manifest_output)))


if __name__ == "__main__":
    main()
