#!/usr/bin/env python3
"""Extract WFA/WFB vector linework from the project-archived official PDF."""

from __future__ import annotations

import argparse
from pathlib import Path

import pdfplumber

from falper_sorgente_linework import EXPECTED, ROOT, SCOPE, relative, sha256, write_json


DEFAULT_INPUT = ROOT / "drawings/evidence/FALPER-official-Sorgente-WFA-WFB.pdf"
DEFAULT_OUTPUT = ROOT / "pipeline/decisions/falper-sorgente-official-pdf-linework.json"


def point(value) -> tuple[float, float]:
    return float(value[0]), float(value[1])


def cubic(a, b, c, d, steps: int = 48):
    values = []
    for index in range(1, steps + 1):
        t = index / steps
        u = 1.0 - t
        values.append(
            (
                u**3 * a[0] + 3 * u**2 * t * b[0] + 3 * u * t**2 * c[0] + t**3 * d[0],
                u**3 * a[1] + 3 * u**2 * t * b[1] + 3 * u * t**2 * c[1] + t**3 * d[1],
            )
        )
    return values


def flatten_curve(curve: dict) -> list[list[list[float]]]:
    paths: list[list[tuple[float, float]]] = []
    current: list[tuple[float, float]] = []
    start = None
    cursor = None
    for command in curve["path"]:
        operation = command[0]
        if operation == "m":
            if current:
                paths.append(current)
            cursor = point(command[1])
            start = cursor
            current = [cursor]
        elif operation == "l":
            cursor = point(command[1])
            current.append(cursor)
        elif operation == "c":
            if cursor is None:
                raise RuntimeError("Bezier curve starts without a move command")
            control_1, control_2, end = point(command[1]), point(command[2]), point(command[3])
            current.extend(cubic(cursor, control_1, control_2, end))
            cursor = end
        elif operation == "h":
            if start is not None and cursor != start:
                current.append(start)
                cursor = start
        else:
            raise RuntimeError(f"unsupported PDF path operation: {operation}")
    if current:
        paths.append(current)
    return [
        [[round(x, 9), round(y, 9)] for x, y in path]
        for path in paths
        if len(path) >= 2
    ]


def is_stroked_unfilled(curve: dict) -> bool:
    return bool(curve.get("stroke")) and not bool(curve.get("fill"))


def select_plan(page) -> list[dict]:
    candidates = [
        curve
        for curve in page.curves
        if is_stroked_unfilled(curve)
        and 90 <= curve["x0"]
        and curve["x1"] <= 192
        and 310 <= curve["top"]
        and curve["bottom"] <= 415
        and curve["width"] >= 10
        and curve["height"] >= 10
    ]
    if len(candidates) != 5:
        raise RuntimeError(f"expected five official PDF plan curves, found {len(candidates)}")
    return sorted(candidates, key=lambda item: item["width"], reverse=True)


def select_elevation(page) -> tuple[list[dict], list[dict]]:
    sides = [
        curve
        for curve in page.curves
        if is_stroked_unfilled(curve)
        and curve["x1"] < 200
        and 68 <= curve["top"] <= 76
        and 150 <= curve["height"] <= 165
        and 30 <= curve["width"] <= 40
    ]
    if len(sides) != 2:
        raise RuntimeError(f"expected two official PDF elevation side curves, found {len(sides)}")
    top = min(curve["top"] for curve in sides)
    bottom = max(curve["bottom"] for curve in sides)
    horizontal = [
        line
        for line in page.lines
        if 90 <= line["x0"]
        and line["x1"] <= 195
        and abs(line["top"] - line["bottom"]) < 0.05
        and (
            (80 <= line["width"] <= 95 and abs(line["top"] - top) <= 0.5)
            or (35 <= line["width"] <= 45 and abs(line["top"] - bottom) <= 0.5)
        )
    ]
    if len(horizontal) != 2:
        raise RuntimeError(f"expected official PDF elevation top/base lines, found {len(horizontal)}")
    return sorted(sides, key=lambda item: item["x0"]), sorted(horizontal, key=lambda item: item["top"])


def normalise_plan(curves: list[dict]) -> tuple[list, dict]:
    outer = curves[0]
    centre_x = (outer["x0"] + outer["x1"]) / 2.0
    centre_y = (outer["top"] + outer["bottom"]) / 2.0
    scale = 520.0 / outer["width"]
    paths = []
    for curve in curves:
        for path in flatten_curve(curve):
            paths.append(
                [
                    [round((x - centre_x) * scale, 6), round((centre_y - y) * scale, 6)]
                    for x, y in path
                ]
            )
    return paths, {
        "published_outer_diameter_mm": 520.0,
        "pdf_outer_bbox_points": [outer["x0"], outer["top"], outer["x1"], outer["bottom"]],
        "points_to_mm_scale": round(scale, 12),
    }


def normalise_elevation(sides: list[dict], horizontal: list[dict]) -> tuple[list, dict]:
    left = min(item["x0"] for item in sides)
    right = max(item["x1"] for item in sides)
    top = min(item["top"] for item in sides)
    bottom = max(item["bottom"] for item in sides)
    centre_x = (left + right) / 2.0
    scale_x = 520.0 / (right - left)
    scale_z = 900.0 / (bottom - top)
    paths = []
    for curve in sides:
        for path in flatten_curve(curve):
            paths.append(
                [
                    [round((x - centre_x) * scale_x, 6), round((bottom - y) * scale_z, 6)]
                    for x, y in path
                ]
            )
    for line in horizontal:
        paths.append(
            [
                [round((line["x0"] - centre_x) * scale_x, 6), round((bottom - line["top"]) * scale_z, 6)],
                [round((line["x1"] - centre_x) * scale_x, 6), round((bottom - line["top"]) * scale_z, 6)],
            ]
        )
    return paths, {
        "published_envelope_mm": [520.0, 900.0],
        "pdf_exterior_bbox_points": [left, top, right, bottom],
        "horizontal_points_to_mm_scale": round(scale_x, 12),
        "vertical_points_to_mm_scale": round(scale_z, 12),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source = args.input.resolve()
    actual_hash = sha256(source)
    if actual_hash != EXPECTED["pdf"]:
        raise RuntimeError(f"official PDF hash mismatch: expected {EXPECTED['pdf']}, got {actual_hash}")
    variants = {}
    with pdfplumber.open(source) as document:
        if len(document.pages) != 2:
            raise RuntimeError(f"expected two official PDF pages, found {len(document.pages)}")
        for index, code in enumerate(("WFA", "WFB")):
            page = document.pages[index]
            plan, plan_normalisation = normalise_plan(select_plan(page))
            sides, horizontal = select_elevation(page)
            elevation, elevation_normalisation = normalise_elevation(sides, horizontal)
            variants[code] = {
                "model_code": code,
                "page": index + 1,
                "source_kind": "official_vector_pdf",
                "scope": SCOPE,
                "views": {
                    "plan": {"paths_mm": plan, "normalisation": plan_normalisation},
                    "elevation": {"paths_mm": elevation, "normalisation": elevation_normalisation},
                },
            }
    write_json(
        args.output.resolve(),
        {
            "schema_version": 2,
            "source_pdf": relative(source),
            "source_pdf_sha256": actual_hash,
            "source_kind": "official_vector_pdf",
            "units": "mm",
            "scope": SCOPE,
            "variants": variants,
        },
    )
    print(f"official PDF linework: {relative(args.output)}; sha256={actual_hash}")


if __name__ == "__main__":
    main()
