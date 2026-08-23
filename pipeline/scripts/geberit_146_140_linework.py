#!/usr/bin/env python3
"""Extract exact Geberit 146.140 three-view linework from official native DWG files."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-146-140"
SOURCE_DIR = PRODUCT_DIR / "official-source"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
ARTICLE = "146.140.11.1"
PRODUCT_PAGE = "https://catalog.international.geberit.com/en-GB/product/PRO_942494"
SCOPE = "official manufacturer product-family reference, not a project shop drawing"
EXPECTED = {
    "A": "8e7f2ec11311e77e8625c908d13586379fb8def19e6f65320ce61b15b2515172",
    "G": "f6253efd2d736ef8fa0830455692f96c0f61b057a597f38ddc0544e3c9de5c0c",
    "L": "e3c29b0f1778de7f3178a6dd6a2261dc3e7d46942e9c81e7a9943308204e9233",
    "P": "229eea526831ff08de65125b3c10cbf3ea9380310a5fd674b74a77634cde7669",
}
VIEW_CODES = {"plan": "G", "front": "A", "side": "L"}
CONTOUR_LAYER = "GEB-KONTUR-A"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def handle_value(handle) -> int:
    return int(handle[-1])


def dwg_json(path: Path) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".json") as target:
        subprocess.run(
            ["dwgread", "-O", "JSON", str(path)],
            check=True,
            stdout=target,
            stderr=subprocess.PIPE,
        )
        target.seek(0)
        return json.load(target)


def de_boor_point(degree: int, knots: list[float], controls: list[tuple[float, float]], value: float):
    count = len(controls)
    if value >= knots[count]:
        span = count - 1
    else:
        span = degree
        while span + 1 < len(knots) and not (knots[span] <= value < knots[span + 1]):
            span += 1
    work = [list(controls[span - degree + index]) for index in range(degree + 1)]
    for level in range(1, degree + 1):
        for index in range(degree, level - 1, -1):
            source = span - degree + index
            denominator = knots[source + degree - level + 1] - knots[source]
            alpha = 0.0 if abs(denominator) < 1e-12 else (value - knots[source]) / denominator
            work[index][0] = (1.0 - alpha) * work[index - 1][0] + alpha * work[index][0]
            work[index][1] = (1.0 - alpha) * work[index - 1][1] + alpha * work[index][1]
    return work[degree][0], work[degree][1]


def sample_spline(entity: dict) -> list[list[float]]:
    degree = int(entity["degree"])
    controls = [(float(point["x"]), float(point["y"])) for point in entity["ctrl_pts"]]
    knots = [float(value) for value in entity["knots"]]
    if entity.get("rational") or len(knots) != len(controls) + degree + 1:
        raise RuntimeError("unsupported Geberit DWG spline identity")
    start, end = knots[degree], knots[len(controls)]
    steps = max(16, min(160, len(controls) * 8))
    return [
        [round(value[0], 6), round(value[1], 6)]
        for value in (
            de_boor_point(degree, knots, controls, start + (end - start) * index / steps)
            for index in range(steps + 1)
        )
    ]


def angular_span(start: float, end: float) -> float:
    span = end - start
    while span <= 0.0:
        span += math.tau
    return span


def sample_arc(entity: dict) -> list[list[float]]:
    center_x, center_y = map(float, entity["center"][:2])
    radius = float(entity["radius"])
    start = float(entity["start_angle"])
    span = angular_span(start, float(entity["end_angle"]))
    steps = max(8, math.ceil(span / (math.pi / 36.0)))
    return [
        [
            round(center_x + radius * math.cos(start + span * index / steps), 6),
            round(center_y + radius * math.sin(start + span * index / steps), 6),
        ]
        for index in range(steps + 1)
    ]


def sample_ellipse(entity: dict) -> list[list[float]]:
    center_x, center_y = map(float, entity["center"][:2])
    major_x, major_y = map(float, entity["sm_axis"][:2])
    ratio = float(entity["axis_ratio"])
    minor_x, minor_y = -major_y * ratio, major_x * ratio
    start = float(entity["start_angle"])
    span = angular_span(start, float(entity["end_angle"]))
    steps = max(8, math.ceil(span / (math.pi / 36.0)))
    return [
        [
            round(center_x + major_x * math.cos(angle) + minor_x * math.sin(angle), 6),
            round(center_y + major_y * math.cos(angle) + minor_y * math.sin(angle), 6),
        ]
        for angle in (start + span * index / steps for index in range(steps + 1))
    ]


def entity_path(entity: dict) -> list[list[float]]:
    kind = entity["entity"]
    if kind == "LINE":
        return [
            [round(float(entity["start"][0]), 6), round(float(entity["start"][1]), 6)],
            [round(float(entity["end"][0]), 6), round(float(entity["end"][1]), 6)],
        ]
    if kind == "SPLINE":
        return sample_spline(entity)
    if kind == "ARC":
        return sample_arc(entity)
    if kind == "ELLIPSE":
        return sample_ellipse(entity)
    raise RuntimeError(f"unsupported Geberit contour entity: {kind}")


def path_bounds(paths: list[list[list[float]]]) -> dict:
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return {
        "minimum": [round(value, 6) for value in minimum],
        "maximum": [round(value, 6) for value in maximum],
        "size": [round(maximum[axis] - minimum[axis], 6) for axis in range(2)],
    }


def extract_view(path: Path, view: str, code: str) -> dict:
    payload = dwg_json(path)
    if payload.get("HEADER", {}).get("INSUNITS") != 4:
        raise RuntimeError(f"{path.name}: native DWG is not millimetres")
    layers = {
        handle_value(item["handle"]): item["name"]
        for item in payload["OBJECTS"]
        if item.get("object") == "LAYER"
    }
    entities = [
        item
        for item in payload["OBJECTS"]
        if item.get("entity") in {"LINE", "SPLINE", "ARC", "ELLIPSE"}
        and layers.get(handle_value(item["layer"])) == CONTOUR_LAYER
    ]
    paths = [entity_path(entity) for entity in entities]
    counts: dict[str, int] = {}
    for entity in entities:
        counts[entity["entity"]] = counts.get(entity["entity"], 0) + 1
    return {
        "view": view,
        "native_dwg_code": code,
        "source_dwg": relative(path),
        "source_dwg_url": f"https://cdn.data.geberit.com/cad/{ARTICLE}_{code}.dwg",
        "source_dwg_sha256": sha256(path),
        "source_kind": "native_dwg",
        "units": "mm",
        "native_contour_layer": CONTOUR_LAYER,
        "native_entity_counts": counts,
        "path_count": len(paths),
        "point_count": sum(len(path) for path in paths),
        "bounds_mm": path_bounds(paths),
        "paths_mm": paths,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=SOURCE_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source_dir = args.source_dir.resolve()
    sources = {code: source_dir / f"{ARTICLE}_{code}.dwg" for code in EXPECTED}
    for code, path in sources.items():
        if not path.is_file() or sha256(path) != EXPECTED[code]:
            raise RuntimeError(f"official Geberit {code} DWG hash mismatch")
    views = {
        view: extract_view(sources[code], view, code)
        for view, code in VIEW_CODES.items()
    }
    expected_sizes = {
        "plan": (385.714, 577.769),
        "front": (385.714, 444.235),
        "side": (577.769, 444.235),
    }
    for view, (expected_width, expected_height) in expected_sizes.items():
        width, height = views[view]["bounds_mm"]["size"]
        if abs(width - expected_width) > 3.0 or abs(height - expected_height) > 3.0:
            raise RuntimeError(f"official Geberit {view} view bounds drifted: {width} x {height}")
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/geberit_146_140_linework.py",
        "manufacturer": "Geberit",
        "family": "AquaClean Sela wall-hung WC",
        "article_number": ARTICLE,
        "ifc_type_name": "Geberit 146.140",
        "product_page": PRODUCT_PAGE,
        "source_kind": "native_dwg",
        "scope": SCOPE,
        "three_view_code_mapping": {
            "G": "Grundriss / plan",
            "A": "Ansicht / front elevation",
            "L": "left side elevation",
            "P": "official 3D model; archived for identity only, never substituted for 2D views",
        },
        "official_sources": {
            code: {
                "path": relative(path),
                "url": f"https://cdn.data.geberit.com/cad/{ARTICLE}_{code}.dwg",
                "sha256": sha256(path),
            }
            for code, path in sources.items()
        },
        "views": views,
        "pass": True,
    }
    write_json(args.output.resolve(), payload)
    print(json.dumps({"output": relative(args.output), "views": {k: v["bounds_mm"] for k, v in views.items()}, "pass": True}, indent=2))


if __name__ == "__main__":
    main()
