#!/usr/bin/env python3
"""Extract Plan/Front/Side linework from the exact official Gessi 54294 DWG."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, relative, sha256, write_json
from geberit_146_140_linework import dwg_json, entity_path, handle_value, path_bounds


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_DWG = SOURCE_DIR / "GPF5429400000G000_3.dwg"
SOURCE_ZIP = SOURCE_DIR / "GPF5429400000G000_arc.zip"
SOURCE_PDF = SOURCE_DIR / "GPF5429400000G000_1.pdf"
REVALIDATION = SOURCE_DIR / "official-source-revalidation.json"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
SOURCE_DWG_URL = "https://gessistorage.blob.core.windows.net/zwa/GPF5429400000G000_arc.zip"
SOURCE_PDF_URL = "https://gessistorage.blob.core.windows.net/zc4/GPF5429400000G000_1.pdf"
EXPECTED = {
    "dwg": "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4",
    "zip": "fad0d98e94a83bb488b5f0703472862d628759f21700c1a87c3af343f4c1bdb9",
    "pdf": "82058bb27752f27e6fc92029783d67c15fd443befcd393f7d572fd0147a581e3",
}
SOURCE_LAYER = "00_COMPONENTI"
ENTITY_TYPES = {"LINE", "ARC", "ELLIPSE", "SPLINE", "CIRCLE", "LWPOLYLINE"}
EXPECTED_PATH_COUNTS = {"plan": 1481, "front": 1099, "side": 745}
SCOPE = "official Gessi exact 54294 family reference, not a project shop drawing"


def circle_path(entity: dict) -> list[list[float]]:
    center_x, center_y = map(float, entity["center"][:2])
    radius = float(entity["radius"])
    return [
        [
            round(center_x + radius * math.cos(math.tau * index / 72), 6),
            round(center_y + radius * math.sin(math.tau * index / 72), 6),
        ]
        for index in range(73)
    ]


def lwpolyline_path(entity: dict) -> list[list[float]]:
    if entity.get("bulges"):
        raise RuntimeError("unsupported bulged Gessi LWPOLYLINE")
    return [[round(float(point[0]), 6), round(float(point[1]), 6)] for point in entity["points"]]


def gessi_entity_path(entity: dict) -> list[list[float]]:
    if entity["entity"] == "CIRCLE":
        return circle_path(entity)
    if entity["entity"] == "LWPOLYLINE":
        return lwpolyline_path(entity)
    return entity_path(entity)


def bounds(path: list[list[float]]) -> tuple[float, float, float, float]:
    return (
        min(point[0] for point in path),
        min(point[1] for point in path),
        max(point[0] for point in path),
        max(point[1] for point in path),
    )


def is_view_path(view: str, value: tuple[float, float, float, float]) -> bool:
    minimum_x, minimum_y, maximum_x, maximum_y = value
    if view == "front":
        return minimum_x >= 1.0 and maximum_x <= 364.6 and minimum_y >= 315.0 and maximum_y < 419.0
    if view == "side":
        return minimum_x >= 364.0 and maximum_x < 593.5 and minimum_y >= 315.0 and maximum_y < 419.0
    if view == "plan":
        in_plan_window = (
            minimum_x >= 1.0
            and maximum_x <= 364.6
            and minimum_y >= 70.0
            and maximum_y < 315.0
        )
        # The clean native DWG also contains isometric details below/right of the plan.
        # Below y=200 only the central 54294 spout belongs to the orthographic plan.
        low_region_is_central_spout = maximum_y >= 200.0 or (
            minimum_x >= 150.0 and maximum_x <= 215.0
        )
        return in_plan_window and low_region_is_central_spout
    raise RuntimeError(f"unknown Gessi view: {view}")


def normalize(view: str, paths: list[list[list[float]]]) -> tuple[list[list[list[float]]], dict]:
    source_bounds = path_bounds(paths)
    minimum_x, minimum_y = source_bounds["minimum"]
    maximum_x, maximum_y = source_bounds["maximum"]
    center_x = (minimum_x + maximum_x) / 2.0
    normalized = []
    for path in paths:
        if view == "plan":
            normalized.append([
                [round(x - center_x, 6), round(y - maximum_y, 6)]
                for x, y in path
            ])
        elif view == "front":
            normalized.append([
                [round(x - center_x, 6), round(y - minimum_y, 6)]
                for x, y in path
            ])
        else:
            normalized.append([
                [round(minimum_x - x, 6), round(y - minimum_y, 6)]
                for x, y in path
            ])
    return normalized, {
        "source_sheet_bounds_mm": source_bounds,
        "normalization_scale": 1.0,
        "normalization_origin": (
            "assembly_center_x_and_wall_y"
            if view == "plan"
            else "assembly_center_x_and_visible_bottom_z"
            if view == "front"
            else "wall_plane_and_visible_bottom_z"
        ),
        "axis_policy": (
            "x preserved; sheet y becomes local y with wall at zero"
            if view == "plan"
            else "x preserved; sheet y becomes local z"
            if view == "front"
            else "sheet x becomes negative local y from wall; sheet y becomes local z"
        ),
    }


def extract_view(payload: dict, view: str, layers: dict[int, str]) -> dict:
    selected = []
    counts: dict[str, int] = {}
    handles = []
    for entity in payload["OBJECTS"]:
        if entity.get("entity") not in ENTITY_TYPES:
            continue
        if layers.get(handle_value(entity["layer"])) != SOURCE_LAYER:
            continue
        path = gessi_entity_path(entity)
        if not is_view_path(view, bounds(path)):
            continue
        selected.append(path)
        counts[entity["entity"]] = counts.get(entity["entity"], 0) + 1
        handles.append(handle_value(entity["handle"]))
    if len(selected) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"Gessi {view} path-count drift: {len(selected)}")
    paths, normalization = normalize(view, selected)
    return {
        "view": view,
        "source_kind": "native_dwg",
        "source_dwg": relative(SOURCE_DWG),
        "source_dwg_sha256": EXPECTED["dwg"],
        "source_dwg_url": SOURCE_DWG_URL,
        "source_layer": SOURCE_LAYER,
        "units": "mm",
        "native_entity_counts": counts,
        "native_entity_handles": handles,
        "path_count": len(paths),
        "point_count": sum(len(path) for path in paths),
        "bounds_mm": path_bounds(paths),
        "paths_mm": paths,
        "normalization": normalization,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    for key, path in (("dwg", SOURCE_DWG), ("zip", SOURCE_ZIP), ("pdf", SOURCE_PDF)):
        if not path.is_file() or sha256(path) != EXPECTED[key]:
            raise RuntimeError(f"Gessi exact official {key} hash mismatch")
    revalidation = json.loads(REVALIDATION.read_text(encoding="utf-8"))
    if (
        revalidation.get("pass") is not True
        or revalidation.get("public_access", {}).get("exact_54294_native_dwg_publicly_downloadable") is not True
        or revalidation.get("drawing_geometry_policy", {}).get("companion_45089_dwg_used_as_54294_geometry") is not False
        or revalidation.get("drawing_geometry_policy", {}).get("adjacent_gessi_product_cad_used") is not False
    ):
        raise RuntimeError("Gessi official source revalidation gate failed")
    payload = dwg_json(SOURCE_DWG)
    if payload.get("HEADER", {}).get("INSUNITS") != 4:
        raise RuntimeError("Gessi native DWG is not declared in millimetres")
    layers = {
        handle_value(item["handle"]): item["name"]
        for item in payload["OBJECTS"]
        if item.get("object") == "LAYER"
    }
    views = {view: extract_view(payload, view, layers) for view in ("plan", "front", "side")}
    plan_size = views["plan"]["bounds_mm"]["size"]
    nominal_plan_size = [362.0, 210.0]
    plan_delta = [abs(plan_size[index] - nominal_plan_size[index]) for index in range(2)]
    if max(plan_delta) > 0.5:
        raise RuntimeError(f"Gessi native DWG plan envelope drifted: {plan_size}")
    output = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/gessi316_54294_linework.py",
        "manufacturer": "Gessi",
        "family": "Gessi316 Meccanica",
        "article_number": "45089_54294",
        "drawing_product_code": "54294",
        "companion_built_in_product_code": "45089",
        "ifc_type_name": "Gessi316 54294",
        "source_kind": "native_dwg",
        "source_label_zh": "Gessi 官方精确型号 54294 原生 DWG 图纸表达",
        "source_label_en": "drawing representation from the exact Gessi 54294 official native DWG",
        "scope": SCOPE,
        "official_sources": {
            "native_dwg_zip": {"path": relative(SOURCE_ZIP), "url": SOURCE_DWG_URL, "sha256": EXPECTED["zip"]},
            "native_dwg": {"path": relative(SOURCE_DWG), "url": SOURCE_DWG_URL, "sha256": EXPECTED["dwg"]},
            "technical_vector_pdf": {"path": relative(SOURCE_PDF), "url": SOURCE_PDF_URL, "sha256": EXPECTED["pdf"]},
            "revalidation": {"path": relative(REVALIDATION), "sha256": sha256(REVALIDATION)},
        },
        "sheet": {
            "native_units": "mm",
            "extents_mm": [594.0, 420.0],
            "source_layer": SOURCE_LAYER,
            "orthographic_view_selection": "mechanical coordinate windows with isometric exclusion gate",
        },
        "nominal_dimension_cross_check": {
            "official_api_width_mm": 362.0,
            "official_api_depth_mm": 210.0,
            "official_api_height_mm": 68.0,
            "native_dwg_plan_envelope_mm": plan_size,
            "absolute_delta_mm": [round(value, 6) for value in plan_delta],
            "tolerance_mm": 0.5,
            "pass": max(plan_delta) <= 0.5,
        },
        "views": views,
        "identity_gates": {
            "54294_native_dwg_used": True,
            "45089_companion_dwg_used_as_54294_geometry": False,
            "adjacent_gessi_product_cad_used": False,
            "third_party_cad_used": False,
        },
        "pass": True,
    }
    write_json(args.output.resolve(), output)
    print(json.dumps({"output": relative(args.output), "path_counts": EXPECTED_PATH_COUNTS, "pass": True}, indent=2))


if __name__ == "__main__":
    main()
