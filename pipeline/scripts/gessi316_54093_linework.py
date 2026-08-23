#!/usr/bin/env python3
"""Extract exact Plan/Front/Side linework from the official Gessi 54093 G000 DWG."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from geberit_146_140_linework import dwg_json, handle_value, path_bounds
from gessi316_54294_linework import bounds, gessi_entity_path


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54093"
SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_DWG = SOURCE_DIR / "GPF5409300000G000_3.dwg"
SOURCE_ZIP = SOURCE_DIR / "GPF5409300000G000_arc.zip"
SOURCE_PDF = SOURCE_DIR / "GPF5409300000G000_1.pdf"
REVALIDATION = SOURCE_DIR / "official-source-revalidation.json"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
SOURCE_DWG_URL = "https://gessistorage.blob.core.windows.net/zwa/GPF5409300000G000_arc.zip"
SOURCE_PDF_URL = "https://gessistorage.blob.core.windows.net/zc4/GPF5409300000G000_1.pdf"
EXPECTED = {
    "dwg": "a9308fc8498d34c8cf2f68fd28aaf90b11e62bfc59424e0a5a3ac7f0d28d5047",
    "zip": "615d4747639b6acfa1daca28a894124fef26e74239de89ef68455b92db31d0f9",
    "pdf": "2eb67af1be75093a12ab8bc3a739e5bdd6e4b1c4c544948795c87cb9598515c6",
}
SOURCE_LAYER_HANDLE = 509
SOURCE_LAYER_NAME = "00_COMPONENTI"
ENTITY_TYPES = {"LINE", "ARC", "ELLIPSE", "SPLINE", "CIRCLE", "LWPOLYLINE"}
VIEW_WINDOWS_MM = {
    "plan": (10.0, -190.0, 90.0, 120.0),
    "front": (10.0, 120.0, 90.0, 420.0),
    "side": (80.0, 120.0, 392.0, 420.0),
}
EXPECTED_PATH_COUNTS = {"plan": 7, "front": 914, "side": 21}
SCOPE = "official Gessi exact 54093 G000 family reference; not a project shop drawing"


def in_window(value, window) -> bool:
    return value[0] >= window[0] and value[1] >= window[1] and value[2] <= window[2] and value[3] <= window[3]


def normalize(view: str, paths: list[list[list[float]]]):
    source_bounds = path_bounds(paths)
    minimum_x, minimum_y = source_bounds["minimum"]
    maximum_x, maximum_y = source_bounds["maximum"]
    center_x = (minimum_x + maximum_x) / 2.0
    if view == "plan":
        normalized = [[[round(x - center_x, 6), round(y - minimum_y, 6)] for x, y in path] for path in paths]
        origin = "assembly_center_x_and_minimum_depth"
    elif view == "front":
        normalized = [[[round(x - center_x, 6), round(y - maximum_y, 6)] for x, y in path] for path in paths]
        origin = "assembly_center_x_and_top_z"
    else:
        normalized = [[[round(x - minimum_x, 6), round(y - maximum_y, 6)] for x, y in path] for path in paths]
        origin = "minimum_depth_and_top_z"
    return normalized, {"source_sheet_bounds_mm": source_bounds, "normalization_scale": 1.0, "normalization_origin": origin}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if sha256(SOURCE_DWG) != EXPECTED["dwg"] or sha256(SOURCE_ZIP) != EXPECTED["zip"] or sha256(SOURCE_PDF) != EXPECTED["pdf"]:
        raise RuntimeError("Gessi 54093 official source hash gate failed")
    revalidation = load_json(REVALIDATION)
    if (
        revalidation.get("pass") is not True
        or revalidation.get("public_access", {}).get("exact_54093_g000_native_dwg_publicly_downloadable") is not True
        or revalidation.get("drawing_geometry_policy", {}).get("g001_or_a004_variant_cad_used") is not False
        or revalidation.get("drawing_geometry_policy", {}).get("adjacent_gessi_product_cad_used") is not False
        or revalidation.get("drawing_geometry_policy", {}).get("third_party_cad_used") is not False
    ):
        raise RuntimeError("Gessi 54093 source revalidation gate failed")
    data = dwg_json(SOURCE_DWG)
    selected = {view: [] for view in VIEW_WINDOWS_MM}
    source_type_counts = {view: {} for view in VIEW_WINDOWS_MM}
    for entity in data["OBJECTS"]:
        kind = entity.get("entity")
        if kind not in ENTITY_TYPES or handle_value(entity.get("layer")) != SOURCE_LAYER_HANDLE:
            continue
        path = gessi_entity_path(entity)
        value = bounds(path)
        for view, window in VIEW_WINDOWS_MM.items():
            if in_window(value, window):
                selected[view].append((handle_value(entity["handle"]), path))
                source_type_counts[view][kind] = source_type_counts[view].get(kind, 0) + 1
                break
    views = {}
    for view, items in selected.items():
        items.sort(key=lambda item: item[0])
        paths, metadata = normalize(view, [path for _, path in items])
        if len(paths) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi 54093 {view} path-count drift: {len(paths)}")
        views[view] = {
            "paths_mm": paths, "path_count": len(paths), "bounds_mm": path_bounds(paths),
            "source_sheet_bounds_mm": metadata["source_sheet_bounds_mm"], "source_window_mm": list(VIEW_WINDOWS_MM[view]),
            "source_layer": SOURCE_LAYER_NAME, "source_layer_handle": SOURCE_LAYER_HANDLE,
            "source_entity_type_counts": source_type_counts[view], "source_kind": "native_dwg",
            "source_dwg_sha256": EXPECTED["dwg"], "normalization_scale": 1.0,
            "normalization_origin": metadata["normalization_origin"],
        }
    plan_size, front_size, side_size = views["plan"]["bounds_mm"]["size"], views["front"]["bounds_mm"]["size"], views["side"]["bounds_mm"]["size"]
    delta = [abs(plan_size[0] - 50.0), abs(plan_size[1] - 190.0), abs(front_size[1] - 273.0), abs(side_size[0] - 190.0)]
    payload = {
        "schema_version": 1, "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/gessi316_54093_linework.py", "manufacturer": "Gessi",
        "family": "Gessi316 Meccanica", "article_number": "54093", "configuration": "G000",
        "ifc_type_name": "Gessi316 54093", "source_kind": "native_dwg",
        "source_label_zh": "Gessi 官方精确型号 54093 G000 原生 DWG 图纸表达",
        "source_label_en": "drawing representation from the exact Gessi 54093 G000 official native DWG",
        "scope": SCOPE,
        "sheet": {"units": "mm", "source_layer": SOURCE_LAYER_NAME},
        "official_sources": {
            "native_dwg_zip": {"path": relative(SOURCE_ZIP), "url": SOURCE_DWG_URL, "sha256": EXPECTED["zip"]},
            "native_dwg": {"path": relative(SOURCE_DWG), "url": SOURCE_DWG_URL, "sha256": EXPECTED["dwg"]},
            "technical_vector_pdf": {"path": relative(SOURCE_PDF), "url": SOURCE_PDF_URL, "sha256": EXPECTED["pdf"]},
            "revalidation": {"path": relative(REVALIDATION), "sha256": sha256(REVALIDATION)},
        },
        "identity_gates": {"54093_g000_native_dwg_used": True, "g001_or_a004_variant_cad_used": False, "adjacent_gessi_product_cad_used": False, "third_party_cad_used": False},
        "nominal_dimension_cross_check": {
            "official_api_width_depth_height_mm": [50.0, 190.0, 273.0],
            "native_dwg_plan_envelope_mm": plan_size, "native_dwg_front_envelope_mm": front_size,
            "native_dwg_side_envelope_mm": side_size, "absolute_delta_mm": [round(value, 6) for value in delta],
            "tolerance_mm": 0.5, "pass": max(delta) <= 0.5,
        },
        "views": views, "pass": max(delta) <= 0.5,
    }
    write_json(args.output.resolve(), payload)
    if not payload["pass"]:
        raise RuntimeError("Gessi 54093 DWG/PDF/API dimension cross-check failed")
    print(relative(args.output))


if __name__ == "__main__":
    main()
