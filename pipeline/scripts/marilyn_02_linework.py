#!/usr/bin/env python3
"""Extract the exact Marilyn pouf Plan/Front/Side paths from Baxter native DWG."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, relative, sha256, write_json
from marilyn_01_linework import bounds, dwg_json, entity_path, handle_value, normalize


PRODUCT_DIR = ROOT / "output/review/highpoly-types/marilyn-02"
SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_DWG = SOURCE_DIR / "Marilyn_Abaco.dwg"
SOURCE_ZIP = SOURCE_DIR / "Baxter_Marilyn_Armchair_2D_3D.zip"
SOURCE_3DS = SOURCE_DIR / "Marilyn_pouf_80x62xh45.3ds"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
PRODUCT_PAGE = "https://www.baxter.it/en/products/marilyn-sofas-and-armchairs"
ZIP_URL = "https://dam.baxter.it/asset/9b02ac3a-7166-40df-9acc-1e8a2ad46705/Baxter_Marilyn_Armchair_2D_3D.zip"
EXPECTED = {
    "zip": "db054a193c9c5b8dc45f0df39199b2d96dc3712a5a7e48821fe94dce1566fedb",
    "dwg": "cb9825eecab7c78ee28e7668e933c2fe413395a63bf75f3afab8cc4959864724",
    "3ds": "53819a0b1321f0177cc162a4b93d62031b4c9aab77d06370e1de90c961cb97f9",
}
REGIONS = {
    "plan": (17900.0, 600.0, 19000.0, 1400.0),
    "front": (17900.0, 1700.0, 19050.0, 2250.0),
    "side": (19100.0, 1700.0, 20000.0, 2250.0),
}
EXPECTED_PATH_COUNTS = {"plan": 10, "front": 52, "side": 56}
PRODUCT_LAYERS = {"Make2D$Visibile$Curve", "0"}
PROJECT_BODY_SIZE_XYZ_MM = [810.61496, 591.948944, 456.581987]
OFFICIAL_NOMINAL_SIZE_XYZ_MM = [800.0, 620.0, 450.0]
DIMENSION_TOLERANCE_MM = 35.0
SCOPE = "exact Baxter Marilyn pouf with swivel base 80 x 62 x 45 cm family CAD reference; not a project shop drawing"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE_DWG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source = args.source.resolve()
    if sha256(SOURCE_ZIP) != EXPECTED["zip"] or sha256(source) != EXPECTED["dwg"] or sha256(SOURCE_3DS) != EXPECTED["3ds"]:
        raise RuntimeError("Baxter Marilyn pouf official native package hash mismatch")
    document, parser_version = dwg_json(source)
    if document.get("HEADER", {}).get("INSUNITS") != 4:
        raise RuntimeError("Marilyn native DWG is not millimetres")
    objects = document["OBJECTS"]
    objects_by_handle = {handle_value(item["handle"]): item for item in objects if "handle" in item}
    layers = {handle_value(item["handle"]): item["name"] for item in objects if item.get("object") == "LAYER"}
    extracted = []
    for entity in objects:
        if entity.get("entmode") != 2 or layers.get(handle_value(entity.get("layer", [0]))) not in PRODUCT_LAYERS:
            continue
        try:
            path = entity_path(entity, objects_by_handle)
        except (RuntimeError, KeyError):
            continue
        if path:
            extracted.append((entity, path, bounds([path]), layers.get(handle_value(entity.get("layer", [0])))))
    views = {}
    for view, region in REGIONS.items():
        selected = [
            (entity, path, layer)
            for entity, path, path_bounds, layer in extracted
            if path_bounds["minimum"][0] >= region[0]
            and path_bounds["minimum"][1] >= region[1]
            and path_bounds["maximum"][0] <= region[2]
            and path_bounds["maximum"][1] <= region[3]
        ]
        paths, normalisation = normalize(view, [path for _, path, _ in selected])
        if len(paths) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Marilyn pouf {view} native path count drifted: {len(paths)}")
        views[view] = {
            "view": view,
            "paths_mm": paths,
            "path_count": len(paths),
            "point_count": sum(len(path) for path in paths),
            "bounds_mm": bounds(paths),
            "native_entity_counts": dict(sorted(Counter(entity["entity"] for entity, _, _ in selected).items())),
            "native_layer_counts": dict(sorted(Counter(layer for _, _, layer in selected).items())),
            "native_entity_indices": [int(entity["index"]) for entity, _, _ in selected],
            "normalisation": normalisation,
        }
    dimension_labels = sorted({
        item.get("text", "").replace("\\A1;", "")
        for item in objects
        if item.get("entity") == "MTEXT"
        and 17500.0 <= float(item.get("ins_pt", [0.0])[0]) <= 20500.0
        and 0.0 <= float(item.get("ins_pt", [0.0, 0.0])[1]) <= 2500.0
    })
    for required in ("800", "620", "450"):
        if required not in dimension_labels:
            raise RuntimeError(f"Marilyn exact pouf DWG dimension {required} missing")
    view_sizes = {view: item["bounds_mm"]["size"] for view, item in views.items()}
    project_sizes = {
        "plan": PROJECT_BODY_SIZE_XYZ_MM[:2],
        "front": [PROJECT_BODY_SIZE_XYZ_MM[0], PROJECT_BODY_SIZE_XYZ_MM[2]],
        "side": [PROJECT_BODY_SIZE_XYZ_MM[1], PROJECT_BODY_SIZE_XYZ_MM[2]],
    }
    maximum_dwg_body_delta = max(
        abs(view_sizes[view][axis] - project_sizes[view][axis])
        for view in views for axis in range(2)
    )
    maximum_nominal_body_delta = max(abs(PROJECT_BODY_SIZE_XYZ_MM[axis] - OFFICIAL_NOMINAL_SIZE_XYZ_MM[axis]) for axis in range(3))
    pass_dimensions = max(maximum_dwg_body_delta, maximum_nominal_body_delta) <= DIMENSION_TOLERANCE_MM
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/marilyn_02_linework.py",
        "manufacturer": "Baxter",
        "family": "Marilyn",
        "model": "pouf with swivel base, 80 x 62 x 45 cm",
        "project_ifc_type_name": "Marilyn 02",
        "project_ifc_type_description": "Pouf with swivel base W80D62H45",
        "product_page": PRODUCT_PAGE,
        "source_kind": "native_dwg",
        "scope": SCOPE,
        "source_zip": relative(SOURCE_ZIP),
        "source_zip_url": ZIP_URL,
        "source_zip_sha256": sha256(SOURCE_ZIP),
        "source_dwg": relative(source),
        "source_dwg_sha256": sha256(source),
        "exact_model_3ds": relative(SOURCE_3DS),
        "exact_model_3ds_sha256": sha256(SOURCE_3DS),
        "dwg_file_version": document["FILEHEADER"]["version"],
        "parser": parser_version,
        "units": "mm",
        "native_dimension_labels_in_exact_cluster": dimension_labels,
        "dimension_cross_check": {
            "official_nominal_width_depth_height_mm": OFFICIAL_NOMINAL_SIZE_XYZ_MM,
            "official_native_dwg_sizes_by_view_mm": view_sizes,
            "project_ifc_body_local_xyz_mm": PROJECT_BODY_SIZE_XYZ_MM,
            "maximum_native_dwg_to_body_delta_mm": round(maximum_dwg_body_delta, 6),
            "maximum_nominal_to_body_delta_mm": round(maximum_nominal_body_delta, 6),
            "tolerance_mm": DIMENSION_TOLERANCE_MM,
            "pass": pass_dimensions,
        },
        "views": views,
        "pass": pass_dimensions,
    }
    write_json(args.output.resolve(), payload)
    print(json.dumps({"output": relative(args.output), "path_counts": EXPECTED_PATH_COUNTS, "view_sizes_mm": view_sizes, "pass": payload["pass"]}, indent=2))


if __name__ == "__main__":
    main()
