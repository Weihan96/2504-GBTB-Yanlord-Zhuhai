#!/usr/bin/env python3
"""Extract the exact Marilyn bergere three-view linework from Baxter native DWG."""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, relative, sha256, write_json
from geberit_146_140_linework import sample_arc, sample_spline


PRODUCT_DIR = ROOT / "output/review/highpoly-types/marilyn-01"
SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_DWG = SOURCE_DIR / "Marilyn_Abaco.dwg"
SOURCE_ZIP = SOURCE_DIR / "Baxter_Marilyn_Armchair_2D_3D.zip"
SOURCE_3DS = SOURCE_DIR / "Marilyn_bergere_86x100xh94.3ds"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
PRODUCT_PAGE = "https://www.baxter.it/en/products/marilyn-sofas-and-armchairs"
ZIP_URL = "https://dam.baxter.it/asset/9b02ac3a-7166-40df-9acc-1e8a2ad46705/Baxter_Marilyn_Armchair_2D_3D.zip"
EXPECTED = {
    "zip": "db054a193c9c5b8dc45f0df39199b2d96dc3712a5a7e48821fe94dce1566fedb",
    "dwg": "cb9825eecab7c78ee28e7668e933c2fe413395a63bf75f3afab8cc4959864724",
    "3ds": "171dbd72b2298f201d834fb54a5d39bca7ca5b86a552b97a5cbc64e23dd5e385",
}
REGIONS = {
    "plan": (9900.0, 300.0, 11100.0, 1400.0),
    "front": (9900.0, 1650.0, 11100.0, 2750.0),
    "side": (11100.0, 1650.0, 12400.0, 2750.0),
}
EXPECTED_PATH_COUNTS = {"plan": 24, "front": 68, "side": 54}
PRODUCT_LAYER = "_ARREDO"
SCOPE = "exact Baxter Marilyn bergere 86 x 100 x 94 cm family CAD reference; not a project shop drawing"


def handle_value(handle) -> int:
    return int(handle[-1])


def dwg_json(path: Path) -> tuple[dict, str]:
    with tempfile.NamedTemporaryFile(suffix=".json") as target:
        subprocess.run(["dwgread", "-O", "JSON", str(path)], check=True, stdout=target, stderr=subprocess.PIPE)
        target.seek(0)
        payload = json.load(target)
    version = subprocess.run(["dwgread", "--version"], check=True, capture_output=True, text=True).stdout.splitlines()[0]
    return payload, version


def entity_path(entity: dict, objects_by_handle: dict[int, dict]) -> list[list[float]] | None:
    kind = entity.get("entity")
    if kind == "SPLINE":
        return sample_spline(entity)
    if kind == "LINE":
        return [[float(value) for value in entity["start"][:2]], [float(value) for value in entity["end"][:2]]]
    if kind == "ARC":
        return sample_arc(entity)
    if kind == "LWPOLYLINE":
        points = [[float(value) for value in point[:2]] for point in entity["points"]]
        return points + ([points[0]] if int(entity.get("flag", 0)) & 1 else [])
    if kind == "POLYLINE_2D":
        points = [
            [float(value) for value in objects_by_handle[handle_value(handle)]["point"][:2]]
            for handle in entity["vertex"]
        ]
        return points + ([points[0]] if int(entity.get("flag", 0)) & 1 else [])
    return None


def bounds(paths: list[list[list[float]]]) -> dict:
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return {
        "minimum": [round(value, 6) for value in minimum],
        "maximum": [round(value, 6) for value in maximum],
        "size": [round(maximum[axis] - minimum[axis], 6) for axis in range(2)],
    }


def normalize(view: str, paths: list[list[list[float]]]) -> tuple[list[list[list[float]]], dict]:
    native = bounds(paths)
    minimum, maximum = native["minimum"], native["maximum"]
    center_x = (minimum[0] + maximum[0]) / 2.0
    if view == "plan":
        center_y = (minimum[1] + maximum[1]) / 2.0
        normalized = [[[round(point[0] - center_x, 6), round(point[1] - center_y, 6)] for point in path] for path in paths]
        origin = [center_x, center_y]
        origin_rule = "native visible-bounds centre"
    else:
        normalized = [[[round(point[0] - center_x, 6), round(point[1] - minimum[1], 6)] for point in path] for path in paths]
        origin = [center_x, minimum[1]]
        origin_rule = "native horizontal centre and floor baseline"
    return normalized, {
        "origin_native_mm": [round(value, 6) for value in origin],
        "origin_rule": origin_rule,
        "scale": 1.0,
        "native_bounds_mm": native,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE_DWG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source = args.source.resolve()
    if sha256(SOURCE_ZIP) != EXPECTED["zip"] or sha256(source) != EXPECTED["dwg"] or sha256(SOURCE_3DS) != EXPECTED["3ds"]:
        raise RuntimeError("Baxter Marilyn official native package hash mismatch")
    document, parser_version = dwg_json(source)
    if document.get("HEADER", {}).get("INSUNITS") != 4:
        raise RuntimeError("Marilyn native DWG is not millimetres")
    objects = document["OBJECTS"]
    objects_by_handle = {handle_value(item["handle"]): item for item in objects if "handle" in item}
    layers = {handle_value(item["handle"]): item["name"] for item in objects if item.get("object") == "LAYER"}
    extracted = []
    for entity in objects:
        if entity.get("entmode") != 2 or layers.get(handle_value(entity.get("layer", [0]))) != PRODUCT_LAYER:
            continue
        try:
            path = entity_path(entity, objects_by_handle)
        except (RuntimeError, KeyError):
            # The source contains rational Make2D splines outside the selected
            # exact bergere clusters. The audited product layer clusters use
            # only the supported non-rational curves below.
            continue
        if path:
            extracted.append((entity, path, bounds([path])))
    views = {}
    for view, region in REGIONS.items():
        selected = [
            (entity, path)
            for entity, path, path_bounds in extracted
            if path_bounds["minimum"][0] >= region[0]
            and path_bounds["minimum"][1] >= region[1]
            and path_bounds["maximum"][0] <= region[2]
            and path_bounds["maximum"][1] <= region[3]
        ]
        paths, normalisation = normalize(view, [path for _, path in selected])
        if len(paths) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Marilyn {view} native path count drifted: {len(paths)}")
        views[view] = {
            "view": view,
            "paths_mm": paths,
            "path_count": len(paths),
            "point_count": sum(len(path) for path in paths),
            "bounds_mm": bounds(paths),
            "native_entity_counts": dict(sorted(Counter(entity["entity"] for entity, _ in selected).items())),
            "native_entity_indices": [int(entity["index"]) for entity, _ in selected],
            "native_product_layer": PRODUCT_LAYER,
            "normalisation": normalisation,
        }
    dimension_labels = sorted(
        {
            item.get("text", "").replace("\\A1;", "")
            for item in objects
            if item.get("entity") == "MTEXT"
            and 9500.0 <= float(item.get("ins_pt", [0.0])[0]) <= 12500.0
            and 0.0 <= float(item.get("ins_pt", [0.0, 0.0])[1]) <= 2800.0
        }
    )
    for required in ("860", "1000", "940", "450"):
        if required not in dimension_labels:
            raise RuntimeError(f"Marilyn exact bergere DWG dimension {required} missing")
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/marilyn_01_linework.py",
        "manufacturer": "Baxter",
        "family": "Marilyn",
        "model": "bergere armchair with swivel base, 86 x 100 x 94 cm",
        "project_ifc_type_name": "Marilyn 01",
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
        "views": views,
        "pass": True,
    }
    write_json(args.output.resolve(), payload)
    print(json.dumps({"output": relative(args.output), "views": {view: item["bounds_mm"] for view, item in views.items()}, "pass": True}, indent=2))


if __name__ == "__main__":
    main()
