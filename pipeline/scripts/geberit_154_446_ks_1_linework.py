#!/usr/bin/env python3
"""Extract exact Geberit CleanLine50 154.446.KS.1 views from official DWGs."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT
from geberit_146_140_linework import (
    CONTOUR_LAYER,
    dwg_json,
    entity_path,
    handle_value,
    path_bounds,
    relative,
    sha256,
    write_json,
)


ARTICLE = "154.446.KS.1"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-154-446-ks-1"
SOURCE_DIR = PRODUCT_DIR / "official-source"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
PRODUCT_PAGE = "https://catalog.geberit-global.com/en-XB/product/PRO_4702328"
ASSORTMENT_OVERVIEW = "https://cdn.data.geberit.com/overviews/INT-en/DAS_331713.pdf"
SCOPE = "official manufacturer exact archived article reference, not a project shop drawing"
EXPECTED = {
    "A": "de09a6bbf99ecbe3c832cbee7ac9d398e95b927246f274643b513bfe960586c3",
    "G": "bc41f8db7c7989de9d9089374e03c3b9deea94e2634bbb0f9418cad982b4b6e4",
    "L": "dedb47964cc69310bd9af3982379c428c8b20494b7588296aedc5d07380e3421",
    "P": "46e1d9fe9b7afda6f4a72970a5e35d7bc0935c2bdc890ccca232f0b0896a8e82",
}
VIEW_CODES = {"plan": "G", "front": "A", "side": "L"}
EXPECTED_PATH_COUNTS = {"plan": 122, "front": 118, "side": 176}


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
    header_minimum = [float(value) for value in payload["HEADER"]["EXTMIN"][:2]]
    header_maximum = [float(value) for value in payload["HEADER"]["EXTMAX"][:2]]
    header_size = [header_maximum[index] - header_minimum[index] for index in range(2)]
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
        "contour_bounds_mm": path_bounds(paths),
        "native_header_extents_mm": {
            "minimum": [round(value, 6) for value in header_minimum],
            "maximum": [round(value, 6) for value in header_maximum],
            "size": [round(value, 6) for value in header_size],
        },
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
            raise RuntimeError(f"official Geberit CleanLine50 {code} DWG hash mismatch")
    views = {
        view: extract_view(sources[code], view, code)
        for view, code in VIEW_CODES.items()
    }
    for view, expected in EXPECTED_PATH_COUNTS.items():
        if views[view]["path_count"] != expected:
            raise RuntimeError(f"official CleanLine50 {view} contour path count drifted")
    expected_contour_sizes = {
        "plan": [900.0, 53.4],
        "front": [900.0, 50.5],
        "side": [53.4, 50.5],
    }
    for view, expected in expected_contour_sizes.items():
        actual = views[view]["contour_bounds_mm"]["size"]
        if any(abs(actual[index] - expected[index]) > 0.01 for index in range(2)):
            raise RuntimeError(f"official CleanLine50 {view} contour bounds drifted: {actual}")
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/geberit_154_446_ks_1_linework.py",
        "manufacturer": "Geberit",
        "family": "CleanLine50 shower channel L90cm",
        "article_number": ARTICLE,
        "replacement_article_number": "154.446.KS.2",
        "ifc_type_name": "Geberit 154.446.KS.1",
        "product_page": PRODUCT_PAGE,
        "assortment_overview": ASSORTMENT_OVERVIEW,
        "source_kind": "native_dwg",
        "scope": SCOPE,
        "three_view_code_mapping": {
            "G": "Grundriss / plan",
            "A": "Ansicht / front elevation",
            "L": "left side elevation",
            "P": "official 3D model; archived for identity only, never substituted for 2D views",
        },
        "retirement_identity": {
            "exact_archived_article": ARTICLE,
            "official_replacement_article": "154.446.KS.2",
            "replacement_cad_used": False,
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
    print(json.dumps({
        "output": relative(args.output.resolve()),
        "article_number": ARTICLE,
        "replacement_cad_used": False,
        "views": {view: {"paths": data["path_count"], "bounds": data["contour_bounds_mm"]} for view, data in views.items()},
        "pass": True,
    }, indent=2))


if __name__ == "__main__":
    main()
