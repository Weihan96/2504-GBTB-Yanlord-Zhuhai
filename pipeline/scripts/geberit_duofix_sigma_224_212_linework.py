#!/usr/bin/env python3
"""Extract exact Geberit 224.212.00.2 views from official native DWGs."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT
from geberit_146_140_linework import (
    CONTOUR_LAYER,
    dwg_json,
    handle_value,
    path_bounds,
    relative,
    sample_arc,
    sample_ellipse,
    sha256,
    write_json,
)


ARTICLE = "224.212.00.2"
IFC_TYPE_NAME = "Geberit Duofix Sigma"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-duofix-sigma-224-212"
SOURCE_DIR = PRODUCT_DIR / "official-source"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
OFFICIAL_CATALOGUE = "https://cdn-geberit-country-sg.prod.web.geberit.com/_assets/local-media/brochures/2025-nsea-sanitary-and-bathroom-catalogue-web.pdf"
LOCAL_CATALOGUE_EXTRACT = SOURCE_DIR / "geberit-224212-catalogue-page11.pdf"
LOCAL_BOX_LABEL = SOURCE_DIR / "project-received-224212-box-label.jpg"
SCOPE = "official manufacturer exact article family reference, not a project shop drawing"
EXPECTED = {
    "A": "f06952ac7bf91f0a0037ae85644df3df0938d549322cbfb59b7dfd9be01939d4",
    "G": "841dd3c72d594cfd2c7f921b84f10a2305ab9957b531771662fb7f03eff97207",
    "L": "e6a26c660ae9cd4fc8cb4b1effadcea38a4794b66523f0fea2c02c5e02b77f25",
    "P": "5c92a663de35357a4bf0ff7eb0838bc5e178030afe7d769ead5632f50821652f",
}
CATALOGUE_EXTRACT_SHA256 = "8c3a17569e6fa53d9a769ab584fed79a865118c4a27018facd619b4a5cc5cdd3"
BOX_LABEL_SHA256 = "522952db10076af64be712e2843033dd69747428c15bf1c5a9085a82db840fa0"
VIEW_CODES = {"plan": "G", "front": "A", "side": "L"}
EXPECTED_PATH_COUNTS = {"plan": 253, "front": 790, "side": 225}
EXPECTED_BOUNDS = {
    "plan": [500.0, 336.5],
    "front": [500.0, 1215.000112],
    "side": [336.5, 1215.0],
}


def de_boor(degree: int, knots: list[float], controls: list[tuple[float, ...]], value: float):
    count = len(controls)
    if value >= knots[count]:
        span = count - 1
    else:
        spans = [
            index
            for index in range(degree, count)
            if knots[index] <= value < knots[index + 1]
        ]
        if not spans:
            raise RuntimeError("unsupported Geberit spline knot span")
        span = spans[-1]
    work = [list(controls[span - degree + index]) for index in range(degree + 1)]
    for level in range(1, degree + 1):
        for index in range(degree, level - 1, -1):
            source = span - degree + index
            denominator = knots[source + degree - level + 1] - knots[source]
            alpha = 0.0 if abs(denominator) < 1e-12 else (value - knots[source]) / denominator
            work[index] = [
                (1.0 - alpha) * work[index - 1][axis] + alpha * work[index][axis]
                for axis in range(len(work[index]))
            ]
    return work[degree]


def sample_spline(entity: dict) -> list[list[float]]:
    degree = int(entity["degree"])
    knots = [float(value) for value in entity["knots"]]
    rational = bool(entity.get("rational"))
    controls = []
    for point in entity["ctrl_pts"]:
        weight = float(point.get("w", 1.0))
        if rational:
            controls.append((float(point["x"]) * weight, float(point["y"]) * weight, weight))
        else:
            controls.append((float(point["x"]), float(point["y"])))
    if len(knots) != len(controls) + degree + 1:
        raise RuntimeError("unsupported Geberit spline knot/control identity")
    start, end = knots[degree], knots[len(controls)]
    steps = max(16, min(200, len(controls) * 10))
    path = []
    for index in range(steps + 1):
        point = de_boor(degree, knots, controls, start + (end - start) * index / steps)
        if rational:
            if abs(point[2]) < 1e-12:
                raise RuntimeError("invalid zero-weight Geberit rational spline")
            point = (point[0] / point[2], point[1] / point[2])
        path.append([round(point[0], 6), round(point[1], 6)])
    return path


def sample_circle(entity: dict) -> list[list[float]]:
    center_x, center_y = map(float, entity["center"][:2])
    radius = float(entity["radius"])
    return [
        [
            round(center_x + radius * math.cos(math.tau * index / 72.0), 6),
            round(center_y + radius * math.sin(math.tau * index / 72.0), 6),
        ]
        for index in range(73)
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
    if kind == "CIRCLE":
        return sample_circle(entity)
    raise RuntimeError(f"unsupported Geberit contour entity: {kind}")


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
        if item.get("entity") in {"LINE", "SPLINE", "ARC", "ELLIPSE", "CIRCLE"}
        and layers.get(handle_value(item["layer"])) == CONTOUR_LAYER
    ]
    paths = [entity_path(entity) for entity in entities]
    return {
        "view": view,
        "native_dwg_code": code,
        "source_dwg": relative(path),
        "source_dwg_url": f"https://cdn.data.geberit.com/cad/{ARTICLE}_{code}.dwg",
        "source_dwg_sha256": sha256(path),
        "source_kind": "native_dwg",
        "units": "mm",
        "native_contour_layer": CONTOUR_LAYER,
        "native_entity_counts": dict(sorted(Counter(entity["entity"] for entity in entities).items())),
        "path_count": len(paths),
        "point_count": sum(len(path) for path in paths),
        "contour_bounds_mm": path_bounds(paths),
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
            raise RuntimeError(f"official Geberit Duofix {code} DWG hash mismatch")
    if sha256(LOCAL_CATALOGUE_EXTRACT) != CATALOGUE_EXTRACT_SHA256:
        raise RuntimeError("official 224.212.00.2 catalogue extract hash mismatch")
    if sha256(LOCAL_BOX_LABEL) != BOX_LABEL_SHA256:
        raise RuntimeError("224.212.00.2 received-box label hash mismatch")
    views = {
        view: extract_view(sources[code], view, code)
        for view, code in VIEW_CODES.items()
    }
    for view, expected in EXPECTED_PATH_COUNTS.items():
        if views[view]["path_count"] != expected:
            raise RuntimeError(f"official Duofix {view} contour path count drifted")
    for view, expected in EXPECTED_BOUNDS.items():
        actual = views[view]["contour_bounds_mm"]["size"]
        if any(abs(actual[index] - expected[index]) > 0.01 for index in range(2)):
            raise RuntimeError(f"official Duofix {view} contour bounds drifted: {actual}")
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/geberit_duofix_sigma_224_212_linework.py",
        "manufacturer": "Geberit",
        "family": "Duofix element for wall-hung WC with Sigma concealed cistern 12 cm",
        "article_number": ARTICLE,
        "ifc_type_name": IFC_TYPE_NAME,
        "official_catalogue": OFFICIAL_CATALOGUE,
        "local_official_catalogue_extract": {
            "path": relative(LOCAL_CATALOGUE_EXTRACT),
            "sha256": CATALOGUE_EXTRACT_SHA256,
        },
        "received_box_label": {
            "path": relative(LOCAL_BOX_LABEL),
            "sha256": BOX_LABEL_SHA256,
            "scope": "project received-product identity evidence for one of two IFC family instances",
        },
        "source_kind": "native_dwg",
        "scope": SCOPE,
        "nominal_dimensions_mm": {"width": 500.0, "height": 1120.0, "cistern_depth": 120.0},
        "three_view_code_mapping": {
            "G": "Grundriss / plan",
            "A": "Ansicht / front elevation",
            "L": "left side elevation",
            "P": "official 3D model; archived for identity only and never substituted for 2D views",
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
        "views": {
            view: {"paths": data["path_count"], "bounds": data["contour_bounds_mm"]}
            for view, data in views.items()
        },
        "pass": True,
    }, indent=2))


if __name__ == "__main__":
    main()
