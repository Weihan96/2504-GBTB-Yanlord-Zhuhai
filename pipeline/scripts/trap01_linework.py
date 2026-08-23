#!/usr/bin/env python3
"""Parse Geberit 151.116.11.1 native DWG family views for TRAP01."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, relative, sha256, write_json
import geberit_146_140_linework as shared


PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
SOURCE_DIR = PRODUCT_DIR / "official-source"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
ARTICLE = "151.116.11.1"
PRODUCT_PAGE = "https://catalog.geberit.us/en-US/product/PRO_185224"
SCOPE = "exact Geberit 151.116.11.1 adjustable family reference; project instance is a shortened installation configuration; not a project shop drawing"
EXPECTED = {
    "A": "01fa50da7cf39c6a6ab2a89339d1e4523270cb0943e8dce148f3d471f1a4fff2",
    "G": "af55992e60aa1d84c9834da995f6ef8269593f43f46a6311cdf3e2e59530dd3f",
    "L": "761eec37c809fc4319e9beaba6eaf472e19eefcbd16adf2815a4644d0a36fdc8",
    "P": "c99cee2a0c9ae94b3a5535db4973982795eba64cedc5da135d5ca8f12264a8e1",
}
VIEW_CODES = {"plan": "G", "front": "L", "side": "A"}
EXPECTED_PATH_COUNTS = {"plan": 92, "front": 77, "side": 117}
IFC_BODY_SIZE_XYZ_MM = [252.502579, 76.540974, 191.93998]


def de_boor_vector(degree: int, knots: list[float], controls: list[tuple[float, ...]], value: float):
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
            work[index] = [
                (1.0 - alpha) * work[index - 1][axis] + alpha * work[index][axis]
                for axis in range(len(work[index]))
            ]
    return work[degree]


def sample_spline(entity: dict) -> list[list[float]]:
    degree = int(entity["degree"])
    knots = [float(value) for value in entity["knots"]]
    controls = []
    for point in entity["ctrl_pts"]:
        weight = float(point.get("w", 1.0))
        controls.append((float(point["x"]) * weight, float(point["y"]) * weight, weight))
    if len(knots) != len(controls) + degree + 1:
        raise RuntimeError("unsupported TRAP01 Geberit spline knot identity")
    start, end = knots[degree], knots[len(controls)]
    steps = max(16, min(160, len(controls) * 8))
    paths = []
    for index in range(steps + 1):
        x, y, weight = de_boor_vector(
            degree,
            knots,
            controls,
            start + (end - start) * index / steps,
        )
        if abs(weight) < 1e-12:
            raise RuntimeError("TRAP01 Geberit rational spline has a zero weight")
        paths.append([round(x / weight, 6), round(y / weight, 6)])
    return paths


def entity_path(entity: dict) -> list[list[float]]:
    if entity["entity"] == "SPLINE":
        return sample_spline(entity)
    return shared.entity_path(entity)


def normalise(paths, view: str):
    if view == "plan":
        oriented = [[[-float(y), float(x)] for x, y in path] for path in paths]
        orientation = "rotate_native_G_quarter_turn_for_ifc_local_xy"
    else:
        oriented = [[[float(x), float(y)] for x, y in path] for path in paths]
        orientation = "native_axes"
    bounds = shared.path_bounds(oriented)
    minimum = bounds["minimum"]
    result = [
        [[round(x - minimum[0], 6), round(y - minimum[1], 6)] for x, y in path]
        for path in oriented
    ]
    return result, orientation


def extract_view(path: Path, view: str, code: str) -> dict:
    payload = shared.dwg_json(path)
    if payload.get("HEADER", {}).get("INSUNITS") != 4:
        raise RuntimeError(f"{path.name}: native DWG is not millimetres")
    layers = {
        shared.handle_value(item["handle"]): item["name"]
        for item in payload["OBJECTS"]
        if item.get("object") == "LAYER"
    }
    entities = [
        item
        for item in payload["OBJECTS"]
        if item.get("entity") in {"LINE", "SPLINE", "ARC", "ELLIPSE"}
        and layers.get(shared.handle_value(item["layer"])) == shared.CONTOUR_LAYER
    ]
    raw_paths = [entity_path(entity) for entity in entities]
    paths, orientation = normalise(raw_paths, view)
    counts = {
        kind: sum(entity["entity"] == kind for entity in entities)
        for kind in sorted({entity["entity"] for entity in entities})
    }
    return {
        "view": view,
        "native_dwg_code": code,
        "source_dwg": relative(path),
        "source_dwg_url": f"https://cdn.data.geberit.com/cad/{ARTICLE}_{code}.dwg",
        "source_dwg_sha256": sha256(path),
        "source_kind": "native_dwg",
        "units": "mm",
        "native_contour_layer": shared.CONTOUR_LAYER,
        "native_entity_counts": counts,
        "path_count": len(paths),
        "point_count": sum(len(path) for path in paths),
        "view_orientation_transform": orientation,
        "bounds_mm": shared.path_bounds(paths),
        "paths_mm": paths,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=SOURCE_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source_dir = args.source_dir.resolve()
    sources = {code: source_dir / f"{ARTICLE}_{code}.dwg" for code in EXPECTED}
    for code, path in sources.items():
        if not path.is_file() or sha256(path) != EXPECTED[code]:
            raise RuntimeError(f"official TRAP01 Geberit {code} DWG hash mismatch")
    views = {view: extract_view(sources[code], view, code) for view, code in VIEW_CODES.items()}
    if {view: item["path_count"] for view, item in views.items()} != EXPECTED_PATH_COUNTS:
        raise RuntimeError("TRAP01 native-DWG path counts drifted")
    family_sizes = {view: item["bounds_mm"]["size"] for view, item in views.items()}
    configured_sizes = {
        "plan": [IFC_BODY_SIZE_XYZ_MM[0], IFC_BODY_SIZE_XYZ_MM[1]],
        "front": [IFC_BODY_SIZE_XYZ_MM[0], IFC_BODY_SIZE_XYZ_MM[2]],
        "side": [IFC_BODY_SIZE_XYZ_MM[1], IFC_BODY_SIZE_XYZ_MM[2]],
    }
    fixed_width_deltas = {
        "plan_width_mm": round(abs(family_sizes["plan"][1] - configured_sizes["plan"][1]), 6),
        "side_width_mm": round(abs(family_sizes["side"][0] - configured_sizes["side"][0]), 6),
    }
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/trap01_linework.py",
        "manufacturer": "Geberit",
        "family": "dip tube trap for washbasin, space-saving model, horizontal outlet",
        "article_number": ARTICLE,
        "ifc_type_name": "TRAP01",
        "ifc_type_description": "Space Saving Dip Tube Trap",
        "product_page": PRODUCT_PAGE,
        "source_kind": "native_dwg",
        "scope": SCOPE,
        "three_view_code_mapping": {
            "G": "plan, rotated to project local X/Y axes",
            "L": "front elevation, project local X/Z axes",
            "A": "side elevation, project local Y/Z axes",
            "P": "official 3D model, archived for exact article identity only",
        },
        "official_sources": {
            code: {"path": relative(path), "url": f"https://cdn.data.geberit.com/cad/{ARTICLE}_{code}.dwg", "sha256": sha256(path)}
            for code, path in sources.items()
        },
        "configuration_cross_check": {
            "official_family_default_full_envelope_by_view_mm": family_sizes,
            "project_ifc_configured_body_by_view_mm": configured_sizes,
            "fixed_diameter_width_absolute_delta_mm": fixed_width_deltas,
            "fixed_diameter_tolerance_mm": 0.02,
            "fixed_diameter_match_pass": max(fixed_width_deltas.values()) <= 0.02,
            "official_product_adjustment_ranges_mm": {"horizontal_extension": [0.0, 252.0], "h": [85.0, 334.0]},
            "project_instance_is_shortened_configuration": True,
            "official_default_family_paths_used_as_project_representation": False,
            "reason": "The exact official article CAD is drawn at a longer adjustable configuration. Its fixed 76.54 mm body width matches the project Body, but forcing the full family envelope into the shortened project instance would be geometrically false.",
            "pass": max(fixed_width_deltas.values()) <= 0.02,
        },
        "views": views,
        "pass": max(fixed_width_deltas.values()) <= 0.02,
    }
    write_json(args.output.resolve(), payload)
    print(json.dumps({"output": relative(args.output), "path_counts": EXPECTED_PATH_COUNTS, "fixed_width_deltas": fixed_width_deltas, "pass": payload["pass"]}, indent=2))


if __name__ == "__main__":
    main()
