#!/usr/bin/env python3
"""Extract exact Plan/Front/Side linework from the official Gessi 54145 G000 DWG."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from geberit_146_140_linework import dwg_json, handle_value, path_bounds
from gessi316_54294_linework import bounds, gessi_entity_path


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54145"
SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_DWG = SOURCE_DIR / "GPF5414500000G000_3.dwg"
SOURCE_ZIP = SOURCE_DIR / "GPF5414500000G000_arc.zip"
SOURCE_PDF = SOURCE_DIR / "GPF5414500000G000_1.pdf"
REVALIDATION = SOURCE_DIR / "official-source-revalidation.json"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
SOURCE_DWG_URL = "https://gessistorage.blob.core.windows.net/zwa/GPF5414500000G000_arc.zip"
SOURCE_PDF_URL = "https://gessistorage.blob.core.windows.net/zc4/GPF5414500000G000_1.pdf"
EXPECTED = {
    "dwg": "9978b68468a61875acb0736aab45a62a98efadfd6be61d347e7ecaa94e08fd09",
    "zip": "9678c2de6a276c0d1f6e3b29e6c39763c1be0763bccf1b8974e68408ede699fd",
    "pdf": "148a16717193cbc066d0314c6d5369492ed582da0f026dbbe6a58d06f8279025",
}
SOURCE_LAYER_HANDLE = 509
SOURCE_LAYER_NAME = "00_COMPONENTI"
ENTITY_TYPES = {"LINE", "ARC", "ELLIPSE", "SPLINE", "CIRCLE", "LWPOLYLINE"}
VIEW_WINDOWS_MM = {
    "plan": (500.0, 400.0, 1200.0, 760.0),
    "front": (500.0, 760.0, 1200.0, 1050.0),
    "side": (100.0, 400.0, 500.0, 760.0),
}
EXPECTED_PATH_COUNTS = {"plan": 22, "front": 615, "side": 653}
IFC_BODY_LOCAL_XYZ_MM = [599.818665, 299.818832, 119.07444]
SCOPE = "official Gessi exact 54145 G000 family reference; not a project shop drawing"


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
        normalized = [[[round(y - minimum_y, 6), round(x - maximum_x, 6)] for x, y in path] for path in paths]
        origin = "rotated_minimum_depth_and_top_z"
    return normalized, {
        "source_sheet_bounds_mm": source_bounds,
        "normalization_scale": 1.0,
        "normalization_origin": origin,
        "source_axes_swapped": view == "side",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if sha256(SOURCE_DWG) != EXPECTED["dwg"] or sha256(SOURCE_ZIP) != EXPECTED["zip"] or sha256(SOURCE_PDF) != EXPECTED["pdf"]:
        raise RuntimeError("Gessi 54145 official source hash gate failed")
    revalidation = load_json(REVALIDATION)
    if (
        revalidation.get("pass") is not True
        or revalidation.get("public_access", {}).get("exact_54145_g000_native_dwg_publicly_downloadable") is not True
        or revalidation.get("drawing_geometry_policy", {}).get("g001_variant_cad_used") is not False
        or revalidation.get("drawing_geometry_policy", {}).get("adjacent_gessi_product_cad_used") is not False
        or revalidation.get("drawing_geometry_policy", {}).get("third_party_cad_used") is not False
    ):
        raise RuntimeError("Gessi 54145 source revalidation gate failed")
    data = dwg_json(SOURCE_DWG)
    selected = {view: [] for view in VIEW_WINDOWS_MM}
    source_type_counts = {view: {} for view in VIEW_WINDOWS_MM}
    selected_handles: set[int] = set()
    eligible_count = 0
    for entity in data["OBJECTS"]:
        kind = entity.get("entity")
        if kind not in ENTITY_TYPES or handle_value(entity.get("layer")) != SOURCE_LAYER_HANDLE:
            continue
        eligible_count += 1
        path = gessi_entity_path(entity)
        value = bounds(path)
        matches = [view for view, window in VIEW_WINDOWS_MM.items() if in_window(value, window)]
        if len(matches) != 1:
            raise RuntimeError(f"Gessi 54145 entity must match exactly one orthographic window: {entity.get('handle')} / {matches}")
        view = matches[0]
        handle = handle_value(entity["handle"])
        selected[view].append((handle, path))
        selected_handles.add(handle)
        source_type_counts[view][kind] = source_type_counts[view].get(kind, 0) + 1
    if len(selected_handles) != eligible_count:
        raise RuntimeError("Gessi 54145 orthographic windows did not consume the complete native drawing layer")
    views = {}
    for view, items in selected.items():
        items.sort(key=lambda item: item[0])
        paths, metadata = normalize(view, [path for _, path in items])
        if len(paths) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi 54145 {view} path-count drift: {len(paths)}")
        views[view] = {
            "paths_mm": paths,
            "path_count": len(paths),
            "bounds_mm": path_bounds(paths),
            "source_sheet_bounds_mm": metadata["source_sheet_bounds_mm"],
            "source_window_mm": list(VIEW_WINDOWS_MM[view]),
            "source_layer": SOURCE_LAYER_NAME,
            "source_layer_handle": SOURCE_LAYER_HANDLE,
            "source_entity_type_counts": source_type_counts[view],
            "source_kind": "native_dwg",
            "source_dwg_sha256": EXPECTED["dwg"],
            "normalization_scale": 1.0,
            "normalization_origin": metadata["normalization_origin"],
            "source_axes_swapped": metadata["source_axes_swapped"],
        }
    plan_size = views["plan"]["bounds_mm"]["size"]
    front_size = views["front"]["bounds_mm"]["size"]
    side_size = views["side"]["bounds_mm"]["size"]
    ifc_projection_sizes = {
        "plan": IFC_BODY_LOCAL_XYZ_MM[:2],
        "front": [IFC_BODY_LOCAL_XYZ_MM[0], IFC_BODY_LOCAL_XYZ_MM[2]],
        "side": [IFC_BODY_LOCAL_XYZ_MM[1], IFC_BODY_LOCAL_XYZ_MM[2]],
    }
    native_sizes = {"plan": plan_size, "front": front_size, "side": side_size}
    deltas = {
        view: [round(abs(native_sizes[view][axis] - ifc_projection_sizes[view][axis]), 6) for axis in range(2)]
        for view in native_sizes
    }
    maximum_delta = max(value for view_delta in deltas.values() for value in view_delta)
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/gessi316_54145_linework.py",
        "manufacturer": "Gessi",
        "family": "Gessi316 Meccanica",
        "article_number": "54145",
        "configuration": "G000",
        "ifc_type_name": "Gessi316 54145",
        "source_kind": "native_dwg",
        "source_label_zh": "Gessi 官方精确型号 54145 G000 原生 DWG 图纸表达",
        "source_label_en": "drawing representation from the exact Gessi 54145 G000 official native DWG",
        "scope": SCOPE,
        "sheet": {"units": "mm", "source_layer": SOURCE_LAYER_NAME, "eligible_entity_count": eligible_count, "all_eligible_entities_partitioned_once": True},
        "official_sources": {
            "native_dwg_zip": {"path": relative(SOURCE_ZIP), "url": SOURCE_DWG_URL, "sha256": EXPECTED["zip"]},
            "native_dwg": {"path": relative(SOURCE_DWG), "url": SOURCE_DWG_URL, "sha256": EXPECTED["dwg"]},
            "technical_vector_pdf": {"path": relative(SOURCE_PDF), "url": SOURCE_PDF_URL, "sha256": EXPECTED["pdf"]},
            "revalidation": {"path": relative(REVALIDATION), "sha256": sha256(REVALIDATION)},
        },
        "identity_gates": {
            "54145_g000_native_dwg_used": True,
            "g001_variant_cad_used": False,
            "adjacent_gessi_product_cad_used": False,
            "third_party_cad_used": False,
        },
        "nominal_dimension_cross_check": {
            "official_api_width_depth_height_mm": [300.0, 600.0, 115.0],
            "official_api_showerhead_diameter_mm": 300.0,
            "official_technical_pdf_arm_reach_mm": 600.0,
            "official_technical_pdf_nominal_height_mm": 115.0,
            "nominal_height_note": "The 115 mm annotation is a product datum; the 118.95 mm orthographic envelope includes the spray-face/nozzle linework.",
            "project_ifc_body_local_xyz_mm": IFC_BODY_LOCAL_XYZ_MM,
            "native_dwg_plan_envelope_mm": plan_size,
            "native_dwg_front_envelope_mm": front_size,
            "native_dwg_side_envelope_mm": side_size,
            "ifc_projection_absolute_delta_mm": deltas,
            "maximum_ifc_projection_delta_mm": round(maximum_delta, 6),
            "tolerance_mm": 0.5,
            "geometry_stretched": False,
            "pass": maximum_delta <= 0.5,
        },
        "views": views,
        "pass": maximum_delta <= 0.5,
    }
    write_json(args.output.resolve(), payload)
    if not payload["pass"]:
        raise RuntimeError("Gessi 54145 native DWG / IFC projection cross-check failed")
    print(relative(args.output))


if __name__ == "__main__":
    main()
