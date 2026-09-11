#!/usr/bin/env python3
"""Extract exact Plan/Front/Side linework from the official Gessi 54146 G000 DWG."""

from __future__ import annotations

import argparse
import math
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from geberit_146_140_linework import dwg_json, handle_value, path_bounds
from gessi316_54294_linework import bounds, gessi_entity_path


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54146"
SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_DWG = SOURCE_DIR / "GPF5414600000G000_3.dwg"
SOURCE_ZIP = SOURCE_DIR / "GPF5414600000G000_arc.zip"
SOURCE_PDF = SOURCE_DIR / "GPF5414600000G000_1.pdf"
REVALIDATION = SOURCE_DIR / "official-source-revalidation.json"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dwg-linework.json"
SOURCE_DWG_URL = "https://gessistorage.blob.core.windows.net/zwa/GPF5414600000G000_arc.zip"
SOURCE_PDF_URL = "https://gessistorage.blob.core.windows.net/zc4/GPF5414600000G000_1.pdf"
EXPECTED = {
    "dwg": "c8ddb90f61565d5273a32574f777812bbb1d9ef33e2b98479f1e7c381e69950d",
    "zip": "aaf9964d57674a03b7e7a8f9eb647df84af74bab1b575c9bad306d83bd40237e",
    "pdf": "ee7909886b1181306902d0a0e146b73f938bf96e3e12e8f1ba27ded996b2b8ce",
}
SOURCE_LAYER_HANDLE = 509
SOURCE_LAYER_NAME = "00_COMPONENTI"
ENTITY_TYPES = {"LINE", "ARC", "ELLIPSE", "SPLINE", "CIRCLE", "LWPOLYLINE"}
VIEW_WINDOWS_MM = {
    "plan": (550.0, 400.0, 1000.0, 800.0),
    "front": (550.0, 800.0, 1000.0, 1250.0),
    "side": (50.0, 400.0, 550.0, 800.0),
}
EXPECTED_PATH_COUNTS = {"plan": 16, "front": 607, "side": 660}
IFC_BODY_LOCAL_XYZ_MM = [300.0, 299.944916, 279.371241]
SCOPE = "official Gessi exact 54146 G000 family reference; not a project shop drawing"


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
        normalized = [[[round(y - minimum_y, 6), round(maximum_x - x, 6)] for x, y in path] for path in paths]
        origin = "rotated_minimum_depth_and_bottom_z_after_vertical_reflection"
    return normalized, {
        "source_sheet_bounds_mm": source_bounds,
        "normalization_origin": origin,
        "source_axes_swapped": view == "side",
        "vertical_reflection_applied": view == "side",
        "normalization_translation_mm": (
            [-minimum_y, maximum_x] if view == "side" else None
        ),
        "normalization_scale": 1.0,
    }


def line_density(paths: list[list[list[float]]]) -> dict:
    node_count = sum(len(path) for path in paths)
    return {
        "entity_count": len(paths),
        "path_count": len(paths),
        "node_count": node_count,
        "segment_count": sum(max(0, len(path) - 1) for path in paths),
        "closed_path_count": sum(
            1 for path in paths
            if len(path) > 2 and math.dist(path[0], path[-1]) <= 0.000001
        ),
    }


def fine_spray_detail_indices(view: str, paths: list[list[list[float]]]) -> list[int]:
    """Identify only sub-8 mm paths in the lower spray-face band.

    The 10.25 mm band is tied to the verified 280.25 mm native-DWG
    envelope, so stem joints and ceiling installation nodes are excluded.
    Plan is intentionally untouched because its small concentric circles are
    installation/connection identity, not elevation spray-nozzle texture.
    """
    if view == "plan":
        return []
    lower = min(point[1] for path in paths for point in path)
    indices = []
    for index, path in enumerate(paths):
        xs = [point[0] for point in path]
        ys = [point[1] for point in path]
        maximum_dimension = max(max(xs) - min(xs), max(ys) - min(ys))
        centre_y = (min(ys) + max(ys)) / 2.0
        if maximum_dimension <= 8.0 and centre_y <= lower + 10.25:
            indices.append(index)
    return indices


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if sha256(SOURCE_DWG) != EXPECTED["dwg"] or sha256(SOURCE_ZIP) != EXPECTED["zip"] or sha256(SOURCE_PDF) != EXPECTED["pdf"]:
        raise RuntimeError("Gessi 54146 official source hash gate failed")
    revalidation = load_json(REVALIDATION)
    policy = revalidation.get("drawing_geometry_policy", {})
    if (
        revalidation.get("pass") is not True
        or revalidation.get("public_access", {}).get("exact_54146_g000_native_dwg_publicly_downloadable") is not True
        or policy.get("g001_variant_cad_used") is not False
        or policy.get("wall_mounted_54145_cad_used") is not False
        or policy.get("adjacent_gessi_product_cad_used") is not False
        or policy.get("third_party_cad_used") is not False
    ):
        raise RuntimeError("Gessi 54146 source revalidation gate failed")
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
            raise RuntimeError(f"Gessi 54146 entity must match exactly one orthographic window: {entity.get('handle')} / {matches}")
        view = matches[0]
        handle = handle_value(entity["handle"])
        selected[view].append((handle, path))
        selected_handles.add(handle)
        source_type_counts[view][kind] = source_type_counts[view].get(kind, 0) + 1
    if len(selected_handles) != eligible_count:
        raise RuntimeError("Gessi 54146 orthographic windows did not consume the complete native drawing layer")
    views = {}
    for view, items in selected.items():
        items.sort(key=lambda item: item[0])
        paths, metadata = normalize(view, [path for _, path in items])
        if len(paths) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi 54146 {view} path-count drift: {len(paths)}")
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
            "vertical_reflection_applied": metadata["vertical_reflection_applied"],
            "normalization_translation_mm": metadata["normalization_translation_mm"],
            "line_density": line_density(paths),
        }
        detail_indices = fine_spray_detail_indices(view, paths)
        detail_paths = [paths[index] for index in detail_indices]
        detail_density = line_density(detail_paths)
        views[view]["fine_spray_nozzle_detail"] = {
            "classification": "sub_8mm_paths_in_lower_10_25mm_spray_face_band",
            **detail_density,
            "path_share_percent": round(100.0 * len(detail_indices) / len(paths), 3),
            "node_share_percent": round(
                100.0 * detail_density["node_count"] / views[view]["line_density"]["node_count"],
                3,
            ),
        }
    native_sizes = {view: views[view]["bounds_mm"]["size"] for view in views}
    ifc_projection_sizes = {
        "plan": IFC_BODY_LOCAL_XYZ_MM[:2],
        "front": [IFC_BODY_LOCAL_XYZ_MM[0], IFC_BODY_LOCAL_XYZ_MM[2]],
        "side": [IFC_BODY_LOCAL_XYZ_MM[1], IFC_BODY_LOCAL_XYZ_MM[2]],
    }
    deltas = {
        view: [round(abs(native_sizes[view][axis] - ifc_projection_sizes[view][axis]), 6) for axis in range(2)]
        for view in native_sizes
    }
    maximum_delta = max(value for view_delta in deltas.values() for value in view_delta)
    payload = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/gessi316_54146_linework.py",
        "manufacturer": "Gessi",
        "family": "Gessi316 Meccanica",
        "article_number": "54146",
        "configuration": "G000",
        "ifc_type_name": "Gessi316 54146",
        "source_kind": "native_dwg",
        "source_label_zh": "Gessi 官方精确型号 54146 G000 原生 DWG 图纸表达",
        "source_label_en": "drawing representation from the exact Gessi 54146 G000 official native DWG",
        "scope": SCOPE,
        "sheet": {"units": "mm", "source_layer": SOURCE_LAYER_NAME, "eligible_entity_count": eligible_count, "all_eligible_entities_partitioned_once": True},
        "official_sources": {
            "native_dwg_zip": {"path": relative(SOURCE_ZIP), "url": SOURCE_DWG_URL, "sha256": EXPECTED["zip"]},
            "native_dwg": {"path": relative(SOURCE_DWG), "url": SOURCE_DWG_URL, "sha256": EXPECTED["dwg"]},
            "technical_vector_pdf": {"path": relative(SOURCE_PDF), "url": SOURCE_PDF_URL, "sha256": EXPECTED["pdf"]},
            "revalidation": {"path": relative(REVALIDATION), "sha256": sha256(REVALIDATION)},
        },
        "identity_gates": {
            "54146_g000_native_dwg_used": True,
            "g001_variant_cad_used": False,
            "wall_mounted_54145_cad_used": False,
            "adjacent_gessi_product_cad_used": False,
            "third_party_cad_used": False,
        },
        "nominal_dimension_cross_check": {
            "official_api_width_depth_height_mm": [300.0, 300.0, 276.0],
            "official_api_showerhead_diameter_mm": 300.0,
            "official_technical_pdf_nominal_height_mm": 276.0,
            "nominal_height_note": "The 276 mm annotation is a product datum; the 280.25 mm orthographic envelope includes the spray-face/nozzle linework.",
            "project_ifc_body_local_xyz_mm": IFC_BODY_LOCAL_XYZ_MM,
            "native_dwg_plan_envelope_mm": native_sizes["plan"],
            "native_dwg_front_envelope_mm": native_sizes["front"],
            "native_dwg_side_envelope_mm": native_sizes["side"],
            "ifc_projection_absolute_delta_mm": deltas,
            "maximum_ifc_projection_delta_mm": round(maximum_delta, 6),
            "tolerance_mm": 1.0,
            "geometry_stretched": False,
            "pass": maximum_delta <= 1.0,
        },
        "views": views,
        "line_density_audit": {
            "three_view_entity_count": sum(view["line_density"]["entity_count"] for view in views.values()),
            "three_view_path_count": sum(view["line_density"]["path_count"] for view in views.values()),
            "three_view_node_count": sum(view["line_density"]["node_count"] for view in views.values()),
            "three_view_segment_count": sum(view["line_density"]["segment_count"] for view in views.values()),
            "fine_spray_nozzle_detail_path_count": sum(view["fine_spray_nozzle_detail"]["path_count"] for view in views.values()),
            "fine_spray_nozzle_detail_node_count": sum(view["fine_spray_nozzle_detail"]["node_count"] for view in views.values()),
        },
        "pass": maximum_delta <= 1.0,
    }
    write_json(args.output.resolve(), payload)
    if not payload["pass"]:
        raise RuntimeError("Gessi 54146 native DWG / IFC projection cross-check failed")
    print(relative(args.output))


if __name__ == "__main__":
    main()
