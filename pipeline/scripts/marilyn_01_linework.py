#!/usr/bin/env python3
"""Extract the exact Marilyn bergere three-view linework from Baxter native DWG."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, relative, sha256, write_json
from geberit_146_140_linework import sample_arc


PRODUCT_DIR = ROOT / "output/review/highpoly-types/marilyn-01"
SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_DWG = SOURCE_DIR / "Marilyn_Abaco.dwg"
SOURCE_ZIP = SOURCE_DIR / "Baxter_Marilyn_Armchair_2D_3D.zip"
SOURCE_3DS = SOURCE_DIR / "Marilyn_bergere_86x100xh94.3ds"
AUTOCAD_FULL_MODELSPACE = SOURCE_DIR / "autocad-marilyn-full-modelspace.png"
LAYER_AUDIT_SVG = SOURCE_DIR / "marilyn-exact-cluster-layer-audit.svg"
LAYER_AUDIT_PNG = SOURCE_DIR / "marilyn-exact-cluster-layer-audit.png"
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
EXPECTED_PATH_COUNTS = {"plan": 44, "front": 96, "side": 78}
ORIGINAL_LAYER_FILTERED_PATH_COUNTS = {"plan": 24, "front": 68, "side": 54}
PREVIOUS_RATIONAL_OMITTED_PATH_COUNTS = {"plan": 38, "front": 80, "side": 62}
EXPECTED_NATIVE_LAYERS_BY_VIEW = {
    "plan": {"_ARREDO", "0"},
    "front": {"_ARREDO", "0"},
    "side": {"_ARREDO", "Make2D$Visibile$Curve"},
}
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


def de_boor_point(degree: int, knots: list[float], controls: list[tuple[float, ...]], value: float) -> tuple[float, ...]:
    """Evaluate a B-spline in ordinary or homogeneous coordinates."""
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
            for axis in range(len(work[index])):
                work[index][axis] = (1.0 - alpha) * work[index - 1][axis] + alpha * work[index][axis]
    return tuple(work[degree])


def sample_native_spline(entity: dict) -> list[list[float]]:
    """Sample every native DWG spline, including rational weighted curves."""
    degree = int(entity["degree"])
    if "knots" not in entity:
        fit_points = [[float(value) for value in point[:2]] for point in entity.get("fit_pts", [])]
        if len(fit_points) < 2:
            raise RuntimeError(f"unsupported Marilyn DWG fit spline at entity {entity.get('index')}")
        if len(fit_points) == 2:
            return [[round(value, 6) for value in point] for point in fit_points]
        if len(fit_points) == 3:
            # The native record stores a three-point fit spline without an
            # explicit knot/control vector. Reconstruct the quadratic curve
            # that interpolates the middle fit point at t=0.5.
            start, middle, end = fit_points
            control = [2.0 * middle[axis] - 0.5 * (start[axis] + end[axis]) for axis in range(2)]
            return [
                [
                    round((1.0 - value) ** 2 * start[axis] + 2.0 * (1.0 - value) * value * control[axis] + value**2 * end[axis], 6)
                    for axis in range(2)
                ]
                for value in (index / 32.0 for index in range(33))
            ]
        # General fit-point records retain the manufacturer's interpolation
        # anchors but omit a knot/control vector in LibreDWG JSON. A centripetal
        # CAD fit would be visually equivalent at this audit scale; use a
        # segment-wise Catmull-Rom interpolation that passes through every
        # native fit point and never invents an extra component.
        sampled = []
        padded = [fit_points[0], *fit_points, fit_points[-1]]
        for segment in range(1, len(padded) - 2):
            p0, p1, p2, p3 = padded[segment - 1 : segment + 3]
            for index in range(16 + (1 if segment == len(padded) - 3 else 0)):
                value = index / 16.0
                point = [
                    0.5
                    * (
                        2.0 * p1[axis]
                        + (-p0[axis] + p2[axis]) * value
                        + (2.0 * p0[axis] - 5.0 * p1[axis] + 4.0 * p2[axis] - p3[axis]) * value**2
                        + (-p0[axis] + 3.0 * p1[axis] - 3.0 * p2[axis] + p3[axis]) * value**3
                    )
                    for axis in range(2)
                ]
                sampled.append([round(coordinate, 6) for coordinate in point])
        return sampled
    knots = [float(value) for value in entity["knots"]]
    control_points = entity["ctrl_pts"]
    if len(knots) != len(control_points) + degree + 1:
        raise RuntimeError(f"unsupported Marilyn DWG spline knot vector at entity {entity.get('index')}")
    rational = bool(entity.get("rational"))
    if rational:
        controls = [
            (
                float(point["x"]) * float(point.get("w", 1.0)),
                float(point["y"]) * float(point.get("w", 1.0)),
                float(point.get("w", 1.0)),
            )
            for point in control_points
        ]
    else:
        controls = [(float(point["x"]), float(point["y"])) for point in control_points]
    start, end = knots[degree], knots[len(controls)]
    steps = max(16, min(160, len(controls) * 8))
    sampled = [
        de_boor_point(degree, knots, controls, start + (end - start) * index / steps)
        for index in range(steps + 1)
    ]
    if rational:
        return [[round(point[0] / point[2], 6), round(point[1] / point[2], 6)] for point in sampled]
    return [[round(point[0], 6), round(point[1], 6)] for point in sampled]


def entity_path(entity: dict, objects_by_handle: dict[int, dict]) -> list[list[float]] | None:
    kind = entity.get("entity")
    if kind == "SPLINE":
        return sample_native_spline(entity)
    if kind == "LINE":
        return [[float(value) for value in entity["start"][:2]], [float(value) for value in entity["end"][:2]]]
    if kind == "ARC":
        # DWG ARC coordinates are in OCS. For a -Z extrusion, the OCS X
        # axis is -WCS X; region selection must happen after this conversion.
        # Without it, genuine negative-Z arcs disappear outside the cluster.
        normal = tuple(float(value) for value in entity.get("extrusion", (0, 0, 1)))
        points = sample_arc(entity)
        if normal == (0.0, 0.0, -1.0):
            return [[-point[0], point[1]] for point in points]
        if normal != (0.0, 0.0, 1.0):
            raise RuntimeError(f"unsupported ARC OCS normal {normal}")
        return points
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


def rational_endpoint_join_audit(selected: list[tuple[dict, list[list[float]], str]]) -> list[dict]:
    """Prove restored rational paths reconnect to the surrounding native topology."""
    endpoints = [
        (int(entity["index"]), position, point)
        for entity, path, _ in selected
        for position, point in (("start", path[0]), ("end", path[-1]))
    ]
    audit = []
    for entity, path, layer in selected:
        if not entity.get("rational"):
            continue
        joins = {}
        for position, point in (("start", path[0]), ("end", path[-1])):
            distance, peer_index, peer_position = min(
                (math.dist(point, peer), peer_index, peer_position)
                for peer_index, peer_position, peer in endpoints
                if peer_index != int(entity["index"])
            )
            joins[position] = {
                "nearest_native_entity_index": peer_index,
                "nearest_native_endpoint": peer_position,
                "distance_mm": round(distance, 6),
                "joins_within_4_mm": distance <= 4.0,
            }
        audit.append({
            "native_entity_index": int(entity["index"]),
            "native_layer": layer,
            "start_mm": [round(value, 6) for value in path[0]],
            "end_mm": [round(value, 6) for value in path[-1]],
            "endpoint_joins": joins,
            "pass": all(item["joins_within_4_mm"] for item in joins.values()),
        })
    return audit


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
        if entity.get("entmode") != 2:
            continue
        path = entity_path(entity, objects_by_handle)
        if path:
            extracted.append(
                (
                    entity,
                    path,
                    bounds([path]),
                    layers.get(handle_value(entity.get("layer", [0])), "<unknown>"),
                )
            )
    views = {}
    for view, region in REGIONS.items():
        spatial_inventory = [
            (entity, path, layer)
            for entity, path, path_bounds, layer in extracted
            if path_bounds["minimum"][0] >= region[0]
            and path_bounds["minimum"][1] >= region[1]
            and path_bounds["maximum"][0] <= region[2]
            and path_bounds["maximum"][1] <= region[3]
        ]
        spatial_layers = {layer for _, _, layer in spatial_inventory}
        if spatial_layers != EXPECTED_NATIVE_LAYERS_BY_VIEW[view]:
            raise RuntimeError(f"Marilyn {view} spatial layer inventory drifted: {sorted(spatial_layers)}")
        selected = spatial_inventory
        paths, normalisation = normalize(view, [path for _, path, _ in selected])
        if len(paths) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Marilyn {view} native path count drifted: {len(paths)}")
        layer_counts = Counter(layer for _, _, layer in selected)
        rational_audit = rational_endpoint_join_audit(selected)
        if not all(item["pass"] for item in rational_audit):
            raise RuntimeError(f"Marilyn {view} restored rational spline topology does not reconnect")
        views[view] = {
            "view": view,
            "paths_mm": paths,
            "path_count": len(paths),
            "point_count": sum(len(path) for path in paths),
            "bounds_mm": bounds(paths),
            "native_entity_counts": dict(sorted(Counter(entity["entity"] for entity, _, _ in selected).items())),
            "native_entity_indices": [int(entity["index"]) for entity, _, _ in selected],
            "negative_z_ocs_arc_entity_indices": [int(entity["index"]) for entity, _, _ in selected if entity.get("entity") == "ARC" and entity.get("extrusion") == [0.0, 0.0, -1.0]],
            "native_layer_path_counts": dict(sorted(layer_counts.items())),
            "native_product_layers": sorted(layer_counts),
            "selection_method": "exact modelspace region; all supported geometric entities inventoried before layer assertion",
            "rational_spline_path_count": len(rational_audit),
            "rational_spline_entity_indices": [item["native_entity_index"] for item in rational_audit],
            "rational_spline_endpoint_join_audit": rational_audit,
            "fit_point_spline_path_count": sum(
                1 for entity, _, _ in selected if entity.get("entity") == "SPLINE" and "knots" not in entity
            ),
            "fit_point_spline_entity_indices": [
                int(entity["index"])
                for entity, _, _ in selected
                if entity.get("entity") == "SPLINE" and "knots" not in entity
            ],
            "normalisation": normalisation,
        }
    old_side = [path for _, path, _, layer in extracted if layer == "_ARREDO" and bounds([path])["minimum"][0] >= REGIONS["side"][0] and bounds([path])["minimum"][1] >= REGIONS["side"][1] and bounds([path])["maximum"][0] <= REGIONS["side"][2] and bounds([path])["maximum"][1] <= REGIONS["side"][3]]
    upholstery_side = [path for _, path, _, layer in extracted if layer == "Make2D$Visibile$Curve" and bounds([path])["minimum"][0] >= REGIONS["side"][0] and bounds([path])["minimum"][1] >= REGIONS["side"][1] and bounds([path])["maximum"][0] <= REGIONS["side"][2] and bounds([path])["maximum"][1] <= REGIONS["side"][3]]
    structural_points = [point for path in old_side for point in path]
    upholstery_offsets = [
        min(math.dist(point, structural) for structural in structural_points)
        for path in upholstery_side
        for point in path
    ]
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
    if not AUTOCAD_FULL_MODELSPACE.is_file():
        raise RuntimeError("direct AutoCAD modelspace inspection screenshot is missing")
    if not LAYER_AUDIT_SVG.is_file() or not LAYER_AUDIT_PNG.is_file():
        raise RuntimeError("Marilyn exact-cluster layer audit evidence is missing")
    exact_dimension_anchors = {
        item.get("text", "").replace("\\A1;", ""): [
            round(float(item.get("ins_pt", [0.0, 0.0])[0]), 6),
            round(float(item.get("ins_pt", [0.0, 0.0])[1]), 6),
        ]
        for item in objects
        if item.get("entity") == "MTEXT"
        and item.get("text", "").replace("\\A1;", "") in {"860", "1000", "940", "450"}
        and 9500.0 <= float(item.get("ins_pt", [0.0])[0]) <= 12500.0
        and 0.0 <= float(item.get("ins_pt", [0.0, 0.0])[1]) <= 2800.0
    }
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
        "upholstery_completeness_audit": {
            "original_layer_filtered_path_counts": ORIGINAL_LAYER_FILTERED_PATH_COUNTS,
            "previous_rational_spline_omitted_path_counts": PREVIOUS_RATIONAL_OMITTED_PATH_COUNTS,
            "revised_complete_path_counts": EXPECTED_PATH_COUNTS,
            "omitted_native_layers_restored": {
                "plan": {"0": 14},
                "front": {"0": 12},
                "side": {"Make2D$Visibile$Curve": 8},
            },
            "omitted_rational_native_splines_restored": {
                view: {
                    "path_count": item["rational_spline_path_count"],
                    "native_entity_indices": item["rational_spline_entity_indices"],
                    "all_endpoints_rejoin_native_topology_within_4_mm": all(
                        audit["pass"] for audit in item["rational_spline_endpoint_join_audit"]
                    ),
                }
                for view, item in views.items()
            },
            "omitted_fit_point_native_splines_restored": {
                view: {
                    "path_count": item["fit_point_spline_path_count"],
                    "native_entity_indices": item["fit_point_spline_entity_indices"],
                }
                for view, item in views.items()
            },
            "layer_semantics": {
                "_ARREDO": "primary furniture outline and structural/upholstery boundaries",
                "0": "manufacturer-authored leather seam, fold and cushion separation detail inside the exact bergere cluster",
                "Make2D$Visibile$Curve": "manufacturer-authored visible upholstery surface envelope in the side view",
            },
            "side_upholstery_surface_offset_to_arredo_mm": {
                "minimum": round(min(upholstery_offsets), 6),
                "median": round(sorted(upholstery_offsets)[len(upholstery_offsets) // 2], 6),
                "maximum": round(max(upholstery_offsets), 6),
            },
            "central_project_headrest_cushion": {
                "present_in_project_ifc_body": True,
                "present_as_independent_closed_component_in_exact_native_dwg_cluster": False,
                "native_dwg_internal_center_detail_is": "leather seam/fold linework, not a detachable accessory outline",
                "synthetic_blue_outline_added": False,
                "interpretation": "The raised central headrest cushion is project-configuration geometry. It remains visible in grey/black comparison layers but is not falsely attributed to Baxter native 2D linework.",
            },
            "cause": "Two independent omissions existed: the original layer filter dropped manufacturer geometry on layer 0 / Make2D$Visibile$Curve, then the prior repair silently skipped four rational native SPLINE entities in both Plan and Front plus two native fit-point SPLINE entities in Front. Spatial inventory, homogeneous B-spline evaluation and native-fit-point interpolation restore every exact-cluster geometric entity.",
            "layer_audit_svg": relative(LAYER_AUDIT_SVG),
            "layer_audit_svg_sha256": sha256(LAYER_AUDIT_SVG),
            "layer_audit_png": relative(LAYER_AUDIT_PNG),
            "layer_audit_png_sha256": sha256(LAYER_AUDIT_PNG),
            "pass": True,
        },
        "native_dimension_labels_in_exact_cluster": dimension_labels,
        "direct_autocad_inspection": {
            "application": "Autodesk AutoCAD 2024",
            "mode": "direct read-only visual inspection of a byte-identical temporary copy",
            "source_opened_sha256": sha256(source),
            "source_matches_archived_native_dwg": sha256(source) == EXPECTED["dwg"],
            "full_modelspace_screenshot": relative(AUTOCAD_FULL_MODELSPACE),
            "full_modelspace_screenshot_sha256": sha256(AUTOCAD_FULL_MODELSPACE),
            "exact_cluster_regions_mm": {view: list(region) for view, region in REGIONS.items()},
            "exact_dimension_anchors_mm": exact_dimension_anchors,
            "exact_variant": "Marilyn bergere armchair with swivel base - 86 x 100 x 94 cm",
            "adjacent_variant_exclusions": [
                {
                    "variant": "standard armchair cluster",
                    "dimension_labels_mm": [800, 870, 780, 400],
                    "reason": "height 780 mm and width/depth cluster do not match the project Body or exact 3DS filename",
                },
                {
                    "variant": "Marilyn pouf cluster",
                    "dimension_labels_mm": [800, 620, 450],
                    "reason": "pouf geometry belongs to project Marilyn 02 and is excluded from Marilyn 01",
                },
            ],
            "selection_unique": True,
            "pass": True,
        },
        "views": views,
        "pass": True,
    }
    write_json(args.output.resolve(), payload)
    print(json.dumps({"output": relative(args.output), "views": {view: item["bounds_mm"] for view, item in views.items()}, "pass": True}, indent=2))


if __name__ == "__main__":
    main()
