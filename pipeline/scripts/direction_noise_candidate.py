#!/usr/bin/env python3
"""Build a candidate that removes near-cardinal IfcDirection ratio noise.

The source IFC is never modified. Only three-dimensional IfcDirection values
used by exactly one IfcAxis2Placement3D Axis or RefDirection are eligible.
Direction-ratio tolerance is dimensionless; world-geometry tolerance is in mm.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import ifcopenshell

from geometry_alignment_audit import (
    geometry_settings,
    placement_deltas,
    placement_matrix,
    sha256,
    world_mesh_mm,
)


CARDINAL_DIRECTIONS_3D: tuple[tuple[float, float, float], ...] = (
    (1.0, 0.0, 0.0),
    (-1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, -1.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.0, 0.0, -1.0),
)


def finite_positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise argparse.ArgumentTypeError("must be a finite number greater than zero")
    return number


def vector_distance(first: Sequence[float], second: Sequence[float]) -> float:
    if len(first) != len(second):
        raise ValueError("vectors must have equal dimensions")
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))


def nearest_cardinal_3d(
    ratios: Sequence[float],
) -> tuple[tuple[float, float, float], float] | None:
    """Return the nearest cardinal to the normalized 3D direction.

    IfcDirection ratios are dimensionless and need not have unit magnitude, so
    classification uses the normalized vector. Two-dimensional and invalid
    directions are outside this candidate's write boundary.
    """

    if len(ratios) != 3:
        return None
    values = tuple(float(value) for value in ratios)
    if not all(math.isfinite(value) for value in values):
        return None
    magnitude = math.sqrt(sum(value * value for value in values))
    if magnitude == 0.0:
        return None
    normalized = tuple(value / magnitude for value in values)
    target = min(
        CARDINAL_DIRECTIONS_3D,
        key=lambda cardinal: vector_distance(normalized, cardinal),
    )
    return target, vector_distance(normalized, target)


def placement_usage(
    model: ifcopenshell.file,
    direction: ifcopenshell.entity_instance,
) -> dict[str, Any]:
    inverses = list(model.get_inverse(direction))
    roles: list[dict[str, Any]] = []
    for inverse in inverses:
        if not inverse.is_a("IfcAxis2Placement3D"):
            continue
        if inverse.Axis == direction:
            roles.append(
                {
                    "placement_id": inverse.id(),
                    "placement_class": inverse.is_a(),
                    "role": "Axis",
                }
            )
        if inverse.RefDirection == direction:
            roles.append(
                {
                    "placement_id": inverse.id(),
                    "placement_class": inverse.is_a(),
                    "role": "RefDirection",
                }
            )
    return {
        "inverse_count": len(inverses),
        "inverse_classes": sorted(inverse.is_a() for inverse in inverses),
        "roles": roles,
        "eligible": len(inverses) == 1 and len(roles) == 1,
    }


def direction_inventory(
    model: ifcopenshell.file,
    ratio_tolerance: float,
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    for direction in model.by_type("IfcDirection"):
        ratios = tuple(float(value) for value in direction.DirectionRatios)
        nearest = nearest_cardinal_3d(ratios)
        if nearest is None:
            category = "outside_3d_write_boundary"
            record = {
                "direction_id": direction.id(),
                "ratios": list(ratios),
                "category": category,
            }
        else:
            target, ratio_distance = nearest
            usage = placement_usage(model, direction)
            if ratios == target:
                category = "exact_cardinal"
            elif ratio_distance <= ratio_tolerance and usage["eligible"]:
                category = "near_cardinal_eligible"
            elif ratio_distance <= ratio_tolerance:
                category = "near_cardinal_ineligible"
            else:
                category = "intentional_noncardinal"
            record = {
                "direction_id": direction.id(),
                "ratios": list(ratios),
                "nearest_cardinal": list(target),
                "normalized_ratio_distance": ratio_distance,
                "category": category,
                "usage": usage,
            }
        counts[category] += 1
        records.append(record)
    return {"counts": dict(counts), "records": records}


def apply_direction_cleanup(
    model: ifcopenshell.file,
    ratio_tolerance: float,
    excluded_direction_ids: set[int] | None = None,
) -> list[dict[str, Any]]:
    excluded_direction_ids = excluded_direction_ids or set()
    inventory = direction_inventory(model, ratio_tolerance)
    changes: list[dict[str, Any]] = []
    for record in inventory["records"]:
        if record["category"] != "near_cardinal_eligible":
            continue
        if record["direction_id"] in excluded_direction_ids:
            continue
        direction = model.by_id(record["direction_id"])
        before = tuple(float(value) for value in direction.DirectionRatios)
        after = tuple(float(value) for value in record["nearest_cardinal"])
        direction.DirectionRatios = after
        changes.append(
            {
                "direction_id": direction.id(),
                "before": list(before),
                "after": list(after),
                "normalized_ratio_distance": record[
                    "normalized_ratio_distance"
                ],
                "usage": record["usage"],
            }
        )
    return changes


def entity_changes(
    source: ifcopenshell.file,
    candidate: ifcopenshell.file,
) -> list[int]:
    changed: list[int] = []
    for source_entity in source:
        try:
            candidate_entity = candidate.by_id(source_entity.id())
        except RuntimeError:
            candidate_entity = None
        if candidate_entity is None or str(source_entity) != str(candidate_entity):
            changed.append(source_entity.id())
    return changed


def maximum_corresponding_vertex_delta_mm(
    first: Sequence[Sequence[float]],
    second: Sequence[Sequence[float]],
) -> float:
    if len(first) != len(second):
        return math.inf
    return max(
        (vector_distance(a, b) for a, b in zip(first, second)),
        default=0.0,
    )


def directed_spatial_hash_distance_mm(
    first: Sequence[Sequence[float]],
    second: Sequence[Sequence[float]],
    tolerance_mm: float,
) -> float:
    """Return a directed nearest-vertex distance within a proven tolerance.

    This avoids an O(n²) scan for complex meshes. A missing neighbour in the
    surrounding 27 tolerance-sized cells proves that the distance exceeds the
    requested gate, so the function returns infinity.
    """

    if not first or not second:
        return 0.0 if not first and not second else math.inf
    buckets: dict[tuple[int, int, int], list[Sequence[float]]] = {}
    for point in second:
        key = tuple(math.floor(value / tolerance_mm) for value in point)
        buckets.setdefault(key, []).append(point)
    maximum = 0.0
    for point in first:
        base = tuple(math.floor(value / tolerance_mm) for value in point)
        nearest = math.inf
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other in buckets.get(
                        (base[0] + dx, base[1] + dy, base[2] + dz), ()
                    ):
                        nearest = min(nearest, vector_distance(point, other))
        if nearest > tolerance_mm:
            return math.inf
        maximum = max(maximum, nearest)
    return maximum


def vertex_hausdorff_within_tolerance_mm(
    first: Sequence[Sequence[float]],
    second: Sequence[Sequence[float]],
    tolerance_mm: float,
) -> float:
    return max(
        directed_spatial_hash_distance_mm(first, second, tolerance_mm),
        directed_spatial_hash_distance_mm(second, first, tolerance_mm),
    )


def all_product_geometry_difference(
    candidate: ifcopenshell.file,
    source: ifcopenshell.file,
    tolerance_mm: float,
) -> dict[str, Any]:
    """Compare every represented IfcProduct without treating index order as topology.

    STEP reload can reorder tessellation indices even when the IFC topology is
    untouched. Product geometry therefore uses a tolerance-bounded spatial
    hash, while the separate entity gate proves that representation entities
    other than the selected IfcDirection values did not change.
    """

    source_products = {
        product.GlobalId: product
        for product in source.by_type("IfcProduct")
        if product.Representation is not None
    }
    candidate_products = {
        product.GlobalId: product
        for product in candidate.by_type("IfcProduct")
        if product.Representation is not None
    }
    source_settings = geometry_settings()
    candidate_settings = geometry_settings()
    records: list[dict[str, Any]] = []
    for global_id in sorted(set(source_products) | set(candidate_products)):
        source_product = source_products.get(global_id)
        candidate_product = candidate_products.get(global_id)
        if source_product is None or candidate_product is None:
            records.append(
                {
                    "global_id": global_id,
                    "status": "added" if source_product is None else "removed",
                    "within_tolerance": False,
                    "mesh_counts_equal": False,
                    "face_index_order_equal": False,
                    "both_unrenderable": False,
                    "asymmetric_shape_failure": True,
                }
            )
            continue
        source_error = None
        candidate_error = None
        try:
            source_vertices, source_faces = world_mesh_mm(
                source_settings, source_product
            )
        except Exception as error:  # pragma: no cover - model-specific evidence
            source_vertices = []
            source_faces = []
            source_error = str(error)
        try:
            candidate_vertices, candidate_faces = world_mesh_mm(
                candidate_settings, candidate_product
            )
        except Exception as error:  # pragma: no cover - model-specific evidence
            candidate_vertices = []
            candidate_faces = []
            candidate_error = str(error)
        both_unrenderable = source_error is not None and candidate_error is not None
        asymmetric_shape_failure = (source_error is None) != (candidate_error is None)
        if both_unrenderable:
            mesh_counts_equal = None
            face_index_order_equal = None
            geometry_delta = 0.0
        elif asymmetric_shape_failure:
            mesh_counts_equal = False
            face_index_order_equal = False
            geometry_delta = math.inf
        else:
            mesh_counts_equal = (
                len(source_vertices) == len(candidate_vertices)
                and len(source_faces) == len(candidate_faces)
            )
            face_index_order_equal = source_faces == candidate_faces
            if not mesh_counts_equal:
                geometry_delta = math.inf
            elif face_index_order_equal:
                geometry_delta = maximum_corresponding_vertex_delta_mm(
                    source_vertices, candidate_vertices
                )
            else:
                geometry_delta = vertex_hausdorff_within_tolerance_mm(
                    source_vertices, candidate_vertices, tolerance_mm
                )
        translation_delta, rotation_delta = placement_deltas(
            placement_matrix(source_product), placement_matrix(candidate_product)
        )
        within_tolerance = (
            both_unrenderable
            or (
                not asymmetric_shape_failure
                and bool(mesh_counts_equal)
                and geometry_delta <= tolerance_mm
            )
        )
        if both_unrenderable:
            status = "unrenderable_in_both_models"
        elif not within_tolerance:
            status = "world_geometry_changed_over_tolerance"
        elif geometry_delta > 0.0 or (rotation_delta or 0.0) > 0.0:
            status = "world_geometry_changed_within_tolerance"
        else:
            status = "unchanged"
        records.append(
            {
                "global_id": global_id,
                "ifc_class": candidate_product.is_a(),
                "name": candidate_product.Name,
                "status": status,
                "within_tolerance": within_tolerance,
                "world_corresponding_vertex_max_delta_mm": geometry_delta,
                "placement_translation_delta_mm": translation_delta,
                "placement_rotation_matrix_max_delta": rotation_delta,
                "source_vertices": len(source_vertices),
                "candidate_vertices": len(candidate_vertices),
                "source_triangles": len(source_faces),
                "candidate_triangles": len(candidate_faces),
                "mesh_counts_equal": mesh_counts_equal,
                "face_index_order_equal": face_index_order_equal,
                "source_shape_error": source_error,
                "candidate_shape_error": candidate_error,
                "both_unrenderable": both_unrenderable,
                "asymmetric_shape_failure": asymmetric_shape_failure,
            }
        )
    status_counts = Counter(record["status"] for record in records)
    finite_deltas = [
        record["world_corresponding_vertex_max_delta_mm"]
        for record in records
        if math.isfinite(record.get("world_corresponding_vertex_max_delta_mm", math.inf))
    ]
    return {
        "tolerance_mm": tolerance_mm,
        "total": len(records),
        "within_tolerance": sum(record["within_tolerance"] for record in records),
        "over_tolerance": sum(not record["within_tolerance"] for record in records),
        "mesh_count_mismatches": sum(
            record["mesh_counts_equal"] is False for record in records
        ),
        "face_index_order_changes": sum(
            record["face_index_order_equal"] is False
            and record["mesh_counts_equal"] is True
            for record in records
        ),
        "unrenderable_in_both_models": sum(
            record["both_unrenderable"] for record in records
        ),
        "asymmetric_shape_failures": sum(
            record["asymmetric_shape_failure"] for record in records
        ),
        "maximum_world_vertex_delta_mm": max(finite_deltas, default=0.0),
        "status_counts": dict(status_counts),
        "records": records,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument(
        "--direction-ratio-tolerance",
        type=finite_positive_float,
        default=1e-6,
        help="Dimensionless distance from a normalized direction to a cardinal axis.",
    )
    parser.add_argument(
        "--geometry-tolerance-mm",
        type=finite_positive_float,
        default=0.1,
        help="Maximum permitted world-geometry movement in millimetres.",
    )
    parser.add_argument(
        "--exclude-direction-id",
        action="append",
        type=int,
        default=[],
        help=(
            "IfcDirection STEP id to preserve as a verified topology-sensitive "
            "exception; repeat for multiple ids."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.direction_ratio_tolerance > 1e-3:
        raise SystemExit("--direction-ratio-tolerance must be at most 0.001")
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    source_inventory = direction_inventory(source, args.direction_ratio_tolerance)
    excluded_direction_ids = set(args.exclude_direction_id)
    eligible_source_records = {
        record["direction_id"]: record
        for record in source_inventory["records"]
        if record["category"] == "near_cardinal_eligible"
    }
    unknown_exclusions = sorted(excluded_direction_ids - set(eligible_source_records))
    if unknown_exclusions:
        raise SystemExit(
            "Excluded direction ids are not eligible near-cardinal directions: "
            + ", ".join(str(value) for value in unknown_exclusions)
        )
    changes = apply_direction_cleanup(
        candidate,
        args.direction_ratio_tolerance,
        excluded_direction_ids=excluded_direction_ids,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    candidate_inventory = direction_inventory(
        candidate, args.direction_ratio_tolerance
    )
    geometry_difference = all_product_geometry_difference(
        candidate, source, args.geometry_tolerance_mm
    )
    changed_entity_ids = entity_changes(source, candidate)
    changed_direction_ids = sorted(change["direction_id"] for change in changes)
    candidate_near_cardinal_ids = sorted(
        record["direction_id"]
        for record in candidate_inventory["records"]
        if record["category"] == "near_cardinal_eligible"
    )
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "changed_directions": len(changes),
        "source_near_cardinal_eligible": source_inventory["counts"].get(
            "near_cardinal_eligible", 0
        ),
        "candidate_near_cardinal_eligible": candidate_inventory["counts"].get(
            "near_cardinal_eligible", 0
        ),
        "excluded_direction_ids_match_candidate_residual": candidate_near_cardinal_ids
        == sorted(excluded_direction_ids),
        "intentional_noncardinal_count_equal": source_inventory["counts"].get(
            "intentional_noncardinal", 0
        )
        == candidate_inventory["counts"].get("intentional_noncardinal", 0),
        "only_selected_direction_entities_changed": set(changed_entity_ids)
        == set(changed_direction_ids),
        "all_product_geometry_over_tolerance": geometry_difference[
            "over_tolerance"
        ],
        "mesh_count_mismatches": geometry_difference["mesh_count_mismatches"],
        "asymmetric_shape_failures": geometry_difference[
            "asymmetric_shape_failures"
        ],
        "schema_equal": source.schema == candidate.schema,
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
        "entity_count_equal": len(list(source)) == len(list(candidate)),
    }
    gates["pass"] = (
        gates["changed_directions"] > 0
        and gates["changed_directions"]
        + gates["candidate_near_cardinal_eligible"]
        == gates["source_near_cardinal_eligible"]
        and gates["candidate_near_cardinal_eligible"]
        == len(excluded_direction_ids)
        and gates["excluded_direction_ids_match_candidate_residual"]
        and gates["intentional_noncardinal_count_equal"]
        and gates["only_selected_direction_entities_changed"]
        and gates["all_product_geometry_over_tolerance"] == 0
        and gates["mesh_count_mismatches"] == 0
        and gates["asymmetric_shape_failures"] == 0
        and gates["schema_equal"]
        and gates["root_global_ids_equal"]
        and gates["entity_count_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-near-cardinal-direction-candidate",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "candidate": {
            "path": str(args.output),
            "sha256": sha256(args.output),
            "schema": candidate.schema,
        },
        "direction_ratio_tolerance": args.direction_ratio_tolerance,
        "world_geometry_tolerance_mm": args.geometry_tolerance_mm,
        "source_inventory_counts": source_inventory["counts"],
        "candidate_inventory_counts": candidate_inventory["counts"],
        "excluded_directions": [
            eligible_source_records[direction_id]
            for direction_id in sorted(excluded_direction_ids)
        ],
        "changes": changes,
        "changed_entity_ids": changed_entity_ids,
        "all_product_geometry_difference": geometry_difference,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "report": str(args.report),
                "source_inventory_counts": source_inventory["counts"],
                "candidate_inventory_counts": candidate_inventory["counts"],
                "maximum_world_vertex_delta_mm": geometry_difference[
                    "maximum_world_vertex_delta_mm"
                ],
                **gates,
            },
            ensure_ascii=False,
        )
    )
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
