#!/usr/bin/env python3
"""Build and verify one contact-preserving fixed-furniture anchor candidate."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import ifcopenshell
import ifcopenshell.api.geometry
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.placement
import numpy as np

from beam_origin_reset_candidate import reset_product_origin
from direction_noise_candidate import all_product_geometry_difference
from furniture_anchor_audit import (
    classify_furniture,
    find_role_decision,
    read_role_decisions,
)
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm


TARGET_FIELDS = {
    "expected_class",
    "global_id",
    "source_anchor_x_mm",
    "source_anchor_y_mm",
    "source_anchor_z_mm",
    "target_x_mm",
    "target_y_mm",
    "target_z_mm",
    "anchor_kind",
    "contact_group",
    "basis",
    "confidence",
    "human_review_required",
    "status",
}
ALLOWANCE_FIELDS = {
    "first_global_id",
    "second_global_id",
    "relation_rule",
    "source_clearance_mm",
    "basis",
    "confidence",
    "human_review_required",
    "status",
}
ALLOWED_ANCHOR_KINDS = {
    "existing_vertex",
    "physical_edge_midpoint",
    "existing_surface_point",
}
ALLOWED_RELATION_RULE = "allow_close_to_contact_without_intersection"
CONTEXT_CLASSES = ("IfcWall", "IfcSlab", "IfcCovering")


def read_targets(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or set(rows[0]) != TARGET_FIELDS:
        raise RuntimeError(f"invalid furniture anchor target schema: {path}")
    result = []
    seen: set[str] = set()
    for row in rows:
        if row["global_id"] in seen:
            raise RuntimeError(f"duplicate furniture anchor target: {row['global_id']}")
        seen.add(row["global_id"])
        if row["expected_class"] != "IfcFurniture":
            raise RuntimeError(f"unsupported furniture target class: {row}")
        if row["anchor_kind"] not in ALLOWED_ANCHOR_KINDS:
            raise RuntimeError(f"unsupported furniture anchor kind: {row}")
        if row["human_review_required"] != "no" or row["status"] != "approved":
            raise RuntimeError(f"unapproved furniture anchor target: {row}")
        source_anchor = np.array(
            [
                float(row["source_anchor_x_mm"]),
                float(row["source_anchor_y_mm"]),
                float(row["source_anchor_z_mm"]),
            ],
            dtype=float,
        )
        target = np.array(
            [
                float(row["target_x_mm"]),
                float(row["target_y_mm"]),
                float(row["target_z_mm"]),
            ],
            dtype=float,
        )
        if np.max(np.abs(target - np.round(target))) > 1e-9:
            raise RuntimeError(f"furniture target is not integer millimetres: {row}")
        confidence = float(row["confidence"])
        if not 0.0 <= confidence <= 1.0:
            raise RuntimeError(f"invalid furniture target confidence: {row}")
        result.append(
            {
                **row,
                "source_anchor_mm": source_anchor,
                "target_mm": target,
                "translation_mm": target - source_anchor,
                "confidence": confidence,
            }
        )
    return result


def pair_key(first: str, second: str) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def read_allowances(path: Path) -> dict[tuple[str, str], dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or set(rows[0]) != ALLOWANCE_FIELDS:
        raise RuntimeError(f"invalid furniture contact allowance schema: {path}")
    result = {}
    for row in rows:
        key = pair_key(row["first_global_id"], row["second_global_id"])
        if key in result:
            raise RuntimeError(f"duplicate furniture contact allowance: {key}")
        if row["relation_rule"] != ALLOWED_RELATION_RULE:
            raise RuntimeError(f"unsupported furniture contact allowance: {row}")
        if row["human_review_required"] != "no" or row["status"] != "approved":
            raise RuntimeError(f"unapproved furniture contact allowance: {row}")
        result[key] = {
            **row,
            "source_clearance_mm": float(row["source_clearance_mm"]),
            "confidence": float(row["confidence"]),
        }
    return result


def fixed_furniture(
    model: ifcopenshell.file,
    role_decisions: list[dict[str, Any]],
) -> list[Any]:
    result = []
    for product in model.by_type("IfcFurniture"):
        assigned_type = ifcopenshell.util.element.get_type(product)
        type_name = getattr(assigned_type, "Name", None)
        decision = find_role_decision(product.GlobalId, type_name, role_decisions)
        classification = classify_furniture(
            getattr(assigned_type, "PredefinedType", None), decision
        )
        if classification["installation_role"] == "fixed_furniture":
            result.append(product)
    return result


def wrapper_global_id(entity: Any) -> str:
    """Read GlobalId from the low-level entity wrappers returned by geom.tree."""
    return str(entity.get_argument(0))


def add_products_to_tree(
    tree: ifcopenshell.geom.tree,
    products: Iterable[Any],
) -> None:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    for product in products:
        tree.add_element(ifcopenshell.geom.create_shape(settings, product))


def clearance_map(
    set_a: list[Any],
    set_b: list[Any],
    clearance_mm: float,
) -> dict[tuple[str, str], float]:
    tree = ifcopenshell.geom.tree()
    unique = {product.GlobalId: product for product in [*set_a, *set_b]}
    add_products_to_tree(tree, unique.values())
    result: dict[tuple[str, str], float] = {}
    for clash in tree.clash_clearance_many(
        set_a, set_b, clearance_mm / 1000.0, True
    ):
        first = wrapper_global_id(clash.a)
        second = wrapper_global_id(clash.b)
        if first == second:
            continue
        key = pair_key(first, second)
        distance = float(clash.distance) * 1000.0
        result[key] = min(result.get(key, math.inf), distance)
    return result


def intersection_pairs(
    set_a: list[Any],
    set_b: list[Any],
    tolerance_mm: float,
) -> set[tuple[str, str]]:
    tree = ifcopenshell.geom.tree()
    unique = {product.GlobalId: product for product in [*set_a, *set_b]}
    add_products_to_tree(tree, unique.values())
    result = set()
    for clash in tree.clash_intersection_many(
        set_a, set_b, tolerance_mm / 1000.0, True
    ):
        first = wrapper_global_id(clash.a)
        second = wrapper_global_id(clash.b)
        if first != second:
            result.add(pair_key(first, second))
    return result


def classify_contact_changes(
    source: dict[tuple[str, str], float],
    candidate: dict[tuple[str, str], float],
    tolerance_mm: float,
    allowed_new_contacts: set[tuple[str, str]],
) -> dict[str, Any]:
    def finite_or_none(value: float) -> float | None:
        return None if math.isinf(value) else value

    regressions = []
    new_contacts = []
    for key, source_distance in source.items():
        candidate_distance = candidate.get(key, math.inf)
        if source_distance <= tolerance_mm and candidate_distance > tolerance_mm:
            regressions.append(
                {
                    "pair": key,
                    "source_clearance_mm": source_distance,
                    "candidate_clearance_mm": finite_or_none(candidate_distance),
                }
            )
    for key, candidate_distance in candidate.items():
        source_distance = source.get(key, math.inf)
        if candidate_distance <= tolerance_mm and source_distance > tolerance_mm:
            new_contacts.append(
                {
                    "pair": key,
                    "source_clearance_mm": finite_or_none(source_distance),
                    "candidate_clearance_mm": candidate_distance,
                    "allowed": key in allowed_new_contacts,
                }
            )
    observed_allowed = {
        tuple(record["pair"]) for record in new_contacts if record["allowed"]
    }
    return {
        "source_contacts": sum(value <= tolerance_mm for value in source.values()),
        "candidate_contacts": sum(
            value <= tolerance_mm for value in candidate.values()
        ),
        "regressions": regressions,
        "new_contacts": new_contacts,
        "unexpected_new_contacts": [
            record for record in new_contacts if not record["allowed"]
        ],
        "missing_allowed_contacts": sorted(allowed_new_contacts - observed_allowed),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--contact-allowances", required=True, type=Path)
    parser.add_argument("--role-decisions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--maximum-axis-shift-mm", type=float, default=0.5)
    parser.add_argument("--clearance-window-mm", type=float, default=1.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= args.maximum_axis_shift_mm <= 1.0:
        raise SystemExit("expected 0 < tolerance <= maximum axis shift <= 1 mm")
    if args.clearance_window_mm < args.tolerance_mm:
        raise SystemExit("clearance window must not be smaller than tolerance")

    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    targets = read_targets(args.targets)
    allowances = read_allowances(args.contact_allowances)
    role_decisions = read_role_decisions(args.role_decisions)
    source_fixed = fixed_furniture(source, role_decisions)
    source_fixed_ids = {product.GlobalId for product in source_fixed}
    target_ids = {row["global_id"] for row in targets}
    if not target_ids.issubset(source_fixed_ids):
        raise RuntimeError(
            f"targets are not confirmed fixed furniture: {sorted(target_ids - source_fixed_ids)}"
        )
    if any(
        max(abs(value) for value in row["translation_mm"])
        > args.maximum_axis_shift_mm
        for row in targets
    ):
        raise RuntimeError("one or more furniture anchor shifts exceed the maximum")

    settings = geometry_settings()
    source_anchor_distances = {}
    for row in targets:
        product = source.by_guid(row["global_id"])
        vertices, faces = world_mesh_mm(settings, product)
        from beam_origin_reset_candidate import point_mesh_distance

        source_anchor_distances[row["global_id"]] = point_mesh_distance(
            row["source_anchor_mm"], vertices, faces
        )
    if max(source_anchor_distances.values()) > args.tolerance_mm:
        raise RuntimeError("one or more source furniture anchors are off geometry")

    results = []
    for row in targets:
        product = candidate.by_guid(row["global_id"])
        matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
            dtype=float,
        )
        matrix[:3, 3] += row["translation_mm"]
        ifcopenshell.api.geometry.edit_object_placement(
            candidate,
            product=product,
            matrix=matrix,
            is_si=False,
            should_transform_children=True,
        )
        reset_result = reset_product_origin(
            candidate,
            row,
            args.tolerance_mm,
            "IfcFurniture",
        )
        results.append(
            {
                **reset_result,
                "source_anchor_mm": row["source_anchor_mm"].tolist(),
                "rigid_geometry_translation_mm": row["translation_mm"].tolist(),
                "contact_group": row["contact_group"],
                "human_review_required": False,
                "status": row["status"],
                "source_anchor_surface_distance_mm": source_anchor_distances[
                    row["global_id"]
                ],
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    candidate_fixed = fixed_furniture(candidate, role_decisions)

    all_geometry = all_product_geometry_difference(
        candidate, source, args.tolerance_mm
    )
    geometry_records = {
        record["global_id"]: record for record in all_geometry["records"]
    }
    non_target_geometry_changes = [
        record
        for record in all_geometry["records"]
        if record["global_id"] not in target_ids and not record["within_tolerance"]
    ]
    target_mesh_count_mismatches = [
        record["global_id"]
        for record in all_geometry["records"]
        if record["global_id"] in target_ids and not record["mesh_counts_equal"]
    ]
    target_shift_mismatches = []
    for row in targets:
        record = geometry_records[row["global_id"]]
        expected = float(np.linalg.norm(row["translation_mm"]))
        actual = record["world_corresponding_vertex_max_delta_mm"]
        if abs(actual - expected) > args.tolerance_mm:
            target_shift_mismatches.append(
                {
                    "global_id": row["global_id"],
                    "expected_mm": expected,
                    "actual_mm": actual,
                }
            )

    source_internal_clearance = clearance_map(
        source_fixed, source_fixed, args.clearance_window_mm
    )
    candidate_internal_clearance = clearance_map(
        candidate_fixed, candidate_fixed, args.clearance_window_mm
    )
    internal_contacts = classify_contact_changes(
        source_internal_clearance,
        candidate_internal_clearance,
        args.tolerance_mm,
        set(allowances),
    )

    source_targets = [source.by_guid(global_id) for global_id in sorted(target_ids)]
    candidate_targets = [
        candidate.by_guid(global_id) for global_id in sorted(target_ids)
    ]
    source_context = [
        product
        for ifc_class in CONTEXT_CLASSES
        for product in source.by_type(ifc_class)
        if getattr(product, "Representation", None)
    ]
    candidate_context = [
        candidate.by_guid(product.GlobalId) for product in source_context
    ]
    context_contacts = classify_contact_changes(
        clearance_map(
            source_targets, source_context, args.clearance_window_mm
        ),
        clearance_map(
            candidate_targets, candidate_context, args.clearance_window_mm
        ),
        args.tolerance_mm,
        set(),
    )

    source_internal_intersections = intersection_pairs(
        source_fixed, source_fixed, args.tolerance_mm
    )
    candidate_internal_intersections = intersection_pairs(
        candidate_fixed, candidate_fixed, args.tolerance_mm
    )
    source_context_intersections = intersection_pairs(
        source_targets, source_context, args.tolerance_mm
    )
    candidate_context_intersections = intersection_pairs(
        candidate_targets, candidate_context, args.tolerance_mm
    )
    new_internal_intersections = sorted(
        candidate_internal_intersections - source_internal_intersections
    )
    new_context_intersections = sorted(
        candidate_context_intersections - source_context_intersections
    )

    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    entity_count_delta = len(list(candidate)) - len(list(source))
    expected_entity_count_delta = sum(
        result["representation_entities_added"] for result in results
    )
    gates = {
        "targets": len(targets),
        "targets_by_contact_group": dict(
            sorted(Counter(row["contact_group"] for row in targets).items())
        ),
        "maximum_axis_shift_mm": max(
            max(abs(value) for value in row["translation_mm"])
            for row in targets
        ),
        "maximum_euclidean_shift_mm": max(
            float(np.linalg.norm(row["translation_mm"])) for row in targets
        ),
        "source_anchors_over_tolerance": sum(
            distance > args.tolerance_mm
            for distance in source_anchor_distances.values()
        ),
        "candidate_anchors_over_tolerance": sum(
            result["anchor_surface_distance_mm"] > args.tolerance_mm
            for result in results
        ),
        "target_placement_mismatches": sum(
            result["placement_residual_mm"] > args.tolerance_mm
            for result in results
        ),
        "target_shift_mismatches": target_shift_mismatches,
        "target_mesh_count_mismatches": target_mesh_count_mismatches,
        "non_target_geometry_changes": [
            record["global_id"] for record in non_target_geometry_changes
        ],
        "internal_contacts": internal_contacts,
        "context_contacts": context_contacts,
        "new_internal_intersections": new_internal_intersections,
        "new_context_intersections": new_context_intersections,
        "schema_equal": source.schema == candidate.schema,
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
        "entity_count_delta": entity_count_delta,
        "expected_entity_count_delta": expected_entity_count_delta,
        "entity_count_delta_matches_expected": entity_count_delta
        == expected_entity_count_delta,
    }
    gates["pass"] = (
        gates["source_anchors_over_tolerance"] == 0
        and gates["candidate_anchors_over_tolerance"] == 0
        and gates["target_placement_mismatches"] == 0
        and not gates["target_shift_mismatches"]
        and not gates["target_mesh_count_mismatches"]
        and not gates["non_target_geometry_changes"]
        and not internal_contacts["regressions"]
        and not internal_contacts["unexpected_new_contacts"]
        and not internal_contacts["missing_allowed_contacts"]
        and not context_contacts["regressions"]
        and not context_contacts["new_contacts"]
        and not new_internal_intersections
        and not new_context_intersections
        and gates["schema_equal"]
        and gates["root_global_ids_equal"]
        and gates["entity_count_delta_matches_expected"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-fixed-furniture-contact-preserving-anchor-candidate",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "candidate": {
            "path": str(args.output.resolve()),
            "sha256": sha256(args.output),
            "schema": candidate.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "maximum_axis_shift_mm": args.maximum_axis_shift_mm,
        "clearance_window_mm": args.clearance_window_mm,
        "targets_path": str(args.targets.resolve()),
        "contact_allowances_path": str(args.contact_allowances.resolve()),
        "allowances": list(allowances.values()),
        "results": results,
        "all_product_geometry": all_geometry,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
