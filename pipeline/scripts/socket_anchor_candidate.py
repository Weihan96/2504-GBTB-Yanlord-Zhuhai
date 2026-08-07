#!/usr/bin/env python3
"""Build and validate one SOC01/SOC04/SOCF04 socket-anchor candidate."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import numpy as np

from direction_noise_candidate import all_product_geometry_difference
from furniture_group_origin_candidate import (
    classify_contact_changes,
    clearance_map,
    intersection_pairs,
)
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm
from light_fixture_anchor_candidate import (
    apply_targets,
    origin_mm,
    renderable_context,
)


SOCKET_TYPES = {"SOC01", "SOC04", "SOCF04"}
CONTROLLED_EXCEPTION_IDS = {
    "0laejMoxn8Lu_X3FZaCXmi",
    "2OOjqQDMHDjRcXQCniWXnp",
    "27MTenki57DQsfMryX_1U0",
    "3KXtmVvejA78j_iAS$kydj",
}
ISLAND_GLOBAL_ID = "1OU00zHknFfQOKVh$gmE1r"
ALLOWED_NEW_CONTEXT_CONTACTS = {
    tuple(sorted((ISLAND_GLOBAL_ID, socket_id)))
    for socket_id in {
        "3Cz0fQK8f54AFcm1SrjsGr",
        "22tO8BAFTDWPLxqXvDzGct",
        "23EqZuJ7zEw8JjKNbjSYrg",
    }
}
CONTEXT_CLASSES = (
    "IfcWall",
    "IfcCovering",
    "IfcSlab",
    "IfcFurniture",
    "IfcBuildingElementProxy",
)


def socket_evidence(
    type_name: str | None,
    origin: np.ndarray,
    bbox_min: np.ndarray,
    bbox_max: np.ndarray,
    evidence_tolerance_mm: float = 0.01,
) -> dict[str, Any]:
    size = bbox_max - bbox_min
    centre_offset = (bbox_min + bbox_max) / 2.0 - origin
    if type_name in {"SOC01", "SOC04"}:
        expected_size = np.array([10.0, 100.0, 100.0])
        mounting_face_centre = (
            abs(abs(float(centre_offset[0])) - 5.0) <= evidence_tolerance_mm
            and max(abs(float(value)) for value in centre_offset[1:])
            <= evidence_tolerance_mm
        )
    elif type_name == "SOCF04":
        expected_size = np.array([100.0, 10.0, 99.952778])
        mounting_face_centre = (
            abs(abs(float(centre_offset[1])) - 5.0) <= evidence_tolerance_mm
            and abs(float(centre_offset[0])) <= evidence_tolerance_mm
            and abs(float(centre_offset[2])) <= evidence_tolerance_mm
        )
    else:
        expected_size = np.zeros(3)
        mounting_face_centre = False
    checks = {
        "supported_type": type_name in SOCKET_TYPES,
        "uniform_bbox_size": bool(
            np.max(np.abs(size - expected_size)) <= evidence_tolerance_mm
        ),
        "origin_is_mounting_face_centre": bool(mounting_face_centre),
    }
    target = np.round(origin)
    return {
        "checks": checks,
        "pass": all(checks.values()),
        "origin_mm": origin.tolist(),
        "target_mm": target.tolist(),
        "translation_mm": (target - origin).tolist(),
        "bbox_min_mm": bbox_min.tolist(),
        "bbox_max_mm": bbox_max.tolist(),
        "bbox_size_mm": size.tolist(),
        "bbox_centre_minus_origin_mm": centre_offset.tolist(),
    }


def select_targets(
    model: ifcopenshell.file,
    tolerance_mm: float,
    maximum_axis_shift_mm: float,
    evidence_tolerance_mm: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    settings = geometry_settings()
    inventory = []
    targets = []
    for appliance in model.by_type("IfcElectricAppliance"):
        assigned_type = ifcopenshell.util.element.get_type(appliance)
        type_name = getattr(assigned_type, "Name", None)
        if type_name not in SOCKET_TYPES:
            continue
        vertices, _ = world_mesh_mm(settings, appliance)
        points = np.array(vertices, dtype=float)
        evidence = socket_evidence(
            type_name,
            origin_mm(appliance),
            points.min(axis=0),
            points.max(axis=0),
            evidence_tolerance_mm,
        )
        record = {
            "global_id": appliance.GlobalId,
            "name": appliance.Name,
            "type_name": type_name,
            "predefined_type": getattr(assigned_type, "PredefinedType", None),
            **evidence,
        }
        if appliance.GlobalId in CONTROLLED_EXCEPTION_IDS:
            record.update(
                {
                    "normalization_disposition": (
                        "controlled_exception_preserve_mounting_contact"
                    ),
                    "automatic_write_allowed": False,
                    "human_review_required": False,
                }
            )
        inventory.append(record)
        maximum_residual = max(
            abs(value) for value in record["translation_mm"]
        )
        if maximum_residual <= tolerance_mm:
            continue
        if appliance.GlobalId in CONTROLLED_EXCEPTION_IDS:
            continue
        if not record["pass"]:
            raise RuntimeError(
                f"socket lacks uniform mounting-face evidence: {appliance.GlobalId}"
            )
        if maximum_residual > maximum_axis_shift_mm:
            raise RuntimeError(
                f"socket shift exceeds maximum: {appliance.GlobalId}"
            )
        targets.append(
            {
                **record,
                "anchor_kind": "existing_mounting_face_centre",
                "basis": (
                    f"assigned IfcElectricApplianceType {type_name}; uniform "
                    "10x100x100 mm wall socket body (SOCF04 rotated in plan); "
                    "ObjectPlacement is the exact mounting-face centre"
                ),
                "confidence": 1.0,
                "human_review_required": False,
                "status": "approved_by_recorded_rule",
            }
        )
    return inventory, targets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--maximum-axis-shift-mm", type=float, default=0.5)
    parser.add_argument("--clearance-window-mm", type=float, default=1.0)
    parser.add_argument("--evidence-tolerance-mm", type=float, default=0.01)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.evidence_tolerance_mm <= args.tolerance_mm:
        raise SystemExit("expected 0 < evidence tolerance <= QA tolerance")
    if not 0.0 < args.tolerance_mm <= args.maximum_axis_shift_mm <= 0.5:
        raise SystemExit("expected 0 < tolerance <= maximum axis shift <= 0.5 mm")
    if args.clearance_window_mm < args.tolerance_mm:
        raise SystemExit("clearance window must not be smaller than tolerance")

    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    inventory, targets = select_targets(
        source,
        args.tolerance_mm,
        args.maximum_axis_shift_mm,
        args.evidence_tolerance_mm,
    )
    if len(inventory) != 11 or len(targets) != 7:
        raise RuntimeError(
            f"expected 11 socket instances and 7 targets, got {len(inventory)} / {len(targets)}"
        )
    inventory_exception_ids = {
        row["global_id"]
        for row in inventory
        if row.get("normalization_disposition")
        == "controlled_exception_preserve_mounting_contact"
    }
    if inventory_exception_ids != CONTROLLED_EXCEPTION_IDS:
        raise RuntimeError(
            "socket controlled-exception inventory does not match the exact rule"
        )
    candidate = ifcopenshell.open(source_path)
    results = apply_targets(candidate, targets)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)

    target_ids = {row["global_id"] for row in targets}
    all_geometry = all_product_geometry_difference(
        candidate, source, args.tolerance_mm
    )
    geometry_records = {
        record["global_id"]: record for record in all_geometry["records"]
    }
    target_shift_mismatches = []
    for row in targets:
        record = geometry_records[row["global_id"]]
        expected = float(
            np.linalg.norm(np.array(row["translation_mm"], dtype=float))
        )
        actual = record["world_corresponding_vertex_max_delta_mm"]
        if not record["mesh_counts_equal"] or abs(actual - expected) > args.tolerance_mm:
            target_shift_mismatches.append(
                {
                    "global_id": row["global_id"],
                    "expected_mm": expected,
                    "actual_mm": actual,
                    "mesh_counts_equal": record["mesh_counts_equal"],
                }
            )
    non_target_geometry_changes = [
        record["global_id"]
        for record in all_geometry["records"]
        if record["global_id"] not in target_ids and not record["within_tolerance"]
    ]

    source_appliances = source.by_type("IfcElectricAppliance")
    candidate_appliances = candidate.by_type("IfcElectricAppliance")
    source_targets = [source.by_guid(global_id) for global_id in sorted(target_ids)]
    candidate_targets = [
        candidate.by_guid(global_id) for global_id in sorted(target_ids)
    ]
    source_context, skipped_source_context = renderable_context(
        source, CONTEXT_CLASSES
    )
    candidate_context = [
        candidate.by_guid(product.GlobalId) for product in source_context
    ]
    candidate_renderable, skipped_candidate_context = renderable_context(
        candidate, CONTEXT_CLASSES
    )
    candidate_context_ids = {product.GlobalId for product in candidate_renderable}
    source_context_ids = {product.GlobalId for product in source_context}

    internal_contacts = classify_contact_changes(
        clearance_map(
            source_appliances, source_appliances, args.clearance_window_mm
        ),
        clearance_map(
            candidate_appliances,
            candidate_appliances,
            args.clearance_window_mm,
        ),
        args.tolerance_mm,
        set(),
    )
    context_contacts = classify_contact_changes(
        clearance_map(source_targets, source_context, args.clearance_window_mm),
        clearance_map(
            candidate_targets, candidate_context, args.clearance_window_mm
        ),
        args.tolerance_mm,
        ALLOWED_NEW_CONTEXT_CONTACTS,
    )
    source_internal_intersections = intersection_pairs(
        source_appliances, source_appliances, args.tolerance_mm
    )
    candidate_internal_intersections = intersection_pairs(
        candidate_appliances, candidate_appliances, args.tolerance_mm
    )
    source_context_intersections = intersection_pairs(
        source_targets, source_context, args.tolerance_mm
    )
    candidate_context_intersections = intersection_pairs(
        candidate_targets, candidate_context, args.tolerance_mm
    )

    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {
        root.GlobalId for root in candidate.by_type("IfcRoot")
    }
    gates = {
        "socket_instances": len(inventory),
        "targets": len(targets),
        "controlled_exceptions": len(inventory_exception_ids),
        "controlled_exception_ids": sorted(inventory_exception_ids),
        "allowed_new_context_contact_pairs": sorted(
            ALLOWED_NEW_CONTEXT_CONTACTS
        ),
        "maximum_axis_shift_mm": max(
            max(abs(value) for value in row["translation_mm"])
            for row in targets
        ),
        "maximum_euclidean_shift_mm": max(
            float(np.linalg.norm(np.array(row["translation_mm"], dtype=float)))
            for row in targets
        ),
        "target_placement_mismatches": sum(
            row["placement_residual_mm"] > args.tolerance_mm
            for row in results
        ),
        "target_shift_mismatches": target_shift_mismatches,
        "non_target_geometry_changes": non_target_geometry_changes,
        "internal_contacts": internal_contacts,
        "context_contacts": context_contacts,
        "new_internal_intersections": sorted(
            candidate_internal_intersections - source_internal_intersections
        ),
        "lost_internal_intersections": sorted(
            source_internal_intersections - candidate_internal_intersections
        ),
        "new_context_intersections": sorted(
            candidate_context_intersections - source_context_intersections
        ),
        "lost_context_intersections": sorted(
            source_context_intersections - candidate_context_intersections
        ),
        "renderable_context_ids_equal": source_context_ids
        == candidate_context_ids,
        "skipped_unrenderable_source_context": sorted(
            record["global_id"] for record in skipped_source_context
        ),
        "skipped_unrenderable_candidate_context": sorted(
            record["global_id"] for record in skipped_candidate_context
        ),
        "schema_equal": source.schema == candidate.schema,
        "entity_count_equal": len(list(source)) == len(list(candidate)),
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
    }
    gates["pass"] = (
        gates["target_placement_mismatches"] == 0
        and not gates["target_shift_mismatches"]
        and not gates["non_target_geometry_changes"]
        and not internal_contacts["regressions"]
        and not internal_contacts["unexpected_new_contacts"]
        and not context_contacts["regressions"]
        and not context_contacts["unexpected_new_contacts"]
        and not context_contacts["missing_allowed_contacts"]
        and not gates["new_internal_intersections"]
        and not gates["lost_internal_intersections"]
        and not gates["new_context_intersections"]
        and not gates["lost_context_intersections"]
        and gates["renderable_context_ids_equal"]
        and gates["schema_equal"]
        and gates["entity_count_equal"]
        and gates["root_global_ids_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-socket-mounting-face-anchor-candidate",
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
        "evidence_tolerance_mm": args.evidence_tolerance_mm,
        "inventory": inventory,
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
