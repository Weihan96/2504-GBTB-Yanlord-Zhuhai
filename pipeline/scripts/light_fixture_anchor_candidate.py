#!/usr/bin/env python3
"""Build and validate one RA.LP light-fixture installation-anchor candidate."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api.geometry
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.placement
import numpy as np

from direction_noise_candidate import all_product_geometry_difference
from furniture_group_origin_candidate import (
    classify_contact_changes,
    clearance_map,
    intersection_pairs,
)
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm


EXPECTED_TYPE_NAME = "RA.LP"
EXPECTED_PREDEFINED_TYPE = "DIRECTIONSOURCE"
EXPECTED_BBOX_MM = np.array([55.0, 55.0, 89.0], dtype=float)
CONTEXT_CLASSES = (
    "IfcCovering",
    "IfcSlab",
    "IfcBuildingElementProxy",
    "IfcBeam",
    "IfcWall",
)


def origin_mm(product: Any) -> np.ndarray:
    matrix = ifcopenshell.util.placement.get_local_placement(
        product.ObjectPlacement
    )
    return np.array(matrix[:3, 3], dtype=float)


def light_fixture_evidence(
    type_name: str | None,
    predefined_type: str | None,
    origin: np.ndarray,
    bbox_min: np.ndarray,
    bbox_max: np.ndarray,
    evidence_tolerance_mm: float = 0.01,
) -> dict[str, Any]:
    size = bbox_max - bbox_min
    bbox_centre = (bbox_min + bbox_max) / 2.0
    residual = np.round(origin) - origin
    target = np.round(origin)
    checks = {
        "type_name": type_name == EXPECTED_TYPE_NAME,
        "predefined_type": predefined_type == EXPECTED_PREDEFINED_TYPE,
        "bbox_size": bool(
            np.max(np.abs(size - EXPECTED_BBOX_MM)) <= evidence_tolerance_mm
        ),
        "xy_centre": bool(
            np.max(np.abs(bbox_centre[:2] - origin[:2]))
            <= evidence_tolerance_mm
        ),
        "integer_z": abs(float(target[2] - origin[2])) <= evidence_tolerance_mm,
    }
    return {
        "checks": checks,
        "pass": all(checks.values()),
        "origin_mm": origin.tolist(),
        "target_mm": target.tolist(),
        "translation_mm": residual.tolist(),
        "bbox_min_mm": bbox_min.tolist(),
        "bbox_max_mm": bbox_max.tolist(),
        "bbox_size_mm": size.tolist(),
        "bbox_centre_minus_origin_mm": (bbox_centre - origin).tolist(),
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
    for fixture in model.by_type("IfcLightFixture"):
        assigned_type = ifcopenshell.util.element.get_type(fixture)
        type_name = getattr(assigned_type, "Name", None)
        predefined_type = getattr(assigned_type, "PredefinedType", None)
        vertices, _ = world_mesh_mm(settings, fixture)
        points = np.array(vertices, dtype=float)
        evidence = light_fixture_evidence(
            type_name,
            predefined_type,
            origin_mm(fixture),
            points.min(axis=0),
            points.max(axis=0),
            evidence_tolerance_mm,
        )
        record = {
            "global_id": fixture.GlobalId,
            "name": fixture.Name,
            "type_name": type_name,
            "predefined_type": predefined_type,
            **evidence,
        }
        inventory.append(record)
        maximum_residual = max(
            abs(value) for value in record["translation_mm"]
        )
        if maximum_residual <= tolerance_mm:
            continue
        if not record["pass"]:
            raise RuntimeError(
                f"light fixture lacks uniform installation evidence: {fixture.GlobalId}"
            )
        if maximum_residual > maximum_axis_shift_mm:
            raise RuntimeError(
                f"light fixture shift exceeds maximum: {fixture.GlobalId}"
            )
        targets.append(
            {
                **record,
                "anchor_kind": "existing_installation_centre",
                "basis": (
                    "assigned IfcLightFixtureType RA.LP / DIRECTIONSOURCE; "
                    "uniform 55x55x89 mm body; ObjectPlacement is the exact XY "
                    "centre and its Z installation datum is already integer"
                ),
                "confidence": 1.0,
                "human_review_required": False,
                "status": "approved_by_recorded_rule",
            }
        )
    return inventory, targets


def apply_targets(
    model: ifcopenshell.file, targets: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    results = []
    for row in targets:
        fixture = model.by_guid(row["global_id"])
        matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(
                fixture.ObjectPlacement
            ),
            dtype=float,
        )
        matrix[:3, 3] = np.array(row["target_mm"], dtype=float)
        ifcopenshell.api.geometry.edit_object_placement(
            model,
            product=fixture,
            matrix=matrix,
            is_si=False,
            should_transform_children=True,
        )
        result = origin_mm(fixture)
        results.append(
            {
                **row,
                "result_mm": result.tolist(),
                "placement_residual_mm": float(
                    np.max(
                        np.abs(result - np.array(row["target_mm"], dtype=float))
                    )
                ),
                "translation_shift_mm": float(
                    np.linalg.norm(np.array(row["translation_mm"], dtype=float))
                ),
            }
        )
    return results


def renderable_context(
    model: ifcopenshell.file,
    context_classes: tuple[str, ...] = CONTEXT_CLASSES,
) -> tuple[list[Any], list[dict[str, Any]]]:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    products = []
    skipped = []
    for ifc_class in context_classes:
        for product in model.by_type(ifc_class):
            if not getattr(product, "Representation", None):
                continue
            try:
                shape = ifcopenshell.geom.create_shape(settings, product)
                probe = ifcopenshell.geom.tree()
                probe.add_element(shape)
            except Exception as error:  # pragma: no cover - model-specific evidence
                skipped.append(
                    {
                        "global_id": product.GlobalId,
                        "ifc_class": product.is_a(),
                        "name": product.Name,
                        "error": str(error),
                    }
                )
                continue
            products.append(product)
    return products, skipped


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
    if not targets:
        raise RuntimeError("light fixture target selection is empty")
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

    source_lights = source.by_type("IfcLightFixture")
    candidate_lights = candidate.by_type("IfcLightFixture")
    source_targets = [source.by_guid(global_id) for global_id in sorted(target_ids)]
    candidate_targets = [
        candidate.by_guid(global_id) for global_id in sorted(target_ids)
    ]
    source_context, skipped_source_context = renderable_context(source)
    candidate_context = [
        candidate.by_guid(product.GlobalId) for product in source_context
    ]
    skipped_candidate_context_ids = {
        record["global_id"] for record in renderable_context(candidate)[1]
    }
    skipped_source_context_ids = {
        record["global_id"] for record in skipped_source_context
    }

    internal_contacts = classify_contact_changes(
        clearance_map(source_lights, source_lights, args.clearance_window_mm),
        clearance_map(
            candidate_lights, candidate_lights, args.clearance_window_mm
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
        set(),
    )
    source_internal_intersections = intersection_pairs(
        source_lights, source_lights, args.tolerance_mm
    )
    candidate_internal_intersections = intersection_pairs(
        candidate_lights, candidate_lights, args.tolerance_mm
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
        "light_fixtures": len(inventory),
        "targets": len(targets),
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
        "skipped_context_ids_equal": skipped_source_context_ids
        == skipped_candidate_context_ids,
        "skipped_unrenderable_context": sorted(skipped_source_context_ids),
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
        and not context_contacts["new_contacts"]
        and not gates["new_internal_intersections"]
        and not gates["lost_internal_intersections"]
        and not gates["new_context_intersections"]
        and not gates["lost_context_intersections"]
        and gates["skipped_context_ids_equal"]
        and gates["schema_equal"]
        and gates["entity_count_equal"]
        and gates["root_global_ids_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-light-fixture-installation-anchor-candidate",
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
        "skipped_source_context": skipped_source_context,
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
