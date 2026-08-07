#!/usr/bin/env python3
"""Build the atomic no-world-movement candidate for eight ready C003 origins."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.placement
import numpy as np

from beam_origin_reset_candidate import reset_product_origin
from direction_noise_candidate import all_product_geometry_difference
from element_assembly_origin_reset_candidate import (
    maximum_descendant_world_placement_residual,
    read_targets as read_assembly_targets,
    relation_signature,
    reset_assembly_origin,
)
from geometry_alignment_audit import sha256
from structural_origin_reset_candidate import read_targets as read_product_targets


def world_matrix(product: Any) -> np.ndarray:
    return np.asarray(
        ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
        dtype=float,
    )


def placement_difference(
    source: ifcopenshell.file,
    candidate: ifcopenshell.file,
    allowed_changed_ids: set[str],
    tolerance_mm: float,
) -> dict[str, Any]:
    records = []
    for product in source.by_type("IfcProduct"):
        if not getattr(product, "ObjectPlacement", None):
            continue
        candidate_product = candidate.by_guid(product.GlobalId)
        if candidate_product is None or not getattr(
            candidate_product, "ObjectPlacement", None
        ):
            records.append(
                {
                    "global_id": product.GlobalId,
                    "missing_candidate_placement": True,
                    "maximum_matrix_residual_mm": None,
                    "allowed_change": product.GlobalId in allowed_changed_ids,
                }
            )
            continue
        residual = float(
            np.max(np.abs(world_matrix(product) - world_matrix(candidate_product)))
        )
        if residual > tolerance_mm or product.GlobalId in allowed_changed_ids:
            records.append(
                {
                    "global_id": product.GlobalId,
                    "missing_candidate_placement": False,
                    "maximum_matrix_residual_mm": residual,
                    "allowed_change": product.GlobalId in allowed_changed_ids,
                }
            )
    return {
        "records": records,
        "non_target_changes": [
            record
            for record in records
            if not record["allowed_change"]
            and (
                record["missing_candidate_placement"]
                or record["maximum_matrix_residual_mm"] > tolerance_mm
            )
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--cladding-targets", required=True, type=Path)
    parser.add_argument("--bay-window-targets", required=True, type=Path)
    parser.add_argument("--cabinet-targets", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= 1.0:
        raise SystemExit("--tolerance-mm must be greater than zero and at most 1 mm")
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    product_rows = read_product_targets(args.cladding_targets)
    assembly_rows = [
        *read_assembly_targets(args.bay_window_targets),
        *read_assembly_targets(args.cabinet_targets),
    ]
    target_ids = {
        row["global_id"] for row in [*product_rows, *assembly_rows]
    }
    if len(target_ids) != 8:
        raise RuntimeError(f"expected eight unique C003 targets, got {len(target_ids)}")
    source_relations = relation_signature(source)
    product_results = [
        reset_product_origin(
            candidate,
            row,
            args.tolerance_mm,
            row["expected_class"],
        )
        for row in product_rows
    ]
    assembly_results = [
        reset_assembly_origin(candidate, row, args.tolerance_mm)
        for row in assembly_rows
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)

    for result in assembly_results:
        result["maximum_descendant_world_placement_residual_mm"] = (
            maximum_descendant_world_placement_residual(
                source, candidate, result["global_id"]
            )
        )
    geometry = all_product_geometry_difference(
        candidate, source, args.tolerance_mm
    )
    geometry_changes = [
        record for record in geometry["records"] if not record["within_tolerance"]
    ]
    placements = placement_difference(
        source, candidate, target_ids, args.tolerance_mm
    )
    target_residuals = []
    for row in [*product_rows, *assembly_rows]:
        actual = world_matrix(candidate.by_guid(row["global_id"]))[:3, 3]
        target_residuals.append(
            {
                "global_id": row["global_id"],
                "maximum_target_residual_mm": float(
                    np.max(np.abs(actual - row["target_mm"]))
                ),
            }
        )
    source_roots = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_roots = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "target_count": len(target_ids),
        "product_targets": len(product_rows),
        "assembly_targets": len(assembly_rows),
        "target_placement_mismatches": [
            record
            for record in target_residuals
            if record["maximum_target_residual_mm"] > args.tolerance_mm
        ],
        "product_anchor_distances_over_tolerance": sum(
            result["anchor_surface_distance_mm"] > args.tolerance_mm
            for result in product_results
        ),
        "assembly_anchor_distances_over_tolerance": sum(
            result["anchor_surface_distance_mm"] > args.tolerance_mm
            for result in assembly_results
        ),
        "assembly_descendant_placement_residuals_over_tolerance": sum(
            result["maximum_descendant_world_placement_residual_mm"]
            > args.tolerance_mm
            for result in assembly_results
        ),
        "product_geometry_changes_over_tolerance": len(geometry_changes),
        "non_target_placement_changes_over_tolerance": len(
            placements["non_target_changes"]
        ),
        "aggregate_relations_equal": source_relations
        == relation_signature(candidate),
        "schema_equal": source.schema == candidate.schema,
        "root_global_ids_equal": source_roots == candidate_roots,
        "entity_count_delta": len(list(candidate)) - len(list(source)),
    }
    gates["pass"] = (
        gates["target_count"] == 8
        and gates["product_targets"] == 5
        and gates["assembly_targets"] == 3
        and not gates["target_placement_mismatches"]
        and gates["product_anchor_distances_over_tolerance"] == 0
        and gates["assembly_anchor_distances_over_tolerance"] == 0
        and gates["assembly_descendant_placement_residuals_over_tolerance"] == 0
        and gates["product_geometry_changes_over_tolerance"] == 0
        and gates["non_target_placement_changes_over_tolerance"] == 0
        and gates["aggregate_relations_equal"]
        and gates["schema_equal"]
        and gates["root_global_ids_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-c003-pending-origin-batch-candidate",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
            "entity_count": len(list(source)),
        },
        "candidate": {
            "path": str(args.output.resolve()),
            "sha256": sha256(args.output),
            "schema": candidate.schema,
            "entity_count": len(list(candidate)),
        },
        "tolerance_mm": args.tolerance_mm,
        "product_results": product_results,
        "assembly_results": assembly_results,
        "placement_difference": placements,
        "all_product_geometry": geometry,
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
