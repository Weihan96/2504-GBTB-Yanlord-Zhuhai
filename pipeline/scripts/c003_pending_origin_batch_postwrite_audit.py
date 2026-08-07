#!/usr/bin/env python3
"""Compare the formal IFC with the approved eight-origin C003 candidate."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.placement
import numpy as np

from direction_noise_candidate import (
    all_product_geometry_difference,
    vertex_hausdorff_within_tolerance_mm,
)
from element_assembly_origin_reset_candidate import relation_signature
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm
from pipe_bundle_split_candidate import TARGET_IDS as PVC110_TARGET_IDS


def world_matrix(product: Any) -> np.ndarray:
    return np.asarray(
        ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
        dtype=float,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formal", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--candidate-report", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    formal_path = args.formal.resolve()
    candidate_path = args.candidate.resolve()
    candidate_report = json.loads(
        args.candidate_report.read_text(encoding="utf-8")
    )
    formal = ifcopenshell.open(formal_path)
    candidate = ifcopenshell.open(candidate_path)
    target_rows = [
        *candidate_report["product_results"],
        *candidate_report["assembly_results"],
    ]
    target_ids = {row["global_id"] for row in target_rows}
    target_placement_records = []
    for row in target_rows:
        formal_matrix = world_matrix(formal.by_guid(row["global_id"]))
        candidate_matrix = world_matrix(candidate.by_guid(row["global_id"]))
        target_placement_records.append(
            {
                "global_id": row["global_id"],
                "formal_to_candidate_matrix_residual_mm": float(
                    np.max(np.abs(formal_matrix - candidate_matrix))
                ),
                "formal_to_approved_target_residual_mm": float(
                    np.max(
                        np.abs(
                            formal_matrix[:3, 3]
                            - np.asarray(row["target_mm"], dtype=float)
                        )
                    )
                ),
            }
        )

    placement_mismatches = []
    for product in candidate.by_type("IfcProduct"):
        if not getattr(product, "ObjectPlacement", None):
            continue
        formal_product = formal.by_guid(product.GlobalId)
        if formal_product is None or not getattr(
            formal_product, "ObjectPlacement", None
        ):
            placement_mismatches.append(
                {"global_id": product.GlobalId, "reason": "missing formal placement"}
            )
            continue
        residual = float(
            np.max(np.abs(world_matrix(product) - world_matrix(formal_product)))
        )
        if residual > args.tolerance_mm:
            placement_mismatches.append(
                {"global_id": product.GlobalId, "residual_mm": residual}
            )

    geometry = all_product_geometry_difference(
        formal, candidate, args.tolerance_mm
    )
    geometry_changes = [
        record for record in geometry["records"] if not record["within_tolerance"]
    ]
    pvc110_records = []
    settings = geometry_settings()
    for global_id in PVC110_TARGET_IDS:
        baseline_product = candidate.by_guid(global_id)
        formal_product = formal.by_guid(global_id)
        baseline_vertices, _ = world_mesh_mm(settings, baseline_product)
        formal_vertices, _ = world_mesh_mm(settings, formal_product)
        body = [
            representation
            for representation in formal_product.Representation.Representations
            if representation.RepresentationIdentifier == "Body"
        ][0]
        pvc110_records.append(
            {
                "global_id": global_id,
                "formal_ifc_class": formal_product.is_a(),
                "formal_body_item_count": len(body.Items),
                "world_geometry_change_from_prewrite_baseline_mm": (
                    vertex_hausdorff_within_tolerance_mm(
                        baseline_vertices,
                        formal_vertices,
                        args.tolerance_mm,
                    )
                ),
            }
        )

    formal_roots = {root.GlobalId for root in formal.by_type("IfcRoot")}
    candidate_roots = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "approved_candidate_hash_matches_report": sha256(candidate_path)
        == candidate_report["candidate"]["sha256"],
        "approved_candidate_source_hash": candidate_report["source"]["sha256"],
        "candidate_prewrite_geometry_gate_passed": candidate_report["gates"][
            "product_geometry_changes_over_tolerance"
        ]
        == 0,
        "target_count": len(target_ids),
        "target_placement_mismatches": [
            record
            for record in target_placement_records
            if record["formal_to_candidate_matrix_residual_mm"]
            > args.tolerance_mm
            or record["formal_to_approved_target_residual_mm"]
            > args.tolerance_mm
        ],
        "all_product_placement_mismatches": placement_mismatches,
        "all_product_geometry_changes_over_tolerance": len(geometry_changes),
        "pvc110_product_count": len(pvc110_records),
        "pvc110_products_preserved": all(
            record["formal_ifc_class"] == "IfcFlowSegment"
            and record["formal_body_item_count"] == 3
            for record in pvc110_records
        ),
        "pvc110_max_world_geometry_change_from_prewrite_baseline_mm": max(
            record["world_geometry_change_from_prewrite_baseline_mm"]
            for record in pvc110_records
        ),
        "aggregate_relations_equal": relation_signature(formal)
        == relation_signature(candidate),
        "schema_equal": formal.schema == candidate.schema,
        "entity_count_equal": len(list(formal)) == len(list(candidate)),
        "root_global_ids_equal": formal_roots == candidate_roots,
    }
    gates["pass"] = (
        gates["approved_candidate_hash_matches_report"]
        and gates["candidate_prewrite_geometry_gate_passed"]
        and gates["target_count"] == 8
        and not gates["target_placement_mismatches"]
        and not gates["all_product_placement_mismatches"]
        and gates["all_product_geometry_changes_over_tolerance"] == 0
        and gates["pvc110_product_count"] == 2
        and gates["pvc110_products_preserved"]
        and gates["pvc110_max_world_geometry_change_from_prewrite_baseline_mm"]
        == 0.0
        and gates["aggregate_relations_equal"]
        and gates["schema_equal"]
        and gates["entity_count_equal"]
        and gates["root_global_ids_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-c003-pending-origin-batch-postwrite-audit",
        "formal": {
            "path": str(formal_path),
            "sha256": sha256(formal_path),
            "entity_count": len(list(formal)),
        },
        "approved_candidate": {
            "path": str(candidate_path),
            "sha256": sha256(candidate_path),
            "entity_count": len(list(candidate)),
        },
        "tolerance_mm": args.tolerance_mm,
        "target_placement_records": target_placement_records,
        "pvc110_records": pvc110_records,
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
