#!/usr/bin/env python3
"""Verify a formal IFC against one approved pure-origin-reset candidate."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
import ifcopenshell.util.placement
import numpy as np

from direction_noise_candidate import all_product_geometry_difference
from geometry_alignment_audit import sha256
from structural_origin_reset_candidate import read_targets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formal", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--candidate-report", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= 1.0:
        raise SystemExit("--tolerance-mm must be greater than zero and at most 1 mm")

    formal_path = args.formal.resolve()
    candidate_path = args.candidate.resolve()
    formal = ifcopenshell.open(formal_path)
    candidate = ifcopenshell.open(candidate_path)
    candidate_report = json.loads(args.candidate_report.read_text(encoding="utf-8"))
    targets = read_targets(args.targets)

    geometry = all_product_geometry_difference(
        formal, candidate, args.tolerance_mm
    )
    geometry_over_tolerance = [
        record for record in geometry["records"] if not record["within_tolerance"]
    ]
    placement_records = []
    for row in targets:
        formal_product = formal.by_guid(row["global_id"])
        candidate_product = candidate.by_guid(row["global_id"])
        formal_origin = np.array(
            ifcopenshell.util.placement.get_local_placement(
                formal_product.ObjectPlacement
            )[:3, 3],
            dtype=float,
        )
        candidate_origin = np.array(
            ifcopenshell.util.placement.get_local_placement(
                candidate_product.ObjectPlacement
            )[:3, 3],
            dtype=float,
        )
        placement_records.append(
            {
                "global_id": row["global_id"],
                "target_mm": row["target_mm"].tolist(),
                "formal_origin_mm": formal_origin.tolist(),
                "candidate_origin_mm": candidate_origin.tolist(),
                "formal_target_residual_mm": float(
                    np.max(np.abs(formal_origin - row["target_mm"]))
                ),
                "formal_candidate_residual_mm": float(
                    np.max(np.abs(formal_origin - candidate_origin))
                ),
            }
        )

    formal_root_ids = {root.GlobalId for root in formal.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    formal_hash = sha256(formal_path)
    candidate_hash = sha256(candidate_path)
    gates = {
        "approved_candidate_hash_matches_disk": candidate_hash
        == candidate_report["candidate"]["sha256"],
        "approved_candidate_prewrite_passed": candidate_report["gates"]["pass"]
        is True,
        "target_products": len(targets),
        "placement_mismatches": [
            record
            for record in placement_records
            if record["formal_target_residual_mm"] > args.tolerance_mm
            or record["formal_candidate_residual_mm"] > args.tolerance_mm
        ],
        "product_geometry_over_tolerance": len(geometry_over_tolerance),
        "schema_equal": formal.schema == candidate.schema,
        "entity_count_equal": len(list(formal)) == len(list(candidate)),
        "root_global_ids_equal": formal_root_ids == candidate_root_ids,
    }
    gates["pass"] = all(
        [
            gates["approved_candidate_hash_matches_disk"],
            gates["approved_candidate_prewrite_passed"],
            not gates["placement_mismatches"],
            gates["product_geometry_over_tolerance"] == 0,
            gates["schema_equal"],
            gates["entity_count_equal"],
            gates["root_global_ids_equal"],
        ]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-structural-origin-reset-postwrite-audit",
        "formal": {
            "path": str(formal_path),
            "sha256": formal_hash,
            "schema": formal.schema,
        },
        "candidate": {
            "path": str(candidate_path),
            "sha256": candidate_hash,
            "schema": candidate.schema,
        },
        "approved_candidate_source_hash": candidate_report["source"]["sha256"],
        "tolerance_mm": args.tolerance_mm,
        "placements": placement_records,
        "all_product_geometry": geometry,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
