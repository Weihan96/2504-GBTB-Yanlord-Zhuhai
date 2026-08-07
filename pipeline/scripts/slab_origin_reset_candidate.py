#!/usr/bin/env python3
"""Build a slab origin-reset candidate without moving world geometry."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
import ifcopenshell.util.placement
import numpy as np

from beam_origin_reset_candidate import read_targets, reset_product_origin
from direction_noise_candidate import all_product_geometry_difference
from geometry_alignment_audit import sha256


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
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
    rows = read_targets(args.targets)
    results = [
        reset_product_origin(candidate, row, args.tolerance_mm, "IfcSlab")
        for row in rows
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    all_geometry = all_product_geometry_difference(
        candidate, source, args.tolerance_mm
    )
    geometry_over_tolerance = [
        record
        for record in all_geometry["records"]
        if not record["within_tolerance"]
    ]
    placement_mismatches = []
    for row in rows:
        slab = candidate.by_guid(row["global_id"])
        matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(slab.ObjectPlacement),
            dtype=float,
        )
        residual = float(np.max(np.abs(matrix[:3, 3] - row["target_mm"])))
        if residual > args.tolerance_mm:
            placement_mismatches.append(
                {"global_id": slab.GlobalId, "residual_mm": residual}
            )
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    entity_count_delta = len(list(candidate)) - len(list(source))
    expected_entity_count_delta = sum(
        result["representation_entities_added"] for result in results
    )
    gates = {
        "target_slabs": len(rows),
        "target_placement_mismatches": placement_mismatches,
        "anchors_over_tolerance": sum(
            result["anchor_surface_distance_mm"] > args.tolerance_mm
            for result in results
        ),
        "product_geometry_over_tolerance": len(geometry_over_tolerance),
        "schema_equal": source.schema == candidate.schema,
        "entity_count_delta": entity_count_delta,
        "expected_entity_count_delta": expected_entity_count_delta,
        "entity_count_delta_matches_expected": entity_count_delta
        == expected_entity_count_delta,
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
    }
    gates["pass"] = (
        not gates["target_placement_mismatches"]
        and gates["anchors_over_tolerance"] == 0
        and gates["product_geometry_over_tolerance"] == 0
        and gates["schema_equal"]
        and gates["entity_count_delta_matches_expected"]
        and gates["root_global_ids_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-slab-origin-reset-candidate",
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
