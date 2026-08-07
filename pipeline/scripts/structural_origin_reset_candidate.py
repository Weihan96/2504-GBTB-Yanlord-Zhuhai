#!/usr/bin/env python3
"""Build one no-world-movement origin-reset candidate for approved fixed products."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.placement
import numpy as np

from beam_origin_reset_candidate import reset_product_origin
from direction_noise_candidate import all_product_geometry_difference
from geometry_alignment_audit import sha256


def read_targets(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    expected = {
        "expected_class",
        "global_id",
        "target_x_mm",
        "target_y_mm",
        "target_z_mm",
        "anchor_kind",
        "basis",
        "confidence",
    }
    if not rows or set(rows[0]) != expected:
        raise RuntimeError(f"invalid structural origin target schema: {path}")
    allowed_classes = {
        "IfcBeam",
        "IfcBuildingElementProxy",
        "IfcCovering",
        "IfcElectricAppliance",
        "IfcFurniture",
        "IfcLightFixture",
        "IfcPipeSegment",
        "IfcSanitaryTerminal",
        "IfcSlab",
    }
    result = []
    for row in rows:
        if row["expected_class"] not in allowed_classes:
            raise RuntimeError(f"unsupported fixed-product class: {row['expected_class']}")
        target = np.array(
            [
                float(row["target_x_mm"]),
                float(row["target_y_mm"]),
                float(row["target_z_mm"]),
            ],
            dtype=float,
        )
        if any(abs(value - round(value)) > 1e-9 for value in target):
            raise RuntimeError(f"target is not integer millimetres: {row}")
        result.append(
            {
                "expected_class": row["expected_class"],
                "global_id": row["global_id"],
                "target_mm": target,
                "anchor_kind": row["anchor_kind"],
                "basis": row["basis"],
                "confidence": float(row["confidence"]),
            }
        )
    if len({row["global_id"] for row in result}) != len(result):
        raise RuntimeError("structural origin target table contains duplicate GlobalIds")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def scoped_owner_history_additions(
    source: ifcopenshell.file,
    candidate: ifcopenshell.file,
    allowed_root_global_ids: set[str],
) -> dict[str, Any]:
    source_ids = {
        history.id()
        for history in source.by_type("IfcOwnerHistory", include_subtypes=False)
    }
    additions = [
        history
        for history in candidate.by_type("IfcOwnerHistory", include_subtypes=False)
        if history.id() not in source_ids
    ]
    records = []
    scoped = True
    for history in additions:
        owners = list(candidate.get_inverse(history))
        owner_ids = sorted(
            owner.GlobalId
            for owner in owners
            if owner.is_a("IfcRoot") and getattr(owner, "GlobalId", None)
        )
        valid = bool(owners) and len(owner_ids) == len(owners) and set(owner_ids) <= allowed_root_global_ids
        scoped = scoped and valid
        records.append(
            {
                "owner_history_id": history.id(),
                "owner_global_ids": owner_ids,
                "scoped_to_target_or_related_opening": valid,
            }
        )
    return {"count": len(additions), "records": records, "scoped": scoped}


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= 1.0:
        raise SystemExit("--tolerance-mm must be greater than zero and at most 1 mm")
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    rows = read_targets(args.targets)
    results = [
        reset_product_origin(
            candidate,
            row,
            args.tolerance_mm,
            row["expected_class"],
        )
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
        product = candidate.by_guid(row["global_id"])
        matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
            dtype=float,
        )
        residual = float(np.max(np.abs(matrix[:3, 3] - row["target_mm"])))
        if residual > args.tolerance_mm:
            placement_mismatches.append(
                {"global_id": product.GlobalId, "residual_mm": residual}
            )
    entity_count_delta = len(list(candidate)) - len(list(source))
    expected_entity_count_delta = sum(
        result["representation_entities_added"] for result in results
    )
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    related_opening_ids = sorted(
        {
            global_id
            for result in results
            for global_id in result["opening_global_ids"]
        }
    )
    allowed_owner_history_roots = {
        row["global_id"] for row in rows
    } | set(related_opening_ids)
    owner_history_additions = scoped_owner_history_additions(
        source, candidate, allowed_owner_history_roots
    )
    expected_entity_count_delta += owner_history_additions["count"]
    gates = {
        "target_products": len(rows),
        "targets_by_class": dict(
            sorted(Counter(row["expected_class"] for row in rows).items())
        ),
        "related_openings": len(related_opening_ids),
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
        "owner_history_additions": owner_history_additions,
    }
    gates["pass"] = (
        not gates["target_placement_mismatches"]
        and gates["anchors_over_tolerance"] == 0
        and gates["product_geometry_over_tolerance"] == 0
        and gates["schema_equal"]
        and gates["entity_count_delta_matches_expected"]
        and gates["root_global_ids_equal"]
        and gates["owner_history_additions"]["scoped"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-structural-origin-reset-candidate",
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
        "related_opening_global_ids": related_opening_ids,
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
