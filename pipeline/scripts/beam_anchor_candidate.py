#!/usr/bin/env python3
"""Build and validate a read-only C003 beam anchor candidate."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api.geometry
import ifcopenshell.util.placement
import numpy as np

from direction_noise_candidate import all_product_geometry_difference
from geometry_alignment_audit import geometry_difference_audit, sha256


def read_targets(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    expected = {
        "global_id",
        "target_x_mm",
        "target_y_mm",
        "target_z_mm",
        "basis",
        "confidence",
        "review_required",
    }
    if not rows or set(rows[0]) != expected:
        raise RuntimeError(f"invalid beam target schema: {path}")
    result = []
    for row in rows:
        if None in row:
            raise RuntimeError(f"extra columns for beam {row['global_id']}")
        target = [
            float(row["target_x_mm"]),
            float(row["target_y_mm"]),
            float(row["target_z_mm"]),
        ]
        if any(abs(value - round(value)) > 1e-9 for value in target):
            raise RuntimeError(f"beam target is not integer millimetres: {row}")
        if row["review_required"] not in {"yes", "no"}:
            raise RuntimeError(f"invalid review flag for {row['global_id']}")
        result.append(
            {
                "global_id": row["global_id"],
                "target_mm": target,
                "basis": row["basis"],
                "confidence": float(row["confidence"]),
                "review_required": row["review_required"] == "yes",
            }
        )
    return result


def select_targets(
    rows: list[dict[str, Any]], review_policy: str
) -> list[dict[str, Any]]:
    if review_policy == "all":
        return rows
    if review_policy == "no-review-only":
        return [row for row in rows if not row["review_required"]]
    raise ValueError(review_policy)


def apply_targets(
    model: ifcopenshell.file,
    rows: list[dict[str, Any]],
    maximum_shift_mm: float,
) -> tuple[list[dict[str, Any]], set[str]]:
    results = []
    allowed_geometry_ids: set[str] = set()
    for row in rows:
        beam = model.by_guid(row["global_id"])
        if beam is None or not beam.is_a("IfcBeam"):
            raise RuntimeError(f"missing target IfcBeam {row['global_id']}")
        source = np.array(
            ifcopenshell.util.placement.get_local_placement(beam.ObjectPlacement),
            dtype=float,
        )
        target = source.copy()
        target[:3, 3] = np.array(row["target_mm"], dtype=float)
        delta = target[:3, 3] - source[:3, 3]
        shift = float(np.linalg.norm(delta))
        if shift > maximum_shift_mm:
            raise RuntimeError(
                f"beam {beam.GlobalId} shift {shift} mm exceeds {maximum_shift_mm} mm"
            )
        child_ids = set()
        for relation in beam.HasOpenings:
            opening = relation.RelatedOpeningElement
            child_ids.add(opening.GlobalId)
            child_ids.update(
                filling.RelatedBuildingElement.GlobalId
                for filling in opening.HasFillings
            )
        ifcopenshell.api.geometry.edit_object_placement(
            model,
            product=beam,
            matrix=target,
            is_si=False,
            should_transform_children=False,
        )
        after = np.array(
            ifcopenshell.util.placement.get_local_placement(beam.ObjectPlacement),
            dtype=float,
        )
        results.append(
            {
                **row,
                "source_mm": source[:3, 3].tolist(),
                "result_mm": after[:3, 3].tolist(),
                "translation_by_axis_mm": delta.tolist(),
                "translation_shift_mm": shift,
                "child_global_ids": sorted(child_ids),
            }
        )
        allowed_geometry_ids.add(beam.GlobalId)
        allowed_geometry_ids.update(child_ids)
    return results, allowed_geometry_ids


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--maximum-shift-mm", type=float, default=1.0)
    parser.add_argument(
        "--review-policy",
        choices=("all", "no-review-only"),
        default="all",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= 1.0:
        raise SystemExit("--tolerance-mm must be greater than zero and at most 1 mm")
    if not 0.0 < args.maximum_shift_mm <= 1.0:
        raise SystemExit("--maximum-shift-mm must be greater than zero and at most 1 mm")
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    all_rows = read_targets(args.targets)
    rows = select_targets(all_rows, args.review_policy)
    if not rows:
        raise RuntimeError("beam target selection is empty")
    candidate = ifcopenshell.open(source_path)
    placement_results, allowed_geometry_ids = apply_targets(
        candidate, rows, args.maximum_shift_mm
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    target_geometry = geometry_difference_audit(
        candidate,
        source,
        str(source_path),
        tolerance_mm=args.maximum_shift_mm + args.tolerance_mm,
        classes=[],
        global_ids=sorted(allowed_geometry_ids),
    )
    all_geometry = all_product_geometry_difference(
        candidate, source, args.tolerance_mm
    )
    unexpected_geometry = [
        record
        for record in all_geometry["records"]
        if not record["within_tolerance"]
        and record.get("global_id") not in allowed_geometry_ids
    ]
    target_mismatches = []
    for row in rows:
        beam = candidate.by_guid(row["global_id"])
        matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(beam.ObjectPlacement),
            dtype=float,
        )
        residual = max(
            abs(float(matrix[index, 3]) - row["target_mm"][index])
            for index in range(3)
        )
        if residual > args.tolerance_mm:
            target_mismatches.append(
                {"global_id": beam.GlobalId, "residual_mm": residual}
            )
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "target_beams": len(rows),
        "allowed_geometry_objects": len(allowed_geometry_ids),
        "target_placement_mismatches": target_mismatches,
        "target_geometry_over_limit": target_geometry["over_tolerance"],
        "unexpected_product_geometry_records": len(unexpected_geometry),
        "schema_equal": source.schema == candidate.schema,
        "entity_count_equal": len(list(source)) == len(list(candidate)),
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
    }
    gates["pass"] = (
        not gates["target_placement_mismatches"]
        and gates["target_geometry_over_limit"] == 0
        and gates["unexpected_product_geometry_records"] == 0
        and gates["schema_equal"]
        and gates["entity_count_equal"]
        and gates["root_global_ids_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-beam-anchor-candidate",
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
        "maximum_shift_mm": args.maximum_shift_mm,
        "review_policy": args.review_policy,
        "available_target_beams": len(all_rows),
        "placement_results": placement_results,
        "allowed_geometry_global_ids": sorted(allowed_geometry_ids),
        "target_geometry": target_geometry,
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
