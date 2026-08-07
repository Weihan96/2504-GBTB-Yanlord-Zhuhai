#!/usr/bin/env python3
"""Classify unfilled Opening origins as host-relative or independent."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.placement
import ifcopenshell.util.unit


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify_opening_relationship(
    host_count: int,
    filling_count: int,
    placement_relative_to_host: bool,
) -> dict[str, Any]:
    if (
        host_count == 1
        and filling_count == 0
        and placement_relative_to_host
    ):
        return {
            "normalization_disposition": "delegated_to_host_placement",
            "basis": (
                "one unfilled Opening is directly relative to its sole host "
                "ObjectPlacement and must follow the host rather than be snapped independently"
            ),
            "review_required": False,
            "automatic_write_allowed": False,
        }
    return {
        "normalization_disposition": "independent_opening_anchor_review",
        "basis": (
            "the Opening placement is not a direct child of one sole unfilled host; "
            "its Boolean intent must be reviewed independently"
        ),
        "review_required": True,
        "automatic_write_allowed": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= 1.0:
        raise SystemExit("--tolerance-mm must be greater than zero and at most 1 mm")
    source_path = args.input.resolve()
    model = ifcopenshell.open(source_path)
    unit_to_mm = ifcopenshell.util.unit.calculate_unit_scale(model) * 1000.0
    records: list[dict[str, Any]] = []
    for opening in model.by_type("IfcOpeningElement"):
        if not opening.ObjectPlacement:
            continue
        matrix = ifcopenshell.util.placement.get_local_placement(
            opening.ObjectPlacement
        )
        current = [float(matrix[index, 3]) * unit_to_mm for index in range(3)]
        delta = [round(value) - value for value in current]
        if max(abs(value) for value in delta) <= args.tolerance_mm:
            continue
        voids = list(opening.VoidsElements or ())
        fillings = list(opening.HasFillings or ())
        hosts = [relation.RelatingBuildingElement for relation in voids]
        placement_relative_to_host = bool(
            len(hosts) == 1
            and opening.ObjectPlacement.PlacementRelTo == hosts[0].ObjectPlacement
        )
        classification = classify_opening_relationship(
            len(hosts), len(fillings), placement_relative_to_host
        )
        records.append(
            {
                "global_id": opening.GlobalId,
                "name": opening.Name,
                "origin_mm": current,
                "nearest_integer_delta_mm": delta,
                "host_global_ids": [host.GlobalId for host in hosts],
                "host_classes": [host.is_a() for host in hosts],
                "filling_global_ids": [
                    relation.RelatedBuildingElement.GlobalId for relation in fillings
                ],
                "placement_relative_to_host": placement_relative_to_host,
                **classification,
            }
        )
    records.sort(key=lambda record: record["global_id"])
    summary = {
        "openings_over_tolerance": len(records),
        "delegated_to_host_placement": sum(
            record["normalization_disposition"] == "delegated_to_host_placement"
            for record in records
        ),
        "independent_opening_anchor_review": sum(
            record["review_required"] for record in records
        ),
        "automatic_write_allowed": sum(
            record["automatic_write_allowed"] for record in records
        ),
        "delegated_by_host_class": dict(
            sorted(
                Counter(
                    record["host_classes"][0]
                    for record in records
                    if record["normalization_disposition"]
                    == "delegated_to_host_placement"
                ).items()
            )
        ),
    }
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-opening-anchor-classification",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": model.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "summary": summary,
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **summary}, ensure_ascii=False))


if __name__ == "__main__":
    main()
