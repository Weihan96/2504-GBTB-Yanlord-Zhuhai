#!/usr/bin/env python3
"""Verify the approved PVC110 bundled-product controlled exception.

The split candidate report freezes the exact source IFC SHA-256 from before
the preserve decision.  Byte identity with that frozen source is stronger than
a tolerance geometry comparison: when it holds, every world-geometry value is
unchanged and the mechanical change is exactly 0.0 mm.
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell

from geometry_alignment_audit import sha256
from pipe_bundle_split_candidate import (
    EXPECTED_BUNDLE_SIZE,
    EXPECTED_ITEM_CLASS,
    EXPECTED_NAME,
    TARGET_IDS,
    body_representation,
)


DECISION_ID = "COORD-FLOW-C003-PVC110-BUNDLE"
DECISION_SCOPE = "flow-segment-bundle-exception"


def parse_ids(value: str) -> set[str]:
    return {token.strip() for token in value.split(";") if token.strip()}


def approved_decision(path: Path) -> dict[str, str] | None:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    matching = [
        row
        for row in rows
        if row["decision_id"] == DECISION_ID
        and row["scope"] == DECISION_SCOPE
        and row["review_required"] == "no"
        and row["status"] == "implemented"
        and parse_ids(row["object_guid"]) == set(TARGET_IDS)
    ]
    return matching[0] if len(matching) == 1 else None


def byte_identity_world_change_mm(
    current_sha256: str, frozen_baseline_sha256: str
) -> float | None:
    return 0.0 if current_sha256 == frozen_baseline_sha256 else None


def formal_inventory(model: ifcopenshell.file) -> list[dict[str, Any]]:
    records = []
    for global_id in TARGET_IDS:
        product = model.by_guid(global_id)
        representation = body_representation(product)
        records.append(
            {
                "global_id": global_id,
                "ifc_class": product.is_a(),
                "name": product.Name,
                "body_item_count": len(representation.Items),
                "body_item_classes": [item.is_a() for item in representation.Items],
                "body_item_ids": [item.id() for item in representation.Items],
            }
        )
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--split-report", required=True, type=Path)
    parser.add_argument("--centerline-report", required=True, type=Path)
    parser.add_argument("--postwrite-report", type=Path)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = args.input.resolve()
    current_sha = sha256(source_path)
    split_report = json.loads(args.split_report.read_text(encoding="utf-8"))
    centerline_report = json.loads(
        args.centerline_report.read_text(encoding="utf-8")
    )
    postwrite_report = (
        json.loads(args.postwrite_report.read_text(encoding="utf-8"))
        if args.postwrite_report
        else None
    )
    decision = approved_decision(args.decisions)
    model = ifcopenshell.open(source_path)
    inventory = formal_inventory(model)
    split_inventory = {
        record["global_id"]: record for record in split_report["source_inventory"]
    }
    centreline_records = {
        record["global_id"]: record for record in centerline_report["records"]
    }
    frozen_sha = split_report["source"]["sha256"]
    world_change = byte_identity_world_change_mm(current_sha, frozen_sha)
    identity_method = "exact IFC SHA-256 byte identity"
    if world_change is None and postwrite_report is not None:
        if (
            postwrite_report["formal"]["sha256"] == current_sha
            and postwrite_report["gates"]["pass"]
        ):
            world_change = postwrite_report["gates"][
                "pvc110_max_world_geometry_change_from_prewrite_baseline_mm"
            ]
            identity_method = "targeted world-mesh comparison in C003 postwrite audit"
    branch_records = [
        {
            "global_id": global_id,
            "read_only_branch_count": centreline_records[global_id][
                "read_only_branch_count"
            ],
            "logical_component_groups": centreline_records[global_id][
                "logical_component_groups"
            ],
            "classification": centreline_records[global_id]["classification"],
            "world_geometry_change_from_frozen_baseline_mm": world_change,
        }
        for global_id in TARGET_IDS
    ]
    gates = {
        "approved_decision_exact": decision is not None,
        "formal_or_postwrite_geometry_baseline_proof": world_change == 0.0,
        "centerline_report_matches_formal": centerline_report["source"]["sha256"]
        == current_sha,
        "formal_product_count": len(inventory),
        "formal_products_preserved": all(
            record["ifc_class"] == "IfcFlowSegment"
            and record["name"] == EXPECTED_NAME
            and record["body_item_count"] == EXPECTED_BUNDLE_SIZE
            and set(record["body_item_classes"]) == {EXPECTED_ITEM_CLASS}
            for record in inventory
        ),
        "body_item_identity_preserved": all(
            record["body_item_ids"]
            == split_inventory[record["global_id"]]["body_item_ids"]
            for record in inventory
        ),
        "six_read_only_branches_recognized": sum(
            record["read_only_branch_count"] for record in branch_records
        )
        == 6
        and all(
            record["classification"] == "controlled_disconnected_bundle"
            for record in branch_records
        ),
        "world_geometry_change_from_frozen_baseline_mm": world_change,
        "formal_write_allowed": False,
    }
    gates["pass"] = (
        gates["approved_decision_exact"]
        and gates["formal_or_postwrite_geometry_baseline_proof"]
        and gates["centerline_report_matches_formal"]
        and gates["formal_product_count"] == 2
        and gates["formal_products_preserved"]
        and gates["body_item_identity_preserved"]
        and gates["six_read_only_branches_recognized"]
        and gates["world_geometry_change_from_frozen_baseline_mm"] == 0.0
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-pvc110-controlled-exception-audit",
        "source": {
            "path": str(source_path),
            "sha256": current_sha,
            "schema": model.schema,
        },
        "frozen_baseline": {
            "source": str(args.split_report.resolve()),
            "sha256": frozen_sha,
            "identity_method": identity_method,
        },
        "decision": {
            "decision_id": DECISION_ID,
            "scope": DECISION_SCOPE,
            "status": decision["status"] if decision else None,
            "review_required": decision["review_required"] if decision else None,
            "global_ids": sorted(parse_ids(decision["object_guid"]))
            if decision
            else [],
        },
        "formal_inventory": inventory,
        "branch_records": branch_records,
        "split_candidate": {
            "path": split_report["candidate"]["path"],
            "sha256": split_report["candidate"]["sha256"],
            "disposition": "verification evidence only; rejected for formal write",
        },
        "postwrite_report": (
            {
                "path": str(args.postwrite_report.resolve()),
                "formal_sha256": postwrite_report["formal"]["sha256"],
                "pass": postwrite_report["gates"]["pass"],
            }
            if postwrite_report is not None
            else None
        ),
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(gates, ensure_ascii=False))
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
