#!/usr/bin/env python3
"""Verify the formal A-102 demolition-wall write against its approved candidate and Git baseline."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element

from a103_wall_plan_candidate import geometry_settings, world_bbox_mm
from geometry_alignment_audit import (
    geometry_difference_audit,
    load_git_ifc,
    sha256,
)


COMPARE_CLASSES = ("IfcWall", "IfcOpeningElement", "IfcDoor", "IfcWindow")
PROTECTED_IDS = ("0hKdvAZkn1TejLgJhK_vDp", "1YxMx6s0r3ZPPohkRKXWbl")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formal", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--register", required=True, type=Path)
    parser.add_argument("--baseline-git-ref", default="424dc17")
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def products_by_id(model: ifcopenshell.file) -> dict[str, Any]:
    return {
        product.GlobalId: product
        for ifc_class in COMPARE_CLASSES
        for product in model.by_type(ifc_class)
        if getattr(product, "GlobalId", None)
    }


def main() -> None:
    args = parse_args()
    formal_path = args.formal.resolve()
    candidate_path = args.candidate.resolve()
    formal = ifcopenshell.open(formal_path)
    candidate = ifcopenshell.open(candidate_path)
    baseline, baseline_source = load_git_ifc(formal_path, args.baseline_git_ref)
    with args.register.open(newline="", encoding="utf-8") as handle:
        register = list(csv.DictReader(handle))

    candidate_ids = [row["candidate_global_id"] for row in register]
    baseline_ids = sorted(products_by_id(baseline))
    candidate_compare_ids = sorted(products_by_id(candidate))
    formal_to_baseline = geometry_difference_audit(
        formal,
        baseline,
        baseline_source,
        tolerance_mm=args.tolerance_mm,
        classes=(),
        global_ids=baseline_ids,
    )
    formal_to_candidate = geometry_difference_audit(
        formal,
        candidate,
        str(candidate_path),
        tolerance_mm=args.tolerance_mm,
        classes=(),
        global_ids=candidate_compare_ids,
    )
    protected_geometry = geometry_difference_audit(
        formal,
        baseline,
        baseline_source,
        tolerance_mm=args.tolerance_mm,
        classes=(),
        global_ids=list(PROTECTED_IDS),
    )

    settings = geometry_settings()
    demolition_records = []
    for row in register:
        wall = formal.by_guid(row["candidate_global_id"])
        psets = ifcopenshell.util.element.get_psets(wall)
        review = psets.get("Pset_A102DemolitionReview", {})
        quantities = psets.get("Qto_WallBaseQuantities", {})
        bbox = world_bbox_mm(settings, wall)
        expected = [
            float(row["candidate_x_min_mm"]),
            float(row["candidate_y_min_mm"]),
            float(row["candidate_z_min_mm"]),
            float(row["candidate_x_max_mm"]),
            float(row["candidate_y_max_mm"]),
            float(row["candidate_z_max_mm"]),
        ]
        actual = bbox["min_mm"] + bbox["max_mm"]
        bbox_delta = max(abs(a - b) for a, b in zip(actual, expected))
        expected_quantities = {
            "Length": float(row["nominal_length_mm"]),
            "Width": float(row["nominal_thickness_mm"]),
            "Height": float(row["candidate_z_max_mm"]) - float(row["candidate_z_min_mm"]),
        }
        actual_quantities = {
            key: quantities.get(key) for key in ("Length", "Width", "Height")
        }
        demolition_records.append(
            {
                "candidate_id": row["candidate_id"],
                "global_id": wall.GlobalId,
                "status": psets.get("Pset_WallCommon", {}).get("Status"),
                "review_status": review.get("ReviewStatus"),
                "formal_ifc_write_allowed": review.get("FormalIfcWriteAllowed"),
                "accuracy_boundary": review.get("AccuracyBoundary"),
                "human_review_basis": review.get("HumanReviewBasis"),
                "maximum_bbox_delta_mm": bbox_delta,
                "quantities_mm": actual_quantities,
                "expected_quantities_mm": expected_quantities,
                "pass": (
                    bbox_delta <= args.tolerance_mm
                    and actual_quantities == expected_quantities
                    and psets.get("Pset_WallCommon", {}).get("Status") == "DEMOLISH"
                    and review.get("ReviewStatus") == "CONFIRMED_APPROXIMATE"
                    and review.get("FormalIfcWriteAllowed") is True
                    and review.get("AccuracyBoundary")
                    == "POSITION_DIRECTION_LENGTH_APPROXIMATE_NOT_SURVEY_GRADE"
                ),
            }
        )

    status_counts = Counter(
        ifcopenshell.util.element.get_psets(wall)
        .get("Pset_WallCommon", {})
        .get("Status")
        for wall in formal.by_type("IfcWall")
    )
    gates = {
        "formal_wall_count": len(formal.by_type("IfcWall")),
        "formal_wall_status_counts": dict(status_counts),
        "demolition_wall_count": len(demolition_records),
        "demolition_global_ids_equal_register": {
            record["global_id"] for record in demolition_records
        }
        == set(candidate_ids),
        "demolition_records_pass": all(record["pass"] for record in demolition_records),
        "maximum_demolition_bbox_delta_mm": max(
            record["maximum_bbox_delta_mm"] for record in demolition_records
        ),
        "baseline_product_count": len(baseline_ids),
        "baseline_products_world_geometry_changes_over_tolerance": formal_to_baseline[
            "over_tolerance"
        ],
        "candidate_products_world_geometry_changes_over_tolerance": formal_to_candidate[
            "over_tolerance"
        ],
        "protected_products_max_world_geometry_change_mm": max(
            record.get("world_vertex_hausdorff_mm", 0.0)
            for record in protected_geometry["records"]
        ),
    }
    gates["pass"] = (
        gates["formal_wall_count"] == 101
        and gates["formal_wall_status_counts"]
        == {"EXISTING": 84, "NEW": 4, "DEMOLISH": 13}
        and gates["demolition_wall_count"] == 13
        and gates["demolition_global_ids_equal_register"]
        and gates["demolition_records_pass"]
        and gates["baseline_product_count"] == 284
        and gates["baseline_products_world_geometry_changes_over_tolerance"] == 0
        and gates["candidate_products_world_geometry_changes_over_tolerance"] == 0
        and gates["protected_products_max_world_geometry_change_mm"] == 0.0
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a102-demolition-postwrite-audit",
        "tolerance_mm": args.tolerance_mm,
        "formal": {"path": str(formal_path), "sha256": sha256(formal_path)},
        "approved_candidate": {
            "path": str(candidate_path),
            "sha256": sha256(candidate_path),
        },
        "baseline": baseline_source,
        "demolition_records": demolition_records,
        "formal_to_git_baseline_geometry": formal_to_baseline,
        "formal_to_approved_candidate_geometry": formal_to_candidate,
        "protected_products_geometry": protected_geometry,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
