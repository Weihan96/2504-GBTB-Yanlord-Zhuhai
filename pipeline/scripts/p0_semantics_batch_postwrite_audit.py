#!/usr/bin/env python3
"""Verify the formal IFC against the approved A-104/A-105 semantic candidate."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
import ifcopenshell.util.element
import numpy as np

import a104_semantics_candidate as a104
import a105_semantics_candidate as a105
from geometry_alignment_audit import geometry_difference_audit, sha256


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formal", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--prewrite-report", required=True, type=Path)
    parser.add_argument("--a104-register", required=True, type=Path)
    parser.add_argument("--a105-register", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    args = parse_args()
    if args.tolerance_mm <= 0:
        raise SystemExit("--tolerance-mm must be positive")

    formal_path = args.formal.resolve()
    candidate_path = args.candidate.resolve()
    prewrite = json.loads(args.prewrite_report.read_text(encoding="utf-8"))
    candidate_sha = sha256(candidate_path)
    if prewrite["candidate"]["sha256"] != candidate_sha or not prewrite["gates"]["mechanical_pass"]:
        raise RuntimeError("approved prewrite candidate evidence does not match the candidate file")

    formal = ifcopenshell.open(formal_path)
    candidate = ifcopenshell.open(candidate_path)
    a104_rows = read_csv(args.a104_register)
    a105_rows = read_csv(args.a105_register)
    protected_ids = list(prewrite["placement_deltas"])

    geometry = geometry_difference_audit(
        formal,
        candidate,
        str(candidate_path),
        tolerance_mm=args.tolerance_mm,
        classes=("IfcDoor", "IfcWindow", "IfcOpeningElement", "IfcCovering", "IfcSlab"),
        global_ids=(),
    )
    placement_deltas = {
        global_id: float(
            np.max(
                np.abs(
                    a105.placement_matrix(formal.by_guid(global_id))
                    - a105.placement_matrix(candidate.by_guid(global_id))
                )
            )
        )
        for global_id in protected_ids
    }
    groups = [a104.group_record(formal, key, spec) for key, spec in a104.GROUPS.items()]
    formal_root_ids = {root.GlobalId for root in formal.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}

    a104_tag_count = sum(
        formal.by_guid(row["global_id"]).Tag == row["candidate_tag"] for row in a104_rows
    )
    a104_name_count = sum(
        formal.by_guid(global_id).Name == name for global_id, name in a104.NAME_UPDATES.items()
    )
    a105_tag_count = sum(
        formal.by_guid(row["global_id"]).Tag == row["candidate_id"] for row in a105_rows
    )
    wet_pset_count = sum(
        "Pset_A105SlopeReview" in ifcopenshell.util.element.get_psets(formal.by_guid(row["global_id"]))
        for row in a105_rows
        if row["kind"] == "sloped_wet_tile"
    )
    finish_pset_count = sum(
        "Pset_A105FinishIntent" in ifcopenshell.util.element.get_psets(formal.by_guid(global_id))
        for global_id in a105.REFERENCE_MATERIALS
    )
    slab_pset_count = sum(
        "Pset_A105DepressedSlabReview" in ifcopenshell.util.element.get_psets(formal.by_guid(global_id))
        for global_id in a105.DEPRESSED_SLABS
    )
    material_count = sum(
        (material := ifcopenshell.util.element.get_material(formal.by_guid(global_id))) is not None
        and material.is_a("IfcMaterial")
        and material.Name == expected
        for global_id, expected in a105.REFERENCE_MATERIALS.items()
    )

    gates = {
        "schema_equal": formal.schema == candidate.schema == "IFC4",
        "entity_count_equal": len(list(formal)) == len(list(candidate)),
        "root_ids_equal": formal_root_ids == candidate_root_ids,
        "fills_relationship_ids_equal": a104.relation_ids(formal, "IfcRelFillsElement")
        == a104.relation_ids(candidate, "IfcRelFillsElement"),
        "voids_relationship_ids_equal": a104.relation_ids(formal, "IfcRelVoidsElement")
        == a104.relation_ids(candidate, "IfcRelVoidsElement"),
        "a104_tag_count": a104_tag_count,
        "a104_confirmed_name_count": a104_name_count,
        "a104_group_pass_count": sum(group["passes"] for group in groups),
        "a105_tag_count": a105_tag_count,
        "a105_wet_slope_pset_count": wet_pset_count,
        "a105_finish_intent_pset_count": finish_pset_count,
        "a105_depressed_slab_pset_count": slab_pset_count,
        "a105_reference_material_count": material_count,
        "maximum_placement_matrix_delta": max(placement_deltas.values(), default=0.0),
        "geometry_over_tolerance": geometry["over_tolerance"],
    }
    gates["mechanical_pass"] = (
        gates["schema_equal"]
        and gates["entity_count_equal"]
        and gates["root_ids_equal"]
        and gates["fills_relationship_ids_equal"]
        and gates["voids_relationship_ids_equal"]
        and gates["a104_tag_count"] == 19
        and gates["a104_confirmed_name_count"] == 4
        and gates["a104_group_pass_count"] == 3
        and gates["a105_tag_count"] == 21
        and gates["a105_wet_slope_pset_count"] == 18
        and gates["a105_finish_intent_pset_count"] == 3
        and gates["a105_depressed_slab_pset_count"] == 2
        and gates["a105_reference_material_count"] == 3
        and gates["maximum_placement_matrix_delta"] <= 1e-12
        and gates["geometry_over_tolerance"] == 0
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "formal-a104-a105-semantics-postwrite-audit",
        "formal": {"path": str(formal_path), "sha256": sha256(formal_path)},
        "approved_candidate": {"path": str(candidate_path), "sha256": candidate_sha},
        "tolerance_mm": args.tolerance_mm,
        "groups": groups,
        "placement_deltas": placement_deltas,
        "geometry_difference": geometry,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["mechanical_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
