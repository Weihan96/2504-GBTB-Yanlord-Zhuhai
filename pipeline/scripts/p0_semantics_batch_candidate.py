#!/usr/bin/env python3
"""Combine the confirmed A-104 and A-105 semantics into one atomic IFC candidate."""

from __future__ import annotations

import argparse
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
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--a104-register", required=True, type=Path)
    parser.add_argument("--a105-register", required=True, type=Path)
    parser.add_argument("--a105-report", required=True, type=Path)
    parser.add_argument("--build-up-report", required=True, type=Path)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.tolerance_mm <= 0:
        raise SystemExit("--tolerance-mm must be positive")
    source_path = args.input.resolve()
    output_path = args.output.resolve()
    if source_path == output_path:
        raise RuntimeError("candidate output must not overwrite the formal IFC")
    source_sha = sha256(source_path)
    a104_rows = a104.read_register(args.a104_register, source_sha)
    a104.require_decisions(args.decisions)
    a105.require_reports(args.a105_report, args.build_up_report, source_sha)
    a105.require_decisions(args.decisions)
    a105_rows = a105.read_floors(args.a105_register)

    source = ifcopenshell.open(source_path)
    protected_ids = (
        [row["global_id"] for row in a104_rows]
        + [row["global_id"] for row in a105_rows]
        + list(a105.DEPRESSED_SLABS)
        + list(a105.DEPRESSED_SLABS.values())
    )
    source_placements = {
        global_id: a105.placement_matrix(source.by_guid(global_id))
        for global_id in protected_ids
    }
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    source_fills = a104.relation_ids(source, "IfcRelFillsElement")
    source_voids = a104.relation_ids(source, "IfcRelVoidsElement")

    candidate = ifcopenshell.open(source_path)
    a104.apply_semantics(candidate, a104_rows)
    a105.apply_semantics(candidate, a105_rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(output_path)
    candidate = ifcopenshell.open(output_path)

    geometry = geometry_difference_audit(
        candidate,
        source,
        str(source_path),
        tolerance_mm=args.tolerance_mm,
        classes=("IfcDoor", "IfcWindow", "IfcOpeningElement", "IfcCovering", "IfcSlab"),
        global_ids=(),
    )
    placement_deltas = {
        global_id: float(np.max(np.abs(a105.placement_matrix(candidate.by_guid(global_id)) - matrix)))
        for global_id, matrix in source_placements.items()
    }
    a104_groups = [a104.group_record(candidate, key, spec) for key, spec in a104.GROUPS.items()]
    a104_tag_count = sum(candidate.by_guid(row["global_id"]).Tag == row["candidate_tag"] for row in a104_rows)
    a104_name_count = sum(candidate.by_guid(global_id).Name == name for global_id, name in a104.NAME_UPDATES.items())
    a105_tag_count = sum(candidate.by_guid(row["global_id"]).Tag == row["candidate_id"] for row in a105_rows)
    a105_wet_psets = sum(
        "Pset_A105SlopeReview" in ifcopenshell.util.element.get_psets(candidate.by_guid(row["global_id"]))
        for row in a105_rows if row["kind"] == "sloped_wet_tile"
    )
    a105_finish_psets = sum(
        "Pset_A105FinishIntent" in ifcopenshell.util.element.get_psets(candidate.by_guid(global_id))
        for global_id in a105.REFERENCE_MATERIALS
    )
    a105_slab_psets = sum(
        "Pset_A105DepressedSlabReview" in ifcopenshell.util.element.get_psets(candidate.by_guid(global_id))
        for global_id in a105.DEPRESSED_SLABS
    )
    a105_material_count = sum(
        (material := ifcopenshell.util.element.get_material(candidate.by_guid(global_id))) is not None
        and material.is_a("IfcMaterial")
        and material.Name == expected
        for global_id, expected in a105.REFERENCE_MATERIALS.items()
    )
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "schema_equal": source.schema == candidate.schema == "IFC4",
        "a104_tag_count": a104_tag_count,
        "a104_confirmed_name_count": a104_name_count,
        "a104_group_pass_count": sum(group["passes"] for group in a104_groups),
        "a105_tag_count": a105_tag_count,
        "a105_wet_slope_pset_count": a105_wet_psets,
        "a105_finish_intent_pset_count": a105_finish_psets,
        "a105_depressed_slab_pset_count": a105_slab_psets,
        "a105_reference_material_count": a105_material_count,
        "geometry_over_tolerance": geometry["over_tolerance"],
        "maximum_placement_matrix_delta": max(placement_deltas.values()),
        "fills_relationship_ids_equal": source_fills == a104.relation_ids(candidate, "IfcRelFillsElement"),
        "voids_relationship_ids_equal": source_voids == a104.relation_ids(candidate, "IfcRelVoidsElement"),
        "source_root_ids_preserved": source_root_ids <= candidate_root_ids,
        "new_root_count": len(candidate_root_ids - source_root_ids),
        "formal_ifc_write_allowed": False,
    }
    gates["mechanical_pass"] = (
        gates["schema_equal"]
        and gates["a104_tag_count"] == 19
        and gates["a104_confirmed_name_count"] == 4
        and gates["a104_group_pass_count"] == 3
        and gates["a105_tag_count"] == 21
        and gates["a105_wet_slope_pset_count"] == 18
        and gates["a105_finish_intent_pset_count"] == 3
        and gates["a105_depressed_slab_pset_count"] == 2
        and gates["a105_reference_material_count"] == 3
        and gates["geometry_over_tolerance"] == 0
        and gates["maximum_placement_matrix_delta"] <= 1e-12
        and gates["fills_relationship_ids_equal"]
        and gates["voids_relationship_ids_equal"]
        and gates["source_root_ids_preserved"]
        and gates["new_root_count"] == 54
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a104-a105-atomic-semantics-candidate",
        "source": {"path": str(source_path), "sha256": source_sha, "schema": source.schema},
        "candidate": {"path": str(output_path), "sha256": sha256(output_path), "schema": candidate.schema},
        "tolerance_mm": args.tolerance_mm,
        "a104_groups": a104_groups,
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
