#!/usr/bin/env python3
"""Build a two-object installation-face anchor candidate for Geberit 115.770 plates."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import numpy as np

from beam_origin_reset_candidate import point_mesh_distance, reset_product_origin
from direction_noise_candidate import all_product_geometry_difference
from electric_appliance_anchor_options import (
    CONTEXT_CLASSES,
    compare_option,
    contact_state,
    json_default,
)
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm
from light_fixture_anchor_candidate import apply_targets, renderable_context


TARGET_RULES = {
    "2gFgcOYEXEaQWAzcKulTFt": {"face_axis": 0, "face_side": "max"},
    "2lDPsdQevFSfeOThtjSlPG": {"face_axis": 0, "face_side": "min"},
}
EXPECTED_TYPE_NAME = "Geberit 115.770"
EXPECTED_TYPE_PREDEFINED = "WCSEAT"
EXPECTED_SIZE_SORTED_MM = np.array([12.445110321, 170.17956543, 254.0])


def plate_evidence(
    model: ifcopenshell.file,
    global_id: str,
    rule: dict[str, Any],
    evidence_tolerance_mm: float,
) -> dict[str, Any]:
    product = model.by_guid(global_id)
    assigned_type = ifcopenshell.util.element.get_type(product)
    container = ifcopenshell.util.element.get_container(product)
    if (
        product is None
        or not product.is_a("IfcSanitaryTerminal")
        or getattr(assigned_type, "Name", None) != EXPECTED_TYPE_NAME
        or getattr(assigned_type, "PredefinedType", None)
        != EXPECTED_TYPE_PREDEFINED
        or getattr(container, "Name", None) != "WC"
    ):
        raise RuntimeError(f"flush-plate identity mismatch: {global_id}")

    vertices, faces = world_mesh_mm(geometry_settings(), product)
    points = np.array(vertices, dtype=float)
    bbox_min, bbox_max = points.min(axis=0), points.max(axis=0)
    size = bbox_max - bbox_min
    if np.max(np.abs(np.sort(size) - EXPECTED_SIZE_SORTED_MM)) > 0.01:
        raise RuntimeError(f"flush-plate size mismatch: {global_id}")

    face_centre = (bbox_min + bbox_max) / 2.0
    axis = rule["face_axis"]
    face_centre[axis] = (
        bbox_min[axis] if rule["face_side"] == "min" else bbox_max[axis]
    )
    surface_distance = point_mesh_distance(face_centre, vertices, faces)
    if surface_distance > evidence_tolerance_mm:
        raise RuntimeError(
            f"installation-face centre is not on rendered mesh: {global_id}"
        )
    target = np.round(face_centre)
    translation = target - face_centre
    return {
        "global_id": global_id,
        "ifc_class": product.is_a(),
        "type_name": assigned_type.Name,
        "stored_predefined_type": assigned_type.PredefinedType,
        "semantic_note": (
            "product code and thin-plate geometry identify a WC actuator plate; "
            "this candidate changes placement only and does not change IFC semantics"
        ),
        "container": container.Name,
        "face_axis": axis,
        "face_side": rule["face_side"],
        "bbox_size_mm": size.tolist(),
        "source_face_centre_mm": face_centre.tolist(),
        "target_mm": target.tolist(),
        "translation_mm": translation.tolist(),
        "maximum_axis_shift_mm": float(np.max(np.abs(translation))),
        "euclidean_shift_mm": float(np.linalg.norm(translation)),
        "source_surface_distance_mm": surface_distance,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--maximum-axis-shift-mm", type=float, default=0.5)
    parser.add_argument("--clearance-window-mm", type=float, default=2.0)
    parser.add_argument("--evidence-tolerance-mm", type=float, default=0.01)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    evidence = [
        plate_evidence(source, global_id, rule, args.evidence_tolerance_mm)
        for global_id, rule in TARGET_RULES.items()
    ]
    if any(
        row["maximum_axis_shift_mm"] > args.maximum_axis_shift_mm
        for row in evidence
    ):
        raise RuntimeError("flush-plate installation-face shift exceeds maximum")

    target_ids = set(TARGET_RULES)
    source_targets = [source.by_guid(global_id) for global_id in target_ids]
    source_context, _ = renderable_context(source, CONTEXT_CLASSES)
    source_context = [p for p in source_context if p.GlobalId not in target_ids]
    source_state = contact_state(
        source_targets,
        source_context,
        args.clearance_window_mm,
        args.tolerance_mm,
    )

    candidate = ifcopenshell.open(source_path)
    reset_results = [
        reset_product_origin(
            candidate,
            {
                "global_id": row["global_id"],
                "target_mm": np.array(row["source_face_centre_mm"]),
                "anchor_kind": "existing_installation_face_centre",
                "basis": "temporary geometry-preserving installation-face reset",
                "confidence": 1.0,
            },
            args.tolerance_mm,
            "IfcSanitaryTerminal",
        )
        for row in evidence
    ]
    rigid_results = apply_targets(
        candidate,
        [
            {
                **row,
                "anchor_kind": "integer_installation_face_centre",
                "basis": (
                    "Geberit 115.770 thin actuator plate; mounting-face centre "
                    "is the shortest directly measurable installation datum"
                ),
                "confidence": 1.0,
            }
            for row in evidence
        ],
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)

    geometry = all_product_geometry_difference(candidate, source, args.tolerance_mm)
    geometry_by_id = {row["global_id"]: row for row in geometry["records"]}
    target_checks = []
    for row, rigid in zip(evidence, rigid_results):
        vertices, faces = world_mesh_mm(
            geometry_settings(), candidate.by_guid(row["global_id"])
        )
        geometry_row = geometry_by_id[row["global_id"]]
        target_checks.append(
            {
                "global_id": row["global_id"],
                "placement_residual_mm": rigid["placement_residual_mm"],
                "surface_distance_mm": point_mesh_distance(
                    row["target_mm"], vertices, faces
                ),
                "expected_shift_mm": row["euclidean_shift_mm"],
                "actual_shift_mm": geometry_row[
                    "world_corresponding_vertex_max_delta_mm"
                ],
                "mesh_counts_equal": geometry_row["mesh_counts_equal"],
            }
        )

    comparison = compare_option(
        source,
        candidate,
        source_state,
        target_ids,
        args.tolerance_mm,
        args.clearance_window_mm,
    )
    source_roots = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_roots = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "targets": len(target_checks),
        "targets_over_placement_tolerance": sum(
            row["placement_residual_mm"] > args.tolerance_mm
            for row in target_checks
        ),
        "targets_over_surface_tolerance": sum(
            row["surface_distance_mm"] > args.tolerance_mm
            for row in target_checks
        ),
        "world_shift_mismatches": sum(
            abs(row["actual_shift_mm"] - row["expected_shift_mm"])
            > args.tolerance_mm
            for row in target_checks
        ),
        "mesh_count_mismatches": sum(
            not row["mesh_counts_equal"] for row in target_checks
        ),
        "non_target_geometry_changes": [
            row["global_id"]
            for row in geometry["records"]
            if row["global_id"] not in target_ids and not row["within_tolerance"]
        ],
        "contact_regressions": comparison["internal_contact_changes"][
            "regressions"
        ]
        + comparison["context_contact_changes"]["regressions"],
        "unexpected_new_contacts": comparison["internal_contact_changes"][
            "unexpected_new_contacts"
        ]
        + comparison["context_contact_changes"]["unexpected_new_contacts"],
        "new_intersections": comparison["new_internal_intersections"]
        + comparison["new_context_intersections"],
        "lost_intersections": comparison["lost_internal_intersections"]
        + comparison["lost_context_intersections"],
        "schema_equal": source.schema == candidate.schema,
        "entity_count_equal": len(list(source)) == len(list(candidate)),
        "root_global_ids_equal": source_roots == candidate_roots,
    }
    gates["pass"] = all(
        [
            gates["targets_over_placement_tolerance"] == 0,
            gates["targets_over_surface_tolerance"] == 0,
            gates["world_shift_mismatches"] == 0,
            gates["mesh_count_mismatches"] == 0,
            not gates["non_target_geometry_changes"],
            not gates["contact_regressions"],
            not gates["unexpected_new_contacts"],
            not gates["new_intersections"],
            not gates["lost_intersections"],
            gates["schema_equal"],
            gates["entity_count_equal"],
            gates["root_global_ids_equal"],
        ]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-geberit-115770-installation-face-candidate",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "candidate": {
            "path": str(args.output.resolve()),
            "sha256": sha256(args.output),
        },
        "tolerance_mm": args.tolerance_mm,
        "maximum_axis_shift_mm": args.maximum_axis_shift_mm,
        "evidence": evidence,
        "reset_results": reset_results,
        "rigid_results": rigid_results,
        "target_checks": target_checks,
        "contact_comparison": comparison,
        "gates": gates,
        "semantic_write_performed": False,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, default=json_default),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "report": str(args.report),
                "candidate": str(args.output),
                "targets": len(evidence),
                "maximum_axis_shift_mm": max(
                    row["maximum_axis_shift_mm"] for row in evidence
                ),
                "pass": gates["pass"],
            },
            ensure_ascii=False,
        )
    )
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
