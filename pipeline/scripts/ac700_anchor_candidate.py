#!/usr/bin/env python3
"""Build mechanically equivalent AC700 upper/lower face-anchor candidates."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.util.placement
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


AC700_GLOBAL_ID = "1PUCikoaP5fgiYt8sJd8$6"
AC700_TYPE_NAME = "AC700"
AC700_DESCRIPTION = "RPIZ-22FSLN5QD/P 700x447x192"


def face_anchor_evidence(
    model: ifcopenshell.file,
    side: str,
    evidence_tolerance_mm: float,
) -> dict[str, Any]:
    if side not in {"lower", "upper"}:
        raise ValueError(f"unsupported AC700 face side: {side}")
    product = model.by_guid(AC700_GLOBAL_ID)
    assigned_type = ifcopenshell.util.element.get_type(product)
    if (
        product is None
        or not product.is_a("IfcElectricAppliance")
        or getattr(assigned_type, "Name", None) != AC700_TYPE_NAME
        or getattr(assigned_type, "Description", None) != AC700_DESCRIPTION
    ):
        raise RuntimeError("AC700 target identity does not match the exact rule")
    settings = geometry_settings()
    vertices, faces = world_mesh_mm(settings, product)
    points = np.array(vertices, dtype=float)
    matrix = np.array(
        ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
        dtype=float,
    )
    homogeneous = np.column_stack((points, np.ones(len(points), dtype=float)))
    local = (np.linalg.inv(matrix) @ homogeneous.T).T[:, :3]
    bbox_min = local.min(axis=0)
    bbox_max = local.max(axis=0)
    local_face_centre = (bbox_min + bbox_max) / 2.0
    local_face_centre[2] = bbox_min[2] if side == "lower" else bbox_max[2]
    source_face_centre = (matrix @ np.append(local_face_centre, 1.0))[:3]
    target = np.round(source_face_centre)
    source_surface_distance = point_mesh_distance(
        source_face_centre, vertices, faces
    )
    if source_surface_distance > evidence_tolerance_mm:
        raise RuntimeError(
            f"AC700 {side} bbox face centre is not on the rendered mesh: "
            f"{source_surface_distance} mm"
        )
    return {
        "side": side,
        "global_id": AC700_GLOBAL_ID,
        "type_name": AC700_TYPE_NAME,
        "type_description": AC700_DESCRIPTION,
        "source_face_centre_mm": source_face_centre.tolist(),
        "target_mm": target.tolist(),
        "translation_mm": (target - source_face_centre).tolist(),
        "maximum_axis_shift_mm": float(
            np.max(np.abs(target - source_face_centre))
        ),
        "euclidean_shift_mm": float(
            np.linalg.norm(target - source_face_centre)
        ),
        "source_surface_distance_mm": source_surface_distance,
    }


def build_option(
    source_path: Path,
    evidence: dict[str, Any],
    output: Path,
    tolerance_mm: float,
    maximum_axis_shift_mm: float,
    clearance_window_mm: float,
    source_state: dict[str, Any],
) -> dict[str, Any]:
    if evidence["maximum_axis_shift_mm"] > maximum_axis_shift_mm:
        raise RuntimeError(
            f"AC700 {evidence['side']} shift exceeds maximum axis shift"
        )
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    reset_result = reset_product_origin(
        candidate,
        {
            "global_id": AC700_GLOBAL_ID,
            "target_mm": np.array(
                evidence["source_face_centre_mm"], dtype=float
            ),
            "anchor_kind": f"existing_{evidence['side']}_face_centre",
            "basis": "temporary geometry-preserving face-centre origin reset",
            "confidence": 1.0,
        },
        tolerance_mm,
        "IfcElectricAppliance",
    )
    rigid_result = apply_targets(
        candidate,
        [
            {
                "global_id": AC700_GLOBAL_ID,
                "target_mm": evidence["target_mm"],
                "translation_mm": evidence["translation_mm"],
                "anchor_kind": f"integer_{evidence['side']}_face_centre",
                "basis": (
                    f"AC700 {evidence['side']} installation face centre; "
                    "semantic choice requires human confirmation"
                ),
                "confidence": 0.95,
            }
        ],
    )[0]
    output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(output)
    candidate = ifcopenshell.open(output)
    all_geometry = all_product_geometry_difference(
        candidate, source, tolerance_mm
    )
    geometry_by_id = {
        row["global_id"]: row for row in all_geometry["records"]
    }
    target_geometry = geometry_by_id[AC700_GLOBAL_ID]
    expected_shift = evidence["euclidean_shift_mm"]
    actual_shift = target_geometry["world_corresponding_vertex_max_delta_mm"]
    candidate_product = candidate.by_guid(AC700_GLOBAL_ID)
    vertices, faces = world_mesh_mm(
        geometry_settings(), candidate_product
    )
    target_surface_distance = point_mesh_distance(
        evidence["target_mm"], vertices, faces
    )
    target_ids = {AC700_GLOBAL_ID}
    comparison = compare_option(
        source,
        candidate,
        source_state,
        target_ids,
        tolerance_mm,
        clearance_window_mm,
    )
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {
        root.GlobalId for root in candidate.by_type("IfcRoot")
    }
    non_target_geometry_changes = [
        row["global_id"]
        for row in all_geometry["records"]
        if row["global_id"] != AC700_GLOBAL_ID and not row["within_tolerance"]
    ]
    gates = {
        "target_placement_residual_mm": rigid_result[
            "placement_residual_mm"
        ],
        "target_surface_distance_mm": target_surface_distance,
        "target_mesh_counts_equal": target_geometry["mesh_counts_equal"],
        "expected_world_shift_mm": expected_shift,
        "actual_world_shift_mm": actual_shift,
        "world_shift_matches_expected": abs(actual_shift - expected_shift)
        <= tolerance_mm,
        "non_target_geometry_changes": non_target_geometry_changes,
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
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
    }
    gates["pass"] = (
        gates["target_placement_residual_mm"] <= tolerance_mm
        and gates["target_surface_distance_mm"] <= tolerance_mm
        and gates["target_mesh_counts_equal"]
        and gates["world_shift_matches_expected"]
        and not gates["non_target_geometry_changes"]
        and not gates["contact_regressions"]
        and not gates["unexpected_new_contacts"]
        and not gates["new_intersections"]
        and not gates["lost_intersections"]
        and gates["schema_equal"]
        and gates["entity_count_equal"]
        and gates["root_global_ids_equal"]
    )
    return {
        "path": str(output.resolve()),
        "sha256": sha256(output),
        "evidence": evidence,
        "reset_result": reset_result,
        "rigid_result": rigid_result,
        "target_geometry": target_geometry,
        "contact_comparison": comparison,
        "gates": gates,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
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
    evidence = {
        side: face_anchor_evidence(
            source, side, args.evidence_tolerance_mm
        )
        for side in ("lower", "upper")
    }
    target = source.by_guid(AC700_GLOBAL_ID)
    source_context, _ = renderable_context(source, CONTEXT_CLASSES)
    source_context = [
        product
        for product in source_context
        if product.GlobalId != AC700_GLOBAL_ID
    ]
    source_state = contact_state(
        [target],
        source_context,
        args.clearance_window_mm,
        args.tolerance_mm,
    )
    options = {}
    for side in ("lower", "upper"):
        option = build_option(
            source_path,
            evidence[side],
            args.output_dir / f"ac700-{side}-face-anchor.ifc",
            args.tolerance_mm,
            args.maximum_axis_shift_mm,
            args.clearance_window_mm,
            source_state,
        )
        options[side] = option
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / f"ac700-{side}-face-anchor.json").write_text(
            json.dumps(
                {
                    "source_sha256": sha256(source_path),
                    "side": side,
                    "option": option,
                },
                ensure_ascii=False,
                indent=2,
                default=json_default,
            ),
            encoding="utf-8",
        )

    lower = ifcopenshell.open(options["lower"]["path"])
    upper = ifcopenshell.open(options["upper"]["path"])
    lower_vertices, lower_faces = world_mesh_mm(
        geometry_settings(), lower.by_guid(AC700_GLOBAL_ID)
    )
    upper_vertices, upper_faces = world_mesh_mm(
        geometry_settings(), upper.by_guid(AC700_GLOBAL_ID)
    )
    cross_option = {
        "mesh_counts_equal": len(lower_vertices) == len(upper_vertices)
        and len(lower_faces) == len(upper_faces),
        "maximum_corresponding_vertex_delta_mm": float(
            np.max(
                np.linalg.norm(
                    np.array(lower_vertices, dtype=float)
                    - np.array(upper_vertices, dtype=float),
                    axis=1,
                )
            )
        ),
    }
    cross_option["geometry_equal"] = (
        cross_option["mesh_counts_equal"]
        and cross_option["maximum_corresponding_vertex_delta_mm"]
        <= args.tolerance_mm
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-ac700-face-anchor-candidates",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "maximum_axis_shift_mm": args.maximum_axis_shift_mm,
        "clearance_window_mm": args.clearance_window_mm,
        "options": options,
        "cross_option": cross_option,
        "human_decision_required": True,
        "human_decision": "choose lower or upper AC700 installation face anchor",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2,
            default=json_default,
        ),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "report": str(args.report),
                "lower_pass": options["lower"]["gates"]["pass"],
                "upper_pass": options["upper"]["gates"]["pass"],
                "cross_option": cross_option,
                "translation_mm": evidence["lower"]["translation_mm"],
                "human_decision_required": True,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
