#!/usr/bin/env python3
"""Build a no-world-movement origin-reset candidate for shared Opening groups."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.placement
import numpy as np

from beam_origin_reset_candidate import point_mesh_distance, transform_body_item
from direction_noise_candidate import all_product_geometry_difference
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm


def parse_ids(value: str) -> set[str]:
    return {part.strip() for part in value.split(";") if part.strip()}


def read_targets(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    expected = {
        "group_id",
        "representative_global_id",
        "member_global_ids",
        "target_x_mm",
        "target_y_mm",
        "target_z_mm",
        "anchor_kind",
        "basis",
        "confidence",
    }
    if not rows or set(rows[0]) != expected:
        raise RuntimeError(f"invalid Opening group target schema: {path}")
    result = []
    for row in rows:
        target = np.array(
            [
                float(row["target_x_mm"]),
                float(row["target_y_mm"]),
                float(row["target_z_mm"]),
            ],
            dtype=float,
        )
        if any(abs(value - round(value)) > 1e-9 for value in target):
            raise RuntimeError(f"Opening group target is not integer millimetres: {row}")
        members = parse_ids(row["member_global_ids"])
        if row["representative_global_id"] not in members:
            raise RuntimeError(f"representative is not a group member: {row['group_id']}")
        result.append(
            {
                "group_id": row["group_id"],
                "representative_global_id": row["representative_global_id"],
                "member_global_ids": members,
                "target_mm": target,
                "anchor_kind": row["anchor_kind"],
                "basis": row["basis"],
                "confidence": float(row["confidence"]),
            }
        )
    if len({row["group_id"] for row in result}) != len(result):
        raise RuntimeError("Opening group target table contains duplicate group IDs")
    all_members = [member for row in result for member in row["member_global_ids"]]
    if len(set(all_members)) != len(all_members):
        raise RuntimeError("Opening group target table contains duplicate members")
    return result


def shared_group_members(
    model: ifcopenshell.file, representative: ifcopenshell.entity_instance
) -> set[str]:
    placement_users = {
        user
        for user in model.get_inverse(representative.ObjectPlacement)
        if user.is_a("IfcOpeningElement")
    }
    representation_users = {
        user
        for user in model.get_inverse(representative.Representation)
        if user.is_a("IfcOpeningElement")
    }
    if placement_users != representation_users:
        raise RuntimeError(
            f"Opening {representative.GlobalId} does not have one exact shared placement/representation group"
        )
    return {member.GlobalId for member in placement_users}


def reset_group_origin(
    model: ifcopenshell.file, row: dict[str, Any], tolerance_mm: float
) -> dict[str, Any]:
    representative = model.by_guid(row["representative_global_id"])
    if representative is None or not representative.is_a("IfcOpeningElement"):
        raise RuntimeError(f"missing representative Opening {row['representative_global_id']}")
    actual_members = shared_group_members(model, representative)
    if actual_members != row["member_global_ids"]:
        raise RuntimeError(
            f"Opening group {row['group_id']} members differ: "
            f"expected {sorted(row['member_global_ids'])}, found {sorted(actual_members)}"
        )
    members = [model.by_guid(global_id) for global_id in sorted(actual_members)]
    source_matrices = [
        np.array(
            ifcopenshell.util.placement.get_local_placement(member.ObjectPlacement),
            dtype=float,
        )
        for member in members
    ]
    source_matrix = source_matrices[0]
    if any(not np.allclose(matrix, source_matrix, atol=1e-9) for matrix in source_matrices):
        raise RuntimeError(f"Opening group {row['group_id']} does not share one world matrix")

    target_matrix = source_matrix.copy()
    target_matrix[:3, 3] = row["target_mm"]
    representation_transform = np.linalg.inv(target_matrix) @ source_matrix
    if not np.allclose(
        representation_transform[:3, :3], np.identity(3), atol=1e-9
    ):
        raise RuntimeError(f"Opening group {row['group_id']} reset unexpectedly rotates geometry")

    placement = representative.ObjectPlacement
    parent_matrix = np.identity(4)
    if placement.PlacementRelTo:
        parent_matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(placement.PlacementRelTo),
            dtype=float,
        )
    target_local_matrix = np.linalg.inv(parent_matrix) @ target_matrix
    current_local_matrix = np.array(
        ifcopenshell.util.placement.get_axis2placement(placement.RelativePlacement),
        dtype=float,
    )
    if not np.allclose(
        target_local_matrix[:3, :3], current_local_matrix[:3, :3], atol=1e-9
    ):
        raise RuntimeError(f"Opening group {row['group_id']} local orientation would change")

    transformed_entity_ids: set[int] = set()
    allowed_entity_ids = {
        entity.id()
        for representation in representative.Representation.Representations
        for entity in model.traverse(representation)
    }
    identifiers = set()
    for representation in representative.Representation.Representations:
        identifier = representation.RepresentationIdentifier
        identifiers.add(identifier)
        if identifier not in {"Body", "Box"}:
            raise RuntimeError(
                f"unsupported shared Opening representation {identifier}"
            )
        for item in representation.Items:
            transform_body_item(
                model,
                item,
                representation_transform,
                transformed_entity_ids,
                allowed_entity_ids,
            )
    if identifiers != {"Body", "Box"}:
        raise RuntimeError(f"Opening group {row['group_id']} requires Body and Box")
    placement.RelativePlacement.Location.Coordinates = tuple(
        float(value) for value in target_local_matrix[:3, 3]
    )

    settings = geometry_settings()
    anchor_distances = {}
    placement_residuals = {}
    hosts = set()
    for member in members:
        matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(member.ObjectPlacement),
            dtype=float,
        )
        placement_residuals[member.GlobalId] = float(
            np.max(np.abs(matrix[:3, 3] - row["target_mm"]))
        )
        vertices, faces = world_mesh_mm(settings, member)
        anchor_distances[member.GlobalId] = point_mesh_distance(
            row["target_mm"], vertices, faces
        )
        hosts.update(
            relation.RelatingBuildingElement.GlobalId
            for relation in member.VoidsElements
        )
    if max(placement_residuals.values()) > tolerance_mm:
        raise RuntimeError(f"Opening group {row['group_id']} placement residual exceeds tolerance")
    if max(anchor_distances.values()) > tolerance_mm:
        raise RuntimeError(f"Opening group {row['group_id']} target is not on its geometry")
    return {
        "group_id": row["group_id"],
        "member_global_ids": sorted(actual_members),
        "host_global_ids": sorted(hosts),
        "source_mm": source_matrix[:3, 3].tolist(),
        "target_mm": row["target_mm"].tolist(),
        "anchor_kind": row["anchor_kind"],
        "basis": row["basis"],
        "confidence": row["confidence"],
        "placement_residuals_mm": placement_residuals,
        "anchor_surface_distances_mm": anchor_distances,
    }


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
    results = [reset_group_origin(candidate, row, args.tolerance_mm) for row in rows]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    geometry = all_product_geometry_difference(candidate, source, args.tolerance_mm)
    geometry_over_tolerance = [
        record for record in geometry["records"] if not record["within_tolerance"]
    ]
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "groups": len(results),
        "openings": sum(len(result["member_global_ids"]) for result in results),
        "hosts": len({host for result in results for host in result["host_global_ids"]}),
        "product_geometry_over_tolerance": len(geometry_over_tolerance),
        "schema_equal": source.schema == candidate.schema,
        "entity_count_equal": len(list(source)) == len(list(candidate)),
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
    }
    gates["pass"] = (
        gates["product_geometry_over_tolerance"] == 0
        and gates["schema_equal"]
        and gates["entity_count_equal"]
        and gates["root_global_ids_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-shared-opening-origin-reset-candidate",
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
        "all_product_geometry": geometry,
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
