#!/usr/bin/env python3
"""Reset shared-Opening Covering hosts without breaking shared placements."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.placement
import numpy as np

from beam_origin_reset_candidate import point_mesh_distance, reset_product_origin
from direction_noise_candidate import all_product_geometry_difference
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm
from structural_origin_reset_candidate import (
    read_targets,
    scoped_owner_history_additions,
)


def related_openings(
    model: ifcopenshell.file, target_global_ids: set[str]
) -> list[ifcopenshell.entity_instance]:
    openings: dict[str, ifcopenshell.entity_instance] = {}
    for global_id in target_global_ids:
        host = model.by_guid(global_id)
        if host is None:
            raise RuntimeError(f"missing shared-Opening host {global_id}")
        for relation in host.HasOpenings:
            opening = relation.RelatedOpeningElement
            openings[opening.GlobalId] = opening
    return [openings[global_id] for global_id in sorted(openings)]


def capture_placement_groups(
    model: ifcopenshell.file,
    openings: list[ifcopenshell.entity_instance],
) -> dict[int, dict[str, Any]]:
    groups: dict[int, dict[str, Any]] = {}
    scoped_ids = {opening.GlobalId for opening in openings}
    for opening in openings:
        placement = opening.ObjectPlacement
        if placement is None or not placement.is_a("IfcLocalPlacement"):
            raise RuntimeError(f"Opening {opening.GlobalId} lacks IfcLocalPlacement")
        if placement.id() in groups:
            continue
        members = sorted(
            entity.GlobalId
            for entity in model.get_inverse(placement)
            if entity.is_a("IfcOpeningElement")
            and getattr(entity, "GlobalId", None) in scoped_ids
        )
        if not members:
            raise RuntimeError(f"Opening placement {placement.id()} has no scoped users")
        world_matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(placement),
            dtype=float,
        )
        local_matrix = np.array(
            ifcopenshell.util.placement.get_axis2placement(
                placement.RelativePlacement
            ),
            dtype=float,
        )
        groups[placement.id()] = {
            "placement": placement,
            "member_global_ids": members,
            "source_world_matrix": world_matrix,
            "source_local_rotation": local_matrix[:3, :3].copy(),
        }
    return groups


def restore_placement_groups(
    groups: dict[int, dict[str, Any]],
    tolerance_mm: float,
) -> list[dict[str, Any]]:
    records = []
    for placement_id, group in sorted(groups.items()):
        placement = group["placement"]
        parent_matrix = np.identity(4)
        if placement.PlacementRelTo:
            parent_matrix = np.array(
                ifcopenshell.util.placement.get_local_placement(
                    placement.PlacementRelTo
                ),
                dtype=float,
            )
        target_local = (
            np.linalg.inv(parent_matrix) @ group["source_world_matrix"]
        )
        if not np.allclose(
            target_local[:3, :3],
            group["source_local_rotation"],
            atol=1e-9,
        ):
            raise RuntimeError(
                f"Opening placement group {placement_id} would rotate"
            )
        placement.RelativePlacement.Location.Coordinates = tuple(
            float(value) for value in target_local[:3, 3]
        )
        restored_world = np.array(
            ifcopenshell.util.placement.get_local_placement(placement),
            dtype=float,
        )
        residual = float(
            np.max(
                np.abs(
                    restored_world[:3, 3]
                    - group["source_world_matrix"][:3, 3]
                )
            )
        )
        if residual > tolerance_mm:
            raise RuntimeError(
                f"Opening placement group {placement_id} residual {residual} mm"
            )
        records.append(
            {
                "source_placement_id": placement_id,
                "member_global_ids": group["member_global_ids"],
                "world_placement_residual_mm": residual,
            }
        )
    return records


def shared_group_signature(
    model: ifcopenshell.file,
    member_global_ids: list[str],
) -> dict[str, Any]:
    members = [model.by_guid(global_id) for global_id in member_global_ids]
    placement_count = len({member.ObjectPlacement.id() for member in members})
    representation_count = len(
        {member.Representation.id() for member in members}
    )
    return {
        "member_global_ids": sorted(member_global_ids),
        "shared_placement_count": placement_count,
        "shared_representation_count": representation_count,
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
    if {row["expected_class"] for row in rows} != {"IfcCovering"}:
        raise RuntimeError("shared surface host batch accepts only IfcCovering")
    target_ids = {row["global_id"] for row in rows}
    openings = related_openings(candidate, target_ids)
    placement_groups = capture_placement_groups(candidate, openings)
    source_group_signatures = [
        shared_group_signature(source, group["member_global_ids"])
        for group in placement_groups.values()
    ]
    results = [
        reset_product_origin(
            candidate,
            row,
            args.tolerance_mm,
            row["expected_class"],
            restore_related_openings=False,
        )
        for row in rows
    ]
    placement_group_results = restore_placement_groups(
        placement_groups, args.tolerance_mm
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    geometry = all_product_geometry_difference(
        candidate, source, args.tolerance_mm
    )
    geometry_over_tolerance = [
        record for record in geometry["records"] if not record["within_tolerance"]
    ]
    placement_mismatches = []
    final_anchor_distances: dict[str, float] = {}
    settings = geometry_settings()
    for row in rows:
        product = candidate.by_guid(row["global_id"])
        matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(
                product.ObjectPlacement
            ),
            dtype=float,
        )
        residual = float(np.max(np.abs(matrix[:3, 3] - row["target_mm"])))
        if residual > args.tolerance_mm:
            placement_mismatches.append(
                {"global_id": product.GlobalId, "residual_mm": residual}
            )
        vertices, faces = world_mesh_mm(settings, product)
        final_anchor_distances[product.GlobalId] = point_mesh_distance(
            row["target_mm"], vertices, faces
        )
    candidate_group_signatures = [
        shared_group_signature(candidate, group["member_global_ids"])
        for group in placement_groups.values()
    ]
    shared_groups_preserved = (
        source_group_signatures == candidate_group_signatures
    )
    related_opening_ids = sorted(
        {opening.GlobalId for opening in openings}
    )
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {
        root.GlobalId for root in candidate.by_type("IfcRoot")
    }
    owner_history_additions = scoped_owner_history_additions(
        source, candidate, target_ids
    )
    representation_additions = sum(
        result["representation_entities_added"] for result in results
    )
    entity_count_delta = len(list(candidate)) - len(list(source))
    expected_entity_count_delta = (
        representation_additions + owner_history_additions["count"]
    )
    gates = {
        "target_coverings": len(rows),
        "related_openings": len(related_opening_ids),
        "opening_placement_groups": len(placement_groups),
        "target_placement_mismatches": placement_mismatches,
        "anchors_over_tolerance": sum(
            distance > args.tolerance_mm
            for distance in final_anchor_distances.values()
        ),
        "opening_group_placement_residual_max_mm": max(
            record["world_placement_residual_mm"]
            for record in placement_group_results
        ),
        "shared_groups_preserved": shared_groups_preserved,
        "product_geometry_over_tolerance": len(geometry_over_tolerance),
        "schema_equal": source.schema == candidate.schema,
        "entity_count_delta": entity_count_delta,
        "expected_entity_count_delta": expected_entity_count_delta,
        "entity_count_delta_matches_expected": (
            entity_count_delta == expected_entity_count_delta
        ),
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
        "owner_history_additions": owner_history_additions,
    }
    gates["pass"] = (
        not gates["target_placement_mismatches"]
        and gates["anchors_over_tolerance"] == 0
        and gates["opening_group_placement_residual_max_mm"]
        <= args.tolerance_mm
        and gates["shared_groups_preserved"]
        and gates["product_geometry_over_tolerance"] == 0
        and gates["schema_equal"]
        and gates["entity_count_delta_matches_expected"]
        and gates["root_global_ids_equal"]
        and gates["owner_history_additions"]["scoped"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-shared-surface-host-origin-reset-candidate",
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
        "final_anchor_surface_distances_mm": final_anchor_distances,
        "placement_group_results": placement_group_results,
        "related_opening_global_ids": related_opening_ids,
        "all_product_geometry": geometry,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps({"report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
