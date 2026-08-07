#!/usr/bin/env python3
"""Generate and mechanically compare isolated wall-normalization candidates.

The source IFC is read-only. Each strategy is written to ``build/`` and is
compared against the source for ObjectPlacement, world geometry, topology,
and wall-alignment candidates. No candidate is promoted to the project IFC.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api.geometry
import ifcopenshell.geom
import ifcopenshell.util.placement
import ifcopenshell.util.representation
import numpy as np

from geometry_alignment_audit import alignment_audit, geometry_difference_audit


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise argparse.ArgumentTypeError("must be a finite number greater than zero")
    return number


def body_solid(
    wall: ifcopenshell.entity_instance,
) -> tuple[ifcopenshell.entity_instance, ifcopenshell.entity_instance]:
    body = ifcopenshell.util.representation.get_representation(wall, "Model", "Body", "MODEL_VIEW")
    if body is None or len(body.Items) != 1 or not body.Items[0].is_a("IfcExtrudedAreaSolid"):
        raise RuntimeError("Target wall must have one Body IfcExtrudedAreaSolid for isolated comparison")
    return body, body.Items[0]


def nearest_proper_cardinal_rotation(rotation: np.ndarray) -> np.ndarray:
    candidates: list[np.ndarray] = []
    identity = np.eye(3)
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1.0, 1.0), repeat=3):
            candidate = np.column_stack(
                [identity[:, permutation[index]] * signs[index] for index in range(3)]
            )
            if np.linalg.det(candidate) > 0.0:
                candidates.append(candidate)
    return min(candidates, key=lambda candidate: float(np.linalg.norm(rotation - candidate)))


def set_axis2placement3d(axis_placement: ifcopenshell.entity_instance, matrix: np.ndarray) -> None:
    axis_placement.Location.Coordinates = tuple(float(value) for value in matrix[:3, 3])
    z_axis = tuple(float(value) for value in matrix[:3, 2])
    x_axis = tuple(float(value) for value in matrix[:3, 0])
    if axis_placement.Axis is None:
        raise RuntimeError("Target solid position has no Axis")
    if axis_placement.RefDirection is None:
        raise RuntimeError("Target solid position has no RefDirection")
    axis_placement.Axis.DirectionRatios = z_axis
    axis_placement.RefDirection.DirectionRatios = x_axis


def placement_target(source_matrix: np.ndarray, cardinal: bool) -> np.ndarray:
    target = source_matrix.copy()
    target[:3, 3] = np.round(source_matrix[:3, 3])
    if cardinal:
        target[:3, :3] = nearest_proper_cardinal_rotation(source_matrix[:3, :3])
    return target


def apply_strategy(
    model: ifcopenshell.file,
    global_id: str,
    strategy: str,
) -> dict[str, Any]:
    wall = model.by_guid(global_id)
    if wall is None or not wall.is_a("IfcWall"):
        raise RuntimeError(f"{global_id} is not an IfcWall")
    _, solid = body_solid(wall)

    source_product_matrix = np.array(
        ifcopenshell.util.placement.get_local_placement(wall.ObjectPlacement), dtype=float
    )
    feature_world_placements = {
        relation.RelatedOpeningElement.GlobalId: np.array(
            ifcopenshell.util.placement.get_local_placement(
                relation.RelatedOpeningElement.ObjectPlacement
            ),
            dtype=float,
        )
        for relation in wall.HasOpenings
    }
    source_solid_matrix = np.array(
        ifcopenshell.util.placement.get_axis2placement(solid.Position), dtype=float
    )

    if strategy == "placement_integer_rebase":
        target_product_matrix = placement_target(source_product_matrix, cardinal=False)
        target_solid_matrix = (
            np.linalg.inv(target_product_matrix) @ source_product_matrix @ source_solid_matrix
        )
    elif strategy == "rigid_integer_translation":
        target_product_matrix = placement_target(source_product_matrix, cardinal=False)
        target_solid_matrix = source_solid_matrix
    elif strategy == "rigid_integer_cardinal":
        target_product_matrix = placement_target(source_product_matrix, cardinal=True)
        target_solid_matrix = source_solid_matrix
    else:
        raise ValueError(strategy)

    ifcopenshell.api.geometry.edit_object_placement(
        model,
        product=wall,
        matrix=target_product_matrix.copy(),
        is_si=False,
        should_transform_children=False,
    )
    # edit_object_placement deliberately moves feature elements with their host.
    # Restore each opening's original world matrix so a large clipping solid does
    # not amplify a sub-millimetre wall normalization into millimetres of drift.
    for opening_id, opening_world_matrix in feature_world_placements.items():
        opening = model.by_guid(opening_id)
        ifcopenshell.api.geometry.edit_object_placement(
            model,
            product=opening,
            matrix=opening_world_matrix.copy(),
            is_si=False,
            should_transform_children=False,
        )
    if strategy == "placement_integer_rebase":
        set_axis2placement3d(solid.Position, target_solid_matrix)

    result_product_matrix = np.array(
        ifcopenshell.util.placement.get_local_placement(wall.ObjectPlacement), dtype=float
    )
    result_solid_matrix = np.array(
        ifcopenshell.util.placement.get_axis2placement(solid.Position), dtype=float
    )
    return {
        "strategy": strategy,
        "source_product_matrix": source_product_matrix.tolist(),
        "target_product_matrix": target_product_matrix.tolist(),
        "result_product_matrix": result_product_matrix.tolist(),
        "source_solid_matrix": source_solid_matrix.tolist(),
        "target_solid_matrix": target_solid_matrix.tolist(),
        "result_solid_matrix": result_solid_matrix.tolist(),
        "placement_translation_delta_mm": float(
            np.linalg.norm(result_product_matrix[:3, 3] - source_product_matrix[:3, 3])
        ),
        "rotation_matrix_max_delta": float(
            np.max(np.abs(result_product_matrix[:3, :3] - source_product_matrix[:3, :3]))
        ),
        "restored_feature_world_placements": sorted(feature_world_placements),
    }


def target_alignment_records(alignment: dict[str, Any], global_id: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key in ("coplanar_edges", "junctions"):
        records = [
            record
            for record in alignment[key]["records"]
            if global_id in (record["wall_a"], record["wall_b"])
        ]
        result[key] = {
            "total": len(records),
            "within_tolerance": sum(record["within_tolerance"] for record in records),
            "over_tolerance": sum(not record["within_tolerance"] for record in records),
            "records": records,
        }
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--global-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=positive_float, default=0.1)
    parser.add_argument("--search-window-mm", type=positive_float, default=1.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = args.input.resolve()
    source_model = ifcopenshell.open(source_path)
    source_wall = source_model.by_guid(args.global_id)
    if source_wall is None:
        raise RuntimeError(f"Target {args.global_id} not found")
    body, solid = body_solid(source_wall)
    profile_inverse_count = len(source_model.get_inverse(solid.SweptArea))
    solid_position_inverse_count = len(source_model.get_inverse(solid.Position))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    strategies = (
        "placement_integer_rebase",
        "rigid_integer_translation",
        "rigid_integer_cardinal",
    )
    candidates = []
    for strategy in strategies:
        candidate_model = ifcopenshell.open(source_path)
        transform = apply_strategy(candidate_model, args.global_id, strategy)
        candidate_path = args.output_dir / f"{args.global_id}-{strategy}.ifc"
        candidate_model.write(candidate_path)
        candidate_model = ifcopenshell.open(candidate_path)
        alignment = alignment_audit(
            candidate_model,
            tolerance_mm=args.tolerance_mm,
            search_window_mm=args.search_window_mm,
            min_segment_length_mm=100.0,
            min_overlap_mm=100.0,
            min_vertical_overlap_mm=100.0,
            angle_tolerance_deg=0.01,
        )
        difference = geometry_difference_audit(
            candidate_model,
            source_model,
            str(source_path),
            tolerance_mm=args.tolerance_mm,
            classes=["IfcWall"],
            global_ids=[],
        )
        target_difference = next(
            record for record in difference["records"] if record["global_id"] == args.global_id
        )
        candidates.append(
            {
                "strategy": strategy,
                "path": str(candidate_path),
                "sha256": sha256(candidate_path),
                "transform": transform,
                "geometry_difference": difference,
                "target_geometry_difference": target_difference,
                "target_alignment": target_alignment_records(alignment, args.global_id),
                "alignment_summary": {
                    key: {
                        summary_key: alignment[key][summary_key]
                        for summary_key in ("total", "within_tolerance", "over_tolerance", "max_gap_mm")
                    }
                    for key in ("coplanar_edges", "junctions")
                },
            }
        )

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-candidate-comparison",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source_model.schema,
        },
        "global_id": args.global_id,
        "tolerance_mm": args.tolerance_mm,
        "search_window_mm": args.search_window_mm,
        "shared_geometry_guard": {
            "profile_inverse_count": profile_inverse_count,
            "solid_position_inverse_count": solid_position_inverse_count,
            "profile_must_not_be_edited_in_place": profile_inverse_count > 1,
        },
        "candidates": candidates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "report": str(args.report),
                "global_id": args.global_id,
                "candidates": [
                    {
                        "strategy": candidate["strategy"],
                        "world_vertex_hausdorff_mm": candidate["target_geometry_difference"].get(
                            "world_vertex_hausdorff_mm"
                        ),
                        "target_coplanar_over": candidate["target_alignment"]["coplanar_edges"][
                            "over_tolerance"
                        ],
                        "target_junction_over": candidate["target_alignment"]["junctions"][
                            "over_tolerance"
                        ],
                    }
                    for candidate in candidates
                ],
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
