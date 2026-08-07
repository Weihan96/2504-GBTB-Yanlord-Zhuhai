#!/usr/bin/env python3
"""Build and verify an approved C003 wall-relationship candidate.

The project IFC is read-only. Approved transforms come from a CSV decision
table. The candidate is written below ``build/`` and compared mechanically
against the source before it can be promoted by the single IFC writer.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api.geometry
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.placement
import numpy as np

from c003_wall_candidate_compare import nearest_proper_cardinal_rotation
from geometry_alignment_audit import alignment_audit, geometry_difference_audit, sha256


AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def finite_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise argparse.ArgumentTypeError("must be finite")
    return number


def positive_float(value: str) -> float:
    number = finite_float(value)
    if number <= 0.0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return number


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "yes", "1"}:
        return True
    if normalized in {"false", "no", "0", ""}:
        return False
    raise ValueError(f"Invalid boolean value: {value}")


def parse_point_indices(value: str) -> list[int]:
    if not value.strip():
        return []
    indices = [int(item.strip()) for item in value.split(";")]
    if any(index < 0 for index in indices):
        raise ValueError("Profile point indices must be non-negative")
    if len(indices) != len(set(indices)):
        raise ValueError("Profile point indices must be unique")
    return indices


def parse_float_list(value: str) -> list[float]:
    if not value.strip():
        return []
    return [finite_float(item.strip()) for item in value.split(";")]


def read_plan(path: Path, group_id: str) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = [row for row in csv.DictReader(handle) if row["group_id"] == group_id]
    if not rows:
        raise RuntimeError(f"No rows found for {group_id} in {path}")
    if any(row["status"] != "approved" for row in rows):
        raise RuntimeError(f"Every {group_id} transform must be approved")
    result = []
    for row in rows:
        result.append(
            {
                **row,
                "cardinalize": parse_bool(row["cardinalize"]),
                "translation_mm": [
                    finite_float(row["translation_x_mm"]),
                    finite_float(row["translation_y_mm"]),
                    finite_float(row["translation_z_mm"]),
                ],
                "snap_bbox_target_mm": (
                    finite_float(row.get("snap_bbox_target_mm", ""))
                    if row.get("snap_bbox_target_mm", "")
                    else None
                ),
                "copy_shared_profile": parse_bool(
                    row.get("copy_shared_profile", "")
                ),
                "profile_coordinate_axis": row.get(
                    "profile_coordinate_axis", ""
                ).strip().lower(),
                "profile_point_indices": parse_point_indices(
                    row.get("profile_point_indices", "")
                ),
                "profile_coordinate_targets_mm": parse_float_list(
                    row.get("profile_coordinate_targets_mm", "")
                ),
                "profile_snap_increment_mm": (
                    positive_float(row.get("profile_snap_increment_mm", ""))
                    if row.get("profile_snap_increment_mm", "")
                    else None
                ),
                "extrusion_depth_target_mm": (
                    finite_float(row.get("extrusion_depth_target_mm", ""))
                    if row.get("extrusion_depth_target_mm", "")
                    else None
                ),
                "confidence": finite_float(row["confidence"]),
            }
        )
    return result


def shape_bbox_mm(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> list[list[float]]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    values = list(shape.geometry.verts)
    points = [
        [float(values[index + axis]) * 1000.0 for axis in range(3)]
        for index in range(0, len(values), 3)
    ]
    if not points:
        raise RuntimeError(f"{product.GlobalId} has no shape vertices")
    return [[min(point[axis] for point in points), max(point[axis] for point in points)] for axis in range(3)]


def giant_opening_world_placements(
    settings: ifcopenshell.geom.settings,
    wall: ifcopenshell.entity_instance,
    threshold_mm: float,
) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for relation in wall.HasOpenings:
        opening = relation.RelatedOpeningElement
        bbox = shape_bbox_mm(settings, opening)
        horizontal_span = max(bbox[0][1] - bbox[0][0], bbox[1][1] - bbox[1][0])
        if not opening.HasFillings and horizontal_span >= threshold_mm:
            result[opening.GlobalId] = np.array(
                ifcopenshell.util.placement.get_local_placement(opening.ObjectPlacement),
                dtype=float,
            )
    return result


def apply_wall_row(
    model: ifcopenshell.file,
    settings: ifcopenshell.geom.settings,
    row: dict[str, Any],
    giant_threshold_mm: float,
) -> dict[str, Any]:
    wall = model.by_guid(row["global_id"])
    if wall is None or not wall.is_a("IfcWall"):
        raise RuntimeError(f'{row["global_id"]} is not an IfcWall')
    source = np.array(
        ifcopenshell.util.placement.get_local_placement(wall.ObjectPlacement), dtype=float
    )
    target = source.copy()
    if row["cardinalize"]:
        target[:3, :3] = nearest_proper_cardinal_rotation(source[:3, :3])
    target[:3, 3] += np.array(row["translation_mm"], dtype=float)
    guarded_openings = giant_opening_world_placements(settings, wall, giant_threshold_mm)
    allowed_moved_descendants: set[str] = set()
    for relation in wall.HasOpenings:
        opening = relation.RelatedOpeningElement
        if opening.GlobalId in guarded_openings:
            continue
        allowed_moved_descendants.add(opening.GlobalId)
        allowed_moved_descendants.update(
            filling.RelatedBuildingElement.GlobalId for filling in opening.HasFillings
        )
    ifcopenshell.api.geometry.edit_object_placement(
        model,
        product=wall,
        matrix=target.copy(),
        is_si=False,
        should_transform_children=False,
    )
    for opening_id, world_matrix in guarded_openings.items():
        ifcopenshell.api.geometry.edit_object_placement(
            model,
            product=model.by_guid(opening_id),
            matrix=world_matrix.copy(),
            is_si=False,
            should_transform_children=False,
        )
    result = np.array(
        ifcopenshell.util.placement.get_local_placement(wall.ObjectPlacement), dtype=float
    )
    return {
        "global_id": wall.GlobalId,
        "source_matrix": source.tolist(),
        "result_matrix": result.tolist(),
        "translation_delta_mm": float(np.linalg.norm(result[:3, 3] - source[:3, 3])),
        "translation_by_axis_mm": [
            float(result[axis, 3] - source[axis, 3]) for axis in range(3)
        ],
        "rotation_matrix_max_delta": float(
            np.max(np.abs(result[:3, :3] - source[:3, :3]))
        ),
        "allowed_moved_descendants": sorted(allowed_moved_descendants),
        "restored_giant_openings": sorted(guarded_openings),
    }


def body_extruded_solid(product: ifcopenshell.entity_instance) -> ifcopenshell.entity_instance:
    representations = (
        product.Representation.Representations if product.Representation else ()
    )
    body_representations = [
        representation
        for representation in representations
        if representation.RepresentationIdentifier == "Body"
    ]
    solids = [
        entity
        for representation in body_representations
        for item in representation.Items
        for entity in product.file.traverse(item)
        if entity.is_a("IfcExtrudedAreaSolid")
    ]
    if len(solids) != 1:
        raise RuntimeError(
            f"{product.GlobalId} must have exactly one Body IfcExtrudedAreaSolid; "
            f"found {len(solids)}"
        )
    return solids[0]


def snap_profile_coordinate(
    value_mm: float, increment_mm: float, maximum_shift_mm: float
) -> tuple[float, float]:
    snapped = round(value_mm / increment_mm) * increment_mm
    shift = snapped - value_mm
    if abs(shift) > maximum_shift_mm:
        raise RuntimeError(
            f"Profile control point shift {shift:.6f} mm exceeds "
            f"{maximum_shift_mm:.6f} mm"
        )
    return snapped, shift


def apply_profile_row(
    model: ifcopenshell.file,
    row: dict[str, Any],
    maximum_shift_mm: float = 0.5,
) -> dict[str, Any]:
    product = model.by_guid(row["global_id"])
    if product is None or not product.is_a(row["ifc_class"]):
        raise RuntimeError(f'{row["global_id"]} is not {row["ifc_class"]}')
    if not row["copy_shared_profile"]:
        raise RuntimeError(f'{product.GlobalId} profile edit must copy the source Profile')
    axis_name = row["profile_coordinate_axis"]
    if axis_name not in {"x", "y"}:
        raise RuntimeError(f"Invalid profile coordinate axis for {product.GlobalId}")
    indices = row["profile_point_indices"]
    targets = row["profile_coordinate_targets_mm"]
    increment_mm = row["profile_snap_increment_mm"]
    if not indices or (increment_mm is None and not targets):
        raise RuntimeError(f"Incomplete profile snap plan for {product.GlobalId}")
    if targets and len(targets) != len(indices):
        raise RuntimeError(
            f"Profile target count does not match point count for {product.GlobalId}"
        )

    solid = body_extruded_solid(product)
    source_profile = solid.SweptArea
    if not source_profile.is_a("IfcArbitraryClosedProfileDef"):
        raise RuntimeError(f"{product.GlobalId} does not use IfcArbitraryClosedProfileDef")
    source_curve = source_profile.OuterCurve
    if not source_curve.is_a("IfcIndexedPolyCurve") or not source_curve.Points.is_a(
        "IfcCartesianPointList2D"
    ):
        raise RuntimeError(f"{product.GlobalId} does not use an indexed 2D Profile")
    source_coordinates = tuple(tuple(point) for point in source_curve.Points.CoordList)
    source_inverse_count = model.get_total_inverses(source_profile)
    copied_profile = ifcopenshell.util.element.copy_deep(model, source_profile)
    solid.SweptArea = copied_profile
    copied_curve = copied_profile.OuterCurve
    coordinates = [list(point) for point in copied_curve.Points.CoordList]
    axis = AXIS_INDEX[axis_name]
    point_results = []
    for position, index in enumerate(indices):
        if index >= len(coordinates):
            raise RuntimeError(
                f"Profile point index {index} out of range for {product.GlobalId}"
            )
        before = float(coordinates[index][axis])
        if targets:
            after = targets[position]
            shift = after - before
            if abs(shift) > maximum_shift_mm:
                raise RuntimeError(
                    f"Profile control point shift {shift:.6f} mm exceeds "
                    f"{maximum_shift_mm:.6f} mm"
                )
        else:
            after, shift = snap_profile_coordinate(
                before, increment_mm, maximum_shift_mm
            )
        coordinates[index][axis] = after
        point_results.append(
            {
                "index": index,
                "axis": axis_name,
                "before_mm": before,
                "after_mm": after,
                "shift_mm": shift,
            }
        )
    copied_curve.Points.CoordList = tuple(tuple(point) for point in coordinates)
    if tuple(tuple(point) for point in source_curve.Points.CoordList) != source_coordinates:
        raise RuntimeError(f"Source Profile changed for {product.GlobalId}")
    return {
        "global_id": product.GlobalId,
        "source_profile_id": source_profile.id(),
        "source_profile_inverse_count": source_inverse_count,
        "copied_profile_id": copied_profile.id(),
        "copied_profile_inverse_count": model.get_total_inverses(copied_profile),
        "coordinate_axis": axis_name,
        "snap_increment_mm": increment_mm,
        "explicit_targets": bool(targets),
        "maximum_absolute_shift_mm": max(
            abs(result["shift_mm"]) for result in point_results
        ),
        "point_results": point_results,
        "source_profile_unchanged": True,
    }


def apply_extrusion_depth_row(
    model: ifcopenshell.file,
    row: dict[str, Any],
    maximum_shift_mm: float = 1.0,
) -> dict[str, Any]:
    product = model.by_guid(row["global_id"])
    if product is None or not product.is_a(row["ifc_class"]):
        raise RuntimeError(f'{row["global_id"]} is not {row["ifc_class"]}')
    solid = body_extruded_solid(product)
    before = float(solid.Depth)
    after = float(row["extrusion_depth_target_mm"])
    shift = after - before
    if abs(shift) > maximum_shift_mm:
        raise RuntimeError(
            f"Extrusion depth shift {shift:.6f} mm exceeds "
            f"{maximum_shift_mm:.6f} mm for {product.GlobalId}"
        )
    solid.Depth = after
    return {
        "global_id": product.GlobalId,
        "ifc_class": product.is_a(),
        "solid_id": solid.id(),
        "before_mm": before,
        "after_mm": after,
        "shift_mm": shift,
    }


def apply_bbox_row(
    model: ifcopenshell.file,
    settings: ifcopenshell.geom.settings,
    row: dict[str, Any],
) -> dict[str, Any]:
    product = model.by_guid(row["global_id"])
    if product is None or not product.is_a(row["ifc_class"]):
        raise RuntimeError(f'{row["global_id"]} is not {row["ifc_class"]}')
    axis_name = row["snap_bbox_axis"].lower()
    side_name = row["snap_bbox_side"].lower()
    if axis_name not in AXIS_INDEX or side_name not in {"min", "max"}:
        raise RuntimeError(f"Invalid bbox snap for {product.GlobalId}")
    axis = AXIS_INDEX[axis_name]
    side = 0 if side_name == "min" else 1
    before_bbox = shape_bbox_mm(settings, product)
    current_value = before_bbox[axis][side]
    delta = float(row["snap_bbox_target_mm"] - current_value)
    source = np.array(
        ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement), dtype=float
    )
    target = source.copy()
    target[axis, 3] += delta
    ifcopenshell.api.geometry.edit_object_placement(
        model,
        product=product,
        matrix=target.copy(),
        is_si=False,
        should_transform_children=False,
    )
    after_bbox = shape_bbox_mm(settings, product)
    return {
        "global_id": product.GlobalId,
        "axis": axis_name,
        "side": side_name,
        "target_mm": row["snap_bbox_target_mm"],
        "before_mm": current_value,
        "after_mm": after_bbox[axis][side],
        "translation_delta_mm": delta,
    }


def pair_results(
    alignment: dict[str, Any], pairs: list[list[str]], tolerance_mm: float
) -> list[dict[str, Any]]:
    result = []
    for raw_pair in pairs:
        pair = tuple(sorted(raw_pair))
        records = [
            {"record_type": record_type, **record}
            for record_type, key in (
                ("coplanar_edge", "coplanar_edges"),
                ("junction", "junctions"),
            )
            for record in alignment[key]["records"]
            if tuple(sorted((record["wall_a"], record["wall_b"]))) == pair
        ]
        result.append(
            {
                "wall_pair": list(pair),
                "record_count": len(records),
                "within_tolerance": bool(records)
                and all(record["gap_mm"] <= tolerance_mm for record in records),
                "maximum_gap_mm": max(
                    (float(record["gap_mm"]) for record in records), default=None
                ),
                "records": records,
            }
        )
    return result


def entity_counts(
    model: ifcopenshell.file, extra_classes: set[str] | None = None
) -> dict[str, int]:
    classes = {
        "IfcRoot",
        "IfcWall",
        "IfcOpeningElement",
        "IfcDoor",
        "IfcWindow",
        "IfcSpace",
        *(extra_classes or set()),
    }
    return {
        ifc_class: len(model.by_type(ifc_class))
        for ifc_class in sorted(classes)
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--alignment-report", required=True, type=Path)
    parser.add_argument("--group-id", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=positive_float, default=0.1)
    parser.add_argument("--search-window-mm", type=positive_float, default=1.0)
    parser.add_argument("--giant-opening-threshold-mm", type=positive_float, default=10000.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    plan = read_plan(args.plan, args.group_id)
    source_alignment_report = json.loads(args.alignment_report.read_text(encoding="utf-8"))
    if source_alignment_report["source"]["sha256"] != sha256(source_path):
        raise RuntimeError("Alignment report does not match the source IFC")
    source_cluster = next(
        cluster
        for cluster in source_alignment_report["alignment"]["review_clusters"]
        if args.group_id in (cluster["review_id"], cluster.get("stable_id"))
    ) if any(
        args.group_id in (cluster["review_id"], cluster.get("stable_id"))
        for cluster in source_alignment_report["alignment"]["review_clusters"]
    ) else {
        "review_id": args.group_id,
        "stable_id": args.group_id,
        "status": "parameter_only_no_wall_pairs",
        "pairs": [],
    }

    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    wall_results = [
        apply_wall_row(candidate, settings, row, args.giant_opening_threshold_mm)
        for row in plan
        if row["ifc_class"] == "IfcWall"
        and (row["cardinalize"] or any(row["translation_mm"]))
    ]
    profile_results = [
        apply_profile_row(candidate, row)
        for row in plan
        if row["copy_shared_profile"]
    ]
    extrusion_depth_results = [
        apply_extrusion_depth_row(candidate, row)
        for row in plan
        if row["extrusion_depth_target_mm"] is not None
    ]
    bbox_results = [
        apply_bbox_row(candidate, settings, row)
        for row in plan
        if row["snap_bbox_axis"]
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    alignment = alignment_audit(
        candidate,
        tolerance_mm=args.tolerance_mm,
        search_window_mm=args.search_window_mm,
        min_segment_length_mm=100.0,
        min_overlap_mm=100.0,
        min_vertical_overlap_mm=100.0,
        angle_tolerance_deg=0.01,
    )
    confirmed_pair_results = pair_results(
        alignment, source_cluster["pairs"], args.tolerance_mm
    )
    comparison_classes = sorted(
        {"IfcWall", "IfcOpeningElement", "IfcDoor", "IfcWindow"}
        | {row["ifc_class"] for row in plan}
    )
    difference = geometry_difference_audit(
        candidate,
        source,
        str(source_path),
        tolerance_mm=args.tolerance_mm,
        classes=comparison_classes,
        global_ids=[],
    )
    topology_mismatches = [
        record
        for record in difference["records"]
        if record.get("topology_counts_equal") is False
    ]
    group_pair_failures = [
        result for result in confirmed_pair_results if not result["within_tolerance"]
    ]
    bbox_snap_failures = [
        result
        for result in bbox_results
        if abs(float(result["after_mm"]) - float(result["target_mm"])) > args.tolerance_mm
    ]
    allowed_changed_ids = {row["global_id"] for row in plan}
    for result in wall_results:
        allowed_changed_ids.update(result["allowed_moved_descendants"])
    unexpected_changes = [
        record
        for record in difference["records"]
        if record["status"] != "unchanged" and record["global_id"] not in allowed_changed_ids
    ]
    plan_classes = {row["ifc_class"] for row in plan}
    source_counts = entity_counts(source, plan_classes)
    candidate_counts = entity_counts(candidate, plan_classes)
    source_root_guids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_guids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    schema_equal = source.schema == candidate.schema
    counts_equal = source_counts == candidate_counts
    root_guids_equal = source_root_guids == candidate_root_guids
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-approved-group-candidate",
        "group_id": args.group_id,
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "candidate": {
            "path": str(args.output),
            "sha256": sha256(args.output),
            "schema": candidate.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "source_cluster": source_cluster,
        "plan": plan,
        "wall_results": wall_results,
        "profile_results": profile_results,
        "extrusion_depth_results": extrusion_depth_results,
        "bbox_results": bbox_results,
        "confirmed_pair_results": confirmed_pair_results,
        "alignment_summary": {
            key: {
                item: alignment[key][item]
                for item in ("total", "within_tolerance", "over_tolerance", "max_gap_mm")
            }
            for key in ("coplanar_edges", "junctions")
        },
        "geometry_difference": difference,
        "allowed_changed_global_ids": sorted(allowed_changed_ids),
        "unexpected_changes": unexpected_changes,
        "integrity": {
            "schema_equal": schema_equal,
            "source_entity_counts": source_counts,
            "candidate_entity_counts": candidate_counts,
            "entity_counts_equal": counts_equal,
            "root_global_ids_equal": root_guids_equal,
        },
        "gates": {
            "confirmed_pairs_total": len(confirmed_pair_results),
            "confirmed_pairs_failed": len(group_pair_failures),
            "bbox_snaps_failed": len(bbox_snap_failures),
            "topology_mismatches": len(topology_mismatches),
            "unexpected_changes": len(unexpected_changes),
            "schema_equal": schema_equal,
            "entity_counts_equal": counts_equal,
            "root_global_ids_equal": root_guids_equal,
            "pass": (
                not group_pair_failures
                and not bbox_snap_failures
                and not topology_mismatches
                and not unexpected_changes
                and schema_equal
                and counts_equal
                and root_guids_equal
            ),
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "report": str(args.report),
                "candidate": str(args.output),
                "candidate_sha256": report["candidate"]["sha256"],
                "confirmed_pairs_total": report["gates"]["confirmed_pairs_total"],
                "confirmed_pairs_failed": report["gates"]["confirmed_pairs_failed"],
                "bbox_snaps_failed": report["gates"]["bbox_snaps_failed"],
                "topology_mismatches": report["gates"]["topology_mismatches"],
                "unexpected_changes": report["gates"]["unexpected_changes"],
                "schema_equal": report["gates"]["schema_equal"],
                "entity_counts_equal": report["gates"]["entity_counts_equal"],
                "root_global_ids_equal": report["gates"]["root_global_ids_equal"],
                "geometry_difference": {
                    key: difference[key]
                    for key in ("total", "within_tolerance", "over_tolerance", "status_counts")
                },
                "pass": report["gates"]["pass"],
            },
            ensure_ascii=False,
        )
    )
    if not report["gates"]["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
