"""Audit an exported Blender HVAC constraint preview against IFC-backed anchors."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cardinal_segment(first: list[float], second: list[float], tolerance_mm: float) -> tuple[bool, list[float]]:
    delta = [second[index] - first[index] for index in range(3)]
    changing_axes = sum(abs(value) > tolerance_mm for value in delta)
    return changing_axes == 1, delta


def cardinal_axis(axis: list[float] | None, tolerance_mm: float) -> bool:
    if axis is None or len(axis) != 3:
        return False
    nonzero = [abs(value) > tolerance_mm for value in axis]
    return sum(nonzero) == 1 and math.isclose(sum(value * value for value in axis), 1.0, abs_tol=tolerance_mm)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--route-report", type=Path, required=True)
    parser.add_argument("--preview", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    args = parser.parse_args()

    formal_sha = sha256(args.input)
    route_report = json.loads(args.route_report.read_text(encoding="utf-8"))
    preview = json.loads(args.preview.read_text(encoding="utf-8"))
    if route_report["source"]["ifc_sha256"] != formal_sha:
        raise RuntimeError("route report does not match the formal IFC")
    if preview["source_ifc_sha256"] != formal_sha:
        raise RuntimeError("Blender preview candidate does not match the formal IFC")
    if preview.get("formal_ifc_write_allowed") is not False:
        raise RuntimeError("Blender preview must not authorize formal IFC writes")

    expected_routes = {
        route["route_id"]: route
        for route in route_report["confirmed_route_graph"]
        if len(route["waypoints"]) >= 2
    }
    preview_routes = {route["route_id"]: route for route in preview["routes"]}
    if set(preview_routes) != set(expected_routes):
        raise RuntimeError("Blender preview route set drift")

    anchor_checks = []
    route_checks = []
    maximum_anchor_deviation = 0.0
    maximum_equipment_service_offset = 0.0
    bend_count = 0
    all_routes_orthogonal = True
    all_route_anchors_on_path = True
    for route_id, expected in expected_routes.items():
        actual = preview_routes[route_id]
        expected_by_id = {row["anchor_id"]: row for row in expected["waypoints"]}
        expected_fixed_order = [row["anchor_id"] for row in expected["waypoints"]]
        actual_fixed_order = [
            row["anchor_id"]
            for row in actual["anchors"]
            if row["anchor_kind"] != "blender_bend"
        ]
        if actual_fixed_order != expected_fixed_order:
            raise RuntimeError(f"fixed waypoint order drift: {route_id}")
        for row in actual["anchors"]:
            if row["anchor_kind"] == "blender_bend":
                bend_count += 1
                anchor_checks.append({
                    "route_id": route_id,
                    "anchor_id": row["anchor_id"],
                    "anchor_kind": "blender_bend",
                    "world_mm": row["world_mm"],
                    "status": "human_review_required_before_register_update",
                })
                continue
            expected_row = expected_by_id[row["anchor_id"]]
            placement_basis = row.get("placement_basis", "ifc_object_bbox_centre")
            is_equipment_service_face = placement_basis == "equipment_local_positive_x_service_face"
            if is_equipment_service_face:
                reference_world = row.get("reference_world_mm")
                if reference_world is None:
                    raise RuntimeError(f"equipment service anchor lacks IFC reference centre: {route_id}")
                deviation = math.dist(reference_world, expected_row["centre_mm"])
                service_offset = math.dist(row["world_mm"], reference_world)
                maximum_equipment_service_offset = max(maximum_equipment_service_offset, service_offset)
                service_axis_is_cardinal = cardinal_axis(row.get("service_axis_world"), args.tolerance_mm)
            else:
                deviation = math.dist(row["world_mm"], expected_row["centre_mm"])
                service_offset = None
                service_axis_is_cardinal = None
            maximum_anchor_deviation = max(maximum_anchor_deviation, deviation)
            anchor_checks.append({
                "route_id": route_id,
                "anchor_id": row["anchor_id"],
                "anchor_kind": row["anchor_kind"],
                "global_id": row["global_id"],
                "deviation_from_current_ifc_anchor_mm": deviation,
                "within_tolerance": deviation <= args.tolerance_mm,
                "placement_basis": placement_basis,
                "equipment_service_offset_from_reference_mm": service_offset,
                "service_axis_is_cardinal": service_axis_is_cardinal,
            })
        orthogonal_points = actual.get("orthogonal_points_world_mm", [])
        if len(orthogonal_points) < 2:
            raise RuntimeError(f"route lacks evaluated orthogonal points: {route_id}")
        segment_rows = []
        for first, second in zip(orthogonal_points, orthogonal_points[1:]):
            is_cardinal, delta = cardinal_segment(first, second, args.tolerance_mm)
            length = math.dist(first, second)
            segment_rows.append({
                "from_world_mm": first,
                "to_world_mm": second,
                "delta_mm": delta,
                "length_mm": length,
                "axis_aligned": is_cardinal,
                "nonzero": length > args.tolerance_mm,
            })
        route_is_orthogonal = all(row["axis_aligned"] and row["nonzero"] for row in segment_rows)
        all_routes_orthogonal = all_routes_orthogonal and route_is_orthogonal
        search_from = 0
        anchor_point_indices = []
        for anchor in actual["anchors"]:
            match = next((
                index
                for index in range(search_from, len(orthogonal_points))
                if math.dist(anchor["world_mm"], orthogonal_points[index]) <= args.tolerance_mm
            ), None)
            if match is None:
                anchor_point_indices.append(None)
                continue
            anchor_point_indices.append(match)
            search_from = match
        route_anchors_on_path = all(index is not None for index in anchor_point_indices)
        all_route_anchors_on_path = all_route_anchors_on_path and route_anchors_on_path
        route_checks.append({
            "route_id": route_id,
            "anchor_order": [row["anchor_id"] for row in actual["anchors"]],
            "orthogonal_points_world_mm": orthogonal_points,
            "segments": segment_rows,
            "total_length_mm": sum(row["length_mm"] for row in segment_rows),
            "all_segments_axis_aligned": route_is_orthogonal,
            "anchor_point_indices": anchor_point_indices,
            "all_anchors_on_path_in_order": route_anchors_on_path,
            "status": "constraint_skeleton_not_fabrication_geometry",
        })

    fixed_anchor_checks = [row for row in anchor_checks if row["anchor_kind"] != "blender_bend"]
    fixed_anchors_pass = all(
        row["within_tolerance"] and row.get("service_axis_is_cardinal") is not False
        for row in fixed_anchor_checks
    )
    output = {
        "mode": "read_only_rcp1_hvac_preview_audit",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source": {
            "ifc": str(args.input.resolve()),
            "ifc_sha256": formal_sha,
            "route_report": str(args.route_report.resolve()),
            "preview_candidate": str(args.preview.resolve()),
        },
        "tolerance_mm": args.tolerance_mm,
        "route_checks": route_checks,
        "anchor_checks": anchor_checks,
        "summary": {
            "route_count": len(route_checks),
            "fixed_anchor_check_count": len(fixed_anchor_checks),
            "temporary_bend_count": bend_count,
            "maximum_fixed_anchor_deviation_mm": maximum_anchor_deviation,
            "maximum_equipment_service_offset_from_reference_mm": maximum_equipment_service_offset,
            "all_routes_orthogonal": all_routes_orthogonal,
            "all_route_anchors_on_path_in_order": all_route_anchors_on_path,
        },
        "gates": {
            "source_hash_matches": True,
            "route_set_matches": True,
            "fixed_anchor_order_matches": True,
            "fixed_anchors_within_tolerance": fixed_anchors_pass,
            "equipment_service_axes_cardinal": all(
                row.get("service_axis_is_cardinal") is not False for row in fixed_anchor_checks
            ),
            "all_preview_segments_axis_aligned": all_routes_orthogonal,
            "all_route_anchors_on_path_in_order": all_route_anchors_on_path,
            "geometry_nodes_fillet_present": preview.get("preview_parameters", {}).get("fillet_radius_m", 0) > 0,
            "temporary_bends_require_human_review": bend_count > 0,
            "fabrication_geometry_ready": False,
            "formal_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output["summary"], ensure_ascii=False))
    return 0 if fixed_anchors_pass and all_routes_orthogonal and all_route_anchors_on_path else 1


if __name__ == "__main__":
    raise SystemExit(main())
