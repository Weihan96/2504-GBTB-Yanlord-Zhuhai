#!/usr/bin/env python3
"""Build and verify the confirmed A-103 wall semantics candidate.

The source IFC is read-only.  The candidate records the user-confirmed wall
phase/load-bearing classification and normalizes the kitchen door-pier
finished thickness from 153.886884 mm to 154 mm.  Mechanical gates protect
the already-adjusted guest-bathroom opening and every non-target product.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api
import ifcopenshell.util.element
import ifcopenshell.util.representation

from a103_wall_plan_candidate import geometry_settings, material_names, world_bbox_mm
from geometry_alignment_audit import alignment_audit, geometry_difference_audit


NEW_WALL_GUIDS = {
    "04rs0EDjn2EvxytEQSxWRB",
    "3jha5L04zBjh0$pMl_tHLy",
    "2Ca2tGerPBj8kvUTjlUtDl",
    "3pAfMJYxPBlwR30CstZ4kK",
}
KITCHEN_PIER_GUID = "2Ca2tGerPBj8kvUTjlUtDl"
PROTECTED_ADJUSTED_WALL_GUID = "0hKdvAZkn1TejLgJhK_vDp"
PROTECTED_ADJUSTED_OPENING_GUID = "1YxMx6s0r3ZPPohkRKXWbl"
KITCHEN_PIER_OPENING_GUID = "3uJdxzrgTCV8XAKZYD$foN"
EXPECTED_WALL_COUNT = 88
EXPECTED_AIRCRETE_EXISTING = 64
EXPECTED_CONCRETE_EXISTING = 20
EXPECTED_KITCHEN_PIER_THICKNESS_MM = 154.0


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


def wall_common_pset(
    model: ifcopenshell.file, wall: ifcopenshell.entity_instance
) -> ifcopenshell.entity_instance:
    psets = ifcopenshell.util.element.get_psets(wall)
    existing = psets.get("Pset_WallCommon")
    if existing:
        return model.by_id(int(existing["id"]))
    return ifcopenshell.api.run("pset.add_pset", model, product=wall, name="Pset_WallCommon")


def expected_wall_semantics(wall: ifcopenshell.entity_instance) -> dict[str, Any]:
    materials = material_names(wall)
    if wall.GlobalId in NEW_WALL_GUIDS:
        return {
            "status": "NEW",
            "load_bearing": None,
            "classification_basis": "user-confirmed NEW wall",
            "materials": materials,
        }
    if "Aircrete" in materials:
        return {
            "status": "EXISTING",
            "load_bearing": False,
            "classification_basis": "user-confirmed orange non-load-bearing EXISTING wall",
            "materials": materials,
        }
    if "Concrete" in materials:
        return {
            "status": "EXISTING",
            "load_bearing": True,
            "classification_basis": "user-confirmed blue-gray load-bearing EXISTING wall",
            "materials": materials,
        }
    raise RuntimeError(
        f"wall is outside the confirmed A-103 classification boundary: {wall.GlobalId} {materials}"
    )


def final_built_walls(model: ifcopenshell.file) -> list[ifcopenshell.entity_instance]:
    """Return only the A-103 final-built wall boundary.

    A-102 demolition walls coexist in the formal IFC after their write.  They
    are intentionally excluded here so rerunning A-103 cannot overwrite their
    DEMOLISH status or count them as final-built walls.
    """
    return sorted(
        (
            wall
            for wall in model.by_type("IfcWall")
            if ifcopenshell.util.element.get_psets(wall)
            .get("Pset_WallCommon", {})
            .get("Status")
            != "DEMOLISH"
        ),
        key=lambda wall: wall.GlobalId,
    )


def apply_wall_semantics(model: ifcopenshell.file) -> list[dict[str, Any]]:
    walls = final_built_walls(model)
    if len(walls) != EXPECTED_WALL_COUNT:
        raise RuntimeError(f"expected {EXPECTED_WALL_COUNT} walls, found {len(walls)}")
    records: list[dict[str, Any]] = []
    for wall in walls:
        expected = expected_wall_semantics(wall)
        properties: dict[str, Any] = {"Status": expected["status"]}
        if expected["load_bearing"] is not None:
            properties["LoadBearing"] = expected["load_bearing"]
        pset = wall_common_pset(model, wall)
        ifcopenshell.api.run("pset.edit_pset", model, pset=pset, properties=properties)
        records.append(
            {
                "global_id": wall.GlobalId,
                "status": expected["status"],
                "load_bearing": expected["load_bearing"],
                "materials": expected["materials"],
                "classification_basis": expected["classification_basis"],
                "pset_id": pset.id(),
            }
        )
    return records


def normalize_kitchen_pier_thickness(model: ifcopenshell.file) -> dict[str, Any]:
    wall = model.by_guid(KITCHEN_PIER_GUID)
    body = ifcopenshell.util.representation.get_representation(
        wall, "Model", "Body", "MODEL_VIEW"
    )
    if body is None or len(body.Items) != 1 or not body.Items[0].is_a("IfcPolygonalFaceSet"):
        raise RuntimeError("kitchen pier must have one Body IfcPolygonalFaceSet")
    face_set = body.Items[0]
    points = [list(map(float, point)) for point in face_set.Coordinates.CoordList]
    y_values = [point[1] for point in points]
    y_min = min(y_values)
    y_max = max(y_values)
    source_thickness = y_max - y_min
    target_thickness = round(source_thickness)
    if target_thickness != EXPECTED_KITCHEN_PIER_THICKNESS_MM:
        raise RuntimeError(f"unexpected kitchen pier target thickness: {target_thickness}")
    target_y_min = y_max - target_thickness
    changed_points = 0
    for point in points:
        if abs(point[1] - y_min) <= 1e-9:
            point[1] = target_y_min
            changed_points += 1
    if changed_points != 4:
        raise RuntimeError(f"expected four kitchen-pier face points, found {changed_points}")
    face_set.Coordinates.CoordList = points
    return {
        "global_id": KITCHEN_PIER_GUID,
        "source_local_y_min_mm": y_min,
        "fixed_local_y_max_mm": y_max,
        "target_local_y_min_mm": target_y_min,
        "source_thickness_mm": source_thickness,
        "target_thickness_mm": target_thickness,
        "changed_coordinate_count": changed_points,
        "maximum_intended_world_shift_mm": abs(target_y_min - y_min),
    }


def status_record(model: ifcopenshell.file, wall: ifcopenshell.entity_instance) -> dict[str, Any]:
    common = ifcopenshell.util.element.get_psets(wall).get("Pset_WallCommon", {})
    return {
        "global_id": wall.GlobalId,
        "status": common.get("Status"),
        "load_bearing": common.get("LoadBearing"),
        "materials": material_names(wall),
    }


def record_by_guid(report: dict[str, Any], global_id: str) -> dict[str, Any]:
    return next(record for record in report["records"] if record["global_id"] == global_id)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=positive_float, default=0.1)
    parser.add_argument("--search-window-mm", type=positive_float, default=1.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = args.input.resolve()
    output_path = args.output.resolve()
    if source_path == output_path:
        raise RuntimeError("candidate output must not overwrite the formal IFC")
    source = ifcopenshell.open(source_path)
    source_hash = sha256(source_path)
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    settings = geometry_settings()
    source_kitchen_bbox = world_bbox_mm(settings, source.by_guid(KITCHEN_PIER_GUID))

    candidate = ifcopenshell.open(source_path)
    semantics = apply_wall_semantics(candidate)
    thickness = normalize_kitchen_pier_thickness(candidate)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(output_path)
    candidate = ifcopenshell.open(output_path)

    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    candidate_statuses = [
        status_record(candidate, wall)
        for wall in final_built_walls(candidate)
    ]
    demolition_statuses = [
        status_record(candidate, wall)
        for wall in sorted(candidate.by_type("IfcWall"), key=lambda wall: wall.GlobalId)
        if ifcopenshell.util.element.get_psets(wall)
        .get("Pset_WallCommon", {})
        .get("Status")
        == "DEMOLISH"
    ]
    phase_counts = Counter(record["status"] for record in candidate_statuses)
    load_bearing_existing_counts = Counter(
        record["load_bearing"]
        for record in candidate_statuses
        if record["status"] == "EXISTING"
    )
    geometry = geometry_difference_audit(
        candidate,
        source,
        str(source_path),
        tolerance_mm=args.tolerance_mm,
        classes=("IfcWall", "IfcOpeningElement", "IfcDoor", "IfcWindow"),
        global_ids=(),
    )
    alignment = alignment_audit(
        candidate,
        tolerance_mm=args.tolerance_mm,
        search_window_mm=args.search_window_mm,
        min_segment_length_mm=100.0,
        min_overlap_mm=100.0,
        min_vertical_overlap_mm=100.0,
        angle_tolerance_deg=0.01,
    )
    kitchen_bbox = world_bbox_mm(settings, candidate.by_guid(KITCHEN_PIER_GUID))
    kitchen_geometry = record_by_guid(geometry, KITCHEN_PIER_GUID)
    protected_wall = record_by_guid(geometry, PROTECTED_ADJUSTED_WALL_GUID)
    protected_opening = record_by_guid(geometry, PROTECTED_ADJUSTED_OPENING_GUID)
    kitchen_opening = record_by_guid(geometry, KITCHEN_PIER_OPENING_GUID)
    unexpected_geometry = [
        record
        for record in geometry["records"]
        if record["global_id"] != KITCHEN_PIER_GUID and not record["within_tolerance"]
    ]
    gates = {
        "wall_count": len(candidate_statuses),
        "demolition_wall_count": len(demolition_statuses),
        "phase_counts": dict(sorted(phase_counts.items())),
        "existing_load_bearing_counts": {
            str(key): value for key, value in sorted(load_bearing_existing_counts.items(), key=lambda item: str(item[0]))
        },
        "kitchen_pier_thickness_mm": kitchen_bbox["dimensions_mm"][1],
        "kitchen_pier_world_geometry_delta_mm": kitchen_geometry["world_vertex_hausdorff_mm"],
        "protected_adjusted_wall_delta_mm": protected_wall["world_vertex_hausdorff_mm"],
        "protected_adjusted_opening_delta_mm": protected_opening["world_vertex_hausdorff_mm"],
        "kitchen_pier_opening_delta_mm": kitchen_opening["world_vertex_hausdorff_mm"],
        "unexpected_geometry_over_tolerance": len(unexpected_geometry),
        "coplanar_edges_over_tolerance": alignment["coplanar_edges"]["over_tolerance"],
        "junctions_over_tolerance": alignment["junctions"]["over_tolerance"],
        "source_root_ids_preserved": source_root_ids <= candidate_root_ids,
        "new_root_count": len(candidate_root_ids - source_root_ids),
    }
    intended_delta = thickness["maximum_intended_world_shift_mm"]
    passed = (
        gates["wall_count"] == EXPECTED_WALL_COUNT
        and gates["demolition_wall_count"] in (0, 13)
        and phase_counts == Counter({"EXISTING": 84, "NEW": 4})
        and load_bearing_existing_counts == Counter({False: EXPECTED_AIRCRETE_EXISTING, True: EXPECTED_CONCRETE_EXISTING})
        and abs(gates["kitchen_pier_thickness_mm"] - EXPECTED_KITCHEN_PIER_THICKNESS_MM) <= 1e-6
        and abs(gates["kitchen_pier_world_geometry_delta_mm"] - intended_delta) <= 1e-6
        and gates["protected_adjusted_wall_delta_mm"] <= args.tolerance_mm
        and gates["protected_adjusted_opening_delta_mm"] <= args.tolerance_mm
        and gates["kitchen_pier_opening_delta_mm"] <= args.tolerance_mm
        and gates["unexpected_geometry_over_tolerance"] == 0
        and gates["coplanar_edges_over_tolerance"] == 0
        and gates["junctions_over_tolerance"] == 0
        and gates["source_root_ids_preserved"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a103-wall-semantics-candidate",
        "source": {
            "path": str(source_path),
            "sha256": source_hash,
            "schema": source.schema,
        },
        "candidate": {
            "path": str(output_path),
            "sha256": sha256(output_path),
            "schema": candidate.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "confirmed_semantics": semantics,
        "kitchen_pier_normalization": thickness,
        "source_kitchen_pier_bbox_mm": source_kitchen_bbox,
        "candidate_kitchen_pier_bbox_mm": kitchen_bbox,
        "wall_statuses": candidate_statuses,
        "excluded_demolition_wall_statuses": demolition_statuses,
        "geometry_difference": geometry,
        "alignment": alignment,
        "gates": gates,
        "pass": passed,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"report": str(args.report), "gates": gates, "pass": passed}, ensure_ascii=False))
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
