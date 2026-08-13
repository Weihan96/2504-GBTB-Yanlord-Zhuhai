#!/usr/bin/env python3
"""Build a geometry-neutral P0 IDS metadata candidate.

The formal IFC is read-only. The candidate adds only mechanically supported
wall base-quantity Length values and broad DOOR/WINDOW predefined types.
Wall tags, unhosted-door sizes, operation types, placements, representations,
and fills/voids relationships remain untouched.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api
import ifcopenshell.util.element
import ifcopenshell.util.placement
import ifcopenshell.util.representation

from a103_wall_plan_candidate import geometry_settings, world_bbox_mm


PROJECT_ROOT = Path(__file__).resolve().parents[2]
FORMAL_IFC = PROJECT_ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
A103_REGISTER = PROJECT_ROOT / "pipeline/decisions/a103-wall-status-review.csv"
A104_REGISTER = PROJECT_ROOT / "pipeline/decisions/a104-door-window-review.csv"
IDS_PATH = PROJECT_ROOT / "pipeline/ids/p0-construction-information.ids"
IDS_VALIDATOR = PROJECT_ROOT / "pipeline/scripts/ids_validate.py"
EXPECTED_SOURCE_IDS_CHECKS = 593
EXPECTED_SOURCE_IDS_PASS = {392, 499, 587}
MINIMUM_CANDIDATE_IDS_PASS = 499
EXPECTED_WALL_COUNT = 88
EXPECTED_DOOR_COUNT = 8
EXPECTED_WINDOW_COUNT = 11

# These five values are already confirmed project dimensions. Their retained
# Axis representations contain legacy sub-millimetre noise and are evidence,
# not a reason to reopen the approved dimensions.
CONFIRMED_LENGTHS_MM = {
    "18Qhf59crDCAy$scVOCK92": 800.0,
    "1aQiUP$CH5uhQBoEX8NjwK": 800.0,
    "21vcsFmEvC2xpOEhH31C0E": 1350.0,
    "2amSzxJIb9ceddTJOslmHk": 1100.0,
    "3OVQygdDn17huGOgJJFTOY": 960.0,
}

# These two final-built walls have no Axis representation. Their unique
# horizontal Body long direction is mechanically checkable at these values.
BODY_LENGTHS_MM = {
    "2Ca2tGerPBj8kvUTjlUtDl": 650.0,
    "3pAfMJYxPBlwR30CstZ4kK": 1360.0,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def forward_entities(roots: tuple[Any, ...]) -> dict[int, str]:
    queue = deque(root for root in roots if root is not None)
    entities: dict[int, str] = {}
    while queue:
        entity = queue.popleft()
        if not hasattr(entity, "id") or entity.id() in entities:
            continue
        entities[entity.id()] = str(entity)
        for value in entity:
            if hasattr(value, "id"):
                queue.append(value)
            elif isinstance(value, (tuple, list)):
                queue.extend(item for item in value if hasattr(item, "id"))
    return entities


def physical_graph(product: Any) -> dict[int, str]:
    representations = tuple(product.Representation.Representations) if product.Representation else ()
    return forward_entities((product.ObjectPlacement, *representations))


def matrix_max_delta(first: Any, second: Any) -> float:
    return max(
        abs(float(first[row][column]) - float(second[row][column]))
        for row in range(4)
        for column in range(4)
    )


def relationship_fingerprint(model: ifcopenshell.file) -> dict[str, list[tuple[str, str, str]]]:
    return {
        "fills": sorted(
            (relation.GlobalId, relation.RelatingOpeningElement.GlobalId, relation.RelatedBuildingElement.GlobalId)
            for relation in model.by_type("IfcRelFillsElement")
        ),
        "voids": sorted(
            (relation.GlobalId, relation.RelatingBuildingElement.GlobalId, relation.RelatedOpeningElement.GlobalId)
            for relation in model.by_type("IfcRelVoidsElement")
        ),
    }


def axis_representation(wall: Any) -> Any | None:
    return (
        ifcopenshell.util.representation.get_representation(wall, "Plan", "Axis", "GRAPH_VIEW")
        or ifcopenshell.util.representation.get_representation(wall, "Model", "Axis", "GRAPH_VIEW")
    )


def axis_length_mm(wall: Any) -> float:
    representation = axis_representation(wall)
    if representation is None or len(representation.Items) != 1:
        raise RuntimeError(f"{wall.GlobalId}: expected exactly one Axis item")
    curve = representation.Items[0]
    if not curve.is_a("IfcIndexedPolyCurve"):
        raise RuntimeError(f"{wall.GlobalId}: unsupported Axis type {curve.is_a()}")
    points = [tuple(map(float, point)) for point in curve.Points.CoordList]
    if len(points) != 2 or len(points[0]) != len(points[1]):
        raise RuntimeError(f"{wall.GlobalId}: expected a two-point Axis")
    return math.dist(points[0], points[1])


def wall_lengths(model: ifcopenshell.file, wall_ids: set[str]) -> list[dict[str, Any]]:
    settings = geometry_settings()
    records: list[dict[str, Any]] = []
    for global_id in sorted(wall_ids):
        wall = model.by_guid(global_id)
        if global_id in CONFIRMED_LENGTHS_MM:
            length = CONFIRMED_LENGTHS_MM[global_id]
            axis_observation = axis_length_mm(wall)
            basis = "confirmed_project_dimension"
        elif global_id in BODY_LENGTHS_MM:
            length = BODY_LENGTHS_MM[global_id]
            if axis_representation(wall) is not None:
                raise RuntimeError(f"{global_id}: Body fallback is invalid because an Axis now exists")
            dimensions = world_bbox_mm(settings, wall)["dimensions_mm"]
            if abs(dimensions[0] - dimensions[1]) <= 0.1:
                raise RuntimeError(f"{global_id}: Body has no unique horizontal long direction")
            horizontal_long = max(dimensions[:2])
            if abs(horizontal_long - length) > 0.1:
                raise RuntimeError(
                    f"{global_id}: Body long direction {horizontal_long} mm does not support {length} mm"
                )
            axis_observation = None
            basis = "unique_body_horizontal_long_direction"
        else:
            length = axis_length_mm(wall)
            axis_observation = length
            basis = "formal_axis_two_point_distance"
        if not math.isfinite(length) or length <= 0:
            raise RuntimeError(f"{global_id}: invalid Length {length}")
        records.append(
            {
                "global_id": global_id,
                "length_mm": length,
                "basis": basis,
                "axis_observation_mm": axis_observation,
            }
        )
    return records


def add_wall_lengths(model: ifcopenshell.file, records: list[dict[str, Any]]) -> None:
    for record in records:
        wall = model.by_guid(record["global_id"])
        existing = ifcopenshell.util.element.get_psets(wall, qtos_only=True).get("Qto_WallBaseQuantities")
        qto = model.by_id(int(existing["id"])) if existing else ifcopenshell.api.run(
            "pset.add_qto", model, product=wall, name="Qto_WallBaseQuantities"
        )
        ifcopenshell.api.run(
            "pset.edit_qto", model, qto=qto, properties={"Length": record["length_mm"]}
        )


def run_ids(ifc_path: Path) -> dict[str, Any]:
    process = subprocess.run(
        [
            sys.executable,
            str(IDS_VALIDATOR),
            "--input",
            str(ifc_path),
            "--ids",
            str(IDS_PATH),
            "--stdout-only",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
    )
    if process.returncode != 0:
        raise RuntimeError(f"IDS validator failed: {process.stderr or process.stdout}")
    return json.loads(process.stdout)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--output", type=Path, default=PROJECT_ROOT / "build/candidates/2504-GBTB-p0-ids-metadata.ifc")
    parser.add_argument("--report", type=Path, default=PROJECT_ROOT / "build/ids/p0-metadata-candidate.json")
    parser.add_argument("--expected-ifc-sha256")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_path = args.input.resolve()
    output_path = args.output.resolve()
    report_path = args.report.resolve()
    if source_path == output_path:
        raise RuntimeError("candidate output must not overwrite the formal IFC")
    source_hash = sha256(source_path)
    if args.expected_ifc_sha256 and source_hash != args.expected_ifc_sha256:
        raise RuntimeError(
            f"formal IFC hash changed: expected {args.expected_ifc_sha256}, got {source_hash}"
        )
    source = ifcopenshell.open(source_path)
    if source.schema != "IFC4":
        raise RuntimeError(f"unexpected schema {source.schema}")

    a103 = read_csv(A103_REGISTER)
    a104 = read_csv(A104_REGISTER)
    wall_ids = {row["global_id"] for row in a103}
    if len(a103) != EXPECTED_WALL_COUNT or len(wall_ids) != EXPECTED_WALL_COUNT:
        raise RuntimeError("A-103 register must contain 88 unique final-built walls")
    final_wall_ids = {
        wall.GlobalId
        for wall in source.by_type("IfcWall")
        if str(
            ifcopenshell.util.element.get_psets(wall)
            .get("Pset_WallCommon", {})
            .get("Status")
            or ""
        ).strip().upper()
        not in {"DEMOLISH", "DEMOLISHED"}
    }
    if wall_ids != final_wall_ids:
        raise RuntimeError("A-103 register does not exactly match the formal final-built wall set")
    if {row["source_ifc_sha256"] for row in a103} != {source_hash}:
        raise RuntimeError("A-103 register is stale against the formal IFC")
    door_rows = [row for row in a104 if row["ifc_class"] == "IfcDoor"]
    window_rows = [row for row in a104 if row["ifc_class"] == "IfcWindow"]
    if len(door_rows) != EXPECTED_DOOR_COUNT or len(window_rows) != EXPECTED_WINDOW_COUNT:
        raise RuntimeError("A-104 register must contain 8 doors and 11 windows")
    a104_ids = [row["global_id"] for row in a104]
    if len(a104_ids) != len(set(a104_ids)):
        raise RuntimeError("A-104 register contains duplicate GlobalIds")
    for row in door_rows + window_rows:
        entity = source.by_guid(row["global_id"])
        if entity is None or not entity.is_a(row["ifc_class"]):
            raise RuntimeError(f"{row['global_id']}: A-104 IFC class does not match the formal IFC")
        expected_predefined_type = "DOOR" if row["ifc_class"] == "IfcDoor" else "WINDOW"
        if row["candidate_predefined_type"] != expected_predefined_type:
            raise RuntimeError(
                f"{row['global_id']}: candidate_predefined_type must be {expected_predefined_type}"
            )
    if {row["source_ifc_sha256"] for row in a104} != {source_hash}:
        raise RuntimeError("A-104 register is stale against the formal IFC")

    lengths = wall_lengths(source, wall_ids)
    source_ids = run_ids(source_path)
    source_relationships = relationship_fingerprint(source)
    source_products = {
        product.GlobalId: product
        for product in source.by_type("IfcProduct")
        if not (product.is_a("IfcAnnotation") and product.ObjectType == "DRAWING")
    }
    source_graphs = {global_id: physical_graph(product) for global_id, product in source_products.items()}
    source_matrices = {
        global_id: ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement)
        for global_id, product in source_products.items()
    }
    source_door_operation_types = {
        row["global_id"]: source.by_guid(row["global_id"]).OperationType for row in door_rows
    }

    candidate = ifcopenshell.open(source_path)
    add_wall_lengths(candidate, lengths)
    for row in door_rows:
        candidate.by_guid(row["global_id"]).PredefinedType = "DOOR"
    for row in window_rows:
        candidate.by_guid(row["global_id"]).PredefinedType = "WINDOW"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(output_path)
    if sha256(source_path) != source_hash:
        raise RuntimeError("formal IFC changed while building the candidate")

    candidate = ifcopenshell.open(output_path)
    candidate_by_guid = {product.GlobalId: product for product in candidate.by_type("IfcProduct")}
    geometry_records = []
    for global_id, source_product in sorted(source_products.items()):
        product = candidate_by_guid.get(global_id)
        graph_equal = product is not None and source_graphs[global_id] == physical_graph(product)
        placement_delta = (
            matrix_max_delta(
                source_matrices[global_id],
                ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
            )
            if product is not None
            else math.inf
        )
        geometry_records.append(
            {
                "global_id": global_id,
                "ifc_class": source_product.is_a(),
                "geometry_graph_exact": graph_equal,
                "placement_matrix_max_delta": placement_delta,
                "world_geometry_delta_mm": 0.0 if graph_equal and placement_delta == 0.0 else math.inf,
            }
        )

    candidate_ids = run_ids(output_path)
    candidate_relationships = relationship_fingerprint(candidate)
    qto_lengths = {
        global_id: ifcopenshell.util.element.get_psets(candidate.by_guid(global_id), qtos_only=True)["Qto_WallBaseQuantities"]["Length"]
        for global_id in wall_ids
    }
    doors_correct = all(candidate.by_guid(row["global_id"]).PredefinedType == "DOOR" for row in door_rows)
    windows_correct = all(candidate.by_guid(row["global_id"]).PredefinedType == "WINDOW" for row in window_rows)
    operation_types_unchanged = all(
        candidate.by_guid(global_id).OperationType == value
        for global_id, value in source_door_operation_types.items()
    )
    length_deltas = [abs(qto_lengths[record["global_id"]] - record["length_mm"]) for record in lengths]
    changed_geometry = [record for record in geometry_records if record["world_geometry_delta_mm"] != 0.0]
    gates = {
        "source_ifc_hash_preserved": sha256(source_path) == source_hash,
        "wall_count": len(lengths),
        "wall_lengths_exact": max(length_deltas, default=0.0) <= 1e-9,
        "length_basis_counts": dict(Counter(record["basis"] for record in lengths)),
        "door_count": len(door_rows),
        "window_count": len(window_rows),
        "door_predefined_types_correct": doors_correct,
        "window_predefined_types_correct": windows_correct,
        "operation_types_unchanged": operation_types_unchanged,
        "fills_voids_relationships_unchanged": candidate_relationships == source_relationships,
        "protected_product_count": len(geometry_records),
        "protected_products_geometry_exact": len(changed_geometry) == 0,
        "maximum_world_geometry_delta_mm": 0.0 if not changed_geometry else math.inf,
        "maximum_placement_matrix_delta": max(
            (record["placement_matrix_max_delta"] for record in geometry_records), default=0.0
        ),
        "source_ids_checks": source_ids["total_checks"],
        "source_ids_pass": source_ids["total_checks_pass"],
        "candidate_ids_checks": candidate_ids["total_checks"],
        "candidate_ids_pass": candidate_ids["total_checks_pass"],
        "candidate_ids_fail": candidate_ids["total_checks_fail"],
        "formal_ifc_write_allowed": False,
    }
    passed = (
        gates["source_ifc_hash_preserved"]
        and gates["wall_count"] == EXPECTED_WALL_COUNT
        and gates["wall_lengths_exact"]
        and gates["length_basis_counts"]
        == {
            "confirmed_project_dimension": len(CONFIRMED_LENGTHS_MM),
            "formal_axis_two_point_distance": EXPECTED_WALL_COUNT - len(CONFIRMED_LENGTHS_MM) - len(BODY_LENGTHS_MM),
            "unique_body_horizontal_long_direction": len(BODY_LENGTHS_MM),
        }
        and gates["door_count"] == EXPECTED_DOOR_COUNT
        and gates["window_count"] == EXPECTED_WINDOW_COUNT
        and gates["door_predefined_types_correct"]
        and gates["window_predefined_types_correct"]
        and gates["operation_types_unchanged"]
        and gates["fills_voids_relationships_unchanged"]
        and gates["protected_products_geometry_exact"]
        and gates["maximum_placement_matrix_delta"] == 0.0
        and gates["source_ids_checks"] == EXPECTED_SOURCE_IDS_CHECKS
        and gates["source_ids_pass"] in EXPECTED_SOURCE_IDS_PASS
        and gates["candidate_ids_checks"] == EXPECTED_SOURCE_IDS_CHECKS
        and gates["candidate_ids_pass"] == max(
            gates["source_ids_pass"], MINIMUM_CANDIDATE_IDS_PASS
        )
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_p0_ids_metadata_candidate",
        "source_ifc_sha256": source_hash,
        "source": {"ifc": str(source_path), "ifc_sha256": source_hash, "schema": source.schema},
        "candidate": {"ifc": str(output_path), "ifc_sha256": sha256(output_path), "schema": candidate.schema},
        "scope": {
            "wall_length_records": lengths,
            "door_global_ids": [row["global_id"] for row in door_rows],
            "window_global_ids": [row["global_id"] for row in window_rows],
            "excluded": ["wall Tag", "door/window OperationType", "M05-M07 OverallWidth/OverallHeight", "host relationships", "placement", "representation"],
        },
        "ids": {"source": source_ids, "candidate": candidate_ids},
        "geometry_records": geometry_records,
        "gates": gates,
        "pass": passed,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(report_path), "candidate": str(output_path), "gates": gates, "pass": passed}, ensure_ascii=False))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
