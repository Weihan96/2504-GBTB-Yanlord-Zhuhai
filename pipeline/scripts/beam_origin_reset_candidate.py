#!/usr/bin/env python3
"""Build a beam origin-reset candidate without moving world geometry."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

import ifcopenshell
import ifcopenshell.api.geometry
import ifcopenshell.util.placement
import numpy as np

from direction_noise_candidate import all_product_geometry_difference
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm


def read_targets(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    expected = {
        "global_id",
        "target_x_mm",
        "target_y_mm",
        "target_z_mm",
        "anchor_kind",
        "basis",
        "confidence",
    }
    if not rows or set(rows[0]) != expected:
        raise RuntimeError(f"invalid beam origin target schema: {path}")
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
            raise RuntimeError(f"beam origin target is not integer millimetres: {row}")
        result.append(
            {
                "global_id": row["global_id"],
                "target_mm": target,
                "anchor_kind": row["anchor_kind"],
                "basis": row["basis"],
                "confidence": float(row["confidence"]),
            }
        )
    return result


def transform_point(point: Sequence[float], matrix: np.ndarray) -> tuple[float, ...]:
    coordinate = np.ones(4, dtype=float)
    coordinate[: len(point)] = np.array(point, dtype=float)
    transformed = matrix @ coordinate
    return tuple(float(value) for value in transformed[: len(point)])


def transform_axis_placement(
    placement: ifcopenshell.entity_instance, matrix: np.ndarray
) -> None:
    placement.Location.Coordinates = transform_point(
        placement.Location.Coordinates, matrix
    )


def transform_curve_points(
    item: ifcopenshell.entity_instance, matrix: np.ndarray
) -> None:
    if item.is_a("IfcIndexedPolyCurve"):
        item.Points.CoordList = tuple(
            transform_point(point, matrix) for point in item.Points.CoordList
        )
    elif item.is_a("IfcPolyline"):
        for point in item.Points:
            point.Coordinates = transform_point(point.Coordinates, matrix)
    else:
        raise RuntimeError(f"unsupported Axis item {item.is_a()} #{item.id()}")


def transform_body_item(
    model: ifcopenshell.file,
    item: ifcopenshell.entity_instance,
    matrix: np.ndarray,
    transformed_entity_ids: Optional[set[int]] = None,
    allowed_entity_ids: Optional[set[int]] = None,
) -> int:
    if transformed_entity_ids is None:
        transformed_entity_ids = set()
    if item.is_a("IfcBooleanResult"):
        return transform_body_item(
            model,
            item.FirstOperand,
            matrix,
            transformed_entity_ids,
            allowed_entity_ids,
        ) + transform_body_item(
            model,
            item.SecondOperand,
            matrix,
            transformed_entity_ids,
            allowed_entity_ids,
        )
    if item.is_a("IfcSweptAreaSolid"):
        if item.Position is None:
            item.Position = model.create_entity(
                "IfcAxis2Placement3D",
                Location=model.create_entity(
                    "IfcCartesianPoint",
                    Coordinates=tuple(float(value) for value in matrix[:3, 3]),
                ),
            )
            return 2
        transform_axis_placement(item.Position, matrix)
        return 0
    if item.is_a("IfcHalfSpaceSolid"):
        surface = item.BaseSurface
        if not surface.is_a("IfcElementarySurface"):
            raise RuntimeError(
                f"unsupported half-space surface {surface.is_a()} #{surface.id()}"
            )
        transform_axis_placement(surface.Position, matrix)
        return 0
    if item.is_a("IfcTessellatedFaceSet"):
        coordinates = item.Coordinates
        if allowed_entity_ids is not None and any(
            inverse.id() not in allowed_entity_ids
            for inverse in model.get_inverse(coordinates)
        ):
            raise RuntimeError(
                f"tessellated coordinates are shared outside the target product "
                f"{coordinates.is_a()} #{coordinates.id()}"
            )
        if coordinates.id() in transformed_entity_ids:
            return 0
        coordinates.CoordList = tuple(
            transform_point(point, matrix) for point in item.Coordinates.CoordList
        )
        transformed_entity_ids.add(coordinates.id())
        return 0
    if item.is_a("IfcBoundingBox"):
        item.Corner.Coordinates = transform_point(
            item.Corner.Coordinates, matrix
        )
        return 0
    if item.is_a("IfcIndexedPolyCurve"):
        if len(model.get_inverse(item.Points)) != 1:
            raise RuntimeError(
                f"shared curve point list cannot be transformed safely "
                f"{item.Points.is_a()} #{item.Points.id()}"
            )
        transform_curve_points(item, matrix)
        return 0
    if item.is_a("IfcPolyline"):
        if any(len(model.get_inverse(point)) != 1 for point in item.Points):
            raise RuntimeError(
                f"shared polyline point cannot be transformed safely "
                f"IfcPolyline #{item.id()}"
            )
        transform_curve_points(item, matrix)
        return 0
    if item.is_a("IfcGeometricCurveSet"):
        return sum(
            transform_body_item(
                model,
                element,
                matrix,
                transformed_entity_ids,
                allowed_entity_ids,
            )
            for element in item.Elements
        )
    if item.is_a("IfcMappedItem"):
        mapping_target = item.MappingTarget
        if len(model.get_inverse(mapping_target)) != 1:
            raise RuntimeError(
                f"shared mapping target cannot be transformed safely "
                f"{mapping_target.is_a()} #{mapping_target.id()}"
            )
        local_origin = mapping_target.LocalOrigin
        if len(model.get_inverse(local_origin)) != 1:
            raise RuntimeError(
                f"shared mapped local origin cannot be transformed safely "
                f"{local_origin.is_a()} #{local_origin.id()}"
            )
        local_origin.Coordinates = transform_point(
            local_origin.Coordinates, matrix
        )
        return 0
    raise RuntimeError(f"unsupported Body item {item.is_a()} #{item.id()}")


def point_segment_distance(
    point: np.ndarray, first: np.ndarray, second: np.ndarray
) -> float:
    segment = second - first
    denominator = float(np.dot(segment, segment))
    if denominator <= 1e-18:
        return float(np.linalg.norm(point - first))
    factor = max(
        0.0,
        min(1.0, float(np.dot(point - first, segment) / denominator)),
    )
    return float(np.linalg.norm(point - (first + factor * segment)))


def point_triangle_distance(
    point: Sequence[float],
    first: Sequence[float],
    second: Sequence[float],
    third: Sequence[float],
) -> float:
    p = np.array(point, dtype=float)
    a = np.array(first, dtype=float)
    b = np.array(second, dtype=float)
    c = np.array(third, dtype=float)
    ab = b - a
    ac = c - a
    normal = np.cross(ab, ac)
    normal_length = float(np.linalg.norm(normal))
    edge_distance = min(
        point_segment_distance(p, a, b),
        point_segment_distance(p, b, c),
        point_segment_distance(p, c, a),
    )
    if normal_length <= 1e-12:
        return edge_distance
    unit_normal = normal / normal_length
    plane_distance = float(np.dot(p - a, unit_normal))
    projected = p - plane_distance * unit_normal
    dot00 = float(np.dot(ab, ab))
    dot01 = float(np.dot(ab, ac))
    dot11 = float(np.dot(ac, ac))
    relative = projected - a
    dot20 = float(np.dot(relative, ab))
    dot21 = float(np.dot(relative, ac))
    denominator = dot00 * dot11 - dot01 * dot01
    if abs(denominator) <= 1e-12:
        return edge_distance
    u = (dot11 * dot20 - dot01 * dot21) / denominator
    v = (dot00 * dot21 - dot01 * dot20) / denominator
    if u >= -1e-9 and v >= -1e-9 and u + v <= 1.0 + 1e-9:
        return abs(plane_distance)
    return edge_distance


def point_mesh_distance(
    point: Sequence[float],
    vertices: Sequence[Sequence[float]],
    faces: Iterable[Sequence[int]],
) -> float:
    return min(
        point_triangle_distance(
            point,
            vertices[face[0]],
            vertices[face[1]],
            vertices[face[2]],
        )
        for face in faces
    )


def related_openings(beam: ifcopenshell.entity_instance) -> list[Any]:
    return [relation.RelatedOpeningElement for relation in beam.HasOpenings]


def reset_product_origin(
    model: ifcopenshell.file,
    row: dict[str, Any],
    tolerance_mm: float,
    expected_class: str,
    restore_related_openings: bool = True,
) -> dict[str, Any]:
    product = model.by_guid(row["global_id"])
    if product is None or not product.is_a(expected_class):
        raise RuntimeError(f"missing target {expected_class} {row['global_id']}")
    source_matrix = np.array(
        ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
        dtype=float,
    )
    target_matrix = source_matrix.copy()
    target_matrix[:3, 3] = row["target_mm"]
    local_transform = np.linalg.inv(target_matrix) @ source_matrix
    if not np.allclose(local_transform[:3, :3], np.identity(3), atol=1e-9):
        raise RuntimeError(
            f"origin reset unexpectedly rotates product {product.GlobalId}"
        )

    openings = related_openings(product)
    opening_world_matrices = {
        opening.GlobalId: np.array(
            ifcopenshell.util.placement.get_local_placement(opening.ObjectPlacement),
            dtype=float,
        )
        for opening in openings
    }
    identifiers = set()
    representation_entities_added = 0
    transformed_entity_ids: set[int] = set()
    allowed_entity_ids = {
        entity.id()
        for representation in product.Representation.Representations
        for entity in model.traverse(representation)
    }
    for representation in product.Representation.Representations:
        identifier = representation.RepresentationIdentifier
        identifiers.add(identifier)
        for item in representation.Items:
            if identifier == "Axis":
                transform_curve_points(item, local_transform)
            elif identifier in {"Body", "Box"}:
                representation_entities_added += transform_body_item(
                    model,
                    item,
                    local_transform,
                    transformed_entity_ids,
                    allowed_entity_ids,
                )
            else:
                raise RuntimeError(
                    f"unsupported representation {identifier} on {product.GlobalId}"
                )
    if "Body" not in identifiers or not identifiers.issubset(
        {"Axis", "Body", "Box"}
    ):
        raise RuntimeError(
            f"product {product.GlobalId} does not have a supported Body/Axis set"
        )

    ifcopenshell.api.geometry.edit_object_placement(
        model,
        product=product,
        matrix=target_matrix,
        is_si=False,
        should_transform_children=False,
    )
    if restore_related_openings:
        for opening in openings:
            ifcopenshell.api.geometry.edit_object_placement(
                model,
                product=opening,
                matrix=opening_world_matrices[opening.GlobalId],
                is_si=False,
                should_transform_children=False,
            )

    result_matrix = np.array(
        ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
        dtype=float,
    )
    settings = geometry_settings()
    vertices, faces = world_mesh_mm(settings, product)
    anchor_distance = point_mesh_distance(row["target_mm"], vertices, faces)
    placement_residual = float(
        np.max(np.abs(result_matrix[:3, 3] - row["target_mm"]))
    )
    if placement_residual > tolerance_mm:
        raise RuntimeError(
            f"product {product.GlobalId} placement residual {placement_residual} mm"
        )
    return {
        "global_id": product.GlobalId,
        "ifc_class": product.is_a(),
        "source_mm": source_matrix[:3, 3].tolist(),
        "target_mm": row["target_mm"].tolist(),
        "translation_by_axis_mm": (
            row["target_mm"] - source_matrix[:3, 3]
        ).tolist(),
        "anchor_kind": row["anchor_kind"],
        "basis": row["basis"],
        "confidence": row["confidence"],
        "placement_residual_mm": placement_residual,
        "anchor_surface_distance_mm": anchor_distance,
        "opening_global_ids": sorted(opening_world_matrices),
        "representation_entities_added": representation_entities_added,
    }


def reset_beam_origin(
    model: ifcopenshell.file, row: dict[str, Any], tolerance_mm: float
) -> dict[str, Any]:
    return reset_product_origin(model, row, tolerance_mm, "IfcBeam")


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
    results = [
        reset_beam_origin(candidate, row, args.tolerance_mm) for row in rows
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    all_geometry = all_product_geometry_difference(
        candidate, source, args.tolerance_mm
    )
    geometry_over_tolerance = [
        record
        for record in all_geometry["records"]
        if not record["within_tolerance"]
    ]
    target_placement_mismatches = []
    for row in rows:
        beam = candidate.by_guid(row["global_id"])
        matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(beam.ObjectPlacement),
            dtype=float,
        )
        residual = float(np.max(np.abs(matrix[:3, 3] - row["target_mm"])))
        if residual > args.tolerance_mm:
            target_placement_mismatches.append(
                {"global_id": beam.GlobalId, "residual_mm": residual}
            )
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "target_beams": len(rows),
        "target_placement_mismatches": target_placement_mismatches,
        "anchors_over_tolerance": sum(
            result["anchor_surface_distance_mm"] > args.tolerance_mm
            for result in results
        ),
        "product_geometry_over_tolerance": len(geometry_over_tolerance),
        "schema_equal": source.schema == candidate.schema,
        "entity_count_equal": len(list(source)) == len(list(candidate)),
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
    }
    gates["pass"] = (
        not gates["target_placement_mismatches"]
        and gates["anchors_over_tolerance"] == 0
        and gates["product_geometry_over_tolerance"] == 0
        and gates["schema_equal"]
        and gates["entity_count_equal"]
        and gates["root_global_ids_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-beam-origin-reset-candidate",
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
        "all_product_geometry": all_geometry,
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
