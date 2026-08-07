#!/usr/bin/env python3
"""Find existing integer geometry anchors for products with noninteger origins."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import ifcopenshell
import numpy as np

from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm


PREFERRED_MODULES_MM = (100.0, 50.0, 30.0, 20.0, 1.0)


def canonical_vertex(point: Sequence[float], precision_mm: float = 1e-6) -> tuple[float, ...]:
    return tuple(round(float(value) / precision_mm) * precision_mm for value in point)


def module_residual(value: float, module_mm: float) -> float:
    return abs(value - round(value / module_mm) * module_mm)


def preferred_module(
    point: Sequence[float], coordinate_tolerance_mm: float
) -> float | None:
    for module in PREFERRED_MODULES_MM:
        if max(module_residual(float(value), module) for value in point) <= coordinate_tolerance_mm:
            return module
    return None


def face_normal(
    vertices: Sequence[Sequence[float]], face: Sequence[int]
) -> np.ndarray:
    first, second, third = (np.array(vertices[index], dtype=float) for index in face)
    normal = np.cross(second - first, third - first)
    length = float(np.linalg.norm(normal))
    return normal / length if length > 1e-12 else np.zeros(3)


def physical_edges(
    vertices: Sequence[Sequence[float]], faces: Iterable[Sequence[int]]
) -> list[tuple[int, int]]:
    face_list = list(faces)
    edge_faces: dict[tuple[int, int], list[int]] = defaultdict(list)
    normals = [face_normal(vertices, face) for face in face_list]
    for face_index, face in enumerate(face_list):
        for start, end in (
            (face[0], face[1]),
            (face[1], face[2]),
            (face[2], face[0]),
        ):
            edge_faces[tuple(sorted((int(start), int(end))))].append(face_index)
    result = []
    for edge, related_faces in edge_faces.items():
        if len(related_faces) != 2:
            result.append(edge)
            continue
        first_normal = normals[related_faces[0]]
        second_normal = normals[related_faces[1]]
        if abs(float(np.dot(first_normal, second_normal))) < 1.0 - 1e-7:
            result.append(edge)
    return sorted(result)


def point_on_axis_aligned_edge(
    first: Sequence[float],
    second: Sequence[float],
    coordinate_tolerance_mm: float,
) -> tuple[tuple[float, ...], float] | None:
    differences = [abs(float(second[index]) - float(first[index])) for index in range(3)]
    varying_axes = [
        index for index, difference in enumerate(differences) if difference > coordinate_tolerance_mm
    ]
    if len(varying_axes) != 1:
        return None
    varying_axis = varying_axes[0]
    constants = [index for index in range(3) if index != varying_axis]
    low = min(float(first[varying_axis]), float(second[varying_axis]))
    high = max(float(first[varying_axis]), float(second[varying_axis]))
    midpoint = (low + high) / 2.0
    for module in PREFERRED_MODULES_MM:
        if any(
            module_residual(float(first[index]), module) > coordinate_tolerance_mm
            or abs(float(first[index]) - float(second[index])) > coordinate_tolerance_mm
            for index in constants
        ):
            continue
        coordinate = round(midpoint / module) * module
        if coordinate < low - coordinate_tolerance_mm or coordinate > high + coordinate_tolerance_mm:
            continue
        point = [
            round((float(first[index]) + float(second[index])) / 2.0 / module) * module
            for index in range(3)
        ]
        point[varying_axis] = coordinate
        return tuple(point), module
    return None


def point_on_non_axis_physical_edge(
    first: Sequence[float],
    second: Sequence[float],
    coordinate_tolerance_mm: float,
) -> tuple[tuple[float, ...], float] | None:
    start = np.array(first, dtype=float)
    end = np.array(second, dtype=float)
    vector = end - start
    varying_axes = [
        index
        for index, difference in enumerate(np.abs(vector))
        if difference > coordinate_tolerance_mm
    ]
    if len(varying_axes) <= 1:
        return None
    dominant_axis = max(varying_axes, key=lambda index: abs(vector[index]))
    low = min(start[dominant_axis], end[dominant_axis])
    high = max(start[dominant_axis], end[dominant_axis])
    midpoint = (low + high) / 2.0
    for module in PREFERRED_MODULES_MM:
        first_multiple = math.ceil((low - coordinate_tolerance_mm) / module)
        last_multiple = math.floor((high + coordinate_tolerance_mm) / module)
        candidates = sorted(
            range(first_multiple, last_multiple + 1),
            key=lambda multiple: abs(multiple * module - midpoint),
        )
        for multiple in candidates:
            coordinate = multiple * module
            factor = (coordinate - start[dominant_axis]) / vector[dominant_axis]
            if factor < -1e-12 or factor > 1.0 + 1e-12:
                continue
            point = start + factor * vector
            rounded = np.round(point / module) * module
            if float(np.max(np.abs(point - rounded))) <= coordinate_tolerance_mm:
                return tuple(float(value) for value in rounded), module
    return None


def point_in_triangle_2d(
    point: np.ndarray,
    triangle: np.ndarray,
    tolerance_mm: float,
) -> bool:
    first, second, third = triangle
    edge_a = second - first
    edge_b = third - first
    relative = point - first
    denominator = edge_a[0] * edge_b[1] - edge_a[1] * edge_b[0]
    if abs(float(denominator)) <= 1e-12:
        return False
    u = (relative[0] * edge_b[1] - relative[1] * edge_b[0]) / denominator
    v = (edge_a[0] * relative[1] - edge_a[1] * relative[0]) / denominator
    scale = max(float(np.linalg.norm(edge_a)), float(np.linalg.norm(edge_b)), 1.0)
    epsilon = tolerance_mm / scale
    return u >= -epsilon and v >= -epsilon and u + v <= 1.0 + epsilon


def point_on_axis_aligned_surface(
    vertices: Sequence[Sequence[float]],
    faces: Sequence[Sequence[int]],
    coordinate_tolerance_mm: float,
) -> tuple[tuple[float, ...], float] | None:
    """Find an integer lattice point inside an axis-aligned rendered face.

    This fallback is intentionally restricted by the caller to construction
    surfaces such as IfcSlab and IfcCovering. It does not turn arbitrary mesh
    triangles on furniture, services, or sanitary fixtures into installation
    semantics.
    """

    candidates: list[tuple[int, float, tuple[float, ...], float]] = []
    for face in faces:
        triangle = np.array([vertices[index] for index in face], dtype=float)
        normal = face_normal(vertices, face)
        dominant_axis = int(np.argmax(np.abs(normal)))
        if abs(float(normal[dominant_axis])) < 1.0 - 1e-7:
            continue
        free_axes = [axis for axis in range(3) if axis != dominant_axis]
        if any(abs(float(normal[axis])) > 1e-7 for axis in free_axes):
            continue
        plane_coordinate = float(np.mean(triangle[:, dominant_axis]))
        if float(np.max(np.abs(triangle[:, dominant_axis] - plane_coordinate))) > coordinate_tolerance_mm:
            continue
        projected = triangle[:, free_axes]
        first_edge = projected[1] - projected[0]
        second_edge = projected[2] - projected[0]
        area = abs(
            float(
                first_edge[0] * second_edge[1]
                - first_edge[1] * second_edge[0]
            )
        ) / 2.0
        if area <= 1e-9:
            continue
        centroid = np.mean(projected, axis=0)
        for module_index, module in enumerate(PREFERRED_MODULES_MM):
            if module_residual(plane_coordinate, module) > coordinate_tolerance_mm:
                continue
            base = np.round(centroid / module) * module
            grid_points = []
            for first_offset in (-2, -1, 0, 1, 2):
                for second_offset in (-2, -1, 0, 1, 2):
                    grid_points.append(
                        base
                        + np.array(
                            [first_offset * module, second_offset * module],
                            dtype=float,
                        )
                    )
            for projected_point in grid_points:
                if not point_in_triangle_2d(
                    projected_point, projected, coordinate_tolerance_mm
                ):
                    continue
                point = [0.0, 0.0, 0.0]
                point[dominant_axis] = round(plane_coordinate / module) * module
                point[free_axes[0]] = float(projected_point[0])
                point[free_axes[1]] = float(projected_point[1])
                candidates.append(
                    (
                        module_index,
                        -area,
                        tuple(point),
                        module,
                    )
                )
                break
    if not candidates:
        return None
    _, _, point, module = min(candidates)
    return point, module


def find_anchor(
    vertices: Sequence[Sequence[float]],
    faces: Sequence[Sequence[int]],
    coordinate_tolerance_mm: float,
    allow_axis_aligned_surface: bool = False,
) -> dict[str, Any] | None:
    unique_vertices = sorted({canonical_vertex(vertex) for vertex in vertices})
    vertex_candidates = []
    for vertex in unique_vertices:
        module = preferred_module(vertex, coordinate_tolerance_mm)
        if module is not None:
            vertex_candidates.append((PREFERRED_MODULES_MM.index(module), vertex, module))
    if vertex_candidates:
        _, source_point, module = min(vertex_candidates)
        point = tuple(
            round(float(value) / module) * module for value in source_point
        )
        snap_distance_mm = max(
            abs(float(source_point[index]) - float(point[index]))
            for index in range(3)
        )
        return {
            "anchor_kind": "integer_point_within_tolerance_of_existing_vertex",
            "point_mm": point,
            "source_geometry_point_mm": source_point,
            "snap_distance_mm": snap_distance_mm,
            "module_mm": module,
            "evidence": "integer lattice point is within the configured tolerance of an existing triangulated world vertex",
        }
    edge_candidates = []
    for start, end in physical_edges(vertices, faces):
        candidate = point_on_axis_aligned_edge(
            vertices[start], vertices[end], coordinate_tolerance_mm
        )
        if candidate is None:
            continue
        point, module = candidate
        edge_candidates.append((PREFERRED_MODULES_MM.index(module), point, module))
    if edge_candidates:
        _, point, module = min(edge_candidates)
        return {
            "anchor_kind": "integer_point_on_physical_axis_edge",
            "point_mm": point,
            "module_mm": module,
            "evidence": "point lies on a non-triangulation axis-aligned mesh edge",
        }
    non_axis_edge_candidates = []
    for start, end in physical_edges(vertices, faces):
        candidate = point_on_non_axis_physical_edge(
            vertices[start], vertices[end], coordinate_tolerance_mm
        )
        if candidate is None:
            continue
        point, module = candidate
        non_axis_edge_candidates.append(
            (PREFERRED_MODULES_MM.index(module), point, module)
        )
    if non_axis_edge_candidates:
        _, point, module = min(non_axis_edge_candidates)
        return {
            "anchor_kind": "integer_point_on_physical_non_axis_edge",
            "point_mm": point,
            "module_mm": module,
            "evidence": "point lies on a non-triangulation sloped or diagonal mesh edge",
        }
    if allow_axis_aligned_surface:
        candidate = point_on_axis_aligned_surface(
            vertices, faces, coordinate_tolerance_mm
        )
        if candidate is not None:
            point, module = candidate
            return {
                "anchor_kind": "integer_point_on_axis_aligned_surface",
                "point_mm": point,
                "module_mm": module,
                "evidence": "point lies inside an axis-aligned rendered construction surface",
            }
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--origin-report", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--origin-tolerance-mm", type=float, default=0.1)
    parser.add_argument("--coordinate-tolerance-mm", type=float, default=0.01)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.coordinate_tolerance_mm <= args.origin_tolerance_mm:
        raise SystemExit(
            "coordinate tolerance must be positive and no larger than origin tolerance"
        )
    input_path = args.input.resolve()
    input_hash = sha256(input_path)
    origin_report = json.loads(args.origin_report.read_text(encoding="utf-8"))
    if origin_report["source"]["sha256"] != input_hash:
        raise RuntimeError("remaining-origin report is stale for the current IFC")
    source_rows = [
        row
        for row in origin_report["records"]
        if not row["within_review_tolerance"]
    ]
    model = ifcopenshell.open(input_path)
    settings = geometry_settings()
    records = []
    for row in source_rows:
        product = model.by_guid(row["global_id"])
        if product is None or not getattr(product, "Representation", None):
            records.append(
                {
                    **row,
                    "shape_status": "no_representation",
                    "anchor": None,
                    "automatic_write_allowed": False,
                }
            )
            continue
        try:
            vertices, faces = world_mesh_mm(settings, product)
        except Exception as error:  # noqa: BLE001 - error is evidence for this object
            records.append(
                {
                    **row,
                    "shape_status": "shape_error",
                    "shape_error": f"{type(error).__name__}: {error}",
                    "anchor": None,
                    "automatic_write_allowed": False,
                }
            )
            continue
        anchor = find_anchor(
            vertices,
            faces,
            args.coordinate_tolerance_mm,
            allow_axis_aligned_surface=(
                product.is_a("IfcSlab") or product.is_a("IfcCovering")
            ),
        )
        records.append(
            {
                **row,
                "shape_status": "ok",
                "mesh_vertices": len(vertices),
                "mesh_triangles": len(faces),
                "anchor": anchor,
                "automatic_write_allowed": False,
                "basis": (
                    "geometry anchor exists; class-specific anchor semantics and a no-world-movement candidate remain required"
                    if anchor
                    else "no existing integer vertex, physical edge point, or approved construction-surface point was proven"
                ),
            }
        )
    anchors = [record for record in records if record.get("anchor")]
    summary = {
        "source_origins_over_tolerance": len(source_rows),
        "records": len(records),
        "anchors_found": len(anchors),
        "anchors_not_found": len(records) - len(anchors),
        "by_anchor_kind": dict(
            sorted(Counter(record["anchor"]["anchor_kind"] for record in anchors).items())
        ),
        "by_module_mm": dict(
            sorted(
                Counter(str(record["anchor"]["module_mm"]) for record in anchors).items(),
                key=lambda item: -float(item[0]),
            )
        ),
        "by_ifc_class": dict(
            sorted(Counter(record["ifc_class"] for record in anchors).items())
        ),
        "by_review_subgroup": dict(
            sorted(Counter(record["review_subgroup"] for record in anchors).items())
        ),
        "automatic_write_allowed": 0,
    }
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-integer-geometry-anchor-audit",
        "source": {
            "path": str(input_path),
            "sha256": input_hash,
            "schema": model.schema,
        },
        "origin_report": str(args.origin_report.resolve()),
        "origin_tolerance_mm": args.origin_tolerance_mm,
        "coordinate_tolerance_mm": args.coordinate_tolerance_mm,
        "preferred_modules_mm": PREFERRED_MODULES_MM,
        "records": records,
        "summary": summary,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **summary}, ensure_ascii=False))


if __name__ == "__main__":
    main()
