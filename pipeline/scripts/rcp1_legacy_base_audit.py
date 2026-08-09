#!/usr/bin/env python3
"""Compare the formal RCP1 geometry with the legacy Blender design base.

Run this script from Blender after opening ``2504_lowpoly.blend``.  It reads
the formal IFC and existing RCP1 coordination report but never writes IFC or
the legacy blend file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any

import bpy
import ifcopenshell
import ifcopenshell.geom
from mathutils import Vector
from mathutils.bvhtree import BVHTree


EXPECTED_IFC_SHA256 = "7521c09991f3d0c7b7d91ca2324fd55ad961d8e32e9e3e9a9777a4cc19b06e81"

FORMAL_TO_LEGACY = {
    "0Ik2RcgGbFOhdYTJPgh5AQ": "Liquid Living Room",
    "0f2ZLauDH8lRnYj6oervDm": "Gas Living Room",
    "0hHnbLj0X4jPDz4o3QJo1l": "Liquid.001",
    "10Wm8ivdX7dAVfz4cV8l5Q": "Drain Pipe",
    "1hZRB0eOX8OA8rcjke67P0": "Gas.001",
    "1QBdVekDnBsOleyo9PM6rT": "BL|RPIZ-22FSN6QD_curve_.009",
    "1yW7DASIz8qA$2j8z9tdl2": "BL|RPIZ-50FSN6QD_curve_.001",
}

PIPE_IDS = {
    "0Ik2RcgGbFOhdYTJPgh5AQ",
    "0f2ZLauDH8lRnYj6oervDm",
    "0hHnbLj0X4jPDz4o3QJo1l",
    "10Wm8ivdX7dAVfz4cV8l5Q",
    "1hZRB0eOX8OA8rcjke67P0",
}

AC_TRANSLATION_ANCHORS = {
    "06GpMzzWj1XQobAhD35cgU": "BL|RPIZ-22FSN6QD_curve_.008",
    "1QBdVekDnBsOleyo9PM6rT": "BL|RPIZ-22FSN6QD_curve_.009",
    "33chLv3TzEOhKalIJJAPNF": "BL|RPIZ-22FSN6QD_curve_.006",
}

LEGACY_EAST_CANDIDATE = "BL|RPIZ-22FSN6QD_curve_.001"
LEGACY_INTERSECTION_PAIRS = {
    tuple(sorted(pair))
    for pair in [
        ("0Ik2RcgGbFOhdYTJPgh5AQ", "0f2ZLauDH8lRnYj6oervDm"),
        ("0Ik2RcgGbFOhdYTJPgh5AQ", "10Wm8ivdX7dAVfz4cV8l5Q"),
        ("0f2ZLauDH8lRnYj6oervDm", "10Wm8ivdX7dAVfz4cV8l5Q"),
        ("0f2ZLauDH8lRnYj6oervDm", "1yW7DASIz8qA$2j8z9tdl2"),
        ("0hHnbLj0X4jPDz4o3QJo1l", "10Wm8ivdX7dAVfz4cV8l5Q"),
        ("0hHnbLj0X4jPDz4o3QJo1l", "1QBdVekDnBsOleyo9PM6rT"),
        ("0hHnbLj0X4jPDz4o3QJo1l", "1hZRB0eOX8OA8rcjke67P0"),
        ("10Wm8ivdX7dAVfz4cV8l5Q", "1hZRB0eOX8OA8rcjke67P0"),
        ("10Wm8ivdX7dAVfz4cV8l5Q", "1yW7DASIz8qA$2j8z9tdl2"),
        ("1QBdVekDnBsOleyo9PM6rT", "1hZRB0eOX8OA8rcjke67P0"),
    ]
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formal-ifc", type=Path, required=True)
    parser.add_argument("--formal-report", type=Path, required=True)
    parser.add_argument("--coordination-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--candidate-tolerance-mm", type=float, default=10.0)
    parser.add_argument(
        "--current-candidate-centre-mm",
        type=float,
        nargs=3,
        default=(-1068.494529, -905.287176, 2620.516300),
    )
    arguments = sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else []
    return parser.parse_args(arguments)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def world_mesh_from_blender(obj: bpy.types.Object) -> tuple[list[Vector], list[tuple[int, ...]]]:
    depsgraph = bpy.context.evaluated_depsgraph_get()
    evaluated = obj.evaluated_get(depsgraph)
    mesh = evaluated.to_mesh(preserve_all_data_layers=False, depsgraph=depsgraph)
    try:
        vertices = [evaluated.matrix_world @ vertex.co for vertex in mesh.vertices]
        polygons = [tuple(polygon.vertices) for polygon in mesh.polygons if len(polygon.vertices) >= 3]
    finally:
        evaluated.to_mesh_clear()
    if not vertices or not polygons:
        raise RuntimeError(f"legacy object has no evaluated mesh: {obj.name}")
    return vertices, polygons


def world_mesh_from_ifc(
    model: ifcopenshell.file,
    settings: ifcopenshell.geom.settings,
    global_id: str,
) -> tuple[list[Vector], list[tuple[int, int, int]]]:
    product = model.by_guid(global_id)
    shape = ifcopenshell.geom.create_shape(settings, product)
    values = shape.geometry.verts
    vertices = [Vector(values[index : index + 3]) for index in range(0, len(values), 3)]
    faces = shape.geometry.faces
    polygons = [tuple(faces[index : index + 3]) for index in range(0, len(faces), 3)]
    if not vertices or not polygons:
        raise RuntimeError(f"formal IFC object has no Body mesh: {global_id}")
    return vertices, polygons


def bbox(vertices: list[Vector]) -> dict[str, list[float]]:
    minimum = [min(vertex[axis] for vertex in vertices) for axis in range(3)]
    maximum = [max(vertex[axis] for vertex in vertices) for axis in range(3)]
    return {
        "min_mm": [value * 1000.0 for value in minimum],
        "max_mm": [value * 1000.0 for value in maximum],
        "centre_mm": [(minimum[axis] + maximum[axis]) * 500.0 for axis in range(3)],
        "dimensions_mm": [(maximum[axis] - minimum[axis]) * 1000.0 for axis in range(3)],
    }


def component_dimensions(
    vertices: list[Vector], polygons: list[tuple[int, ...]]
) -> list[list[float]]:
    parents = list(range(len(vertices)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(first: int, second: int) -> None:
        first_root = find(first)
        second_root = find(second)
        if first_root != second_root:
            parents[second_root] = first_root

    used = set()
    for polygon in polygons:
        used.update(polygon)
        for second in polygon[1:]:
            union(polygon[0], second)

    groups: dict[int, list[Vector]] = {}
    for index in used:
        groups.setdefault(find(index), []).append(vertices[index])

    dimensions = []
    for group in groups.values():
        bounds = bbox(group)
        dimensions.append(sorted(bounds["dimensions_mm"]))
    return sorted(dimensions)


def best_component_difference_mm(
    formal: list[list[float]], legacy: list[list[float]]
) -> float:
    if len(formal) != len(legacy):
        return math.inf
    remaining = list(legacy)
    maximum_difference = 0.0
    for formal_component in sorted(formal):
        best_index, best_difference = min(
            (
                index,
                max(
                    abs(formal_component[axis] - legacy_component[axis])
                    for axis in range(3)
                ),
            )
            for index, legacy_component in enumerate(remaining)
        )
        maximum_difference = max(maximum_difference, best_difference)
        remaining.pop(best_index)
    return maximum_difference


def object_centre_mm(name: str) -> list[float]:
    obj = bpy.data.objects.get(name)
    if obj is None:
        raise RuntimeError(f"missing legacy object: {name}")
    vertices, _ = world_mesh_from_blender(obj)
    return bbox(vertices)["centre_mm"]


def main() -> int:
    args = parse_args()
    if args.tolerance_mm <= 0 or args.candidate_tolerance_mm <= 0:
        raise RuntimeError("tolerances must be positive")

    formal_sha = sha256(args.formal_ifc)
    if formal_sha != EXPECTED_IFC_SHA256:
        raise RuntimeError(f"formal IFC SHA drift: {formal_sha}")

    formal_report = json.loads(args.formal_report.read_text(encoding="utf-8"))
    coordination = json.loads(args.coordination_report.read_text(encoding="utf-8"))
    if formal_report["source"]["sha256"] != formal_sha:
        raise RuntimeError("formal RCP1 report is stale")
    if coordination["source"]["sha256"] != formal_sha:
        raise RuntimeError("formal coordination report is stale")

    model = ifcopenshell.open(args.formal_ifc)
    geometry_settings = ifcopenshell.geom.settings()
    geometry_settings.set(geometry_settings.USE_WORLD_COORDS, True)

    required_names = set(FORMAL_TO_LEGACY.values()) | set(AC_TRANSLATION_ANCHORS.values()) | {
        LEGACY_EAST_CANDIDATE
    }
    missing_names = sorted(name for name in required_names if bpy.data.objects.get(name) is None)

    pipe_comparisons = []
    legacy_meshes: dict[str, tuple[list[Vector], list[tuple[int, ...]]]] = {}
    for global_id in sorted(PIPE_IDS):
        legacy_name = FORMAL_TO_LEGACY[global_id]
        legacy_mesh = world_mesh_from_blender(bpy.data.objects[legacy_name])
        formal_mesh = world_mesh_from_ifc(model, geometry_settings, global_id)
        legacy_meshes[global_id] = legacy_mesh
        formal_components = component_dimensions(*formal_mesh)
        legacy_components = component_dimensions(*legacy_mesh)
        pipe_comparisons.append(
            {
                "global_id": global_id,
                "legacy_object": legacy_name,
                "formal_component_count": len(formal_components),
                "legacy_component_count": len(legacy_components),
                "formal_component_dimensions_mm": formal_components,
                "legacy_component_dimensions_mm": legacy_components,
                "maximum_component_dimension_difference_mm": best_component_difference_mm(
                    formal_components, legacy_components
                ),
            }
        )

    for global_id, legacy_name in FORMAL_TO_LEGACY.items():
        if global_id in legacy_meshes:
            continue
        legacy_meshes[global_id] = world_mesh_from_blender(bpy.data.objects[legacy_name])

    current_review_pairs = [
        relation
        for relation in coordination["pairs"]
        if tuple(relation["pair"]) in LEGACY_INTERSECTION_PAIRS
    ]
    legacy_pair_comparisons = []
    for relation in current_review_pairs:
        first, second = relation["pair"]
        first_tree = BVHTree.FromPolygons(*legacy_meshes[first], all_triangles=False)
        second_tree = BVHTree.FromPolygons(*legacy_meshes[second], all_triangles=False)
        legacy_intersects = bool(first_tree.overlap(second_tree))
        legacy_pair_comparisons.append(
            {
                "pair": [first, second],
                "formal_geometry_state": relation["geometry_state"],
                "legacy_objects": [FORMAL_TO_LEGACY[first], FORMAL_TO_LEGACY[second]],
                "legacy_geometry_state": "intersecting" if legacy_intersects else "not_intersecting",
                "same_intersection_state": relation["geometry_state"] == "intersecting" and legacy_intersects,
            }
        )

    equipment_records = {
        record["global_id"]: record
        for record in formal_report["inventory"]["typed_high_equipment"]
    }
    translations = []
    for global_id, legacy_name in AC_TRANSLATION_ANCHORS.items():
        formal_centre = equipment_records[global_id]["bbox"]["centre_mm"]
        legacy_centre = object_centre_mm(legacy_name)
        translations.append([formal_centre[axis] - legacy_centre[axis] for axis in range(3)])
    translation = [statistics.median(item[axis] for item in translations) for axis in range(3)]

    legacy_candidate_centre = object_centre_mm(LEGACY_EAST_CANDIDATE)
    predicted_candidate_centre = [
        legacy_candidate_centre[axis] + translation[axis] for axis in range(3)
    ]
    current_candidate_centre = list(args.current_candidate_centre_mm)
    candidate_difference = [
        current_candidate_centre[axis] - predicted_candidate_centre[axis] for axis in range(3)
    ]

    pipe_topology_pass = all(
        item["formal_component_count"] == item["legacy_component_count"]
        and item["maximum_component_dimension_difference_mm"] <= args.tolerance_mm
        for item in pipe_comparisons
    )
    pair_intersections_pass = (
        len(legacy_pair_comparisons) == 10
        and all(item["same_intersection_state"] for item in legacy_pair_comparisons)
    )
    candidate_position_pass = max(abs(value) for value in candidate_difference) <= args.candidate_tolerance_mm

    report: dict[str, Any] = {
        "mode": "read_only_rcp1_legacy_design_base_audit",
        "source": {
            "formal_ifc": str(args.formal_ifc.resolve()),
            "formal_ifc_sha256": formal_sha,
            "legacy_blend": str(Path(bpy.data.filepath).resolve()),
            "legacy_blend_sha256": sha256(Path(bpy.data.filepath)),
        },
        "scope": {
            "legacy_geometry_role": "post-demolition approximate design base only",
            "formal_ifc_write_allowed": False,
            "legacy_blend_write_allowed": False,
            "final_remodel_hvac_design_inferred": False,
        },
        "pipe_comparisons": pipe_comparisons,
        "legacy_pair_comparisons": legacy_pair_comparisons,
        "legacy_east_ac_candidate": {
            "legacy_object": LEGACY_EAST_CANDIDATE,
            "legacy_centre_mm": legacy_candidate_centre,
            "legacy_to_formal_translation_mm": translation,
            "predicted_formal_centre_mm": predicted_candidate_centre,
            "current_blender_candidate_centre_mm": current_candidate_centre,
            "difference_mm": candidate_difference,
            "maximum_axis_difference_mm": max(abs(value) for value in candidate_difference),
            "candidate_has_formal_ifc_global_id": False,
        },
        "gates": {
            "required_legacy_objects_present": not missing_names,
            "missing_legacy_objects": missing_names,
            "five_pipe_topology_and_component_bounds_within_tolerance": pipe_topology_pass,
            "ten_formal_intersections_also_exist_in_legacy_base": pair_intersections_pass,
            "east_ac_candidate_matches_legacy_base_within_tolerance": candidate_position_pass,
            "formal_ifc_write_allowed": False,
            "legacy_base_consistency_pass": not missing_names
            and pipe_topology_pass
            and pair_intersections_pass
            and candidate_position_pass,
        },
        "remaining_design_work": [
            "redesign all remodel refrigerant and condensate routes",
            "select the final AC equipment count and positions",
            "add system and port connectivity",
            "verify diameter, insulation, condensate slope, access and supports",
            "retain verified developer AC openings and do not add new openings",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
