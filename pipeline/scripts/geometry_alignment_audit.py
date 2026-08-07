#!/usr/bin/env python3
"""Read-only wall alignment and IFC geometry-difference audit.

This script deliberately keeps three different concepts separate:

1. ObjectPlacement: the IFC placement transform / object origin.
2. World geometry: the rendered vertices after every placement is applied.
3. Wall alignment: near-coplanar wall footprint edges and wall junction gaps.

It never writes IFC. Reports are derived artifacts intended for ``build/``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.placement


DEFAULT_TOLERANCE_MM = 0.1
DEFAULT_SEARCH_WINDOW_MM = 1.0
DEFAULT_MIN_SEGMENT_LENGTH_MM = 100.0
DEFAULT_MIN_OVERLAP_MM = 100.0
DEFAULT_MIN_VERTICAL_OVERLAP_MM = 100.0
DEFAULT_ANGLE_TOLERANCE_DEG = 0.01
PLANE_Z_TOLERANCE_MM = 0.01
POINT_KEY_PRECISION_MM = 1e-6

Point2D = tuple[float, float]
Point3D = tuple[float, float, float]


@dataclass(frozen=True)
class Segment:
    wall_guid: str
    wall_name: str | None
    start: Point2D
    end: Point2D
    length_mm: float
    z_min_mm: float = -math.inf
    z_max_mm: float = math.inf


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vector_length(values: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def subtract_2d(a: Point2D, b: Point2D) -> Point2D:
    return a[0] - b[0], a[1] - b[1]


def dot_2d(a: Point2D, b: Point2D) -> float:
    return a[0] * b[0] + a[1] * b[1]


def cross_2d(a: Point2D, b: Point2D) -> float:
    return a[0] * b[1] - a[1] * b[0]


def distance_2d(a: Point2D, b: Point2D) -> float:
    return vector_length(subtract_2d(a, b))


def canonical_point(point: Sequence[float]) -> tuple[float, ...]:
    return tuple(round(float(value) / POINT_KEY_PRECISION_MM) * POINT_KEY_PRECISION_MM for value in point)


def canonical_edge(a: Sequence[float], b: Sequence[float]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    first = canonical_point(a)
    second = canonical_point(b)
    return (first, second) if first <= second else (second, first)


def geometry_settings() -> ifcopenshell.geom.settings:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    return settings


def world_mesh_mm(
    settings: ifcopenshell.geom.settings, product: ifcopenshell.entity_instance
) -> tuple[list[Point3D], list[tuple[int, int, int]]]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    raw_vertices = list(shape.geometry.verts)
    raw_faces = list(shape.geometry.faces)
    vertices = [
        tuple(float(raw_vertices[index + offset]) * 1000.0 for offset in range(3))
        for index in range(0, len(raw_vertices), 3)
    ]
    faces = [
        (int(raw_faces[index]), int(raw_faces[index + 1]), int(raw_faces[index + 2]))
        for index in range(0, len(raw_faces), 3)
    ]
    return vertices, faces


def horizontal_outline_segments(
    settings: ifcopenshell.geom.settings,
    wall: ifcopenshell.entity_instance,
    min_segment_length_mm: float,
) -> list[Segment]:
    vertices, faces = world_mesh_mm(settings, wall)
    wall_z_min = min((vertex[2] for vertex in vertices), default=0.0)
    wall_z_max = max((vertex[2] for vertex in vertices), default=0.0)
    horizontal_faces_by_z: dict[float, list[tuple[int, int, int]]] = defaultdict(list)
    for face in faces:
        z_values = [vertices[index][2] for index in face]
        if max(z_values) - min(z_values) <= PLANE_Z_TOLERANCE_MM:
            z_key = round(sum(z_values) / 3.0 / PLANE_Z_TOLERANCE_MM) * PLANE_Z_TOLERANCE_MM
            horizontal_faces_by_z[z_key].append(face)

    projected: dict[tuple[tuple[float, ...], tuple[float, ...]], Segment] = {}
    for plane_faces in horizontal_faces_by_z.values():
        edge_counts: Counter[tuple[tuple[float, ...], tuple[float, ...]]] = Counter()
        edge_points: dict[tuple[tuple[float, ...], tuple[float, ...]], tuple[Point3D, Point3D]] = {}
        for face in plane_faces:
            for first_index, second_index in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
                first = vertices[first_index]
                second = vertices[second_index]
                key = canonical_edge(first, second)
                edge_counts[key] += 1
                edge_points[key] = (first, second)
        for key, count in edge_counts.items():
            if count != 1:
                continue
            first, second = edge_points[key]
            start = (first[0], first[1])
            end = (second[0], second[1])
            length = distance_2d(start, end)
            if length < min_segment_length_mm:
                continue
            projected_key = canonical_edge(start, end)
            projected[projected_key] = Segment(
                wall_guid=wall.GlobalId,
                wall_name=wall.Name,
                start=start,
                end=end,
                length_mm=length,
                z_min_mm=wall_z_min,
                z_max_mm=wall_z_max,
            )
    return list(projected.values())


def segment_unit(segment: Segment) -> Point2D:
    delta = subtract_2d(segment.end, segment.start)
    length = vector_length(delta)
    return delta[0] / length, delta[1] / length


def segments_are_parallel(
    first_unit: Point2D,
    second_unit: Point2D,
    sine_tolerance: float,
) -> bool:
    return abs(cross_2d(first_unit, second_unit)) <= sine_tolerance


def segment_overlap_mm(first: Segment, second: Segment, unit: Point2D) -> float:
    origin = first.start
    first_interval = sorted((0.0, dot_2d(subtract_2d(first.end, origin), unit)))
    second_interval = sorted(
        (
            dot_2d(subtract_2d(second.start, origin), unit),
            dot_2d(subtract_2d(second.end, origin), unit),
        )
    )
    return max(0.0, min(first_interval[1], second_interval[1]) - max(first_interval[0], second_interval[0]))


def parallel_line_gap_mm(first: Segment, second: Segment, unit: Point2D) -> float:
    offsets = (
        abs(cross_2d(subtract_2d(second.start, first.start), unit)),
        abs(cross_2d(subtract_2d(second.end, first.start), unit)),
    )
    return max(offsets)


def point_to_segment_distance_mm(point: Point2D, segment: Segment) -> tuple[float, float]:
    delta = subtract_2d(segment.end, segment.start)
    length_squared = dot_2d(delta, delta)
    if length_squared == 0.0:
        return distance_2d(point, segment.start), 0.0
    position = dot_2d(subtract_2d(point, segment.start), delta) / length_squared
    clamped = min(1.0, max(0.0, position))
    closest = (segment.start[0] + delta[0] * clamped, segment.start[1] + delta[1] * clamped)
    return distance_2d(point, closest), position


def vertical_overlap_mm(first: Segment, second: Segment) -> float:
    return max(0.0, min(first.z_max_mm, second.z_max_mm) - max(first.z_min_mm, second.z_min_mm))


def record_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        record["wall_a"],
        record["wall_b"],
        round(record["gap_mm"], 6),
        tuple(round(value, 3) for value in record.get("point_mm", [])),
    )


def build_review_clusters(
    coplanar_records: Sequence[dict[str, Any]],
    junction_records: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Group unresolved wall-pair candidates into deterministic review clusters."""
    typed_records = [
        (record_type, record)
        for record_type, records in (
            ("coplanar_edge", coplanar_records),
            ("junction", junction_records),
        )
        for record in records
        if not record["within_tolerance"]
    ]
    adjacency: dict[str, set[str]] = defaultdict(set)
    for _, record in typed_records:
        wall_a, wall_b = record["wall_a"], record["wall_b"]
        adjacency[wall_a].add(wall_b)
        adjacency[wall_b].add(wall_a)

    components: list[list[str]] = []
    visited: set[str] = set()
    for start in sorted(adjacency):
        if start in visited:
            continue
        stack = [start]
        visited.add(start)
        component: list[str] = []
        while stack:
            wall = stack.pop()
            component.append(wall)
            for neighbour in sorted(adjacency[wall], reverse=True):
                if neighbour not in visited:
                    visited.add(neighbour)
                    stack.append(neighbour)
        components.append(sorted(component))
    components.sort(key=lambda component: (-len(component), component))

    clusters: list[dict[str, Any]] = []
    for index, component in enumerate(components, start=1):
        component_set = set(component)
        records = [
            (record_type, record)
            for record_type, record in typed_records
            if record["wall_a"] in component_set and record["wall_b"] in component_set
        ]
        pairs = sorted({tuple(sorted((record["wall_a"], record["wall_b"]))) for _, record in records})
        gaps = [float(record["gap_mm"]) for _, record in records]
        clusters.append(
            {
                "review_id": f"C003-G{index}",
                "stable_id": f"C003-S{hashlib.sha256(component[0].encode('utf-8')).hexdigest()[:8].upper()}",
                "display_order": index,
                "status": "requires_design_intent_confirmation",
                "wall_count": len(component),
                "wall_global_ids": component,
                "pair_count": len(pairs),
                "pairs": [list(pair) for pair in pairs],
                "coplanar_record_count": sum(
                    record_type == "coplanar_edge" for record_type, _ in records
                ),
                "junction_record_count": sum(
                    record_type == "junction" for record_type, _ in records
                ),
                "minimum_gap_mm": min(gaps),
                "maximum_gap_mm": max(gaps),
            }
        )
    return clusters


def retain_minimum_gap_record(
    records: dict[Any, dict[str, Any]],
    key: Any,
    record: dict[str, Any],
) -> None:
    previous = records.get(key)
    if previous is None or record["gap_mm"] < previous["gap_mm"]:
        records[key] = record


def alignment_audit(
    ifc: ifcopenshell.file,
    tolerance_mm: float,
    search_window_mm: float,
    min_segment_length_mm: float,
    min_overlap_mm: float,
    min_vertical_overlap_mm: float,
    angle_tolerance_deg: float,
) -> dict[str, Any]:
    settings = geometry_settings()
    segments: list[Segment] = []
    shape_failures: list[dict[str, str]] = []
    walls_without_horizontal_outlines: list[str] = []
    walls = list(ifc.by_type("IfcWall"))
    for wall in walls:
        try:
            wall_segments = horizontal_outline_segments(settings, wall, min_segment_length_mm)
        except Exception as error:  # pragma: no cover - model-specific geometry failures
            shape_failures.append({"global_id": wall.GlobalId, "error": str(error)})
            continue
        if not wall_segments:
            walls_without_horizontal_outlines.append(wall.GlobalId)
        segments.extend(wall_segments)

    sine_tolerance = math.sin(math.radians(angle_tolerance_deg))
    coplanar_records: list[dict[str, Any]] = []
    # A non-parallel wall pair has one construction connection even though its
    # two thickness faces can produce several endpoint distances. Keep the
    # minimum pair distance so a harmless face overhang is not reported as a
    # gap after another face already proves contact. Parallel relationships
    # remain in the separate coplanar-edge audit.
    junction_records_by_key: dict[tuple[str, str], dict[str, Any]] = {}

    for first_index, first in enumerate(segments):
        first_unit = segment_unit(first)
        for second in segments[first_index + 1 :]:
            if first.wall_guid == second.wall_guid:
                continue
            if vertical_overlap_mm(first, second) < min_vertical_overlap_mm:
                continue
            wall_a, wall_b = sorted((first.wall_guid, second.wall_guid))
            second_unit = segment_unit(second)
            parallel = segments_are_parallel(
                first_unit, second_unit, sine_tolerance
            )
            if parallel:
                overlap = segment_overlap_mm(first, second, first_unit)
                if overlap >= min_overlap_mm:
                    gap = parallel_line_gap_mm(first, second, first_unit)
                    if gap <= search_window_mm:
                        coplanar_records.append(
                            {
                                "wall_a": wall_a,
                                "wall_b": wall_b,
                                "gap_mm": gap,
                                "overlap_mm": overlap,
                                "within_tolerance": gap <= tolerance_mm,
                                "segment_a_mm": [first.start, first.end],
                                "segment_b_mm": [second.start, second.end],
                            }
                        )

            if parallel:
                continue

            for endpoint, source, target in (
                (first.start, first, second),
                (first.end, first, second),
                (second.start, second, first),
                (second.end, second, first),
            ):
                gap, position = point_to_segment_distance_mm(endpoint, target)
                if gap > search_window_mm or position < -1e-9 or position > 1.0 + 1e-9:
                    continue
                key = (wall_a, wall_b)
                record = {
                    "wall_a": wall_a,
                    "wall_b": wall_b,
                    "gap_mm": gap,
                    "point_mm": endpoint,
                    "source_wall": source.wall_guid,
                    "target_wall": target.wall_guid,
                    "within_tolerance": gap <= tolerance_mm,
                }
                retain_minimum_gap_record(junction_records_by_key, key, record)

    # Collapse duplicate top/bottom or segmented records.
    coplanar_unique: dict[tuple[Any, ...], dict[str, Any]] = {}
    for record in sorted(coplanar_records, key=lambda item: item["gap_mm"]):
        direction = subtract_2d(record["segment_a_mm"][1], record["segment_a_mm"][0])
        angle = math.atan2(direction[1], direction[0]) % math.pi
        pair_key = (record["wall_a"], record["wall_b"], round(angle, 5), round(record["gap_mm"], 4))
        coplanar_unique.setdefault(pair_key, record)
    coplanar = sorted(coplanar_unique.values(), key=record_key)
    junctions = sorted(junction_records_by_key.values(), key=record_key)
    review_clusters = build_review_clusters(coplanar, junctions)

    def summary(records: Sequence[dict[str, Any]]) -> dict[str, Any]:
        return {
            "total": len(records),
            "within_tolerance": sum(record["within_tolerance"] for record in records),
            "over_tolerance": sum(not record["within_tolerance"] for record in records),
            "max_gap_mm": max((record["gap_mm"] for record in records), default=0.0),
            "records": records,
        }

    return {
        "interpretation": (
            "Near-coplanar edges and junctions are geometric candidates, not confirmed construction "
            "relationships. They become gate conditions only after the relationship is explicitly confirmed."
        ),
        "tolerance_mm": tolerance_mm,
        "search_window_mm": search_window_mm,
        "angle_tolerance_deg": angle_tolerance_deg,
        "minimum_segment_length_mm": min_segment_length_mm,
        "minimum_overlap_mm": min_overlap_mm,
        "minimum_vertical_overlap_mm": min_vertical_overlap_mm,
        "walls_total": len(walls),
        "wall_segments_total": len(segments),
        "walls_without_horizontal_outlines": walls_without_horizontal_outlines,
        "shape_failures": shape_failures,
        "coplanar_edges": summary(coplanar),
        "junctions": summary(junctions),
        "review_clusters": review_clusters,
    }


def load_git_ifc(source: Path, git_ref: str) -> tuple[ifcopenshell.file, str]:
    repository = subprocess.check_output(
        ["git", "-C", str(source.parent), "rev-parse", "--show-toplevel"], text=True
    ).strip()
    repository_path = Path(repository)
    relative_path = source.resolve().relative_to(repository_path.resolve())
    blob = subprocess.check_output(
        ["git", "-C", repository, "show", f"{git_ref}:{relative_path.as_posix()}"]
    )
    return ifcopenshell.file.from_string(blob.decode("utf-8")), f"git:{git_ref}:{relative_path.as_posix()}"


def placement_matrix(product: ifcopenshell.entity_instance) -> list[list[float]] | None:
    placement = getattr(product, "ObjectPlacement", None)
    if not placement or not placement.is_a("IfcLocalPlacement"):
        return None
    matrix = ifcopenshell.util.placement.get_local_placement(placement)
    # The project uses millimetres, and get_local_placement returns project-unit
    # translation values. Rotation remains dimensionless.
    return [[float(matrix[row][column]) for column in range(4)] for row in range(4)]


def placement_deltas(
    first: list[list[float]] | None, second: list[list[float]] | None
) -> tuple[float | None, float | None]:
    if first is None or second is None:
        return None, None
    translation_delta = vector_length(
        tuple(first[row][3] - second[row][3] for row in range(3))
    )
    rotation_delta = max(
        abs(first[row][column] - second[row][column])
        for row in range(3)
        for column in range(3)
    )
    return translation_delta, rotation_delta


def unique_vertices(vertices: Iterable[Point3D]) -> list[Point3D]:
    return sorted({canonical_point(vertex) for vertex in vertices})  # type: ignore[return-value]


def directed_vertex_distance_mm(first: Sequence[Point3D], second: Sequence[Point3D]) -> float:
    if not first or not second:
        return math.inf if first != second else 0.0
    maximum = 0.0
    for point in first:
        nearest = min(vector_length(tuple(point[index] - other[index] for index in range(3))) for other in second)
        maximum = max(maximum, nearest)
    return maximum


def vertex_hausdorff_mm(first: Sequence[Point3D], second: Sequence[Point3D]) -> float:
    return max(directed_vertex_distance_mm(first, second), directed_vertex_distance_mm(second, first))


def products_for_comparison(
    ifc: ifcopenshell.file, classes: Sequence[str], global_ids: Sequence[str]
) -> dict[str, ifcopenshell.entity_instance]:
    if global_ids:
        products: dict[str, ifcopenshell.entity_instance] = {}
        for global_id in global_ids:
            try:
                product = ifc.by_guid(global_id)
            except RuntimeError:
                continue
            if product and product.is_a("IfcProduct"):
                products[global_id] = product
        return products
    products = {}
    for ifc_class in classes:
        for product in ifc.by_type(ifc_class):
            global_id = getattr(product, "GlobalId", None)
            if global_id:
                products[global_id] = product
    return products


def geometry_difference_audit(
    current: ifcopenshell.file,
    baseline: ifcopenshell.file,
    baseline_source: str,
    tolerance_mm: float,
    classes: Sequence[str],
    global_ids: Sequence[str],
) -> dict[str, Any]:
    current_products = products_for_comparison(current, classes, global_ids)
    baseline_products = products_for_comparison(baseline, classes, global_ids)
    settings = geometry_settings()
    records: list[dict[str, Any]] = []
    for global_id in sorted(set(current_products) | set(baseline_products)):
        current_product = current_products.get(global_id)
        baseline_product = baseline_products.get(global_id)
        if current_product is None or baseline_product is None:
            records.append(
                {
                    "global_id": global_id,
                    "status": "added" if baseline_product is None else "removed",
                    "within_tolerance": False,
                }
            )
            continue
        try:
            current_vertices, current_faces = world_mesh_mm(settings, current_product)
            baseline_vertices, baseline_faces = world_mesh_mm(settings, baseline_product)
            current_unique = unique_vertices(current_vertices)
            baseline_unique = unique_vertices(baseline_vertices)
            geometry_delta = vertex_hausdorff_mm(current_unique, baseline_unique)
            shape_error = None
        except Exception as error:  # pragma: no cover - model-specific geometry failures
            current_unique = []
            baseline_unique = []
            current_faces = []
            baseline_faces = []
            geometry_delta = math.inf
            shape_error = str(error)
        current_matrix = placement_matrix(current_product)
        baseline_matrix = placement_matrix(baseline_product)
        placement_translation_delta, placement_rotation_delta = placement_deltas(
            current_matrix, baseline_matrix
        )
        placement_changed = (
            placement_translation_delta is not None
            and placement_translation_delta > tolerance_mm
        ) or (
            placement_rotation_delta is not None
            and placement_rotation_delta > 1e-9
        )
        geometry_changed = geometry_delta > tolerance_mm
        if not placement_changed and not geometry_changed:
            status = "unchanged"
        elif placement_changed and not geometry_changed:
            status = "placement_changed_world_geometry_preserved"
        else:
            status = "world_geometry_changed"
        records.append(
            {
                "global_id": global_id,
                "ifc_class": current_product.is_a(),
                "name": current_product.Name,
                "status": status,
                "within_tolerance": not geometry_changed,
                "placement_translation_delta_mm": placement_translation_delta,
                "placement_rotation_matrix_max_delta": placement_rotation_delta,
                "world_vertex_hausdorff_mm": geometry_delta,
                "current_unique_vertices": len(current_unique),
                "baseline_unique_vertices": len(baseline_unique),
                "current_triangles": len(current_faces),
                "baseline_triangles": len(baseline_faces),
                "topology_counts_equal": len(current_unique) == len(baseline_unique)
                and len(current_faces) == len(baseline_faces),
                "shape_error": shape_error,
            }
        )
    status_counts = Counter(record["status"] for record in records)
    return {
        "baseline": baseline_source,
        "tolerance_mm": tolerance_mm,
        "classes": list(classes),
        "requested_global_ids": list(global_ids),
        "total": len(records),
        "within_tolerance": sum(record["within_tolerance"] for record in records),
        "over_tolerance": sum(not record["within_tolerance"] for record in records),
        "status_counts": dict(status_counts),
        "records": records,
    }


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise argparse.ArgumentTypeError("must be a finite number greater than zero")
    return number


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument(
        "--tolerance-mm",
        type=positive_float,
        default=DEFAULT_TOLERANCE_MM,
        help="Pass/fail tolerance for wall alignment and geometry differences (default: 0.1 mm).",
    )
    parser.add_argument(
        "--search-window-mm",
        type=positive_float,
        default=DEFAULT_SEARCH_WINDOW_MM,
        help="Maximum gap considered a candidate wall relationship (default: 1.0 mm).",
    )
    parser.add_argument(
        "--minimum-segment-length-mm",
        type=positive_float,
        default=DEFAULT_MIN_SEGMENT_LENGTH_MM,
    )
    parser.add_argument(
        "--minimum-overlap-mm",
        type=positive_float,
        default=DEFAULT_MIN_OVERLAP_MM,
    )
    parser.add_argument(
        "--minimum-vertical-overlap-mm",
        type=positive_float,
        default=DEFAULT_MIN_VERTICAL_OVERLAP_MM,
        help="Minimum shared wall height required before plan gaps are compared (default: 100 mm).",
    )
    parser.add_argument(
        "--angle-tolerance-deg",
        type=positive_float,
        default=DEFAULT_ANGLE_TOLERANCE_DEG,
    )
    baseline_group = parser.add_mutually_exclusive_group()
    baseline_group.add_argument("--baseline-ifc", type=Path)
    baseline_group.add_argument("--baseline-git-ref")
    parser.add_argument(
        "--compare-class",
        action="append",
        default=[],
        help="IFC class to compare against the baseline; repeatable (default: IfcWall).",
    )
    parser.add_argument(
        "--global-id",
        action="append",
        default=[],
        help="Limit baseline comparison to one GlobalId; repeatable.",
    )
    parser.add_argument(
        "--fail-on-geometry-difference",
        action="store_true",
        help="Exit 2 only when a requested baseline geometry comparison exceeds tolerance.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = args.input.resolve()
    current = ifcopenshell.open(source)
    if args.search_window_mm < args.tolerance_mm:
        raise SystemExit("--search-window-mm must be greater than or equal to --tolerance-mm")
    report: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-geometry-alignment-audit",
        "source": {"path": str(source), "sha256": sha256(source), "schema": current.schema},
        "alignment": alignment_audit(
            current,
            tolerance_mm=args.tolerance_mm,
            search_window_mm=args.search_window_mm,
            min_segment_length_mm=args.minimum_segment_length_mm,
            min_overlap_mm=args.minimum_overlap_mm,
            min_vertical_overlap_mm=args.minimum_vertical_overlap_mm,
            angle_tolerance_deg=args.angle_tolerance_deg,
        ),
    }
    if args.baseline_ifc:
        baseline_path = args.baseline_ifc.resolve()
        baseline = ifcopenshell.open(baseline_path)
        baseline_source = str(baseline_path)
    elif args.baseline_git_ref:
        baseline, baseline_source = load_git_ifc(source, args.baseline_git_ref)
    else:
        baseline = None
        baseline_source = ""
    if baseline is not None:
        report["geometry_difference"] = geometry_difference_audit(
            current,
            baseline,
            baseline_source,
            tolerance_mm=args.tolerance_mm,
            classes=args.compare_class or ["IfcWall"],
            global_ids=args.global_id,
        )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = {
        "report": str(args.report),
        "tolerance_mm": args.tolerance_mm,
        "coplanar_edges": {
            key: report["alignment"]["coplanar_edges"][key]
            for key in ("total", "within_tolerance", "over_tolerance", "max_gap_mm")
        },
        "junctions": {
            key: report["alignment"]["junctions"][key]
            for key in ("total", "within_tolerance", "over_tolerance", "max_gap_mm")
        },
    }
    if "geometry_difference" in report:
        summary["geometry_difference"] = {
            key: report["geometry_difference"][key]
            for key in ("total", "within_tolerance", "over_tolerance", "status_counts")
        }
    print(json.dumps(summary, ensure_ascii=False))
    geometry_over_tolerance = report.get("geometry_difference", {}).get("over_tolerance", 0)
    if args.fail_on_geometry_difference and geometry_over_tolerance:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
