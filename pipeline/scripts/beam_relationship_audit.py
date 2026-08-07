#!/usr/bin/env python3
"""Compare structural face relationships before and after a beam candidate."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import ifcopenshell

from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm


AXES = ("x", "y", "z")
STRUCTURAL_CLASSES = ("IfcBeam", "IfcWall", "IfcSlab")


def read_review_targets(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {"global_id", "review_required"}
    if not rows or not required.issubset(rows[0]):
        raise RuntimeError(f"invalid beam target table: {path}")
    result = [
        row["global_id"] for row in rows if row["review_required"].strip() == "yes"
    ]
    if not result:
        raise RuntimeError("beam target table has no review-required rows")
    return result


def aabb(vertices: Iterable[Sequence[float]]) -> tuple[float, ...]:
    points = list(vertices)
    if not points:
        raise ValueError("cannot derive AABB from an empty vertex set")
    return tuple(
        [min(point[index] for point in points) for index in range(3)]
        + [max(point[index] for point in points) for index in range(3)]
    )


def interval_overlap(first: Sequence[float], second: Sequence[float], axis: int) -> float:
    return max(
        0.0,
        min(first[axis + 3], second[axis + 3])
        - max(first[axis], second[axis]),
    )


def aabb_distance(first: Sequence[float], second: Sequence[float]) -> float:
    deltas = [
        max(
            second[index] - first[index + 3],
            first[index] - second[index + 3],
            0.0,
        )
        for index in range(3)
    ]
    return math.sqrt(sum(delta * delta for delta in deltas))


def face_metrics(
    target: Sequence[float],
    neighbor: Sequence[float],
    axis: int,
    minimum_cross_overlap_mm: float,
) -> dict[str, Any] | None:
    cross_axes = [index for index in range(3) if index != axis]
    overlaps = [interval_overlap(target, neighbor, index) for index in cross_axes]
    if min(overlaps) < minimum_cross_overlap_mm:
        return None
    opposing = {
        "target_min_to_neighbor_max": abs(target[axis] - neighbor[axis + 3]),
        "target_max_to_neighbor_min": abs(target[axis + 3] - neighbor[axis]),
    }
    aligned = {
        "target_min_to_neighbor_min": abs(target[axis] - neighbor[axis]),
        "target_max_to_neighbor_max": abs(
            target[axis + 3] - neighbor[axis + 3]
        ),
    }
    opposing_label, opposing_residual = min(opposing.items(), key=lambda item: item[1])
    aligned_label, aligned_residual = min(aligned.items(), key=lambda item: item[1])
    return {
        "axis": AXES[axis],
        "cross_overlap_mm": {
            AXES[cross_axes[index]]: overlaps[index] for index in range(2)
        },
        "opposing_face": opposing_label,
        "opposing_residual_mm": opposing_residual,
        "aligned_face": aligned_label,
        "aligned_residual_mm": aligned_residual,
    }


def transition(source_mm: float, candidate_mm: float, tolerance_mm: float) -> str:
    source_within = source_mm <= tolerance_mm
    candidate_within = candidate_mm <= tolerance_mm
    if source_within and not candidate_within:
        return "regression"
    if not source_within and candidate_within:
        return "improvement"
    if source_within and candidate_within:
        return "retained_within_tolerance"
    if abs(candidate_mm - source_mm) <= tolerance_mm:
        return "stable_outside_tolerance"
    return "changed_outside_tolerance"


def products_by_gid(model: ifcopenshell.file) -> dict[str, Any]:
    products: dict[str, Any] = {}
    for ifc_class in STRUCTURAL_CLASSES:
        for product in model.by_type(ifc_class):
            global_id = getattr(product, "GlobalId", None)
            if global_id and getattr(product, "Representation", None):
                products[global_id] = product
    return products


def product_aabbs(model: ifcopenshell.file) -> tuple[dict[str, Any], dict[str, tuple[float, ...]]]:
    products = products_by_gid(model)
    settings = geometry_settings()
    bounds: dict[str, tuple[float, ...]] = {}
    for global_id, product in products.items():
        try:
            vertices, _ = world_mesh_mm(settings, product)
        except Exception:  # noqa: BLE001 - individual invalid shapes remain outside the relation set
            continue
        if vertices:
            bounds[global_id] = aabb(vertices)
    return products, bounds


def compare_target(
    global_id: str,
    source_products: dict[str, Any],
    source_bounds: dict[str, tuple[float, ...]],
    candidate_products: dict[str, Any],
    candidate_bounds: dict[str, tuple[float, ...]],
    tolerance_mm: float,
    search_window_mm: float,
    minimum_cross_overlap_mm: float,
) -> dict[str, Any]:
    if global_id not in source_bounds or global_id not in candidate_bounds:
        raise RuntimeError(f"target beam is missing usable geometry: {global_id}")
    source_target = source_bounds[global_id]
    candidate_target = candidate_bounds[global_id]
    neighbor_ids = sorted(
        (set(source_bounds) | set(candidate_bounds)) - {global_id}
    )
    relations = []
    for neighbor_id in neighbor_ids:
        source_neighbor = source_bounds.get(neighbor_id)
        candidate_neighbor = candidate_bounds.get(neighbor_id)
        if source_neighbor is None or candidate_neighbor is None:
            continue
        if (
            aabb_distance(source_target, source_neighbor) > search_window_mm
            and aabb_distance(candidate_target, candidate_neighbor) > search_window_mm
        ):
            continue
        metrics = []
        for axis in range(3):
            source_metric = face_metrics(
                source_target,
                source_neighbor,
                axis,
                minimum_cross_overlap_mm,
            )
            candidate_metric = face_metrics(
                candidate_target,
                candidate_neighbor,
                axis,
                minimum_cross_overlap_mm,
            )
            if source_metric is None or candidate_metric is None:
                continue
            metrics.append(
                {
                    "axis": AXES[axis],
                    "source": source_metric,
                    "candidate": candidate_metric,
                    "opposing_transition": transition(
                        source_metric["opposing_residual_mm"],
                        candidate_metric["opposing_residual_mm"],
                        tolerance_mm,
                    ),
                    "aligned_transition": transition(
                        source_metric["aligned_residual_mm"],
                        candidate_metric["aligned_residual_mm"],
                        tolerance_mm,
                    ),
                }
            )
        if not metrics:
            continue
        neighbor = candidate_products.get(neighbor_id) or source_products[neighbor_id]
        relations.append(
            {
                "neighbor_global_id": neighbor_id,
                "neighbor_class": neighbor.is_a(),
                "neighbor_name": getattr(neighbor, "Name", None),
                "source_aabb_distance_mm": aabb_distance(source_target, source_neighbor),
                "candidate_aabb_distance_mm": aabb_distance(
                    candidate_target, candidate_neighbor
                ),
                "metrics": metrics,
            }
        )
    regressions = [
        {
            "neighbor_global_id": relation["neighbor_global_id"],
            "neighbor_class": relation["neighbor_class"],
            "axis": metric["axis"],
            "metric": metric_name,
            "source_mm": metric["source"][f"{metric_name}_residual_mm"],
            "candidate_mm": metric["candidate"][f"{metric_name}_residual_mm"],
        }
        for relation in relations
        for metric in relation["metrics"]
        for metric_name, transition_name in (
            ("opposing", metric["opposing_transition"]),
            ("aligned", metric["aligned_transition"]),
        )
        if transition_name == "regression"
    ]
    improvements = [
        {
            "neighbor_global_id": relation["neighbor_global_id"],
            "neighbor_class": relation["neighbor_class"],
            "axis": metric["axis"],
            "metric": metric_name,
            "source_mm": metric["source"][f"{metric_name}_residual_mm"],
            "candidate_mm": metric["candidate"][f"{metric_name}_residual_mm"],
        }
        for relation in relations
        for metric in relation["metrics"]
        for metric_name, transition_name in (
            ("opposing", metric["opposing_transition"]),
            ("aligned", metric["aligned_transition"]),
        )
        if transition_name == "improvement"
    ]
    return {
        "global_id": global_id,
        "source_aabb_mm": source_target,
        "candidate_aabb_mm": candidate_target,
        "relations": relations,
        "regressions": regressions,
        "improvements": improvements,
        "requires_linked_resolution": bool(regressions),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--search-window-mm", type=float, default=2.0)
    parser.add_argument("--minimum-cross-overlap-mm", type=float, default=10.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= args.search_window_mm:
        raise SystemExit("tolerance must be positive and no larger than search window")
    if args.minimum_cross_overlap_mm <= 0.0:
        raise SystemExit("minimum cross overlap must be positive")
    source_path = args.source.resolve()
    candidate_path = args.candidate.resolve()
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(candidate_path)
    target_ids = read_review_targets(args.targets)
    source_products, source_bounds = product_aabbs(source)
    candidate_products, candidate_bounds = product_aabbs(candidate)
    targets = [
        compare_target(
            global_id,
            source_products,
            source_bounds,
            candidate_products,
            candidate_bounds,
            args.tolerance_mm,
            args.search_window_mm,
            args.minimum_cross_overlap_mm,
        )
        for global_id in target_ids
    ]
    regression_count = sum(len(target["regressions"]) for target in targets)
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-beam-relationship-audit",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "candidate": {
            "path": str(candidate_path),
            "sha256": sha256(candidate_path),
            "schema": candidate.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "search_window_mm": args.search_window_mm,
        "minimum_cross_overlap_mm": args.minimum_cross_overlap_mm,
        "structural_classes": STRUCTURAL_CLASSES,
        "targets": targets,
        "gates": {
            "target_count": len(targets),
            "regression_count": regression_count,
            "linked_resolution_target_count": sum(
                bool(target["requires_linked_resolution"]) for target in targets
            ),
            "schema_equal": source.schema == candidate.schema,
            "root_global_ids_equal": {
                root.GlobalId for root in source.by_type("IfcRoot")
            }
            == {root.GlobalId for root in candidate.by_type("IfcRoot")},
            "pass_without_linked_resolution": regression_count == 0,
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **report["gates"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
