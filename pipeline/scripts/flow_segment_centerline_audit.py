#!/usr/bin/env python3
"""Audit centreline evidence for remaining IfcFlowSegment-family origins.

This audit is read-only.  It distinguishes disconnected geometry that should
not share one IFC product from multiple representation components that meet at
an existing connector point.  Derived PCA axes are evidence only; they never
authorize placement or geometry writes.
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.util.system
import numpy as np

from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm


PVC110_BUNDLE_DECISION_ID = "COORD-FLOW-C003-PVC110-BUNDLE"


def mesh_components(
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, int, int]],
) -> list[list[int]]:
    """Return triangle-connected vertex indices in stable order."""

    parent = list(range(len(vertices)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(first: int, second: int) -> None:
        first_root, second_root = find(first), find(second)
        if first_root != second_root:
            parent[second_root] = first_root

    for face in faces:
        for first, second in zip(face, (*face[1:], face[0])):
            union(int(first), int(second))

    groups: dict[int, list[int]] = {}
    for index in range(len(vertices)):
        groups.setdefault(find(index), []).append(index)
    return sorted(groups.values(), key=lambda group: min(group))


def minimum_vertex_distance_mm(
    first: np.ndarray, second: np.ndarray, chunk_size: int = 512
) -> float:
    """Return the exact sampled mesh-vertex distance without a SciPy dependency."""

    closest = float("inf")
    for offset in range(0, len(first), chunk_size):
        chunk = first[offset : offset + chunk_size]
        distances = np.linalg.norm(chunk[:, None, :] - second[None, :, :], axis=2)
        closest = min(closest, float(distances.min()))
    return closest


def linked_component_groups(
    component_points: list[np.ndarray], connection_tolerance_mm: float
) -> tuple[list[list[int]], list[dict[str, Any]]]:
    """Group mesh components only when sampled vertices actually meet."""

    parent = list(range(len(component_points)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(first: int, second: int) -> None:
        first_root, second_root = find(first), find(second)
        if first_root != second_root:
            parent[second_root] = first_root

    pairs = []
    for first in range(len(component_points)):
        for second in range(first + 1, len(component_points)):
            distance = minimum_vertex_distance_mm(
                component_points[first], component_points[second]
            )
            linked = distance <= connection_tolerance_mm
            if linked:
                union(first, second)
            pairs.append(
                {
                    "component_indices": [first, second],
                    "minimum_vertex_distance_mm": distance,
                    "linked_within_tolerance": linked,
                }
            )

    groups: dict[int, list[int]] = {}
    for index in range(len(component_points)):
        groups.setdefault(find(index), []).append(index)
    return sorted(groups.values(), key=lambda group: min(group)), pairs


def integer_residual_mm(point: np.ndarray) -> list[float]:
    return [abs(float(value) - round(float(value))) for value in point]


def principal_axis_record(points: np.ndarray) -> dict[str, Any]:
    """Describe a possible straight circular-prism centreline from mesh points."""

    centre = points.mean(axis=0)
    covariance = np.cov((points - centre).T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    axis = eigenvectors[:, order[0]]
    dominant = int(np.argmax(np.abs(axis)))
    if axis[dominant] < 0:
        axis = -axis
    projections = (points - centre) @ axis
    lower_value, upper_value = float(projections.min()), float(projections.max())
    cap_tolerance = max(1e-6, (upper_value - lower_value) * 1e-9)
    lower = points[np.abs(projections - lower_value) <= cap_tolerance].mean(axis=0)
    upper = points[np.abs(projections - upper_value) <= cap_tolerance].mean(axis=0)
    axial_ratio = float(eigenvalues[0] / max(eigenvalues[1], 1e-12))
    cross_ratio = float(eigenvalues[1] / max(eigenvalues[2], 1e-12))
    straight_circular_prism_candidate = axial_ratio >= 4.0 and cross_ratio <= 1.25
    return {
        "vertex_count": len(points),
        "bbox_size_mm": (points.max(axis=0) - points.min(axis=0)).tolist(),
        "eigenvalues": eigenvalues.tolist(),
        "principal_axis": axis.tolist(),
        "axis_extent_mm": upper_value - lower_value,
        "axial_to_cross_variance_ratio": axial_ratio,
        "cross_variance_ratio": cross_ratio,
        "straight_circular_prism_candidate": straight_circular_prism_candidate,
        "derived_end_centres_mm": [lower.tolist(), upper.tolist()],
        "end_centre_integer_residual_mm": [
            integer_residual_mm(lower),
            integer_residual_mm(upper),
        ],
        "derived_centreline_is_write_authority": False,
    }


def body_item_count(product: Any) -> int:
    if not product.Representation:
        return 0
    return sum(
        len(representation.Items)
        for representation in product.Representation.Representations
        if representation.RepresentationIdentifier == "Body"
    )


def body_item_classes(product: Any) -> list[str]:
    if not product.Representation:
        return []
    return [
        item.is_a()
        for representation in product.Representation.Representations
        if representation.RepresentationIdentifier == "Body"
        for item in representation.Items
    ]


def classify_product(
    logical_run_count: int,
    mesh_component_count: int,
    has_assigned_type: bool,
    controlled_disconnected_bundle: bool = False,
) -> tuple[str, str, bool]:
    """Classify topology without treating typed product parts as separate runs."""

    if logical_run_count > 1 and controlled_disconnected_bundle:
        return (
            "controlled_disconnected_bundle",
            "preserve_existing_ifc_product_and_audit_each_geometry_branch_individually",
            False,
        )
    if logical_run_count > 1 and has_assigned_type:
        return (
            "typed_multibody_product",
            "preserve_as_one_typed_product_and_audit_its_installation_or_connector_datum",
            False,
        )
    if logical_run_count > 1:
        return (
            "disconnected_geometry_bundle",
            "split_into_one_IfcFlowSegment_per_disconnected_run_after_human_confirmation",
            True,
        )
    if mesh_component_count > 1:
        return (
            "multi_component_connected_run",
            "preserve_as_one_logical_run_and_audit_connected_centreline",
            False,
        )
    return (
        "single_component_run",
        "audit_derived_centreline_against_connected_fittings_or_installation_datum",
        False,
    )


def audit_product(
    product: Any,
    origin_record: dict[str, Any],
    connection_tolerance_mm: float,
    controlled_bundle_ids: set[str],
) -> dict[str, Any]:
    vertices, faces = world_mesh_mm(geometry_settings(), product)
    points = np.asarray(vertices, dtype=float)
    components = mesh_components(vertices, faces)
    component_points = [points[indices] for indices in components]
    logical_groups, pair_checks = linked_component_groups(
        component_points, connection_tolerance_mm
    )
    ports = ifcopenshell.util.system.get_ports(product)
    assigned_type = ifcopenshell.util.element.get_type(product)
    classification, next_action, human_decision_required = classify_product(
        len(logical_groups),
        len(components),
        assigned_type is not None,
        product.GlobalId in controlled_bundle_ids,
    )

    return {
        "global_id": product.GlobalId,
        "name": product.Name,
        "current_origin_mm": origin_record["current_mm"],
        "current_origin_max_integer_residual_mm": origin_record["max_residual_mm"],
        "body_item_count": body_item_count(product),
        "body_item_classes": body_item_classes(product),
        "assigned_type": (
            {
                "global_id": assigned_type.GlobalId,
                "ifc_class": assigned_type.is_a(),
                "name": assigned_type.Name,
            }
            if assigned_type is not None
            else None
        ),
        "distribution_port_count": len(ports),
        "distribution_port_global_ids": [port.GlobalId for port in ports],
        "mesh_component_count": len(components),
        "logical_run_count": len(logical_groups),
        "read_only_branch_count": len(logical_groups),
        "logical_component_groups": logical_groups,
        "component_pair_checks": pair_checks,
        "components": [principal_axis_record(component) for component in component_points],
        "classification": classification,
        "next_action": next_action,
        "human_decision_required": human_decision_required,
        "automatic_write_allowed": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--origin-report", required=True, type=Path)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--connection-tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def controlled_bundle_ids(path: Path) -> set[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    matching = [
        row
        for row in rows
        if row["decision_id"] == PVC110_BUNDLE_DECISION_ID
        and row["scope"] == "flow-segment-bundle-exception"
        and row["review_required"] == "no"
        and row["status"] == "implemented"
    ]
    if len(matching) != 1:
        return set()
    return {
        value.strip()
        for value in matching[0]["object_guid"].split(";")
        if value.strip()
    }


def main() -> None:
    args = parse_args()
    if args.connection_tolerance_mm <= 0:
        raise SystemExit("--connection-tolerance-mm must be positive")
    source_path = args.input.resolve()
    model = ifcopenshell.open(source_path)
    origin_report = json.loads(args.origin_report.read_text(encoding="utf-8"))
    origin_records = {
        record["global_id"]: record
        for record in origin_report["records"]
        if record["ifc_class"] in {"IfcFlowSegment", "IfcPipeSegment"}
    }
    controlled_ids = controlled_bundle_ids(args.decisions)
    records = [
        audit_product(
            model.by_guid(global_id),
            record,
            args.connection_tolerance_mm,
            controlled_ids,
        )
        for global_id, record in origin_records.items()
    ]
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-flow-segment-centreline-audit",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": model.schema,
        },
        "origin_report": {
            "path": str(args.origin_report.resolve()),
            "source_sha256": origin_report["source"]["sha256"],
            "matches_input": origin_report["source"]["sha256"] == sha256(source_path),
        },
        "decisions": {
            "path": str(args.decisions.resolve()),
            "controlled_bundle_global_ids": sorted(controlled_ids),
        },
        "connection_tolerance_mm": args.connection_tolerance_mm,
        "summary": {
            "remaining_flow_segment_family_products": len(records),
            "without_distribution_ports": sum(
                record["distribution_port_count"] == 0 for record in records
            ),
            "disconnected_geometry_bundles": sum(
                record["classification"] == "disconnected_geometry_bundle"
                for record in records
            ),
            "controlled_disconnected_bundles": sum(
                record["classification"] == "controlled_disconnected_bundle"
                for record in records
            ),
            "controlled_read_only_branches": sum(
                record["read_only_branch_count"]
                for record in records
                if record["classification"] == "controlled_disconnected_bundle"
            ),
            "typed_multibody_products": sum(
                record["classification"] == "typed_multibody_product"
                for record in records
            ),
            "multi_component_connected_runs": sum(
                record["classification"] == "multi_component_connected_run"
                for record in records
            ),
            "single_component_runs": sum(
                record["classification"] == "single_component_run"
                for record in records
            ),
            "human_decisions_required": sum(
                record["human_decision_required"] for record in records
            ),
        },
        "records": records,
        "automatic_write_allowed": False,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report["summary"], ensure_ascii=False))


if __name__ == "__main__":
    main()
