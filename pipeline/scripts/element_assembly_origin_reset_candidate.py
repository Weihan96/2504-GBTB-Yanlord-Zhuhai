#!/usr/bin/env python3
"""Build a no-world-movement origin-reset candidate for typed assemblies."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api
import ifcopenshell.util.placement
import numpy as np

from beam_origin_reset_candidate import point_mesh_distance
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
        raise RuntimeError(f"invalid assembly origin target schema: {path}")
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
            raise RuntimeError(f"assembly target is not integer millimetres: {row}")
        result.append(
            {
                "global_id": row["global_id"],
                "target_mm": target,
                "anchor_kind": row["anchor_kind"],
                "basis": row["basis"],
                "confidence": float(row["confidence"]),
            }
        )
    if len({row["global_id"] for row in result}) != len(result):
        raise RuntimeError("assembly origin targets contain duplicate GlobalIds")
    return result


def decomposed_products(
    assembly: ifcopenshell.entity_instance,
) -> list[tuple[int, ifcopenshell.entity_instance]]:
    result: list[tuple[int, ifcopenshell.entity_instance]] = []
    queue = [
        (1, child)
        for relation in assembly.IsDecomposedBy
        for child in relation.RelatedObjects
    ]
    seen: set[int] = set()
    while queue:
        depth, product = queue.pop(0)
        if product.id() in seen:
            continue
        seen.add(product.id())
        result.append((depth, product))
        queue.extend(
            (depth + 1, child)
            for relation in product.IsDecomposedBy
            for child in relation.RelatedObjects
        )
    return result


def descendant_mesh_mm(
    settings: Any,
    descendants: list[tuple[int, ifcopenshell.entity_instance]],
) -> tuple[list[tuple[float, ...]], list[tuple[int, int, int]]]:
    vertices: list[tuple[float, ...]] = []
    faces: list[tuple[int, int, int]] = []
    for _, product in descendants:
        if not getattr(product, "Representation", None):
            continue
        product_vertices, product_faces = world_mesh_mm(settings, product)
        offset = len(vertices)
        vertices.extend(product_vertices)
        faces.extend(
            tuple(int(index) + offset for index in face)
            for face in product_faces
        )
    if not vertices or not faces:
        raise RuntimeError("assembly has no represented descendant geometry")
    return vertices, faces


def world_matrix(product: ifcopenshell.entity_instance) -> np.ndarray:
    return np.array(
        ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
        dtype=float,
    )


def relation_signature(model: ifcopenshell.file) -> dict[str, list[str]]:
    return {
        relation.GlobalId: sorted(
            [relation.RelatingObject.GlobalId]
            + [product.GlobalId for product in relation.RelatedObjects]
        )
        for relation in model.by_type("IfcRelAggregates", include_subtypes=False)
    }


def maximum_descendant_world_placement_residual(
    source: ifcopenshell.file,
    candidate: ifcopenshell.file,
    assembly_global_id: str,
) -> float:
    source_assembly = source.by_guid(assembly_global_id)
    candidate_assembly = candidate.by_guid(assembly_global_id)
    source_descendants = {
        product.GlobalId: world_matrix(product)
        for _, product in decomposed_products(source_assembly)
        if getattr(product, "ObjectPlacement", None)
    }
    candidate_descendants = {
        product.GlobalId: world_matrix(product)
        for _, product in decomposed_products(candidate_assembly)
        if getattr(product, "ObjectPlacement", None)
    }
    if source_descendants.keys() != candidate_descendants.keys():
        raise RuntimeError(
            f"assembly descendant set changed: {assembly_global_id}"
        )
    return max(
        (
            float(
                np.max(
                    np.abs(
                        source_descendants[global_id]
                        - candidate_descendants[global_id]
                    )
                )
            )
            for global_id in source_descendants
        ),
        default=0.0,
    )


def reset_assembly_origin(
    model: ifcopenshell.file,
    row: dict[str, Any],
    tolerance_mm: float,
) -> dict[str, Any]:
    assembly = model.by_guid(row["global_id"])
    if assembly is None or not assembly.is_a("IfcElementAssembly"):
        raise RuntimeError(f"missing IfcElementAssembly {row['global_id']}")
    if getattr(assembly, "Representation", None):
        raise RuntimeError(f"assembly unexpectedly owns Body geometry: {row['global_id']}")
    descendants = decomposed_products(assembly)
    descendant_world = {
        product.GlobalId: world_matrix(product)
        for _, product in descendants
        if getattr(product, "ObjectPlacement", None)
    }
    source_matrix = world_matrix(assembly)
    target_matrix = source_matrix.copy()
    target_matrix[:3, 3] = row["target_mm"]
    ifcopenshell.api.geometry.edit_object_placement(
        model,
        product=assembly,
        matrix=target_matrix,
        is_si=False,
        should_transform_children=False,
    )
    for _, product in sorted(descendants, key=lambda item: item[0]):
        if product.GlobalId not in descendant_world:
            continue
        ifcopenshell.api.geometry.edit_object_placement(
            model,
            product=product,
            matrix=descendant_world[product.GlobalId],
            is_si=False,
            should_transform_children=False,
        )
    placement_residual = float(
        np.max(np.abs(world_matrix(assembly)[:3, 3] - row["target_mm"]))
    )
    settings = geometry_settings()
    vertices, faces = descendant_mesh_mm(
        settings, decomposed_products(assembly)
    )
    anchor_distance = point_mesh_distance(row["target_mm"], vertices, faces)
    return {
        "global_id": assembly.GlobalId,
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
        "descendant_global_ids": sorted(descendant_world),
        "maximum_descendant_world_placement_residual_mm": None,
    }


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
    source_relations = relation_signature(source)
    results = [
        reset_assembly_origin(candidate, row, args.tolerance_mm)
        for row in rows
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    for result in results:
        result["maximum_descendant_world_placement_residual_mm"] = (
            maximum_descendant_world_placement_residual(
                source, candidate, result["global_id"]
            )
        )
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
        assembly = candidate.by_guid(row["global_id"])
        residual = float(
            np.max(np.abs(world_matrix(assembly)[:3, 3] - row["target_mm"]))
        )
        if residual > args.tolerance_mm:
            target_placement_mismatches.append(
                {"global_id": row["global_id"], "residual_mm": residual}
            )
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    entity_count_delta = len(list(candidate)) - len(list(source))
    gates = {
        "target_assemblies": len(rows),
        "target_placement_mismatches": target_placement_mismatches,
        "anchors_over_tolerance": sum(
            result["anchor_surface_distance_mm"] > args.tolerance_mm
            for result in results
        ),
        "descendant_placement_regressions": sum(
            result["maximum_descendant_world_placement_residual_mm"]
            > args.tolerance_mm
            for result in results
        ),
        "product_geometry_over_tolerance": len(geometry_over_tolerance),
        "schema_equal": source.schema == candidate.schema,
        "entity_count_delta": entity_count_delta,
        "entity_count_equal": entity_count_delta == 0,
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
        "aggregate_relations_equal": source_relations
        == relation_signature(candidate),
    }
    gates["pass"] = (
        not gates["target_placement_mismatches"]
        and gates["anchors_over_tolerance"] == 0
        and gates["descendant_placement_regressions"] == 0
        and gates["product_geometry_over_tolerance"] == 0
        and gates["schema_equal"]
        and gates["entity_count_equal"]
        and gates["root_global_ids_equal"]
        and gates["aggregate_relations_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-element-assembly-origin-reset-candidate",
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
