#!/usr/bin/env python3
"""Build verification-only evidence for the rejected PVC110 split option."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.guid
import ifcopenshell.util.element
import numpy as np

from direction_noise_candidate import (
    all_product_geometry_difference,
    vertex_hausdorff_within_tolerance_mm,
)
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm


TARGET_IDS = (
    "178mqyyzzFowLcbXcH6prO",
    "0bfVg4Ys1CevZs$qxhkXTo",
)
EXPECTED_NAME = "Sewage Pipe PVC110"
EXPECTED_ITEM_CLASS = "IfcPolygonalFaceSet"
EXPECTED_BUNDLE_SIZE = 3


def body_representation(product: Any) -> Any:
    representations = [
        representation
        for representation in product.Representation.Representations
        if representation.RepresentationIdentifier == "Body"
    ]
    if len(representations) != 1:
        raise RuntimeError(f"expected one Body: {product.GlobalId}")
    return representations[0]


def clone_local_placement(model: ifcopenshell.file, placement: Any) -> Any:
    relative = placement.RelativePlacement
    location = model.create_entity(
        "IfcCartesianPoint", tuple(relative.Location.Coordinates)
    )
    axis = (
        model.create_entity("IfcDirection", tuple(relative.Axis.DirectionRatios))
        if getattr(relative, "Axis", None)
        else None
    )
    ref_direction = (
        model.create_entity(
            "IfcDirection", tuple(relative.RefDirection.DirectionRatios)
        )
        if getattr(relative, "RefDirection", None)
        else None
    )
    relative_copy = model.create_entity(
        relative.is_a(), Location=location, Axis=axis, RefDirection=ref_direction
    )
    return model.create_entity(
        "IfcLocalPlacement",
        PlacementRelTo=placement.PlacementRelTo,
        RelativePlacement=relative_copy,
    )


def item_bbox_distance(first: Any, second: Any) -> float:
    first_points = np.array(first.Coordinates.CoordList, dtype=float)
    second_points = np.array(second.Coordinates.CoordList, dtype=float)
    first_min, first_max = first_points.min(axis=0), first_points.max(axis=0)
    second_min, second_max = second_points.min(axis=0), second_points.max(axis=0)
    separation = np.maximum(
        np.maximum(first_min - second_max, second_min - first_max), 0.0
    )
    return float(np.linalg.norm(separation))


def validate_source(model: ifcopenshell.file) -> list[dict[str, Any]]:
    inventory = []
    for global_id in TARGET_IDS:
        product = model.by_guid(global_id)
        representation = body_representation(product)
        container = ifcopenshell.util.element.get_container(product)
        assigned_type = ifcopenshell.util.element.get_type(product)
        material = ifcopenshell.util.element.get_material(product)
        psets = ifcopenshell.util.element.get_psets(product)
        items = tuple(representation.Items)
        if (
            not product.is_a("IfcFlowSegment")
            or product.Name != EXPECTED_NAME
            or len(items) != EXPECTED_BUNDLE_SIZE
            or any(not item.is_a(EXPECTED_ITEM_CLASS) for item in items)
            or getattr(container, "Name", None) != "WC"
            or assigned_type is not None
            or material is not None
            or psets
        ):
            raise RuntimeError(f"PVC110 bundle identity mismatch: {global_id}")
        pair_distances = [
            item_bbox_distance(items[first], items[second])
            for first in range(len(items))
            for second in range(first + 1, len(items))
        ]
        if min(pair_distances) <= 0.1:
            raise RuntimeError(f"PVC110 bundle items are not disjoint: {global_id}")
        inventory.append(
            {
                "global_id": global_id,
                "ifc_class": product.is_a(),
                "name": product.Name,
                "container": container.Name,
                "body_item_ids": [item.id() for item in items],
                "pair_bbox_distances_mm": pair_distances,
                "minimum_pair_bbox_distance_mm": min(pair_distances),
                "assigned_type": None,
                "material": None,
                "psets": {},
            }
        )
    return inventory


def split_bundle(
    model: ifcopenshell.file, global_id: str
) -> tuple[list[Any], list[int]]:
    source = model.by_guid(global_id)
    representation = body_representation(source)
    items = tuple(representation.Items)
    containment = [
        inverse
        for inverse in model.get_inverse(source)
        if inverse.is_a("IfcRelContainedInSpatialStructure")
        and source in inverse.RelatedElements
    ]
    if len(containment) != 1:
        raise RuntimeError(f"expected one containment relation: {global_id}")

    source.Name = f"{EXPECTED_NAME}-01"
    representation.Items = (items[0],)
    products = [source]
    added_entity_ids: list[int] = []
    for index, item in enumerate(items[1:], start=2):
        new_representation = model.create_entity(
            "IfcShapeRepresentation",
            ContextOfItems=representation.ContextOfItems,
            RepresentationIdentifier=representation.RepresentationIdentifier,
            RepresentationType=representation.RepresentationType,
            Items=(item,),
        )
        new_shape = model.create_entity(
            "IfcProductDefinitionShape",
            Name=source.Representation.Name,
            Description=source.Representation.Description,
            Representations=(new_representation,),
        )
        new_placement = clone_local_placement(model, source.ObjectPlacement)
        new_product = model.create_entity(
            "IfcFlowSegment",
            GlobalId=ifcopenshell.guid.new(),
            OwnerHistory=source.OwnerHistory,
            Name=f"{EXPECTED_NAME}-{index:02d}",
            Description=source.Description,
            ObjectType=source.ObjectType,
            ObjectPlacement=new_placement,
            Representation=new_shape,
            Tag=source.Tag,
        )
        products.append(new_product)
        added_entity_ids.extend(
            [
                new_product.id(),
                new_representation.id(),
                new_shape.id(),
                new_placement.id(),
            ]
        )
    relation = containment[0]
    relation.RelatedElements = tuple(relation.RelatedElements) + tuple(products[1:])
    return products, added_entity_ids


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    inventory = validate_source(source)
    source_roots = {root.GlobalId for root in source.by_type("IfcRoot")}
    source_entity_count = len(list(source))
    source_meshes = {
        global_id: world_mesh_mm(geometry_settings(), source.by_guid(global_id))[0]
        for global_id in TARGET_IDS
    }

    candidate = ifcopenshell.open(source_path)
    split_groups = {}
    for global_id in TARGET_IDS:
        products, added_entity_ids = split_bundle(candidate, global_id)
        split_groups[global_id] = {
            "products": products,
            "added_entity_ids": added_entity_ids,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)

    group_checks = []
    new_global_ids = []
    for global_id, group in split_groups.items():
        group_ids = [product.GlobalId for product in group["products"]]
        new_global_ids.extend(group_ids[1:])
        candidate_vertices = []
        item_counts = []
        names = []
        for product_id in group_ids:
            product = candidate.by_guid(product_id)
            vertices, _ = world_mesh_mm(geometry_settings(), product)
            candidate_vertices.extend(vertices)
            item_counts.append(len(body_representation(product).Items))
            names.append(product.Name)
        group_checks.append(
            {
                "source_global_id": global_id,
                "candidate_global_ids": group_ids,
                "candidate_names": names,
                "body_item_counts": item_counts,
                "union_vertex_hausdorff_mm": vertex_hausdorff_within_tolerance_mm(
                    source_meshes[global_id], candidate_vertices, args.tolerance_mm
                ),
                "same_container": all(
                    getattr(ifcopenshell.util.element.get_container(candidate.by_guid(gid)), "Name", None)
                    == "WC"
                    for gid in group_ids
                ),
            }
        )

    all_geometry = all_product_geometry_difference(
        candidate, source, args.tolerance_mm
    )
    allowed_ids = set(TARGET_IDS) | set(new_global_ids)
    non_target_changes = [
        row["global_id"]
        for row in all_geometry["records"]
        if row["global_id"] not in allowed_ids and not row["within_tolerance"]
    ]
    candidate_roots = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    added_roots = candidate_roots - source_roots
    removed_roots = source_roots - candidate_roots
    gates = {
        "source_bundles": len(TARGET_IDS),
        "candidate_pipe_segments": sum(
            len(row["candidate_global_ids"]) for row in group_checks
        ),
        "all_products_have_one_body_item": all(
            counts == [1, 1, 1]
            for counts in (row["body_item_counts"] for row in group_checks)
        ),
        "union_geometry_over_tolerance": sum(
            row["union_vertex_hausdorff_mm"] > args.tolerance_mm
            for row in group_checks
        ),
        "all_containers_preserved": all(row["same_container"] for row in group_checks),
        "non_target_geometry_changes": non_target_changes,
        "schema_equal": source.schema == candidate.schema,
        "root_global_ids_removed": sorted(removed_roots),
        "root_global_ids_added": sorted(added_roots),
        "expected_new_global_ids": sorted(new_global_ids),
        "root_delta_matches_expected": added_roots == set(new_global_ids)
        and not removed_roots,
        "entity_count_delta": len(list(candidate)) - source_entity_count,
    }
    gates["pass"] = (
        gates["candidate_pipe_segments"] == 6
        and gates["all_products_have_one_body_item"]
        and gates["union_geometry_over_tolerance"] == 0
        and gates["all_containers_preserved"]
        and not gates["non_target_geometry_changes"]
        and gates["schema_equal"]
        and gates["root_delta_matches_expected"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-pvc110-bundle-split-candidate",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
            "entity_count": source_entity_count,
        },
        "candidate": {
            "path": str(args.output.resolve()),
            "sha256": sha256(args.output),
            "entity_count": len(list(candidate)),
        },
        "tolerance_mm": args.tolerance_mm,
        "source_inventory": inventory,
        "group_checks": group_checks,
        "all_product_geometry": all_geometry,
        "gates": gates,
        "formal_write_allowed": False,
        "decision_status": "rejected_for_formal_write",
        "disposition": "verification_evidence_only_preserve_two_source_products",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "report": str(args.report),
                "candidate": str(args.output),
                "candidate_pipe_segments": gates["candidate_pipe_segments"],
                "entity_count_delta": gates["entity_count_delta"],
                "root_global_ids_added": len(gates["root_global_ids_added"]),
                "pass": gates["pass"],
                "formal_write_allowed": False,
                "decision_status": "rejected_for_formal_write",
            },
            ensure_ascii=False,
        )
    )
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
