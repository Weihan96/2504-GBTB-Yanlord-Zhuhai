#!/usr/bin/env python3
"""Build and verify the IFC4 Space relationship and gross-area candidate."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api.aggregate
import ifcopenshell.api.pset
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.shape

from geometry_alignment_audit import geometry_difference_audit, sha256


def footprint_areas(model: ifcopenshell.file) -> dict[str, float]:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    result = {}
    for space in model.by_type("IfcSpace"):
        shape = ifcopenshell.geom.create_shape(settings, space)
        result[space.GlobalId] = float(
            ifcopenshell.util.shape.get_footprint_area(shape.geometry)
        )
    return result


def apply_space_semantics(
    model: ifcopenshell.file,
    storey_global_id: str,
    areas: dict[str, float],
) -> dict[str, Any]:
    storey = model.by_guid(storey_global_id)
    if storey is None or not storey.is_a("IfcBuildingStorey"):
        raise RuntimeError(f"{storey_global_id} is not an IfcBuildingStorey")
    spaces = model.by_type("IfcSpace")
    if {space.GlobalId for space in spaces} != set(areas):
        raise RuntimeError("Area evidence does not match the current Space set")

    mixed_relations = [
        relation
        for relation in model.by_type("IfcRelContainedInSpatialStructure")
        if relation.RelatingStructure == storey
        and any(element.is_a("IfcSpace") for element in relation.RelatedElements)
    ]
    if len(mixed_relations) != 1:
        raise RuntimeError(
            f"Expected one mixed Space containment relation; found {len(mixed_relations)}"
        )
    mixed_relation = mixed_relations[0]
    source_related_count = len(mixed_relation.RelatedElements)
    non_spaces = tuple(
        element
        for element in mixed_relation.RelatedElements
        if not element.is_a("IfcSpace")
    )
    mixed_relation.RelatedElements = non_spaces
    aggregation = ifcopenshell.api.aggregate.assign_object(
        model, products=spaces, relating_object=storey
    )

    qto_results = []
    for space in spaces:
        existing = ifcopenshell.util.element.get_psets(space, qtos_only=True)
        if "Qto_SpaceBaseQuantities" in existing:
            raise RuntimeError(
                f"{space.GlobalId} already has Qto_SpaceBaseQuantities"
            )
        qto = ifcopenshell.api.pset.add_qto(
            model, product=space, name="Qto_SpaceBaseQuantities"
        )
        ifcopenshell.api.pset.edit_qto(
            model,
            qto=qto,
            properties={"GrossFloorArea": areas[space.GlobalId]},
        )
        qto_results.append(
            {
                "global_id": space.GlobalId,
                "long_name": space.LongName,
                "gross_floor_area_m2": areas[space.GlobalId],
                "qto_global_id": qto.GlobalId,
            }
        )
    return {
        "storey_global_id": storey.GlobalId,
        "source_mixed_containment_global_id": mixed_relation.GlobalId,
        "source_mixed_related_count": source_related_count,
        "result_non_space_containment_count": len(non_spaces),
        "aggregation_global_id": aggregation.GlobalId,
        "aggregated_space_count": len(aggregation.RelatedObjects),
        "qto_results": qto_results,
    }


def validate_candidate(
    source: ifcopenshell.file,
    candidate: ifcopenshell.file,
    storey_global_id: str,
    areas: dict[str, float],
    tolerance_mm: float,
) -> dict[str, Any]:
    space_ids = set(areas)
    remaining_contained_spaces = [
        element.GlobalId
        for relation in candidate.by_type("IfcRelContainedInSpatialStructure")
        for element in relation.RelatedElements
        if element.is_a("IfcSpace")
    ]
    aggregations = [
        relation
        for relation in candidate.by_type("IfcRelAggregates")
        if relation.RelatingObject.GlobalId == storey_global_id
        and {item.GlobalId for item in relation.RelatedObjects if item.is_a("IfcSpace")}
        == space_ids
    ]
    qto_values = {}
    missing_or_extra_qto = []
    net_floor_area_present = []
    for space in candidate.by_type("IfcSpace"):
        qtos = ifcopenshell.util.element.get_psets(space, qtos_only=True)
        base = qtos.get("Qto_SpaceBaseQuantities", {})
        if "GrossFloorArea" not in base:
            missing_or_extra_qto.append(space.GlobalId)
            continue
        qto_values[space.GlobalId] = float(base["GrossFloorArea"])
        if "NetFloorArea" in base:
            net_floor_area_present.append(space.GlobalId)
    qto_values_match = set(qto_values) == space_ids and all(
        abs(qto_values[global_id] - areas[global_id]) <= 1e-9
        for global_id in space_ids
    )

    difference = geometry_difference_audit(
        candidate,
        source,
        "source",
        tolerance_mm=tolerance_mm,
        classes=["IfcSpace"],
        global_ids=[],
    )
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    existing_root_ids_preserved = source_root_ids <= candidate_root_ids
    return {
        "remaining_contained_spaces": remaining_contained_spaces,
        "matching_aggregation_count": len(aggregations),
        "qto_values_match": qto_values_match,
        "missing_or_extra_qto": missing_or_extra_qto,
        "net_floor_area_present": net_floor_area_present,
        "gross_floor_area_total_m2": sum(qto_values.values()),
        "space_geometry_difference": difference,
        "existing_root_global_ids_preserved": existing_root_ids_preserved,
        "added_root_global_ids": sorted(candidate_root_ids - source_root_ids),
        "schema_equal": source.schema == candidate.schema,
        "space_count_equal": len(source.by_type("IfcSpace"))
        == len(candidate.by_type("IfcSpace")),
        "pass": (
            not remaining_contained_spaces
            and len(aggregations) == 1
            and qto_values_match
            and not missing_or_extra_qto
            and not net_floor_area_present
            and difference["over_tolerance"] == 0
            and existing_root_ids_preserved
            and source.schema == candidate.schema
            and len(source.by_type("IfcSpace"))
            == len(candidate.by_type("IfcSpace"))
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--storey-global-id", required=True)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    areas = footprint_areas(source)
    application = apply_space_semantics(
        candidate, args.storey_global_id, areas
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    validation = validate_candidate(
        source,
        candidate,
        args.storey_global_id,
        areas,
        args.tolerance_mm,
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-space-semantics-candidate",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "candidate": {
            "path": str(args.output),
            "sha256": sha256(args.output),
            "schema": candidate.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "area_basis": "IfcSpace world-geometry footprint; confirmed Grid semantic boundary",
        "application": application,
        "validation": validation,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "report": str(args.report),
                "candidate": str(args.output),
                "candidate_sha256": report["candidate"]["sha256"],
                "space_count": len(areas),
                "gross_floor_area_total_m2": validation[
                    "gross_floor_area_total_m2"
                ],
                "space_geometry_over_tolerance": validation[
                    "space_geometry_difference"
                ]["over_tolerance"],
                "pass": validation["pass"],
            },
            ensure_ascii=False,
        )
    )
    if not validation["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
