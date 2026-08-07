#!/usr/bin/env python3
"""Classify product origins that remain over the project review tolerance.

This audit is read-only.  It separates an origin residual from construction
geometry and from the object's semantic/install role.  No queue authorizes an
IFC write; each class still requires its own confirmed anchor candidate.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element

from coordinate_normalize import product_origin_records, sha256


SURFACE_CLASSES = {"IfcCovering", "IfcSlab"}
SERVICE_CLASSES = {
    "IfcElectricAppliance",
    "IfcLightFixture",
    "IfcFlowSegment",
    "IfcPipeSegment",
    "IfcWasteTerminal",
    "IfcSanitaryTerminal",
}
SECONDARY_BUILDING_CLASSES = {
    "IfcBeam",
    "IfcBuildingElementProxy",
    "IfcElementAssembly",
    "IfcRailing",
}


def review_queue(product: ifcopenshell.entity_instance) -> tuple[str, str]:
    if product.is_a("IfcDoor"):
        if not product.FillsVoids:
            return (
                "unhosted_door_exception",
                "IfcDoor has no IfcRelFillsElement; host and threshold anchor are unknown",
            )
        return (
            "hosted_door_anchor_review",
            "IfcDoor fills an Opening but still requires a verified threshold anchor",
        )
    if product.is_a("IfcOpeningElement"):
        host_ids = [relation.RelatingBuildingElement.GlobalId for relation in product.VoidsElements]
        fill_ids = [relation.RelatedBuildingElement.GlobalId for relation in product.HasFillings]
        if fill_ids:
            return (
                "filled_opening_anchor_review",
                f"Opening has host(s) {host_ids} and filling(s) {fill_ids}; the whole relationship must be checked together",
            )
        return (
            "unfilled_opening_anchor_review",
            f"Opening has host(s) {host_ids} but no filling; cut intent and Boolean safety must be preserved",
        )
    if product.is_a("IfcFurniture"):
        return (
            "furniture_installation_role_review",
            "complex furniture geometry is preserved until fixed/loose role and installation anchor are confirmed",
        )
    if product.is_a() in SURFACE_CLASSES:
        return (
            "surface_anchor_review",
            "surface role and completed-face anchor must be confirmed before placement or geometry changes",
        )
    if product.is_a() in SERVICE_CLASSES:
        return (
            "service_installation_anchor_review",
            "service object requires a confirmed centreline, face, drain, or mounting datum",
        )
    if product.is_a() in SECONDARY_BUILDING_CLASSES:
        return (
            "secondary_building_element_anchor_review",
            "secondary building element requires a class-specific construction datum",
        )
    return (
        "other_origin_review",
        "no class-specific construction anchor has been confirmed",
    )


def review_subgroup(
    product: ifcopenshell.entity_instance, queue: str
) -> tuple[str, str]:
    """Create discipline-sized queues without inferring an anchor or role."""

    if product.is_a("IfcOpeningElement"):
        host_classes = sorted(
            {
                relation.RelatingBuildingElement.is_a()
                for relation in product.VoidsElements
            }
        )
        label = "+".join(host_classes) if host_classes else "unhosted"
        return (
            f"opening_in_{label}_review",
            f"existing VoidsElements host class(es)={host_classes}",
        )
    if product.is_a("IfcFurniture"):
        container = ifcopenshell.util.element.get_container(product)
        container_name = getattr(container, "Name", None)
        if container_name in {"KITCHEN", "VVD", "BATHM", "NBW"}:
            return (
                "furniture_service_zone_review",
                f"existing container={container_name}; fixed/loose role remains unconfirmed",
            )
        if container_name in {"LIVING", "BEDM", "BEDG"}:
            return (
                "furniture_room_review",
                f"existing container={container_name}; fixed/loose role remains unconfirmed",
            )
        return (
            "furniture_other_context_review",
            f"existing container={container_name}; role remains unconfirmed",
        )
    if product.is_a("IfcCovering"):
        return "covering_anchor_review", "IFC class is IfcCovering"
    if product.is_a("IfcSlab"):
        container = ifcopenshell.util.element.get_container(product)
        return (
            f"slab_{str(getattr(container, 'Name', 'unassigned')).lower()}_anchor_review",
            f"existing container={getattr(container, 'Name', None)}; completed-face role remains unconfirmed",
        )
    service_subgroups = {
        "IfcLightFixture": "light_fixture_anchor_review",
        "IfcElectricAppliance": "electric_appliance_anchor_review",
        "IfcFlowSegment": "flow_segment_centreline_review",
        "IfcPipeSegment": "pipe_segment_centreline_review",
        "IfcSanitaryTerminal": "sanitary_terminal_anchor_review",
        "IfcWasteTerminal": "waste_terminal_anchor_review",
    }
    if product.is_a() in service_subgroups:
        return service_subgroups[product.is_a()], f"IFC class is {product.is_a()}"
    if product.is_a("IfcBeam"):
        return "beam_anchor_review", "IFC class is IfcBeam"
    if product.is_a("IfcBuildingElementProxy"):
        container = ifcopenshell.util.element.get_container(product)
        container_name = getattr(container, "Name", None)
        name = (product.Name or "").lower()
        if "ceiling" in name or container_name == "DCL":
            return (
                "ceiling_proxy_anchor_review",
                f"name={product.Name}; existing container={container_name}",
            )
        if any(
            token in name
            for token in ("furniture", "worktop", "cabinet", "countertop")
        ) or container_name in {"KITCHEN", "VVD", "NBW"}:
            return (
                "joinery_proxy_anchor_review",
                f"name={product.Name}; existing container={container_name}",
            )
        if container_name == "ELEC":
            return (
                "electrical_proxy_anchor_review",
                f"name={product.Name}; existing container=ELEC",
            )
        return (
            "other_proxy_anchor_review",
            f"name={product.Name}; existing container={container_name}",
        )
    if product.is_a("IfcElementAssembly"):
        container = ifcopenshell.util.element.get_container(product)
        container_name = getattr(container, "Name", None)
        child_classes = sorted(
            {
                child.is_a()
                for relation in product.IsDecomposedBy
                for child in relation.RelatedObjects
            }
        )
        if child_classes and set(child_classes) == {"IfcWall"}:
            return (
                "wall_assembly_anchor_review",
                "all direct children are IfcWall; use a represented child-wall construction boundary",
            )
        if child_classes and set(child_classes) == {"IfcFurniture"}:
            return (
                "joinery_assembly_anchor_review",
                "all direct children are fixed IfcFurniture; reuse a confirmed child installation datum",
            )
        if container_name == "WC" or (product.Name or "").startswith("WC"):
            return (
                "sanitary_assembly_anchor_review",
                f"name={product.Name}; container={container_name}; require a rough-in, connector-centre, or fixing datum",
            )
        return (
            "element_assembly_anchor_review",
            f"direct child classes={child_classes}; existing container={container_name}",
        )
    if product.is_a("IfcRailing"):
        return "railing_anchor_review", "IFC class is IfcRailing"
    if queue == "unhosted_door_exception":
        return queue, "door has no filling relationship"
    return queue, "no narrower evidence-only subgroup applies"


def audit(model: ifcopenshell.file, source: Path, tolerance_mm: float) -> dict[str, Any]:
    origin_records = product_origin_records(model, tolerance_mm)
    records = []
    for origin in origin_records:
        if origin["within_review_tolerance"]:
            continue
        product = model.by_id(origin["step_id"])
        queue, basis = review_queue(product)
        subgroup, subgroup_basis = review_subgroup(product, queue)
        relationship_evidence: dict[str, Any] = {}
        if product.is_a("IfcOpeningElement"):
            hosts = [
                relation.RelatingBuildingElement for relation in product.VoidsElements
            ]
            relationship_evidence = {
                "host_global_ids": [host.GlobalId for host in hosts],
                "host_classes": [host.is_a() for host in hosts],
                "filling_global_ids": [
                    relation.RelatedBuildingElement.GlobalId
                    for relation in product.HasFillings
                ],
                "placement_relative_to_host": any(
                    product.ObjectPlacement.PlacementRelTo == host.ObjectPlacement
                    for host in hosts
                ),
            }
        records.append(
            {
                **origin,
                "review_queue": queue,
                "basis": basis,
                "review_subgroup": subgroup,
                "review_subgroup_basis": subgroup_basis,
                "relationship_evidence": relationship_evidence,
                "automatic_write_allowed": False,
                "review_required": True,
            }
        )
    opening_records = [
        record for record in records if record["ifc_class"] == "IfcOpeningElement"
    ]
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-remaining-origin-review-audit",
        "source": {
            "path": str(source),
            "sha256": sha256(source),
            "schema": model.schema,
        },
        "tolerance_mm": tolerance_mm,
        "summary": {
            "all_product_origins": len(origin_records),
            "within_tolerance": sum(
                record["within_review_tolerance"] for record in origin_records
            ),
            "remaining_over_tolerance": len(records),
            "by_class": dict(Counter(record["ifc_class"] for record in records)),
            "by_review_queue": dict(
                Counter(record["review_queue"] for record in records)
            ),
            "by_review_subgroup": dict(
                Counter(record["review_subgroup"] for record in records)
            ),
            "automatic_write_allowed": 0,
            "opening_relationships": {
                "remaining_openings": len(opening_records),
                "relative_to_host_placement": sum(
                    record["relationship_evidence"].get(
                        "placement_relative_to_host", False
                    )
                    for record in opening_records
                ),
                "not_relative_to_host_placement": sum(
                    not record["relationship_evidence"].get(
                        "placement_relative_to_host", False
                    )
                    for record in opening_records
                ),
                "with_fillings": sum(
                    bool(record["relationship_evidence"].get("filling_global_ids"))
                    for record in opening_records
                ),
                "host_class_counts": dict(
                    Counter(
                        host_class
                        for record in opening_records
                        for host_class in record["relationship_evidence"].get(
                            "host_classes", []
                        )
                    )
                ),
            },
        },
        "records": records,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.tolerance_mm <= 0.0:
        raise SystemExit("--tolerance-mm must be positive")
    source = args.input.resolve()
    report = audit(ifcopenshell.open(source), source, args.tolerance_mm)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **report["summary"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
