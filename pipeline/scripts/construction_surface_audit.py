#!/usr/bin/env python3
"""Read-only semantic and anchor audit for construction surfaces and fixtures."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.placement
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union

from geometry_alignment_audit import sha256


AUDITED_CLASSES = (
    "IfcCovering",
    "IfcSlab",
    "IfcFurniture",
    "IfcSanitaryTerminal",
)


def geometry_orientation(dimensions_mm: list[float], thin_mm: float) -> str:
    x, y, z = dimensions_mm
    if z <= thin_mm and max(x, y) >= 500.0:
        return "horizontal_thin"
    if z >= 100.0 and min(x, y) <= thin_mm:
        return "vertical_thin"
    return "other"


def semantic_status(
    product: ifcopenshell.entity_instance, orientation: str
) -> tuple[str, str]:
    occurrence_predefined = getattr(product, "PredefinedType", None)
    type_object = ifcopenshell.util.element.get_type(product)
    type_predefined = getattr(type_object, "PredefinedType", None)
    predefined = (
        occurrence_predefined
        if occurrence_predefined not in {None, "NOTDEFINED"}
        else type_predefined
    )
    name = product.Name or ""
    if predefined and predefined not in {"NOTDEFINED", "USERDEFINED"}:
        source = (
            "occurrence"
            if occurrence_predefined not in {None, "NOTDEFINED"}
            else "type"
        )
        return (
            "explicit_ifc_semantics",
            f"effective PredefinedType={predefined} from {source}",
        )
    if product.is_a("IfcCovering") and name.startswith("Baseboard"):
        return (
            "name_candidate_requires_review",
            "Name suggests baseboard; IFC PredefinedType is missing",
        )
    if product.is_a("IfcCovering"):
        return (
            "orientation_candidate_requires_review",
            f"{orientation} geometry is not sufficient to assign covering semantics",
        )
    if product.is_a("IfcSlab"):
        return (
            "missing_predefined_type_requires_review",
            "IfcSlab PredefinedType is missing",
        )
    if product.is_a("IfcFurniture"):
        return (
            "preserve_complex_geometry_until_installation_review",
            "Furniture requires a confirmed fixed-installation role and anchor",
        )
    return (
        "fixture_candidate_requires_review",
        "Sanitary fixture location is meaningful but installation anchor is unconfirmed",
    )


def automatic_anchor_write_allowed(
    semantic_status_value: str,
    origin_within_tolerance: bool,
    origin_at_vertex_within_tolerance: bool,
) -> tuple[bool, str]:
    """Return the anchor write gate for this generic read-only audit.

    Existing IFC semantics, an integer origin, or an origin that happens to
    coincide with a mesh vertex do not prove which construction anchor the
    object should use.  A class-specific candidate builder must supply that
    proof before any placement or geometry write is allowed.
    """

    evidence = []
    if semantic_status_value == "explicit_ifc_semantics":
        evidence.append("IFC semantics are already explicit")
    if origin_within_tolerance:
        evidence.append("origin is already within integer tolerance")
    if origin_at_vertex_within_tolerance:
        evidence.append("origin coincides with a geometry vertex")
    prefix = "; ".join(evidence) or "no generic anchor evidence"
    return (
        False,
        f"{prefix}; a class-specific confirmed construction anchor is still required",
    )


def review_batch(
    ifc_class: str,
    name: str | None,
    predefined_type: str | None,
    material_name: str | None,
    orientation: str,
    bbox_min_mm: list[float],
    bbox_max_mm: list[float],
) -> str:
    """Group records for Blender review without assigning IFC semantics."""

    object_name = name or ""
    if ifc_class == "IfcCovering":
        if predefined_type and predefined_type not in {"NOTDEFINED", "USERDEFINED"}:
            return "covering_existing_explicit_semantics"
        if object_name.startswith("Baseboard"):
            return "covering_baseboard_name_candidates"
        if material_name == "TerrazzoMosaicTile" and orientation == "horizontal_thin":
            if bbox_max_mm[2] <= -900.0:
                return "covering_tile_instances_below_ffl"
            return "covering_tile_instances_near_ffl"
        if orientation == "vertical_thin":
            return "covering_vertical_candidates"
        if orientation == "horizontal_thin":
            return "covering_horizontal_candidates"
        return "covering_other_candidates"
    if ifc_class == "IfcSlab":
        if bbox_max_mm[2] <= 0.1:
            return "slab_below_or_at_ffl_candidates"
        if bbox_min_mm[2] >= 2400.0:
            return "slab_above_ceiling_candidates"
        return "slab_other_candidates"
    if ifc_class == "IfcFurniture":
        return "furniture_fixed_or_loose_review"
    return "sanitary_installation_anchor_review"


def review_subgroup(
    ifc_class: str,
    container_name: str | None,
    batch: str,
    bbox_min_mm: list[float],
    bbox_max_mm: list[float],
    dimensions_mm: list[float],
) -> tuple[str, str]:
    """Split broad review batches using only observable model evidence.

    These labels describe position, dimensions, and the existing IFC
    container. They deliberately avoid assigning FLOORING, CLADDING, ROOF,
    fixed-furniture, or installation-anchor semantics.
    """

    container = container_name or "UNASSIGNED"
    height_mm = dimensions_mm[2]
    if ifc_class == "IfcCovering":
        if batch == "covering_vertical_candidates":
            if bbox_min_mm[2] >= 2300.0 and height_mm <= 500.0:
                return (
                    "covering_high_level_vertical_strips",
                    "vertical thin geometry starts at or above 2300 mm and is at most 500 mm high",
                )
            return (
                "covering_vertical_full_height_surfaces",
                "vertical thin geometry is taller than the high-level strip gate",
            )
        if batch == "covering_horizontal_candidates":
            if bbox_min_mm[2] >= 2400.0:
                return (
                    "covering_high_level_horizontal_strips",
                    "horizontal thin geometry is entirely at or above 2400 mm",
                )
            if 500.0 <= bbox_min_mm[2] <= 700.0:
                return (
                    "covering_raised_horizontal_surfaces",
                    "horizontal thin geometry starts between 500 and 700 mm",
                )
        return batch, "no narrower evidence-only covering subgroup applies"
    if ifc_class == "IfcSlab":
        return (
            f"slab_{container.lower()}_geometry_review",
            f"existing container={container}; no slab PredefinedType inferred",
        )
    if ifc_class == "IfcFurniture":
        if container in {"KITCHEN", "VVD", "BATHM", "NBW"}:
            return (
                "furniture_service_zone_fixed_or_loose_review",
                f"existing container={container}; fixed-installation role remains unconfirmed",
            )
        if container in {"LIVING", "BEDM", "BEDG"}:
            return (
                "furniture_room_fixed_or_loose_review",
                f"existing container={container}; loose/fixed role remains unconfirmed",
            )
        return (
            "furniture_other_context_review",
            f"existing container={container}; role remains unconfirmed",
        )
    if container in {"BATHM", "BATHG", "WC"}:
        return (
            "sanitary_bathroom_anchor_review",
            f"existing container={container}; installation anchor remains unconfirmed",
        )
    return (
        "sanitary_other_context_anchor_review",
        f"existing container={container}; installation role and anchor remain unconfirmed",
    )


def projected_footprint(
    vertices_mm: np.ndarray,
    faces: np.ndarray,
):
    polygons = []
    for face in faces:
        polygon = Polygon(vertices_mm[face, :2])
        if polygon.is_valid and polygon.area > 1e-6:
            polygons.append(polygon)
    return unary_union(polygons) if polygons else Polygon()


def footprint_space_overlaps(footprint, space_footprints) -> list[dict[str, Any]]:
    if footprint.is_empty or footprint.area <= 1e-6:
        return []
    result = []
    for space in space_footprints:
        intersection_area = footprint.intersection(space["geometry"]).area
        if intersection_area <= 1e-6:
            continue
        result.append(
            {
                "global_id": space["global_id"],
                "long_name": space["long_name"],
                "overlap_area_mm2": intersection_area,
                "covering_overlap_ratio": intersection_area / footprint.area,
            }
        )
    return sorted(
        result,
        key=lambda item: (-item["covering_overlap_ratio"], item["global_id"]),
    )


def build_space_footprints(
    model: ifcopenshell.file,
    settings: ifcopenshell.geom.settings,
) -> list[dict[str, Any]]:
    result = []
    for space in model.by_type("IfcSpace"):
        shape = ifcopenshell.geom.create_shape(settings, space)
        vertices = np.array(shape.geometry.verts, dtype=float).reshape(-1, 3) * 1000.0
        faces = np.array(shape.geometry.faces, dtype=int).reshape(-1, 3)
        result.append(
            {
                "global_id": space.GlobalId,
                "long_name": space.LongName,
                "geometry": projected_footprint(vertices, faces),
            }
        )
    return result


def audit_product(
    model: ifcopenshell.file,
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
    tolerance_mm: float,
    thin_mm: float,
    space_footprints: list[dict[str, Any]],
) -> dict[str, Any]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.array(shape.geometry.verts, dtype=float).reshape(-1, 3) * 1000.0
    faces = np.array(shape.geometry.faces, dtype=int).reshape(-1, 3)
    minimum = vertices.min(axis=0)
    maximum = vertices.max(axis=0)
    dimensions = maximum - minimum
    orientation = geometry_orientation(dimensions.tolist(), thin_mm)
    placement = np.array(
        ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
        dtype=float,
    )
    origin = placement[:3, 3]
    integer_delta = np.round(origin) - origin
    nearest_vertex_distance = float(
        np.min(np.linalg.norm(vertices - origin, axis=1))
    )
    container = ifcopenshell.util.element.get_container(product)
    material = ifcopenshell.util.element.get_material(
        product, should_skip_usage=False
    )
    type_object = ifcopenshell.util.element.get_type(product)
    occurrence_predefined = getattr(product, "PredefinedType", None)
    type_predefined = getattr(type_object, "PredefinedType", None)
    effective_predefined = (
        occurrence_predefined
        if occurrence_predefined not in {None, "NOTDEFINED"}
        else type_predefined
    )
    status, basis = semantic_status(product, orientation)
    origin_within_tolerance = bool(
        np.max(np.abs(integer_delta)) <= tolerance_mm
    )
    origin_at_vertex_within_tolerance = (
        nearest_vertex_distance <= tolerance_mm
    )
    anchor_write_allowed, anchor_write_basis = automatic_anchor_write_allowed(
        status,
        origin_within_tolerance,
        origin_at_vertex_within_tolerance,
    )
    batch = review_batch(
        product.is_a(),
        product.Name,
        effective_predefined,
        getattr(material, "Name", None) if material else None,
        orientation,
        minimum.tolist(),
        maximum.tolist(),
    )
    subgroup, subgroup_basis = review_subgroup(
        product.is_a(),
        container.Name if container else None,
        batch,
        minimum.tolist(),
        maximum.tolist(),
        dimensions.tolist(),
    )
    overlaps = (
        footprint_space_overlaps(
            projected_footprint(vertices, faces), space_footprints
        )
        if orientation == "horizontal_thin"
        else []
    )
    primary_space = overlaps[0] if overlaps else None
    return {
        "ifc_class": product.is_a(),
        "global_id": product.GlobalId,
        "name": product.Name,
        "object_type": getattr(product, "ObjectType", None),
        "predefined_type": occurrence_predefined,
        "type_predefined_type": type_predefined,
        "effective_predefined_type": effective_predefined,
        "container": container.Name if container else None,
        "material": getattr(material, "Name", None) if material else None,
        "orientation": orientation,
        "bbox_min_mm": minimum.tolist(),
        "bbox_max_mm": maximum.tolist(),
        "dimensions_mm": dimensions.tolist(),
        "origin_mm": origin.tolist(),
        "origin_integer_delta_mm": integer_delta.tolist(),
        "origin_within_tolerance": origin_within_tolerance,
        "nearest_vertex_distance_mm": nearest_vertex_distance,
        "origin_at_vertex_within_tolerance": origin_at_vertex_within_tolerance,
        "semantic_status": status,
        "basis": basis,
        "semantic_already_explicit": status == "explicit_ifc_semantics",
        "automatic_write_allowed": anchor_write_allowed,
        "automatic_write_basis": anchor_write_basis,
        "review_batch": batch,
        "review_subgroup": subgroup,
        "review_subgroup_basis": subgroup_basis,
        "space_overlaps": overlaps,
        "primary_space_global_id": (
            primary_space["global_id"] if primary_space else None
        ),
        "primary_space_long_name": (
            primary_space["long_name"] if primary_space else None
        ),
        "primary_space_overlap_ratio": (
            primary_space["covering_overlap_ratio"] if primary_space else None
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--thin-threshold-mm", type=float, default=100.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.tolerance_mm <= 0.0:
        raise SystemExit("--tolerance-mm must be positive")
    if args.thin_threshold_mm <= 0.0:
        raise SystemExit("--thin-threshold-mm must be positive")
    source_path = args.input.resolve()
    model = ifcopenshell.open(source_path)
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    space_footprints = build_space_footprints(model, settings)
    records = [
        audit_product(
            model,
            settings,
            product,
            args.tolerance_mm,
            args.thin_threshold_mm,
            space_footprints,
        )
        for ifc_class in AUDITED_CLASSES
        for product in model.by_type(ifc_class)
        if product.ObjectPlacement and product.Representation
    ]
    summary = {
        "records": len(records),
        "by_class": dict(Counter(record["ifc_class"] for record in records)),
        "by_orientation": dict(
            Counter(record["orientation"] for record in records)
        ),
        "by_semantic_status": dict(
            Counter(record["semantic_status"] for record in records)
        ),
        "by_review_batch": dict(
            Counter(record["review_batch"] for record in records)
        ),
        "by_review_subgroup": dict(
            Counter(record["review_subgroup"] for record in records)
        ),
        "horizontal_by_primary_space": dict(
            Counter(
                record["primary_space_long_name"]
                for record in records
                if record["primary_space_long_name"]
            )
        ),
        "origins_over_tolerance": sum(
            not record["origin_within_tolerance"] for record in records
        ),
        "origins_at_vertex_within_tolerance": sum(
            record["origin_at_vertex_within_tolerance"] for record in records
        ),
        "semantic_already_explicit": sum(
            record["semantic_already_explicit"] for record in records
        ),
        "automatic_write_allowed": sum(
            record["automatic_write_allowed"] for record in records
        ),
        "human_review_required": sum(
            not record["automatic_write_allowed"] for record in records
        ),
    }
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-construction-surface-anchor-audit",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": model.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "thin_threshold_mm": args.thin_threshold_mm,
        "summary": summary,
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **summary}, ensure_ascii=False))


if __name__ == "__main__":
    main()
