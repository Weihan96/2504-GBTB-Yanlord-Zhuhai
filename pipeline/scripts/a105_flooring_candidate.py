#!/usr/bin/env python3
"""Build a read-only A-105 flooring candidate for the confirmed tile geometry.

Eighteen near-FFL TerrazzoMosaicTile IfcCovering objects are assigned
PredefinedType=FLOORING. Four complete tiles at -1020..-1000 mm can be kept or
deleted as an explicit human-reviewed action. The source IFC is never changed.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api.root
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.placement

from direction_noise_candidate import all_product_geometry_difference, entity_changes
from geometry_alignment_audit import geometry_difference_audit, sha256


LOWER_TILE_GLOBAL_IDS = {
    "23YqpP7lXCT978tPUoh3YX",
    "0SbDT5whT0SPIUZKEmMh8x",
    "3_AXo$1Tn7MOj7mUw5rkr5",
    "27RE3Wo$r52gPaI1T2c_af",
}

# Each lower reference tile owns one unfilled drain Opening. Removing the host
# through ifcopenshell.api.root.remove_product intentionally removes that child
# as well, so the deletion boundary must name and validate all eight IfcRoot
# objects instead of treating the four dependent Openings as unexpected loss.
LOWER_DEPENDENT_OPENING_HOSTS = {
    "25FS_Tf$97FAs6LZDDTfW1": "0SbDT5whT0SPIUZKEmMh8x",
    "0Q2QlDeXj8xvVW7jQMDapI": "3_AXo$1Tn7MOj7mUw5rkr5",
    "0cIV3cLCnAeQHv95vrwRz0": "23YqpP7lXCT978tPUoh3YX",
    "0kYDUK7VXFWwK2GdO681O0": "27RE3Wo$r52gPaI1T2c_af",
}

# The visible linear-drain slot is cut into eight current-flooring tiles by
# eight host-relative IfcOpeningElement objects.  Each bathroom has two tiles
# on each side of the slot.  The opening bodies, rather than tile placement,
# own the slot width.
DRAIN_GAP_GROUPS = {
    "guest_bathroom": {
        "low_y_tile_ids": {
            "1CcaP1Bo53phkUAss7qLPT",
            "38mTGSlgvAK8S$ghNqax$m",
        },
        "high_y_tile_ids": {
            "0AN5MRpZ16dRCAg1l5TOyz",
            "3k9cT8nxj3o8yzLjoQc_Yt",
        },
        "swept_opening_ids": {
            "1Qca7wJGL0Q9BV8PW2gwGq",
            "04DzBdgQr99hTwEiZ4Zo0q",
        },
        "tessellated_opening_ids": {
            "2KV8sWa1DC9RDyHYIRPBOP",
            "3qBJcTlbjDOAe_s0cVkL7N",
        },
    },
    "main_bathroom": {
        "low_y_tile_ids": {
            "1XEnNJwc177xHjBydgoQFl",
            "3bR3wts7vEBQHsVMX7kiaY",
        },
        "high_y_tile_ids": {
            "18Hiw4PZP4rQrmw23BQQys",
            "2dEG2g8YT6te1$zV2kUnlv",
        },
        "swept_opening_ids": {
            "1dnB5Jx9b64ft$22nM2as1",
            "0pITNrup52k9_jELXkC_F8",
        },
        "tessellated_opening_ids": {
            "0A58600NfA7QWCD97tbUdf",
            "08Sns8O7n0j9zl5hNpovWh",
        },
    },
}


def expected_removed_global_ids(lower_action: str) -> set[str]:
    if lower_action == "keep":
        return set()
    if lower_action == "delete":
        return LOWER_TILE_GLOBAL_IDS | set(LOWER_DEPENDENT_OPENING_HOSTS)
    raise ValueError(f"unsupported lower action {lower_action}")


def expected_removed_root_global_ids(
    lower_action: str,
    dependency_records: list[dict[str, Any]],
) -> set[str]:
    result = expected_removed_global_ids(lower_action)
    if lower_action == "delete":
        result |= {
            global_id
            for record in dependency_records
            for global_id in record["void_relationship_global_ids"]
        }
    return result


def lower_dependency_inventory(model: ifcopenshell.file) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for opening_global_id, expected_host_global_id in sorted(
        LOWER_DEPENDENT_OPENING_HOSTS.items()
    ):
        opening = model.by_guid(opening_global_id)
        expected_host = model.by_guid(expected_host_global_id)
        if opening is None or not opening.is_a("IfcOpeningElement"):
            raise RuntimeError(f"missing lower-tile Opening {opening_global_id}")
        if expected_host is None or not expected_host.is_a("IfcCovering"):
            raise RuntimeError(f"missing lower-tile host {expected_host_global_id}")
        voids = list(opening.VoidsElements or ())
        fillings = list(opening.HasFillings or ())
        actual_host_global_ids = sorted(
            relation.RelatingBuildingElement.GlobalId for relation in voids
        )
        placement_relative_to_host = bool(
            opening.ObjectPlacement
            and expected_host.ObjectPlacement
            and opening.ObjectPlacement.PlacementRelTo == expected_host.ObjectPlacement
        )
        record = {
            "opening_global_id": opening_global_id,
            "opening_name": opening.Name,
            "expected_host_global_id": expected_host_global_id,
            "actual_host_global_ids": actual_host_global_ids,
            "void_relationship_count": len(voids),
            "void_relationship_global_ids": sorted(
                relation.GlobalId for relation in voids
            ),
            "filling_relationship_count": len(fillings),
            "placement_relative_to_host": placement_relative_to_host,
        }
        record["valid"] = (
            actual_host_global_ids == [expected_host_global_id]
            and len(fillings) == 0
            and placement_relative_to_host
        )
        records.append(record)
    return records


def tile_level(z_min_mm: float, z_max_mm: float) -> str:
    if -1050.0 <= z_min_mm <= -950.0 and -1050.0 <= z_max_mm <= -950.0:
        return "below_ffl_reference"
    if -100.0 <= z_min_mm <= 100.0 and -100.0 <= z_max_mm <= 100.0:
        return "near_ffl_flooring"
    return "outside_a105_scope"


def material_label(product: ifcopenshell.entity_instance) -> str:
    material = ifcopenshell.util.element.get_material(product)
    return str(material) if material else ""


def shape_z_bounds_mm(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> tuple[float, float]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    values = list(shape.geometry.verts)[2::3]
    if not values:
        raise RuntimeError(f"{product.GlobalId} has no shape vertices")
    return min(values) * 1000.0, max(values) * 1000.0


def shape_y_bounds_mm(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> tuple[float, float]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    values = list(shape.geometry.verts)[1::3]
    if not values:
        raise RuntimeError(f"{product.GlobalId} has no shape vertices")
    return min(values) * 1000.0, max(values) * 1000.0


def drain_gap_inventory(
    model: ifcopenshell.file,
) -> list[dict[str, Any]]:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    records = []
    for name, group in DRAIN_GAP_GROUPS.items():
        low_y_bounds = [
            shape_y_bounds_mm(settings, model.by_guid(global_id))
            for global_id in sorted(group["low_y_tile_ids"])
        ]
        high_y_bounds = [
            shape_y_bounds_mm(settings, model.by_guid(global_id))
            for global_id in sorted(group["high_y_tile_ids"])
        ]
        low_y_edge_mm = max(bounds[1] for bounds in low_y_bounds)
        high_y_edge_mm = min(bounds[0] for bounds in high_y_bounds)
        records.append(
            {
                "group": name,
                "low_y_tile_ids": sorted(group["low_y_tile_ids"]),
                "high_y_tile_ids": sorted(group["high_y_tile_ids"]),
                "swept_opening_ids": sorted(group["swept_opening_ids"]),
                "tessellated_opening_ids": sorted(
                    group["tessellated_opening_ids"]
                ),
                "low_y_edge_mm": low_y_edge_mm,
                "high_y_edge_mm": high_y_edge_mm,
                "centreline_y_mm": (low_y_edge_mm + high_y_edge_mm) / 2.0,
                "gap_mm": high_y_edge_mm - low_y_edge_mm,
            }
        )
    return records


def _body_items(
    opening: ifcopenshell.entity_instance,
) -> list[ifcopenshell.entity_instance]:
    if opening is None or not opening.is_a("IfcOpeningElement"):
        raise RuntimeError("missing drain tile Opening")
    return [
        item
        for representation in opening.Representation.Representations
        if representation.RepresentationIdentifier == "Body"
        for item in representation.Items
    ]


def _grow_swept_opening(
    model: ifcopenshell.file,
    opening: ifcopenshell.entity_instance,
    world_growth_mm: float,
) -> None:
    items = _body_items(opening)
    if len(items) != 1 or not items[0].is_a("IfcExtrudedAreaSolid"):
        raise RuntimeError(
            f"unexpected swept drain Opening body {opening.GlobalId}"
        )
    solid = items[0]
    profile = solid.SweptArea
    curve = getattr(profile, "OuterCurve", None)
    if curve is None or not curve.is_a("IfcIndexedPolyCurve"):
        raise RuntimeError(
            f"unexpected drain Opening profile {opening.GlobalId}"
        )
    if len(model.get_inverse(curve.Points)) != 1:
        raise RuntimeError(
            f"shared drain Opening point list {opening.GlobalId}"
        )
    placement = ifcopenshell.util.placement.get_local_placement(
        opening.ObjectPlacement
    )
    axis_scale = abs(float(placement[1, 1]))
    if axis_scale < 0.99:
        raise RuntimeError(
            f"drain Opening local Y is not world-Y aligned {opening.GlobalId}"
        )
    local_growth = world_growth_mm / axis_scale
    coordinates = [list(point) for point in curve.Points.CoordList]
    maximum_y = max(point[1] for point in coordinates)
    changed = 0
    for point in coordinates:
        if abs(point[1] - maximum_y) <= 1e-9:
            point[1] += local_growth
            changed += 1
    if changed < 2:
        raise RuntimeError(
            f"cannot identify swept cutter edge {opening.GlobalId}"
        )
    curve.Points.CoordList = tuple(tuple(point) for point in coordinates)


def _grow_tessellated_opening(
    model: ifcopenshell.file,
    opening: ifcopenshell.entity_instance,
    world_growth_mm: float,
) -> None:
    items = _body_items(opening)
    if len(items) != 1 or not items[0].is_a("IfcPolygonalFaceSet"):
        raise RuntimeError(
            f"unexpected tessellated drain Opening body {opening.GlobalId}"
        )
    face_set = items[0]
    if len(model.get_inverse(face_set.Coordinates)) != 1:
        raise RuntimeError(
            f"shared drain Opening coordinate list {opening.GlobalId}"
        )
    placement = ifcopenshell.util.placement.get_local_placement(
        opening.ObjectPlacement
    )
    axis_scale = abs(float(placement[1, 1]))
    if axis_scale < 0.99:
        raise RuntimeError(
            f"drain Opening local Y is not world-Y aligned {opening.GlobalId}"
        )
    local_growth = world_growth_mm / axis_scale
    coordinates = [list(point) for point in face_set.Coordinates.CoordList]
    minimum_y = min(point[1] for point in coordinates)
    changed = 0
    for point in coordinates:
        if abs(point[1] - minimum_y) <= 1e-9:
            point[1] -= local_growth
            changed += 1
    if changed != 4:
        raise RuntimeError(
            f"cannot identify tessellated cutter edge {opening.GlobalId}"
        )
    face_set.Coordinates.CoordList = tuple(
        tuple(point) for point in coordinates
    )

    box_items = [
        item
        for representation in opening.Representation.Representations
        if representation.RepresentationIdentifier == "Box"
        for item in representation.Items
    ]
    if len(box_items) != 1 or not box_items[0].is_a("IfcBoundingBox"):
        raise RuntimeError(
            f"unexpected drain Opening Box {opening.GlobalId}"
        )
    box = box_items[0]
    if len(model.get_inverse(box.Corner)) != 1:
        raise RuntimeError(
            f"shared drain Opening Box corner {opening.GlobalId}"
        )
    corner = list(box.Corner.Coordinates)
    corner[1] -= local_growth
    box.Corner.Coordinates = tuple(corner)
    box.YDim += local_growth


def apply_drain_gap(
    model: ifcopenshell.file,
    target_gap_mm: float,
    result_tolerance_mm: float = 0.0001,
) -> list[dict[str, Any]]:
    before = {record["group"]: record for record in drain_gap_inventory(model)}
    iterations = {name: 0 for name in DRAIN_GAP_GROUPS}
    for _ in range(4):
        current = {
            record["group"]: record for record in drain_gap_inventory(model)
        }
        unresolved = False
        for name, group in DRAIN_GAP_GROUPS.items():
            residual_mm = target_gap_mm - current[name]["gap_mm"]
            if abs(residual_mm) <= result_tolerance_mm:
                continue
            unresolved = True
            half_growth_mm = residual_mm / 2.0
            for global_id in sorted(group["swept_opening_ids"]):
                _grow_swept_opening(
                    model, model.by_guid(global_id), half_growth_mm
                )
            for global_id in sorted(group["tessellated_opening_ids"]):
                _grow_tessellated_opening(
                    model, model.by_guid(global_id), half_growth_mm
                )
            iterations[name] += 1
        if not unresolved:
            break
    after = {record["group"]: record for record in drain_gap_inventory(model)}
    return [
        {
            "group": name,
            "target_gap_mm": target_gap_mm,
            "source_gap_mm": before[name]["gap_mm"],
            "candidate_gap_mm": after[name]["gap_mm"],
            "source_centreline_y_mm": before[name]["centreline_y_mm"],
            "candidate_centreline_y_mm": after[name]["centreline_y_mm"],
            "iterations": iterations[name],
            "swept_opening_ids": sorted(group["swept_opening_ids"]),
            "tessellated_opening_ids": sorted(
                group["tessellated_opening_ids"]
            ),
        }
        for name, group in DRAIN_GAP_GROUPS.items()
    ]


def mapped_representation_ids(product: ifcopenshell.entity_instance) -> list[int]:
    result: list[int] = []
    representations = (
        product.Representation.Representations if product.Representation else ()
    )
    for representation in representations:
        for item in representation.Items:
            if item.is_a("IfcMappedItem"):
                result.append(item.MappingSource.id())
    return result


def tile_inventory(model: ifcopenshell.file) -> list[dict[str, Any]]:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    result: list[dict[str, Any]] = []
    for covering in model.by_type("IfcCovering"):
        if "TerrazzoMosaicTile" not in material_label(covering):
            continue
        z_min_mm, z_max_mm = shape_z_bounds_mm(settings, covering)
        result.append(
            {
                "global_id": covering.GlobalId,
                "step_id": covering.id(),
                "name": covering.Name,
                "predefined_type": covering.PredefinedType,
                "z_min_mm": z_min_mm,
                "z_max_mm": z_max_mm,
                "level": tile_level(z_min_mm, z_max_mm),
                "representation_map_ids": mapped_representation_ids(covering),
            }
        )
    return sorted(result, key=lambda record: record["global_id"])


def apply_candidate(
    model: ifcopenshell.file,
    current_global_ids: set[str],
    lower_action: str,
    drain_gap_mm: float,
) -> list[dict[str, Any]]:
    for global_id in sorted(current_global_ids):
        covering = model.by_guid(global_id)
        if covering is None or not covering.is_a("IfcCovering"):
            raise RuntimeError(f"missing current tile {global_id}")
        covering.PredefinedType = "FLOORING"
    if lower_action == "delete":
        for global_id in sorted(LOWER_TILE_GLOBAL_IDS):
            covering = model.by_guid(global_id)
            if covering is None or not covering.is_a("IfcCovering"):
                raise RuntimeError(f"missing lower tile {global_id}")
            ifcopenshell.api.root.remove_product(model, product=covering)
    return apply_drain_gap(model, drain_gap_mm)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--lower-action", choices=("keep", "delete"), required=True)
    parser.add_argument("--drain-gap-mm", type=float, required=True)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= 1.0:
        raise SystemExit("--tolerance-mm must be greater than zero and at most 1 mm")
    if not 0.0 < args.drain_gap_mm <= 100.0:
        raise SystemExit("--drain-gap-mm must be greater than zero and at most 100 mm")
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    source_inventory = tile_inventory(source)
    source_current_ids = {
        record["global_id"]
        for record in source_inventory
        if record["level"] == "near_ffl_flooring"
    }
    source_lower_ids = {
        record["global_id"]
        for record in source_inventory
        if record["level"] == "below_ffl_reference"
    }
    source_map_ids = {
        map_id
        for record in source_inventory
        for map_id in record["representation_map_ids"]
    }
    if len(source_inventory) != 22:
        raise RuntimeError(f"expected 22 TerrazzoMosaicTile coverings, found {len(source_inventory)}")
    if len(source_current_ids) != 18:
        raise RuntimeError(f"expected 18 near-FFL tiles, found {len(source_current_ids)}")
    if source_lower_ids != LOWER_TILE_GLOBAL_IDS:
        raise RuntimeError(
            f"lower tile identity drift: {sorted(source_lower_ids)}"
        )
    if len(source_map_ids) != 1:
        raise RuntimeError(f"expected one shared representation map, found {source_map_ids}")
    source_lower_dependencies = lower_dependency_inventory(source)
    if not all(record["valid"] for record in source_lower_dependencies):
        raise RuntimeError(
            "lower-tile dependent Opening relationship drift: "
            f"{source_lower_dependencies}"
        )
    expected_removed_ids = expected_removed_global_ids(args.lower_action)
    expected_removed_root_ids = expected_removed_root_global_ids(
        args.lower_action, source_lower_dependencies
    )

    candidate = ifcopenshell.open(source_path)
    drain_gap_adjustments = apply_candidate(
        candidate,
        source_current_ids,
        args.lower_action,
        args.drain_gap_mm,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    candidate_inventory = tile_inventory(candidate)
    candidate_current_ids = {
        record["global_id"]
        for record in candidate_inventory
        if record["level"] == "near_ffl_flooring"
    }
    candidate_lower_ids = {
        record["global_id"]
        for record in candidate_inventory
        if record["level"] == "below_ffl_reference"
    }
    current_geometry = geometry_difference_audit(
        candidate,
        source,
        str(source_path),
        tolerance_mm=args.tolerance_mm,
        classes=[],
        global_ids=sorted(source_current_ids),
    )
    all_geometry = all_product_geometry_difference(
        candidate, source, args.tolerance_mm
    )
    removed_records = [
        record
        for record in all_geometry["records"]
        if record["status"] == "removed"
    ]
    unexpected_geometry_records = [
        record
        for record in all_geometry["records"]
        if not record["within_tolerance"]
        and record.get("global_id") not in expected_removed_ids
    ]
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    expected_root_ids = source_root_ids - expected_removed_root_ids
    removed_root_ids = source_root_ids - candidate_root_ids
    changed_retained_root_ids = {
        global_id
        for global_id in source_root_ids & candidate_root_ids
        if str(source.by_guid(global_id)) != str(candidate.by_guid(global_id))
    }
    expected_changed_retained_root_ids = set(source_current_ids)
    if args.lower_action == "delete":
        for global_id in expected_removed_ids:
            product = source.by_guid(global_id)
            for inverse in source.get_inverse(product):
                inverse_global_id = getattr(inverse, "GlobalId", None)
                if (
                    inverse_global_id
                    and inverse_global_id in candidate_root_ids
                ):
                    expected_changed_retained_root_ids.add(inverse_global_id)
    unexpected_changed_retained_root_ids = (
        changed_retained_root_ids - expected_changed_retained_root_ids
    )
    changed_entity_ids = entity_changes(source, candidate)
    current_step_ids = sorted(
        source.by_guid(global_id).id() for global_id in source_current_ids
    )
    drain_gap_records = drain_gap_inventory(candidate)
    drain_gap_residual_mm = max(
        abs(record["gap_mm"] - args.drain_gap_mm)
        for record in drain_gap_records
    )
    drain_gap_centreline_delta_mm = max(
        abs(
            record["candidate_centreline_y_mm"]
            - record["source_centreline_y_mm"]
        )
        for record in drain_gap_adjustments
    )
    maximum_current_tile_hausdorff_mm = max(
        (
            record.get("world_vertex_hausdorff_mm") or 0.0
            for record in current_geometry["records"]
        ),
        default=0.0,
    )
    candidate_map_ids = {
        map_id
        for record in candidate_inventory
        if record["global_id"] in candidate_current_ids
        for map_id in record["representation_map_ids"]
    }
    gates = {
        "source_tiles": len(source_inventory),
        "source_current_tiles": len(source_current_ids),
        "source_lower_tiles": len(source_lower_ids),
        "candidate_current_tiles": len(candidate_current_ids),
        "candidate_lower_tiles": len(candidate_lower_ids),
        "current_tiles_are_flooring": all(
            candidate.by_guid(global_id).PredefinedType == "FLOORING"
            for global_id in candidate_current_ids
        ),
        "current_tile_geometry_over_tolerance": current_geometry["over_tolerance"],
        "maximum_current_tile_hausdorff_mm": maximum_current_tile_hausdorff_mm,
        "drain_gap_target_mm": args.drain_gap_mm,
        "drain_gap_residual_mm": drain_gap_residual_mm,
        "drain_gap_centreline_delta_mm": drain_gap_centreline_delta_mm,
        "unexpected_product_geometry_records": len(unexpected_geometry_records),
        "removed_product_global_ids": sorted(
            record["global_id"] for record in removed_records
        ),
        "expected_removed_global_ids": sorted(expected_removed_ids),
        "removed_root_global_ids": sorted(removed_root_ids),
        "expected_removed_root_global_ids": sorted(expected_removed_root_ids),
        "changed_retained_root_global_ids": sorted(changed_retained_root_ids),
        "expected_changed_retained_root_global_ids": sorted(
            expected_changed_retained_root_ids
        ),
        "unexpected_changed_retained_root_global_ids": sorted(
            unexpected_changed_retained_root_ids
        ),
        "lower_dependent_openings_valid": all(
            record["valid"] for record in source_lower_dependencies
        ),
        "shared_representation_map_preserved": candidate_map_ids == source_map_ids,
        "schema_equal": source.schema == candidate.schema,
        "root_global_ids_match_expected": candidate_root_ids == expected_root_ids,
        "semantic_only_entity_changes": None,
    }
    lower_gate = (
        candidate_lower_ids == LOWER_TILE_GLOBAL_IDS
        and not gates["removed_product_global_ids"]
        if args.lower_action == "keep"
        else not candidate_lower_ids
        and set(gates["removed_product_global_ids"]) == expected_removed_ids
    )
    gates["lower_action_matches"] = lower_gate
    gates["pass"] = (
        gates["source_tiles"] == 22
        and gates["source_current_tiles"] == 18
        and gates["source_lower_tiles"] == 4
        and gates["candidate_current_tiles"] == 18
        and gates["current_tiles_are_flooring"]
        and gates["current_tile_geometry_over_tolerance"] == 0
        and gates["maximum_current_tile_hausdorff_mm"] <= args.tolerance_mm
        and gates["drain_gap_residual_mm"] <= 0.0001
        and gates["drain_gap_centreline_delta_mm"] <= 0.0001
        and gates["unexpected_product_geometry_records"] == 0
        and gates["shared_representation_map_preserved"]
        and gates["lower_dependent_openings_valid"]
        and gates["schema_equal"]
        and gates["root_global_ids_match_expected"]
        and not gates["unexpected_changed_retained_root_global_ids"]
        and lower_gate
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a105-flooring-candidate",
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
        "lower_action": args.lower_action,
        "tolerance_mm": args.tolerance_mm,
        "drain_gap_mm": args.drain_gap_mm,
        "drain_gap_adjustments": drain_gap_adjustments,
        "drain_gap_inventory": drain_gap_records,
        "source_inventory": source_inventory,
        "source_lower_dependencies": source_lower_dependencies,
        "candidate_inventory": candidate_inventory,
        "current_tile_geometry": current_geometry,
        "all_product_geometry": all_geometry,
        "changed_entity_ids": changed_entity_ids,
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
