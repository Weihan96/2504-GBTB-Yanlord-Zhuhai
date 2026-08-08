#!/usr/bin/env python3
"""Generate a read-only ELEC existing-condition candidate from the formal IFC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import numpy as np
from ifcopenshell.util.element import get_container, get_psets, get_type
from ifcopenshell.util.placement import get_local_placement
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union


EXPECTED_SOURCE_SHA256 = "c7295688003f3f36775a25f6adc2c9878e203c52c980f8e46faed66a3537c4a8"
LIGHT_TYPE_GLOBAL_ID = "26DmeZ15D7ZeHYavba3FdK"
EXPECTED_PROXY_IDS = {
    "1faflkXXH6M9cnYPE9Liir",
    "1WMWifzpjFOxatyFb9stFb",
    "3ufPf53o1CbQlrIxF$Oogk",
    "2qvK0LPkPAr8uEI8Oa_U1h",
    "3Irt7GGfb5MP7Q6g$zFVaX",
    "2xvG7uDof7OwbTH_w0edo0",
    "2SdKxc_q96SBGrAUSYsWM9",
    "1SNXhKeZb7r9Y0MmzcQnc3",
    "2d2Vw3ZSn0exH1seMBMiVf",
}
CONTROLLED_SOCKET_IDS = {
    "0laejMoxn8Lu_X3FZaCXmi",
    "27MTenki57DQsfMryX_1U0",
    "2OOjqQDMHDjRcXQCniWXnp",
    "3KXtmVvejA78j_iAS$kydj",
}
NETWORK_TYPE_CATEGORIES = {"LAN", "LANFLUSH", "LANSOCKET"}
CONTROL_TYPE_CATEGORIES = {"SWITCHPANEL", "CONTROLPANEL"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument(
        "--review",
        type=Path,
        default=Path("pipeline/decisions/elec-existing-review.csv"),
    )
    parser.add_argument("--output", type=Path, default=Path("build/elec/elec-existing-candidate.json"))
    parser.add_argument("--expected-sha256", default=EXPECTED_SOURCE_SHA256)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def shape_data(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> tuple[np.ndarray, np.ndarray]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3))
    faces = np.asarray(shape.geometry.faces, dtype=int).reshape((-1, 3))
    if not len(vertices):
        raise RuntimeError(f"{product.GlobalId} has no renderable world geometry")
    return vertices, faces


def bbox_mm(vertices_m: np.ndarray) -> dict[str, list[float]]:
    vertices = vertices_m * 1000.0
    minimum = vertices.min(axis=0)
    maximum = vertices.max(axis=0)
    return {
        "min_mm": minimum.tolist(),
        "max_mm": maximum.tolist(),
        "dimensions_mm": (maximum - minimum).tolist(),
        "centre_mm": ((minimum + maximum) / 2.0).tolist(),
    }


def placement_record(product: ifcopenshell.entity_instance) -> dict[str, Any]:
    if product.ObjectPlacement is None:
        raise RuntimeError(f"{product.GlobalId} has no ObjectPlacement")
    origin = np.asarray(get_local_placement(product.ObjectPlacement)[:3, 3], dtype=float)
    residual = np.abs(origin - np.round(origin))
    return {
        "origin_mm": origin.tolist(),
        "maximum_integer_residual_mm": float(residual.max()),
    }


def read_review(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        records = list(csv.DictReader(handle))
    ids = [row["review_id"] for row in records]
    if len(ids) != len(set(ids)):
        raise RuntimeError("ELEC review IDs are not unique")
    return records


def review_index(records: list[dict[str, str]], scope: str) -> dict[str, dict[str, str]]:
    return {row["global_id"]: row for row in records if row["scope"] == scope}


def space_footprints(
    model: ifcopenshell.file,
    settings: ifcopenshell.geom.settings,
) -> list[tuple[ifcopenshell.entity_instance, Any]]:
    records: list[tuple[ifcopenshell.entity_instance, Any]] = []
    for space in model.by_type("IfcSpace"):
        vertices, faces = shape_data(settings, space)
        polygons = []
        for face in faces:
            polygon = Polygon(vertices[face, :2])
            if polygon.area > 1e-10:
                polygons.append(polygon)
        if not polygons:
            raise RuntimeError(f"Space {space.GlobalId} has no plan footprint")
        records.append((space, unary_union(polygons).buffer(1e-8)))
    return records


def candidate_space(
    centre_mm: list[float],
    footprints: list[tuple[ifcopenshell.entity_instance, Any]],
) -> dict[str, Any]:
    point = Point(float(centre_mm[0]) / 1000.0, float(centre_mm[1]) / 1000.0)
    matches = [space for space, footprint in footprints if footprint.covers(point)]
    if len(matches) != 1:
        raise RuntimeError(
            f"electrical point at {centre_mm[:2]} has {len(matches)} candidate Spaces"
        )
    space = matches[0]
    return {
        "global_id": space.GlobalId,
        "name": str(space.Name or ""),
        "long_name": str(space.LongName or ""),
        "basis": "world Body bbox centre contained by current IfcSpace plan footprint",
        "confidence": 0.95,
        "review_required": "yes",
    }


def instance_record(
    product: ifcopenshell.entity_instance,
    settings: ifcopenshell.geom.settings,
    footprints: list[tuple[ifcopenshell.entity_instance, Any]],
) -> dict[str, Any]:
    vertices, _ = shape_data(settings, product)
    bbox = bbox_mm(vertices)
    assigned_type = get_type(product)
    return {
        "global_id": product.GlobalId,
        "ifc_class": product.is_a(),
        "name": str(product.Name or ""),
        "object_type": str(product.ObjectType or "") if hasattr(product, "ObjectType") else "",
        "container": str(getattr(get_container(product), "Name", "") or ""),
        "assigned_type": {
            "global_id": str(getattr(assigned_type, "GlobalId", "") or ""),
            "name": str(getattr(assigned_type, "Name", "") or ""),
            "element_type": str(getattr(assigned_type, "ElementType", "") or ""),
            "predefined_type": str(getattr(assigned_type, "PredefinedType", "") or ""),
        },
        "bbox": bbox,
        "placement": placement_record(product),
        "candidate_space": candidate_space(bbox["centre_mm"], footprints),
        "has_psets": bool(get_psets(product)),
        "formal_ifc_write_allowed": "no",
    }


def plan_sort_key(record: dict[str, Any]) -> tuple[float, float, str]:
    centre = record["bbox"]["centre_mm"]
    return (-round(float(centre[1]), 3), round(float(centre[0]), 3), record["global_id"])


def assign_candidate_ids(records: list[dict[str, Any]], prefix: str) -> None:
    records.sort(key=plan_sort_key)
    for index, record in enumerate(records, 1):
        record["candidate_id"] = f"{prefix}{index:03d}"


def duplicate_centres(records: list[dict[str, Any]], precision_mm: int = 3) -> list[dict[str, Any]]:
    groups: defaultdict[tuple[float, float, float], list[str]] = defaultdict(list)
    for record in records:
        centre = record["bbox"]["centre_mm"]
        key = tuple(round(float(value), precision_mm) for value in centre)
        groups[key].append(record["global_id"])
    return [
        {"centre_mm": list(centre), "global_ids": ids}
        for centre, ids in sorted(groups.items())
        if len(ids) > 1
    ]


def type_record(type_object: ifcopenshell.entity_instance, instance_count: int) -> dict[str, Any]:
    return {
        "global_id": type_object.GlobalId,
        "name": str(type_object.Name or ""),
        "element_type": str(type_object.ElementType or ""),
        "predefined_type": str(type_object.PredefinedType or ""),
        "instance_count": instance_count,
        "is_instance_point": False,
    }


def opening_record(opening: ifcopenshell.entity_instance) -> dict[str, Any]:
    hosts = [relation.RelatingBuildingElement for relation in opening.VoidsElements]
    return {
        "global_id": opening.GlobalId,
        "name": str(opening.Name or ""),
        "hosts": [
            {"ifc_class": host.is_a(), "global_id": host.GlobalId, "name": str(host.Name or "")}
            for host in hosts
        ],
        "formal_ifc_write_allowed": "no",
    }


def main() -> int:
    args = parse_args()
    source_sha = sha256(args.input)
    if source_sha != args.expected_sha256:
        raise RuntimeError(
            f"formal IFC SHA drift: {source_sha} != {args.expected_sha256}"
        )

    model = ifcopenshell.open(args.input)
    if model.schema != "IFC4":
        raise RuntimeError(f"unexpected schema {model.schema}")
    reviews = read_review(args.review)
    proxy_reviews = review_index(reviews, "proxy_handoff")
    socket_exception_reviews = review_index(reviews, "socket_exception")
    if set(proxy_reviews) != EXPECTED_PROXY_IDS:
        raise RuntimeError("proxy handoff review set drift")
    if set(socket_exception_reviews) != CONTROLLED_SOCKET_IDS:
        raise RuntimeError("controlled socket review set drift")

    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    footprints = space_footprints(model, settings)

    lights = [instance_record(product, settings, footprints) for product in model.by_type("IfcLightFixture")]
    for record in lights:
        record["review_basis"] = "existing RA.LP world installation centre; room is a geometric candidate"
        record["review_required"] = "yes"
    assign_candidate_ids(lights, "L")

    appliances = model.by_type("IfcElectricAppliance")
    sockets = [
        instance_record(product, settings, footprints)
        for product in appliances
        if str(getattr(get_type(product), "ElementType", "") or "") == "SOCKET"
    ]
    for record in sockets:
        record["controlled_installation_exception"] = record["global_id"] in CONTROLLED_SOCKET_IDS
        record["review_id"] = socket_exception_reviews.get(record["global_id"], {}).get("review_id", "")
    assign_candidate_ids(sockets, "P")

    equipment = [
        instance_record(product, settings, footprints)
        for product in appliances
        if str(getattr(get_type(product), "ElementType", "") or "") != "SOCKET"
    ]
    for record in equipment:
        record["power_and_interface_status"] = "missing"
        record["review_required"] = "yes"
    assign_candidate_ids(equipment, "EQ")

    proxies = []
    for global_id in EXPECTED_PROXY_IDS:
        product = model.by_guid(global_id)
        if product is None or not product.is_a("IfcBuildingElementProxy"):
            raise RuntimeError(f"ELEC proxy handoff missing or class drift: {global_id}")
        record = instance_record(product, settings, footprints)
        review = proxy_reviews[global_id]
        record.update(
            {
                "candidate_role": review["candidate_role"],
                "review_id": review["review_id"],
                "review_required": review["review_required"],
                "review_status": review["status"],
                "protected_action": review["protected_action"],
                "required_confirmation": review["required_confirmation"],
                "related_openings": [
                    relation.RelatedOpeningElement.GlobalId for relation in product.HasOpenings
                ],
            }
        )
        proxies.append(record)
    assign_candidate_ids(proxies, "PX")

    openings = []
    for product in proxies:
        for global_id in product["related_openings"]:
            openings.append(opening_record(model.by_guid(global_id)))
    openings.sort(key=lambda record: record["global_id"])

    type_use_counts = Counter(
        get_type(product).GlobalId for product in appliances if get_type(product) is not None
    )
    appliance_types = model.by_type("IfcElectricApplianceType")
    used_types = [type_record(item, type_use_counts[item.GlobalId]) for item in appliance_types if type_use_counts[item.GlobalId]]
    unused_types = [type_record(item, 0) for item in appliance_types if not type_use_counts[item.GlobalId]]
    used_types.sort(key=lambda record: (record["element_type"], record["name"], record["global_id"]))
    unused_types.sort(key=lambda record: (record["element_type"], record["name"], record["global_id"]))
    control_types = [record for record in unused_types if record["element_type"] in CONTROL_TYPE_CATEGORIES]
    network_types = [record for record in unused_types if record["element_type"] in NETWORK_TYPE_CATEGORIES]

    root_ids = [root.GlobalId for root in model.by_type("IfcRoot")]
    light_type_ids = {record["assigned_type"]["global_id"] for record in lights}
    socket_ids = {record["global_id"] for record in sockets}
    proxy_ids = {record["global_id"] for record in proxies}
    light_room_counts = Counter(record["candidate_space"]["long_name"] for record in lights)
    light_residuals = [record["placement"]["maximum_integer_residual_mm"] for record in lights]
    socket_residuals = [record["placement"]["maximum_integer_residual_mm"] for record in sockets]
    proxy_residuals = [record["placement"]["maximum_integer_residual_mm"] for record in proxies]

    topology = {
        "ifc_systems": len(model.by_type("IfcSystem")),
        "ifc_distribution_systems": len(model.by_type("IfcDistributionSystem")),
        "distribution_ports": len(model.by_type("IfcDistributionPort")),
        "port_connections": len(model.by_type("IfcRelConnectsPorts")),
        "port_to_element_connections": len(model.by_type("IfcRelConnectsPortToElement")),
        "control_assignments": len(model.by_type("IfcRelAssignsToControl")),
    }
    e302_instances = len(model.by_type("IfcSwitchingDevice"))
    e304_instances = len(model.by_type("IfcCommunicationsAppliance"))

    gates = {
        "root_global_ids_unique": len(root_ids) == len(set(root_ids)),
        "light_count": len(lights),
        "light_type_ids": sorted(light_type_ids),
        "light_geometry_and_type_complete": all(record["assigned_type"]["global_id"] for record in lights),
        "light_space_candidates_single": all(record["candidate_space"]["global_id"] for record in lights),
        "light_duplicate_centres": duplicate_centres(lights),
        "light_origins_within_0_1_mm": sum(value <= args.tolerance_mm + 1e-9 for value in light_residuals),
        "light_maximum_integer_residual_mm": max(light_residuals),
        "socket_count": len(sockets),
        "socket_duplicate_centres": duplicate_centres(sockets),
        "socket_origins_within_0_1_mm": sum(value <= args.tolerance_mm + 1e-9 for value in socket_residuals),
        "controlled_socket_ids": sorted(record["global_id"] for record in sockets if record["controlled_installation_exception"]),
        "typed_equipment_count": len(equipment),
        "typed_equipment_integer_origins": sum(
            record["placement"]["maximum_integer_residual_mm"] <= 1e-9 for record in equipment
        ),
        "proxy_handoff_count": len(proxies),
        "proxy_handoff_ids": sorted(proxy_ids),
        "proxy_handoff_over_0_1_mm": sum(value > args.tolerance_mm + 1e-9 for value in proxy_residuals),
        "proxy_handoff_maximum_integer_residual_mm": max(proxy_residuals),
        "related_opening_count": len(openings),
        "related_openings_have_one_host": all(len(record["hosts"]) == 1 for record in openings),
        "appliance_type_count": len(appliance_types),
        "used_appliance_type_count": len(used_types),
        "unused_appliance_type_count": len(unused_types),
        "type_library_separated_from_instances": not ({record["global_id"] for record in used_types} & {record["global_id"] for record in unused_types}),
        "e302_switch_instances": e302_instances,
        "e302_unused_control_type_definitions": len(control_types),
        "e304_network_instances": e304_instances,
        "e304_unused_network_type_definitions": len(network_types),
        "topology": topology,
        "automatic_ifc_write_allowed": False,
    }
    gates["candidate_pass"] = all(
        [
            gates["root_global_ids_unique"],
            gates["light_count"] == 79,
            gates["light_type_ids"] == [LIGHT_TYPE_GLOBAL_ID],
            gates["light_geometry_and_type_complete"],
            gates["light_space_candidates_single"],
            gates["light_duplicate_centres"] == [],
            gates["light_origins_within_0_1_mm"] == 79,
            gates["socket_count"] == 11,
            gates["socket_duplicate_centres"] == [],
            gates["socket_origins_within_0_1_mm"] == 7,
            set(gates["controlled_socket_ids"]) == CONTROLLED_SOCKET_IDS,
            gates["typed_equipment_count"] == 8,
            gates["typed_equipment_integer_origins"] == 8,
            gates["proxy_handoff_count"] == 9,
            set(gates["proxy_handoff_ids"]) == EXPECTED_PROXY_IDS,
            gates["proxy_handoff_over_0_1_mm"] == 9,
            gates["related_opening_count"] == 4,
            gates["related_openings_have_one_host"],
            gates["appliance_type_count"] == 42,
            gates["used_appliance_type_count"] == 9,
            gates["unused_appliance_type_count"] == 33,
            gates["type_library_separated_from_instances"],
            gates["e302_switch_instances"] == 0,
            gates["e302_unused_control_type_definitions"] == 13,
            gates["e304_network_instances"] == 0,
            gates["e304_unused_network_type_definitions"] == 7,
            all(value == 0 for value in topology.values()),
        ]
    )
    gates["construction_release_ready"] = False
    gates["release_blocks"] = [
        "E-301 fixture product, mounting method, circuit, and control group are not confirmed",
        "E-302 has no switch instances or control relations",
        "E-303 has no equipment power, circuit, waterproofing, or interface properties",
        "E-304 has no network instances, ports, systems, or topology",
        "nine ELEC proxy handoff objects have no confirmed installation anchors",
    ]

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_existing_condition_candidate",
        "source": {"path": str(args.input.resolve()), "sha256": source_sha, "schema": model.schema},
        "tolerance_mm": args.tolerance_mm,
        "sheets": {
            "E-301": {
                "status": "existing_points_candidate",
                "lights": lights,
                "room_counts": dict(sorted(light_room_counts.items())),
                "missing": ["fixture product", "mounting method", "circuit", "control group"],
            },
            "E-302": {
                "status": "pending_layout_no_instances",
                "instances": [],
                "available_uninstantiated_type_definitions": control_types,
                "confirmed_requirements_only": [
                    "entry lighting master control",
                    "entry-to-master-bedroom corridor/living dual control",
                    "text or icon label for every key",
                ],
            },
            "E-303": {
                "status": "existing_points_and_footprints_candidate",
                "sockets": sockets,
                "typed_equipment": equipment,
                "proxy_handoffs": proxies,
                "related_openings": openings,
                "missing": ["power", "voltage", "circuit", "waterproofing", "interface"],
            },
            "E-304": {
                "status": "pending_layout_no_instances_or_topology",
                "instances": [],
                "available_uninstantiated_type_definitions": network_types,
                "topology": topology,
            },
        },
        "type_library": {"used": used_types, "unused": unused_types},
        "review_register": {"path": str(args.review), "records": len(reviews)},
        "gates": gates,
    }
    if not gates["candidate_pass"]:
        raise RuntimeError("ELEC existing-condition candidate gates failed")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "gates": gates}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
