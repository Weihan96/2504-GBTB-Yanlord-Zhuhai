#!/usr/bin/env python3
"""Audit observable M-401/RCP1 objects without writing or inferring IFC data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element


AC_TYPE_NAMES = {"AC1180", "AC700", "AC700F"}
HVAC_FLOW_NAMES = {"Liquid", "Drain Pipe", "Liquid Living Room", "Gas", "Gas Living Room"}
HVAC_FLOW_LEGACY_ROLES = {
    "Liquid": {
        "role": "legacy_base_refrigerant_liquid_geometry",
        "basis": "用户确认该对象在旧方案中表达冷媒液管，但紫色几何不是装修后最终管线；旧 lowpoly 与正式 IFC 的 3 个独立分支及组件包围尺寸在 0.0153 mm 内一致.",
        "missing": "diameter, insulation, system membership, ports, connections, route ownership and installer clearances",
        "stop": "Treat this as a legacy design base only; redesign the remodel route before claiming a connected or installable system.",
    },
    "Liquid Living Room": {
        "role": "legacy_base_refrigerant_liquid_geometry",
        "basis": "用户确认该对象在旧方案中表达冷媒液管，但不是装修后最终管线；旧 lowpoly 与正式 IFC 的单一分支组件包围尺寸在 0.001 mm 内一致.",
        "missing": "diameter, insulation, system membership, ports, connections, route ownership and installer clearances",
        "stop": "Treat this as a legacy design base only; redesign the remodel route before claiming a connected or installable system.",
    },
    "Gas": {
        "role": "legacy_base_refrigerant_gas_geometry",
        "basis": "用户确认该对象在旧方案中表达冷媒气管，但紫色几何不是装修后最终管线；旧 lowpoly 与正式 IFC 的 3 个独立分支及组件包围尺寸在 0.0077 mm 内一致.",
        "missing": "diameter, insulation, system membership, ports, connections, route ownership and installer clearances",
        "stop": "Treat this as a legacy design base only; redesign the remodel route before claiming a connected or installable system.",
    },
    "Gas Living Room": {
        "role": "legacy_base_refrigerant_gas_geometry",
        "basis": "用户确认该对象在旧方案中表达冷媒气管，但不是装修后最终管线；旧 lowpoly 与正式 IFC 的单一分支组件包围尺寸在 0.001 mm 内一致.",
        "missing": "diameter, insulation, system membership, ports, connections, route ownership and installer clearances",
        "stop": "Treat this as a legacy design base only; redesign the remodel route before claiming a connected or installable system.",
    },
    "Drain Pipe": {
        "role": "legacy_base_condensate_geometry",
        "basis": "用户确认该对象在旧方案中表达冷凝水管，但紫色几何不是装修后最终管线；旧 lowpoly 与正式 IFC 的 4 个独立分支及组件包围尺寸在 0.0008 mm 内一致.",
        "missing": "diameter, slope, discharge destination, system membership, ports, connections and access",
        "stop": "Treat this as a legacy design base only; redesign the remodel route before claiming drainage continuity or an installable route.",
    },
}
CHECK_VALVE_ID = "1faflkXXH6M9cnYPE9Liir"
DIFFUSER_PROXY_IDS = {"16Ey9Flj9BK9VRun$ozzjH", "3Bv_Kl3jDC5RvaUMdyge1U"}
DIFFUSER_PROXY_DECISIONS = {
    "16Ey9Flj9BK9VRun$ozzjH": {
        "basis": "用户在 Blender 集中审图中确认该 Proxy 是转角空调风口；正式 IFC 仍无 IfcAirTerminal 类型、系统、端口或风量.",
        "confidence": 1.0,
        "status": "IDENTITY_CONFIRMED_REMODEL_DESIGN_PENDING",
        "missing": "supply/return/exhaust role, grille selection, airflow, connection, installation and access",
        "stop": "Retain the confirmed outlet identity, but do not convert the proxy to a terminal or publish it until the remodel air-side design is complete.",
    },
    "3Bv_Kl3jDC5RvaUMdyge1U": {
        "basis": "用户确认该 Proxy 是旧方案空调出风口，同时明确该方案有缺陷；正式 IFC 仍无 IfcAirTerminal、系统、端口或风量.",
        "confidence": 1.0,
        "status": "LEGACY_SCHEME_CONFIRMED_REDESIGN_REQUIRED",
        "missing": "remodel supply-air design, airflow, grille size, centreline, connection, installation and access",
        "stop": "Retain only as a legacy design base; do not publish the current proxy as the remodel supply outlet.",
    },
}
CEILING_CONTEXT_IDS = {
    "0w_j3LiLHFX8uK$fI3H5IP",
    "1_eTV3TaD0ewHXs5A2gCk0",
    "0goB22sEf0W9caxUqAkTwJ",
    "3$lrsbcgT3Jg4KFoQ9xsL_",
    "2HNrakrjbCTuX3CgB2ORBO",
    "24GcPrcEPB0QT3HxbcvPtq",
    "3Kxq1aQl9EyfzEhA_7IoLA",
    "2B9gGJ5197jvvwqUj7h8_U",
    "1wtOLxfnP3ZfqIIkLZaTFU",
    "30YXVJrQjBEgTfH5xSTneR",
    "26M2TYddDFzeD98_XKwBgZ",
    "0JcqcpduD9YRHOQJ2Q92gO",
    "2PT3O6ohvAxPR570xBmdX6",
    "0jI$MQCH9A_w5aml5rEXbg",
    "3BNIc86L19ZBCt40wmF$KZ",
}

CSV_FIELDS = [
    "record_kind",
    "queue_id",
    "global_id",
    "ifc_class",
    "name",
    "description",
    "object_type",
    "predefined_type",
    "type_global_id",
    "type_name",
    "type_description",
    "type_occurrence_count",
    "container",
    "observable_role",
    "bbox_min_mm",
    "bbox_max_mm",
    "dimensions_mm",
    "basis",
    "confidence",
    "human_review_required",
    "review_status",
    "missing_or_unverified",
    "stop_condition",
    "source_ifc_sha256",
]

MISSING_INPUTS = [
    {
        "queue_id": "M401-MISS-001",
        "observable_role": "typed_air_terminal_and_system_connectivity_missing",
        "basis": "No IfcAirTerminal, IfcAirTerminalType, IfcDuctSegment, IfcFan, IfcDistributionPort, IfcSystem or IfcRelConnectsPorts is present for the two diffuser proxies.",
        "missing_or_unverified": "air-terminal identity, supply/return/exhaust role, system membership, ports and connections",
        "stop_condition": "Do not infer terminal role, duct route, connection or airflow from proxy geometry.",
    },
    {
        "queue_id": "M401-MISS-002",
        "observable_role": "kitchen_gas_alarm_missing",
        "basis": "No observable alarm or sensor instance/type matches the confirmed kitchen gas-alarm requirement. A Xiaomi gas alarm and recessed mount are an aesthetic candidate only; the user states the gas authority may require another model.",
        "missing_or_unverified": "gas-authority-approved device selection, mounting location, power/control and required interlock",
        "stop_condition": "Keep the kitchen gas alarm required, but do not lock the candidate 94 mm opening, 120 mm counterbore, 43 mm embed depth or mark it coordinated until the gas authority and final equipment requirements are confirmed.",
    },
    {
        "queue_id": "M401-MISS-003",
        "observable_role": "bedroom_living_fire_sensor_missing",
        "basis": "No observable fire/smoke detector or sensor instance/type exists in the formal IFC. Xiaomi smoke alarms with recessed mounts are accepted as a visual candidate for the master bedroom, guest bedroom and living room.",
        "missing_or_unverified": "approved device type, final quantity and coverage, exact location, power/communication and ceiling clearance",
        "stop_condition": "Do not publish final sensor positions or coverage until the fire/safety equipment basis is confirmed; the candidate mount may hide only the base and must leave smoke entry, indicator and downward removal unobstructed.",
    },
    {
        "queue_id": "M401-MISS-004",
        "observable_role": "bathroom_warm_air_unit_missing",
        "basis": "No observable bathroom heater, warm-air unit or ventilation fan instance/type exists in the formal IFC.",
        "missing_or_unverified": "equipment size, power, airflow path, grille, mounting and access conditions",
        "stop_condition": "Do not reserve an opening, circuit, grille or access zone from an assumed equipment size.",
    },
    {
        "queue_id": "M401-MISS-005",
        "observable_role": "operational_and_access_data_missing",
        "basis": "Existing AC appliances, flow-segment meshes and proxies contain observable geometry but no verified operational schedule or installation/access record.",
        "missing_or_unverified": "airflow, pipe/duct diameter, circuit, controls, condensate/refrigerant identity, access and manufacturer installation clearances",
        "stop_condition": "Do not infer performance, circuit, diameter, access or connection relations from names, bounding boxes or proximity.",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rounded(values: list[float]) -> list[float]:
    return [round(value, 6) for value in values]


def world_bbox(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> tuple[list[float], list[float], list[float]]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    values = list(shape.geometry.verts)
    if not values:
        raise RuntimeError(f"M-401 instance has no observable geometry: {product.GlobalId}")
    vertices = [
        [values[index] * 1000.0, values[index + 1] * 1000.0, values[index + 2] * 1000.0]
        for index in range(0, len(values), 3)
    ]
    minimum = [min(vertex[axis] for vertex in vertices) for axis in range(3)]
    maximum = [max(vertex[axis] for vertex in vertices) for axis in range(3)]
    dimensions = [maximum[axis] - minimum[axis] for axis in range(3)]
    return rounded(minimum), rounded(maximum), rounded(dimensions)


def container_name(product: ifcopenshell.entity_instance) -> str:
    container = ifcopenshell.util.element.get_container(product)
    return str(getattr(container, "Name", "") or "")


def type_occurrence_count(type_object: ifcopenshell.entity_instance) -> int:
    return sum(len(rel.RelatedObjects) for rel in getattr(type_object, "Types", ()) or ())


def instance_record(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
    source_hash: str,
    queue_id: str,
    observable_role: str,
    basis: str,
    confidence: float,
    human_review_required: bool,
    review_status: str,
    missing_or_unverified: str,
    stop_condition: str,
) -> dict[str, Any]:
    assigned_type = ifcopenshell.util.element.get_type(product)
    minimum, maximum, dimensions = world_bbox(settings, product)
    return {
        "record_kind": "actual_instance",
        "queue_id": queue_id,
        "global_id": product.GlobalId,
        "ifc_class": product.is_a(),
        "name": product.Name or "",
        "description": getattr(product, "Description", None) or "",
        "object_type": getattr(product, "ObjectType", None) or "",
        "predefined_type": getattr(product, "PredefinedType", None) or "",
        "type_global_id": getattr(assigned_type, "GlobalId", None) or "",
        "type_name": getattr(assigned_type, "Name", None) or "",
        "type_description": getattr(assigned_type, "Description", None) or "",
        "type_occurrence_count": type_occurrence_count(assigned_type) if assigned_type else "",
        "container": container_name(product),
        "observable_role": observable_role,
        "bbox_min_mm": minimum,
        "bbox_max_mm": maximum,
        "dimensions_mm": dimensions,
        "basis": basis,
        "confidence": confidence,
        "human_review_required": human_review_required,
        "review_status": review_status,
        "missing_or_unverified": missing_or_unverified,
        "stop_condition": stop_condition,
        "source_ifc_sha256": source_hash,
    }


def type_record(type_object: ifcopenshell.entity_instance, source_hash: str) -> dict[str, Any]:
    count = type_occurrence_count(type_object)
    return {
        "record_kind": "type_definition",
        "queue_id": "M401-TYPE-001",
        "global_id": "",
        "ifc_class": type_object.is_a(),
        "name": "",
        "description": "",
        "object_type": "",
        "predefined_type": getattr(type_object, "PredefinedType", None) or "",
        "type_global_id": type_object.GlobalId,
        "type_name": type_object.Name or "",
        "type_description": type_object.Description or "",
        "type_occurrence_count": count,
        "container": "",
        "observable_role": "assigned_ac_equipment_type_definition",
        "bbox_min_mm": "",
        "bbox_max_mm": "",
        "dimensions_mm": "",
        "basis": "IfcElectricApplianceType ElementType=AC and exact occurrence assignment.",
        "confidence": 1.0,
        "human_review_required": False,
        "review_status": "OBSERVED",
        "missing_or_unverified": "manufacturer installation data, capacity, airflow, circuit, controls and connections",
        "stop_condition": "Type identity does not authorize performance, service or installation assumptions.",
        "source_ifc_sha256": source_hash,
    }


def missing_record(item: dict[str, str], source_hash: str) -> dict[str, Any]:
    return {
        "record_kind": "missing_input",
        "queue_id": item["queue_id"],
        "global_id": "",
        "ifc_class": "",
        "name": "",
        "description": "",
        "object_type": "",
        "predefined_type": "",
        "type_global_id": "",
        "type_name": "",
        "type_description": "",
        "type_occurrence_count": "",
        "container": "",
        "observable_role": item["observable_role"],
        "bbox_min_mm": "",
        "bbox_max_mm": "",
        "dimensions_mm": "",
        "basis": item["basis"],
        "confidence": 1.0,
        "human_review_required": True,
        "review_status": "BLOCK",
        "missing_or_unverified": item["missing_or_unverified"],
        "stop_condition": item["stop_condition"],
        "source_ifc_sha256": source_hash,
    }


def encode(value: Any) -> Any:
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, bool):
        return "yes" if value else "no"
    return value


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({field: encode(record.get(field, "")) for field in CSV_FIELDS})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--decision-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = args.input.resolve()
    source_hash = sha256(source)
    model = ifcopenshell.open(source)
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    records: list[dict[str, Any]] = []

    ac_instances = []
    ac_types: dict[str, ifcopenshell.entity_instance] = {}
    for product in model.by_type("IfcElectricAppliance"):
        assigned_type = ifcopenshell.util.element.get_type(product)
        if not assigned_type or assigned_type.Name not in AC_TYPE_NAMES:
            continue
        ac_instances.append(product)
        ac_types[assigned_type.GlobalId] = assigned_type
        records.append(
            instance_record(
                settings, product, source_hash, "M401-EQUIP-001",
                "assigned_ac_equipment_instance",
                "Assigned IfcElectricApplianceType has ElementType=AC and preserves the exact model/type identity.",
                1.0, True, "HUMAN_REVIEW_REQUIRED",
                "final equipment selection, manufacturer installation data, performance, services, controls and access",
                "Do not infer airflow, circuit, pipe/duct connection, access or final equipment selection from the current shell.",
            )
        )

    hvac_flows = []
    for product in model.by_type("IfcFlowSegment"):
        if product.Name not in HVAC_FLOW_NAMES or product.ObjectType != "HVAC":
            continue
        hvac_flows.append(product)
        design_role = HVAC_FLOW_LEGACY_ROLES[product.Name]
        records.append(
            instance_record(
                settings, product, source_hash, "M401-FLOW-001",
                design_role["role"],
                design_role["basis"],
                1.0, True, "LEGACY_BASE_CONFIRMED_REMODEL_DESIGN_PENDING",
                design_role["missing"],
                design_role["stop"],
            )
        )

    check_valve = model.by_guid(CHECK_VALVE_ID)
    records.append(
        instance_record(
            settings, check_valve, source_hash, "M401-VALVE-001",
            "named_flue_check_valve_proxy",
            "Exact IfcBuildingElementProxy name is Electric flue check valve; existing Boolean-sensitive geometry is observable.",
            0.9, True, "HUMAN_REVIEW_REQUIRED",
            "served exhaust, product specification, power/control, connection and access conditions",
            "Do not infer the served appliance, duct connection, powered action or access requirement from proximity or proxy shape.",
        )
    )

    for global_id in sorted(DIFFUSER_PROXY_IDS):
        product = model.by_guid(global_id)
        decision = DIFFUSER_PROXY_DECISIONS[global_id]
        records.append(
            instance_record(
                settings, product, source_hash, "M401-DIFFUSER-001",
                "named_embedded_ac_diffuser_proxy",
                decision["basis"],
                decision["confidence"], True, decision["status"],
                decision["missing"],
                decision["stop"],
            )
        )

    for global_id in sorted(CEILING_CONTEXT_IDS):
        product = model.by_guid(global_id)
        records.append(
            instance_record(
                settings, product, source_hash, "M401-RCP-CONTEXT-001",
                "high_level_ceiling_or_led_coordination_context",
                "Exact IfcBuildingElementProxy name and elevated world geometry are observable; the object is retained only as RCP coordination context.",
                1.0, False, "OBSERVED_CONTEXT",
                "ceiling build-up, support, access zones and relation to unmodelled MEP equipment",
                "Do not treat a ceiling/LED proxy as HVAC equipment, a service opening or an access panel.",
            )
        )

    for type_object in sorted(ac_types.values(), key=lambda item: item.Name):
        records.append(type_record(type_object, source_hash))
    records.extend(missing_record(item, source_hash) for item in MISSING_INPUTS)

    instance_records = [record for record in records if record["record_kind"] == "actual_instance"]
    instance_ids = [record["global_id"] for record in instance_records]
    if len(instance_ids) != len(set(instance_ids)):
        raise RuntimeError("duplicate M-401 actual-instance assignment")
    counts = Counter(record["observable_role"] for record in instance_records)
    expected_counts = Counter(
        {
            "assigned_ac_equipment_instance": 5,
            "legacy_base_refrigerant_liquid_geometry": 2,
            "legacy_base_refrigerant_gas_geometry": 2,
            "legacy_base_condensate_geometry": 1,
            "named_flue_check_valve_proxy": 1,
            "named_embedded_ac_diffuser_proxy": 2,
            "high_level_ceiling_or_led_coordination_context": 15,
        }
    )
    if counts != expected_counts:
        raise RuntimeError(f"unexpected M-401 instance inventory: {counts}")
    if Counter(item.Name for item in ac_types.values()) != Counter({"AC1180": 1, "AC700": 1, "AC700F": 1}):
        raise RuntimeError("unexpected M-401 AC type inventory")
    if len(MISSING_INPUTS) != 5:
        raise RuntimeError("expected five explicit M-401 missing-input gates")

    write_csv(args.decision_csv, records)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    review_queue = [
        {
            "queue_id": record["queue_id"],
            "global_id": record["global_id"],
            "observable_role": record["observable_role"],
            "name": record["name"],
            "type_name": record["type_name"],
            "missing_or_unverified": record["missing_or_unverified"],
            "stop_condition": record["stop_condition"],
        }
        for record in records
        if record["human_review_required"] is True
    ]
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_m401_existing_candidate",
        "source": {"ifc": str(source), "ifc_sha256": source_hash},
        "summary": {
            "actual_instances": len(instance_records),
            "instance_role_counts": dict(sorted(counts.items())),
            "type_definitions": len(ac_types),
            "missing_input_blocks": len(MISSING_INPUTS),
            "human_review_queue": len(review_queue),
            "ifc_air_terminal_instances": len(model.by_type("IfcAirTerminal")),
            "ifc_fan_instances": len(model.by_type("IfcFan")),
            "ifc_sensor_instances": len(model.by_type("IfcSensor")),
            "ifc_alarm_instances": len(model.by_type("IfcAlarm")),
            "ifc_distribution_ports": len(model.by_type("IfcDistributionPort")),
            "ifc_systems": len(model.by_type("IfcSystem")),
            "ifc_port_connections": len(model.by_type("IfcRelConnectsPorts")),
        },
        "gates": {
            "inventory_pass": True,
            "formal_ifc_write_allowed": False,
            "m401_design_ready": False,
            "rcp1_completion_pass": False,
        },
        "review_queue": review_queue,
        "records": records,
    }
    (args.output_dir / "m401-existing-report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# M-401 / RCP1 Existing Candidate",
        "",
        f"- IFC SHA-256: `{source_hash}`",
        f"- Actual instances: {len(instance_records)}",
        f"- AC type definitions: {len(ac_types)}",
        f"- Human review queue: {len(review_queue)}",
        f"- Explicit missing-input BLOCKs: {len(MISSING_INPUTS)}",
        "",
        "## Concentrated human review",
        "",
    ]
    for item in review_queue:
        identity = item["global_id"] or item["queue_id"]
        label = item["type_name"] or item["name"] or item["observable_role"]
        lines.append(f"- `{identity}` — {label}: {item['missing_or_unverified']}")
    lines.extend(
        [
            "",
            "> Bounding boxes are observable coordination envelopes only. No airflow, diameter, circuit, access, connection or system relation is inferred.",
            "",
        ]
    )
    (args.output_dir / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
