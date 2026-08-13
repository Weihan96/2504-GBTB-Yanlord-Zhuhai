"""Build read-only RCP1C route constraints from confirmed HVAC decisions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell


EQUIPMENT = {
    "A01": {"global_id": "33chLv3TzEOhKalIJJAPNF", "label": "主卧入口"},
    "A02": {"global_id": "1QBdVekDnBsOleyo9PM6rT", "label": "走廊"},
    "A03": {"global_id": "06GpMzzWj1XQobAhD35cgU", "label": "次卧"},
    "A04": {"global_id": "1PUCikoaP5fgiYt8sJd8$6", "label": "西厨"},
    "A05": {"global_id": "1yW7DASIz8qA$2j8z9tdl2", "label": "书房"},
    "A06": {"global_id": "0YzEUom7522RIg1TonOQXn", "label": "东侧固定机位"},
}

OPENINGS = {
    "H01": "0DOeKdT3DE9f2LMR_x$G7q",
    "H02": "1JVGu2xtb1ZhGElhvSmyd$",
    "H03": "2oHWzdjkr8X8Yt0Gd2e3RQ",
    "H04": "3ERXx822H9jOPo6CetKX9r",
    "H05": "3k2aFe_p5Ah8HJEqWUBB7F",
    "H06": "3qgu$TepT0J8Yat7ZQ90Mf",
    "H07": "1_EX1UWfL8ShcGm5BI1UPc",
}

DIRECT_SERVICE = {
    "A01": ("R09", "主卧", "geometry_probe"),
    "A02": ("R07", "餐厅", "user_confirmed"),
    "A03": ("R14", "次卧", "geometry_probe"),
    "A04": ("R03", "西厨", "geometry_probe"),
    "A05": ("R22", "书房", "geometry_probe"),
    "A06": ("R07", "餐厅", "user_confirmed"),
}

EQUIPMENT_OPENING_MAPPING = {
    "A01": {
        "opening_id": "H05",
        "status": "geometry_candidate",
        "basis": "nearest existing opening by current IFC world-AABB clearance",
    },
    "A02": {
        "opening_id": "H03",
        "status": "user_confirmed",
        "basis": "confirmed RCP1-SERVICE-A02 waypoint path A02 → H03",
    },
    "A03": {
        "opening_id": "H04",
        "status": "user_confirmed",
        "basis": "confirmed RCP1-SERVICE-A03 first equipment-side waypoint A03 → H04",
    },
    "A04": {
        "opening_id": "H06",
        "status": "geometry_candidate",
        "basis": "nearest existing opening by current IFC world-AABB clearance",
    },
    "A05": {
        "opening_id": "H07",
        "status": "geometry_candidate",
        "basis": "nearest existing opening by current IFC world-AABB clearance; H07 endpoint role is separately user-confirmed",
    },
    "A06": {
        "opening_id": "H03",
        "status": "geometry_candidate",
        "basis": "nearest existing opening by current IFC world-AABB clearance; shared-opening use remains a candidate",
    },
}

INTERFACE_ROLES = {
    "H01": "user_confirmed_outdoor_unit_interface",
    "H02": "user_confirmed_multihop_downstream_anchor",
    "H03": "user_confirmed_A02_service_opening_with_A06_candidate",
    "H04": "user_confirmed_A03_equipment_side_opening",
    "H05": "geometry_candidate_equipment_service_opening",
    "H06": "geometry_candidate_equipment_service_opening",
    "H07": "user_confirmed_condensate_endpoint_with_A05_candidate",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def requirement_value(row: dict[str, str]) -> str | float | None:
    text_value = row.get("value_text", "").strip()
    if text_value:
        return text_value
    number_value = row.get("value_number", "").strip()
    if not number_value:
        return None
    value = float(number_value)
    return int(value) if value.is_integer() else value


def manufacturer_interface_inputs(
    equipment_rows: list[dict[str, str]], requirement_rows: list[dict[str, str]]
) -> dict:
    candidate_ids_by_equipment: dict[str, list[str]] = {}
    for row in equipment_rows:
        equipment_id = row.get("equipment_id", "")
        if not equipment_id.startswith("HVAC-"):
            continue
        candidate_ids_by_equipment[equipment_id] = sorted(
            set(re.findall(r"A0[1-6]", row.get("item_name", "")))
        )

    by_equipment = {
        candidate_id: {"requirements": [], "release_blocking_requirements": []}
        for candidate_id in EQUIPMENT
    }
    consumed = 0
    blocking = 0
    for row in requirement_rows:
        candidate_ids = candidate_ids_by_equipment.get(row.get("equipment_id", ""), [])
        if not candidate_ids or "HVAC" not in row.get("discipline", ""):
            continue
        requirement = {
            "requirement_id": row["requirement_id"],
            "requirement_name": row["parameter_key"],
            "value": requirement_value(row),
            "unit": row.get("unit", "") or None,
            "status": row.get("status", ""),
            "basis_kind": row.get("value_origin", ""),
            "evidence_id": row.get("source_id", "") or None,
            "blocks_release": row.get("blocks_release", "").lower() == "yes",
            "notes": row.get("notes", "") or None,
        }
        for candidate_id in candidate_ids:
            by_equipment[candidate_id]["requirements"].append(requirement)
            consumed += 1
            if requirement["blocks_release"]:
                by_equipment[candidate_id]["release_blocking_requirements"].append(
                    requirement["requirement_name"]
                )
                blocking += 1

    for payload in by_equipment.values():
        payload["requirements"].sort(key=lambda row: row["requirement_id"])
        payload["release_blocking_requirements"] = sorted(
            set(payload["release_blocking_requirements"])
        )
    return {
        "source_of_truth": "equipment-register.csv + equipment-installation-requirements.csv",
        "consumed_requirement_count": consumed,
        "release_blocking_requirement_count": blocking,
        "by_equipment": by_equipment,
    }


def bbox_clearance_mm(first: dict, second: dict) -> float:
    deltas = [
        max(0.0, second["min_mm"][axis] - first["max_mm"][axis], first["min_mm"][axis] - second["max_mm"][axis])
        for axis in range(3)
    ]
    return math.sqrt(sum(delta * delta for delta in deltas))


def relation_distance(item: dict, opening_global_id: str) -> float:
    relation = next(
        relation
        for relation in item["opening_relations"]
        if relation["global_id"] == opening_global_id
    )
    return float(relation["minimum_clearance_candidate_mm"])


def equipment_row(candidate_id: str, report: dict, opening_by_id: dict[str, dict]) -> dict:
    definition = EQUIPMENT[candidate_id]
    if candidate_id == "A06":
        legacy = report["legacy_east_ac_candidate"]
        bbox = legacy["predicted_formal_bbox"]
        ranked = sorted(
            (
                {
                    "opening_id": opening_id,
                    "global_id": global_id,
                    "clearance_mm": bbox_clearance_mm(bbox, opening_by_id[global_id]["bbox"]),
                }
                for opening_id, global_id in OPENINGS.items()
            ),
            key=lambda row: (row["clearance_mm"], row["opening_id"]),
        )
        centre = legacy["predicted_formal_centre_mm"]
        formal_identity_status = "present_placement_only"
    else:
        item = next(
            item
            for item in report["formal_equipment_pairing"]
            if item["global_id"] == definition["global_id"]
        )
        ranked = sorted(
            (
                {
                    "opening_id": opening_id,
                    "global_id": global_id,
                    "clearance_mm": relation_distance(item, global_id),
                }
                for opening_id, global_id in OPENINGS.items()
            ),
            key=lambda row: (row["clearance_mm"], row["opening_id"]),
        )
        centre = item["bbox"]["centre_mm"]
        formal_identity_status = "present"
    service = DIRECT_SERVICE[candidate_id]
    return {
        "equipment_id": candidate_id,
        "global_id": definition["global_id"],
        "label": definition["label"],
        "position_status": "confirmed_fixed",
        "formal_identity_status": formal_identity_status,
        "centre_mm": centre,
        "direct_service": {
            "space_reference": service[0],
            "space_long_name": service[1],
            "status": service[2],
        },
        "opening_ranking": ranked,
    }


def validate_decisions(route_rows: list[dict[str, str]], waypoint_rows: list[dict[str, str]]) -> None:
    routes = {row["route_id"]: row for row in route_rows}
    expected_service = {
        "RCP1-AIRSIDE-A02": ("A02", "R07"),
        "RCP1-AIRSIDE-A06": ("A06", "R07"),
    }
    for route_id, expected in expected_service.items():
        row = routes.get(route_id)
        if row is None or (row["equipment_id"], row["served_space_reference"]) != expected:
            raise RuntimeError(f"confirmed airside service drift: {route_id}")
    ordered = {}
    for row in waypoint_rows:
        ordered.setdefault(row["route_id"], []).append(row)
    for rows in ordered.values():
        rows.sort(key=lambda row: int(row["sequence"]))
    expected_paths = {
        "RCP1-SERVICE-A02": ["A02", "H03"],
        "RCP1-SERVICE-A03": ["A03", "H04", "H02"],
        "RCP1-CONDENSATE-ENDPOINT": ["H07"],
        "RCP1-OUTDOOR-ENDPOINT": ["H01"],
    }
    for route_id, expected in expected_paths.items():
        actual = [row["anchor_id"] for row in ordered.get(route_id, [])]
        if actual != expected:
            raise RuntimeError(f"confirmed waypoint path drift: {route_id}: {actual}")
    if routes["RCP1-CONDENSATE-ENDPOINT"]["terminal_anchor_id"] != "H07":
        raise RuntimeError("condensate endpoint must remain H07")
    if routes["RCP1-OUTDOOR-ENDPOINT"]["terminal_anchor_id"] != "H01":
        raise RuntimeError("outdoor-unit interface must remain H01")


def resolve_waypoint(
    row: dict[str, str], equipment_by_id: dict[str, dict], opening_by_code: dict[str, dict]
) -> dict:
    anchor_id = row["anchor_id"]
    if row["anchor_kind"] == "equipment":
        centre = equipment_by_id[anchor_id]["centre_mm"]
    elif row["anchor_kind"] == "existing_opening":
        centre = opening_by_code[anchor_id]["centre_mm"]
    elif row["anchor_kind"] == "blender_bend":
        centre = [float(row[key]) for key in ("x_mm", "y_mm", "z_mm")]
    else:
        raise RuntimeError(f"unsupported route anchor kind: {row['anchor_kind']}")
    return {
        **row,
        "sequence": int(row["sequence"]),
        "centre_mm": centre,
        "confidence": float(row["confidence"]),
        "review_required": row["review_required"].lower() == "yes",
    }


def route_graph(
    route_rows: list[dict[str, str]],
    waypoint_rows: list[dict[str, str]],
    equipment_by_id: dict[str, dict],
    opening_by_code: dict[str, dict],
) -> list[dict]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in waypoint_rows:
        grouped.setdefault(row["route_id"], []).append(row)
    routes = []
    for route in route_rows:
        points = [
            resolve_waypoint(row, equipment_by_id, opening_by_code)
            for row in sorted(grouped.get(route["route_id"], []), key=lambda row: int(row["sequence"]))
        ]
        segments = []
        for first, second in zip(points, points[1:]):
            length = math.dist(first["centre_mm"], second["centre_mm"])
            segments.append({
                "from_anchor_id": first["anchor_id"],
                "to_anchor_id": second["anchor_id"],
                "centre_to_centre_length_mm": length,
                "status": "constraint_skeleton_not_fabrication_geometry",
            })
        routes.append({
            **route,
            "confidence": float(route["confidence"]),
            "review_required": route["review_required"].lower() == "yes",
            "formal_ifc_write_allowed": route["formal_ifc_write_allowed"].lower() == "yes",
            "waypoints": points,
            "segments": segments,
        })
    return routes


def mapping_coverage(equipment: list[dict], routes: list[dict]) -> tuple[list[dict], list[dict]]:
    equipment_by_id = {row["equipment_id"]: row for row in equipment}
    mappings: list[dict] = []
    for equipment_id in sorted(EQUIPMENT_OPENING_MAPPING):
        definition = EQUIPMENT_OPENING_MAPPING[equipment_id]
        opening_id = definition["opening_id"]
        ranked = equipment_by_id[equipment_id]["opening_ranking"]
        relation = next(row for row in ranked if row["opening_id"] == opening_id)
        rank = next(index for index, row in enumerate(ranked, 1) if row["opening_id"] == opening_id)
        if definition["status"] == "geometry_candidate" and rank != 1:
            raise RuntimeError(
                f"geometry-derived mapping is no longer nearest: {equipment_id} → {opening_id} rank={rank}"
            )
        mappings.append({
            "equipment_id": equipment_id,
            "equipment_global_id": equipment_by_id[equipment_id]["global_id"],
            "opening_id": opening_id,
            "opening_global_id": relation["global_id"],
            "minimum_clearance_candidate_mm": relation["clearance_mm"],
            "opening_rank_by_clearance": rank,
            "status": definition["status"],
            "basis": definition["basis"],
            "direct_service": equipment_by_id[equipment_id]["direct_service"],
            "formal_ifc_write_allowed": False,
        })

    route_ids_by_anchor: dict[str, set[str]] = {opening_id: set() for opening_id in OPENINGS}
    for route in routes:
        for waypoint in route["waypoints"]:
            anchor_id = waypoint["anchor_id"]
            if anchor_id in route_ids_by_anchor:
                route_ids_by_anchor[anchor_id].add(route["route_id"])
        terminal_anchor_id = route.get("terminal_anchor_id", "")
        if terminal_anchor_id in route_ids_by_anchor:
            route_ids_by_anchor[terminal_anchor_id].add(route["route_id"])

    coverage = []
    for opening_id in sorted(OPENINGS):
        assigned = [row for row in mappings if row["opening_id"] == opening_id]
        coverage.append({
            "opening_id": opening_id,
            "opening_global_id": OPENINGS[opening_id],
            "role": INTERFACE_ROLES[opening_id],
            "route_ids": sorted(route_ids_by_anchor[opening_id]),
            "equipment_mappings": [
                {"equipment_id": row["equipment_id"], "status": row["status"]}
                for row in assigned
            ],
            "coverage_status": (
                "user_confirmed_route_or_endpoint"
                if route_ids_by_anchor[opening_id]
                else "geometry_candidate"
            ),
            "formal_ifc_write_allowed": False,
        })
    return mappings, coverage


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--hvac-report", type=Path, required=True)
    parser.add_argument("--legacy-report", type=Path, required=True)
    parser.add_argument("--delivery-dwg", type=Path, required=True)
    parser.add_argument("--routes", type=Path, required=True)
    parser.add_argument("--waypoints", type=Path, required=True)
    parser.add_argument("--equipment-register", type=Path, required=True)
    parser.add_argument("--requirements", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    formal_sha = sha256(args.input)
    hvac = json.loads(args.hvac_report.read_text(encoding="utf-8"))
    legacy = json.loads(args.legacy_report.read_text(encoding="utf-8"))
    if hvac["source"]["ifc_sha256"] != formal_sha:
        raise RuntimeError("HVAC report does not match the formal IFC")
    if legacy["source"]["formal_ifc_sha256"] != formal_sha:
        raise RuntimeError("legacy blend audit does not match the formal IFC")
    route_rows = read_csv(args.routes)
    waypoint_rows = read_csv(args.waypoints)
    equipment_rows = read_csv(args.equipment_register)
    requirement_rows = read_csv(args.requirements)
    validate_decisions(route_rows, waypoint_rows)
    manufacturer_inputs = manufacturer_interface_inputs(equipment_rows, requirement_rows)

    model = ifcopenshell.open(args.input)
    opening_by_global_id = {item["global_id"]: item for item in hvac["developer_opening_pairing"]}
    equipment = [equipment_row(candidate_id, hvac, opening_by_global_id) for candidate_id in sorted(EQUIPMENT)]
    equipment_by_id = {row["equipment_id"]: row for row in equipment}
    openings = [
        {
            "opening_id": opening_id,
            "global_id": global_id,
            "ifc_class": "IfcOpeningElement",
            "name": opening_by_global_id[global_id]["name"],
            "hosts": opening_by_global_id[global_id]["hosts"],
            "centre_mm": opening_by_global_id[global_id]["bbox"]["centre_mm"],
        }
        for opening_id, global_id in OPENINGS.items()
    ]
    opening_by_code = {row["opening_id"]: row for row in openings}
    routes = route_graph(route_rows, waypoint_rows, equipment_by_id, opening_by_code)
    equipment_opening_mappings, interface_coverage = mapping_coverage(equipment, routes)

    port_count = len(model.by_type("IfcDistributionPort"))
    system_count = len(model.by_type("IfcSystem"))
    distribution_system_count = len(model.by_type("IfcDistributionSystem"))
    outdoor_equipment_count = sum(len(model.by_type(name)) for name in ("IfcCondenser", "IfcCompressor"))
    fixed_diagnostics = hvac["fixed_equipment_airside_diagnostics"]
    a05_diagnostic = next(
        item for item in fixed_diagnostics if item["equipment_candidate"] == EQUIPMENT["A05"]["global_id"]
    )
    legacy_topology_ok = all(
        row["formal_component_count"] == row["legacy_component_count"]
        for row in legacy["pipe_comparisons"]
    )
    maximum_legacy_difference = max(
        row["maximum_component_dimension_difference_mm"]
        for row in legacy["pipe_comparisons"]
    )

    output = {
        "mode": "read_only_rcp1_hvac_constraint_graph",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source": {
            "ifc": str(args.input.resolve()),
            "ifc_sha256": formal_sha,
            "hvac_report": str(args.hvac_report.resolve()),
            "legacy_blend": legacy["source"]["legacy_blend"],
            "legacy_blend_sha256": legacy["source"]["legacy_blend_sha256"],
            "delivery_dwg": str(args.delivery_dwg.resolve()),
            "delivery_dwg_sha256": sha256(args.delivery_dwg),
            "route_register": str(args.routes.resolve()),
            "waypoint_register": str(args.waypoints.resolve()),
            "equipment_register": str(args.equipment_register.resolve()),
            "equipment_register_sha256": sha256(args.equipment_register),
            "installation_requirements": str(args.requirements.resolve()),
            "installation_requirements_sha256": sha256(args.requirements),
        },
        "source_dependencies": [
            {"path": str(path.resolve()), "sha256": sha256(path)}
            for path in (
                args.hvac_report,
                args.routes,
                args.waypoints,
                args.equipment_register,
                args.requirements,
            )
        ],
        "scope": {
            "formal_ifc_write_allowed": False,
            "equipment_positions_may_move": False,
            "new_openings_allowed": False,
            "demolition_walls_are_permanent_obstacles": False,
            "legacy_pipe_geometry_is_final_route": False,
            "blender_preview_is_source_of_truth": False,
        },
        "equipment": equipment,
        "openings": openings,
        "confirmed_route_graph": routes,
        "equipment_opening_mapping": equipment_opening_mappings,
        "interface_coverage": interface_coverage,
        "interface_roles": {
            "H01": {
                "role": "outdoor_unit_interface_at_existing_opening",
                "ifc_entity_remains": "IfcOpeningElement",
                "status": "user_confirmed_location_semantic_pending",
            },
            "H07": {
                "role": "condensate_discharge_interface_at_existing_opening",
                "ifc_entity_remains": "IfcOpeningElement",
                "status": "user_confirmed_location_semantic_pending",
            },
        },
        "legacy_reference": {
            "component_topology_matches": legacy_topology_ok,
            "maximum_component_dimension_difference_mm": maximum_legacy_difference,
            "role": "post_demolition_approximate_design_base_not_final_route",
        },
        "delivery_dwg_reference": {
            "role": "developer_delivery_hvac_and_existing_opening_reference",
            "coordinate_copy_allowed": False,
            "reason": "the drawing documents developer equipment, systems and openings but is not the remodel routing model",
        },
        "airside_readiness": {
            "direct_service": {item["equipment_id"]: item["direct_service"] for item in equipment},
            "a02_and_a06_both_serve_dining": True,
            "living_room_R20_direct_candidate_count": 0,
            "a05_to_R20_demolition_wall_crossings": a05_diagnostic["route"]["demolition_wall_crossings"],
            "a05_to_R20_permanent_wall_crossings": a05_diagnostic["route"]["permanent_wall_crossings"],
            "status": "service_rooms_partly_confirmed_supply_return_geometry_pending",
        },
        "refrigerant_and_condensate_readiness": {
            "distribution_port_count": port_count,
            "ifc_system_count": system_count,
            "ifc_distribution_system_count": distribution_system_count,
            "outdoor_condenser_or_compressor_count": outdoor_equipment_count,
            "legacy_pipe_products": hvac["summary"]["legacy_pipe_products"],
            "legacy_pipe_independent_components": hvac["summary"]["legacy_pipe_independent_components"],
            "status": "constraint_skeleton_started_real_ports_sections_slopes_and_fittings_pending",
        },
        "manufacturer_interface_inputs": manufacturer_inputs,
        "authoring_contract": {
            "source_of_truth": "IFC plus route and waypoint decision registers",
            "blender_role": "rebuildable constraint editing and Geometry Nodes preview only",
            "update_trigger": "explicit rebuild preview or compile candidate command; never an automatic IFC write from depsgraph",
            "future_ifc_topology": "typed duct/pipe segments and fittings with nested IfcDistributionPort connections",
        },
        "minimum_required_inputs": [
            {
                "input_id": "RCP1C-I04",
                "question": "逐台确认真实送回风口、冷媒液/气管口和冷凝水口的厂家坐标与朝向。",
                "why_required": "当前约束点使用设备包围盒中心，不能冒充厂家接口。",
            },
            {
                "input_id": "RCP1C-I05",
                "question": "确认尚未被厂家证据覆盖的风量与风管截面、保温外径、弯曲半径、最终锚固/吊架和检修净距；补齐 A05/A06 冷媒与排水参数。",
                "why_required": "A01–A04 已有名义冷媒管径、排水外径与坡度，A05 已有吊杆/风口/检修候选；当前仍缺施工路线所需的完整截面、端口和现场构造。",
            },
        ],
        "gates": {
            "source_hash_matches": True,
            "six_fixed_equipment_positions_registered": len(equipment) == 6,
            "seven_existing_openings_registered": len(openings) == 7,
            "six_equipment_opening_mappings_registered": len(equipment_opening_mappings) == 6,
            "seven_interfaces_have_controlled_roles": (
                len(interface_coverage) == 7
                and all(row["role"] for row in interface_coverage)
            ),
            "confirmed_and_candidate_mapping_split_preserved": (
                sum(row["status"] == "user_confirmed" for row in equipment_opening_mappings) == 2
                and sum(row["status"] == "geometry_candidate" for row in equipment_opening_mappings) == 4
                and all(not row["formal_ifc_write_allowed"] for row in equipment_opening_mappings)
            ),
            "confirmed_shared_and_multihop_graph_ready": True,
            "manufacturer_constraints_consumed_from_ssot": manufacturer_inputs["consumed_requirement_count"] > 0,
            "a01_a04_refrigerant_and_condensate_nominals_available": all(
                {
                    "gas_pipe_od",
                    "liquid_pipe_od",
                    "drain_pipe_od",
                    "drain_slope_min",
                    "drain_slope_max",
                }.issubset({
                    row["requirement_name"]
                    for row in manufacturer_inputs["by_equipment"][candidate_id]["requirements"]
                    if row["value"] is not None
                })
                for candidate_id in ("A01", "A02", "A03", "A04")
            ),
            "legacy_pipe_topology_matches": legacy_topology_ok,
            "demolition_walls_excluded_from_permanent_obstacles": (
                len(a05_diagnostic["route"]["demolition_wall_crossings"]) >= 1
                and len(a05_diagnostic["route"]["permanent_wall_crossings"]) == 0
            ),
            "final_airside_routes_ready": False,
            "final_refrigerant_routes_ready": False,
            "final_condensate_routes_ready": False,
            "formal_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "fixed_equipment": len(equipment),
        "existing_openings": len(openings),
        "equipment_opening_mappings": len(equipment_opening_mappings),
        "controlled_interface_roles": len(interface_coverage),
        "confirmed_routes": len(routes),
        "confirmed_route_segments": sum(len(route["segments"]) for route in routes),
        "legacy_topology_matches": legacy_topology_ok,
        "required_inputs": len(output["minimum_required_inputs"]),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
