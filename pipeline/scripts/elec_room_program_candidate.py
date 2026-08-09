#!/usr/bin/env python3
"""Compile a read-only room-by-room electrical function program and equipment demand register."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from collections import Counter
from pathlib import Path
from typing import Any


EXPECTED_IFC_SHA256 = "7521c09991f3d0c7b7d91ca2324fd55ad961d8e32e9e3e9a9777a4cc19b06e81"
ROOM_REFERENCE_BY_NAME = {
    "玄关": "R01", "走廊": "R02", "西厨": "R03", "中厨": "R04", "中厨飘窗": "R05",
    "餐厅飘窗": "R06", "餐厅": "R07", "主卧飘窗": "R08", "主卧": "R09", "主卧入口": "R10",
    "主卫干区": "R11", "主卫湿区": "R12", "主卫湿区飘窗": "R13", "次卧": "R14",
    "次卧飘窗": "R15", "客卫干区": "R16", "客卫": "R17", "客卫飘窗": "R18", "阳台": "R19",
    "客厅": "R20", "书房飘窗": "R21", "书房": "R22",
}
POWERED_PROXY_ROLES = {
    "refrigerator_volume_candidate",
    "food_waste_disposer",
    "dishwasher_model_candidate",
    "island_dishwasher",
    "electric_flue_check_valve",
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument("--program", type=Path, default=root / "pipeline/decisions/elec-room-program.csv")
    parser.add_argument("--design-rules", type=Path, default=root / "pipeline/decisions/elec-design-rules.csv")
    parser.add_argument("--delta", type=Path, default=root / "build/mep-positioning/mep-renovation-delta-candidate.json")
    parser.add_argument("--elec-existing", type=Path, default=root / "build/elec/elec-existing-candidate.json")
    parser.add_argument("--elec-positioning", type=Path, default=root / "build/elec/elec-positioning-candidate.json")
    parser.add_argument("--round1", type=Path, default=root / "build/elec/elec-renovation-round1-candidate.json")
    parser.add_argument("--int1", type=Path, default=root / "build/int1/int1-existing-report.json")
    parser.add_argument("--route-readiness", type=Path, default=root / "build/rcp1/route-readiness-candidate.json")
    parser.add_argument("--output", type=Path, default=root / "build/elec/elec-room-program-candidate.json")
    parser.add_argument("--output-svg", type=Path, default=root / "drawings/E302-E304-room-program-candidate.svg")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_program(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key in (
            "lighting_control_groups_min", "general_power_groups_min", "network_data_groups_min",
            "manual_dedicated_power_min", "cabinet_light_feed_zones_min",
        ):
            row[key] = int(row[key])
        row["confidence"] = float(row["confidence"])
    refs = [row["room_reference"] for row in rows]
    if refs != [f"R{index:02d}" for index in range(1, 23)]:
        raise RuntimeError("electrical room program must contain R01-R22 in order")
    if len({row["room_name"] for row in rows}) != 22:
        raise RuntimeError("electrical room names are not unique")
    return rows


def read_design_rules(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ids = [row["rule_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise RuntimeError("electrical design rule IDs are not unique")
    for row in rows:
        row["confidence"] = float(row["confidence"])
    return rows


def compile_design_rule_summary(rules: list[dict[str, Any]]) -> dict[str, Any]:
    confirmed = [row for row in rules if row["status"] == "confirmed"]
    datums = {
        row["controlled_or_served_scope"]: int(row["value"])
        for row in confirmed
        if row["rule_kind"] == "installation_datum"
    }
    network_roles = Counter(
        row["device_or_datum"]
        for row in confirmed
        if row["rule_kind"] == "network_role"
    )
    return {
        "physical_wired_switch_only": any(
            row["rule_kind"] == "control_method" and row["value"] == "physical_wired_only"
            for row in confirmed
        ),
        "panel_bottom_AFF_mm": datums,
        "bedroom_AP_count": network_roles["wireless_access_point"],
        "living_study_router_count": network_roles["router_no_AP"],
        "paired_two_way_control_groups": sum(row["rule_kind"] == "two_way_control" for row in confirmed),
        "entrance_master_lighting_switches": sum(row["rule_kind"] == "master_control" for row in confirmed),
        "confirmed_island_dishwashers": sum(row["rule_kind"] == "equipment_identity" for row in confirmed),
        "demolition_walls_hidden_in_general_reviews": any(
            row["rule_kind"] == "review_visibility" and row["value"] == "hidden"
            for row in confirmed
        ),
    }


def centre(bbox: dict[str, list[float]]) -> list[float]:
    return [round(float(value), 6) for value in bbox["centre_mm"]]


def typed_equipment_demands(elec: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for record in elec["sheets"]["E-303"]["typed_equipment"]:
        name = record["assigned_type"]["name"] or record["name"]
        room_name = record["candidate_space"]["long_name"]
        rows.append({
            "source_kind": "formal_typed_equipment",
            "source_id": record["candidate_id"],
            "global_id": record["global_id"],
            "equipment_role": name,
            "placement_room_reference": ROOM_REFERENCE_BY_NAME[room_name],
            "placement_room_name": room_name,
            "coordination_centre_mm": centre(record["bbox"]),
            "basis": "formal IFC equipment instance and assigned type; centre is coordination evidence, not a socket position",
            "confidence": 1.0,
            "review_required": True,
            "position_status": "equipment_envelope_only",
            "automatic_ifc_write_allowed": False,
        })
    return rows


def proxy_equipment_demands(elec: dict[str, Any], positioning: dict[str, Any]) -> list[dict[str, Any]]:
    proxy_by_id = {row["global_id"]: row for row in elec["sheets"]["E-303"]["proxy_handoffs"]}
    rows = []
    for identity in positioning["proxy_identity_candidates"]:
        if identity["candidate_role"] not in POWERED_PROXY_ROLES:
            continue
        record = proxy_by_id[identity["global_id"]]
        room_name = record["candidate_space"]["long_name"]
        rows.append({
            "source_kind": "formal_proxy_equipment_candidate",
            "source_id": record["candidate_id"],
            "global_id": record["global_id"],
            "equipment_role": identity["candidate_role"],
            "placement_room_reference": ROOM_REFERENCE_BY_NAME[room_name],
            "placement_room_name": room_name,
            "coordination_centre_mm": centre(record["bbox"]),
            "basis": identity["basis"] + "; current proxy centre is not a final connection point",
            "confidence": identity["confidence"],
            "review_required": True,
            "position_status": "proxy_identity_and_envelope_only",
            "automatic_ifc_write_allowed": False,
        })
    return rows


def supplemental_equipment_demands(int1: dict[str, Any], route: dict[str, Any]) -> list[dict[str, Any]]:
    wash_tower = next(row for row in int1["records"] if row["global_id"] == "3JAkt8PsX7vPfGKWLK5EKp")
    wash_centre = [
        round((float(wash_tower["bbox_min_mm"][axis]) + float(wash_tower["bbox_max_mm"][axis])) / 2.0, 6)
        for axis in range(3)
    ]
    a06 = next(row for row in route["equipment"] if row["equipment_id"] == "A06")
    return [
        {
            "source_kind": "formal_named_equipment_proxy",
            "source_id": "INT1-WASHTOWER",
            "global_id": wash_tower["global_id"],
            "equipment_role": "LG WashTower",
            "placement_room_reference": "R04",
            "placement_room_name": "中厨",
            "coordination_centre_mm": wash_centre,
            "basis": "formal IFC named LG WashTower envelope contained by R04; manufacturer connection point remains unknown",
            "confidence": 1.0,
            "review_required": True,
            "position_status": "equipment_envelope_only",
            "automatic_ifc_write_allowed": False,
        },
        {
            "source_kind": "hash_fixed_legacy_equipment_candidate",
            "source_id": "A06",
            "global_id": "",
            "equipment_role": "east AC serving dining room",
            "placement_room_reference": "R02",
            "placement_room_name": "走廊",
            "served_room_reference": "R07",
            "coordination_centre_mm": [round(float(value), 6) for value in a06["centre_mm"]],
            "basis": "user-confirmed fixed A06 position recovered from the hash-fixed legacy blend; formal IFC identity is still missing",
            "confidence": 1.0,
            "review_required": True,
            "position_status": "equipment_identity_missing_connection_pending",
            "automatic_ifc_write_allowed": False,
        },
    ]


def round1_counts(round1: dict[str, Any]) -> dict[str, Counter[str]]:
    counts: dict[str, Counter[str]] = {
        reference: Counter() for reference in ROOM_REFERENCE_BY_NAME.values()
    }
    for row in round1["bedside_light_candidates"]:
        counts[row["room_reference"]]["bedside_lights"] += 1
    for row in round1["new_socket_candidates"]:
        counts[row["room_reference"]]["new_socket_candidates"] += 1
    for row in round1["cabinet_power_zones"]:
        counts[row["room_reference"]]["cabinet_power_zones"] += 1
    for row in round1["kitchen_socket_rechecks"]:
        counts[row["room_reference"]]["socket_rechecks"] += 1
    return counts


def compile_rooms(
    program: list[dict[str, Any]], delta: dict[str, Any], equipment: list[dict[str, Any]], round1: dict[str, Any]
) -> list[dict[str, Any]]:
    observed = {row["reference"]: row for row in delta["room_program"]}
    equipment_counts = Counter(row["placement_room_reference"] for row in equipment)
    r1_counts = round1_counts(round1)
    rooms = []
    for rule in program:
        reference = rule["room_reference"]
        current = observed[reference]
        rooms.append({
            **rule,
            "current_light_instances": current["current_light_instances"],
            "current_socket_instances": current["current_socket_instances"],
            "developer_power_references": current["developer_power_or_other_electrical_references"],
            "developer_control_references": current["developer_switch_or_control_references"],
            "developer_weak_references": current["developer_weak_point_references"],
            "modelled_equipment_power_demands": equipment_counts[reference],
            "round1_bedside_lights": r1_counts[reference]["bedside_lights"],
            "round1_new_socket_candidates": r1_counts[reference]["new_socket_candidates"],
            "round1_cabinet_power_zones": r1_counts[reference]["cabinet_power_zones"],
            "round1_socket_rechecks": r1_counts[reference]["socket_rechecks"],
            "unresolved_control_groups": rule["lighting_control_groups_min"],
            "unresolved_general_power_groups": rule["general_power_groups_min"],
            "unresolved_network_data_groups": rule["network_data_groups_min"],
            "unresolved_manual_dedicated_power": rule["manual_dedicated_power_min"],
            "program_complete": False,
        })
    return rooms


def svg_text(x: float, y: float, value: Any, css: str = "cell", anchor: str = "start") -> str:
    return f'<text class="{css}" x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}">{html.escape(str(value))}</text>'


def render_svg(report: dict[str, Any]) -> str:
    columns = [
        ("Ref", 10, 18), ("房间", 28, 46), ("控制", 79, 21), ("普通用电", 104, 27),
        ("模型设备", 135, 27), ("补充专用", 166, 27), ("弱电", 197, 21), ("柜体灯", 222, 24),
        ("现有灯", 250, 24), ("现有插座", 278, 28), ("交付参考", 310, 32), ("本轮候选", 346, 40),
    ]
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="500mm" height="400mm" viewBox="0 0 500 400">',
        '<style>@page{size:500mm 400mm;margin:0}text{font-family:Arial,"Noto Sans CJK SC",sans-serif}.title{font-size:5px;font-weight:700;fill:#102f43}.subtitle{font-size:2.6px;fill:#526777}.head{font-size:2.45px;font-weight:700;fill:#fff}.cell{font-size:2.35px;fill:#102f43}.num{font-size:2.5px;font-weight:700;fill:#102f43}.warn{font-size:2.5px;font-weight:700;fill:#c92a2a}.small{font-size:2.15px;fill:#526777}.grid{stroke:#b7c5ce;stroke-width:.35}.rowalt{fill:#f5f8fa}.header{fill:#284b63}.panel{fill:#fbfcfd;stroke:#102f43;stroke-width:.5}</style>',
        '<rect width="500" height="400" fill="#fff"/>',
        svg_text(10, 14, "E-302/E-304 逐房间用电与弱电功能程序候选", "title"),
        svg_text(10, 21, "最小功能组，不是面板/回路数｜只读候选｜非施工发布｜不写 IFC", "subtitle"),
        '<rect class="header" x="8" y="30" width="386" height="11"/>',
    ]
    for title, x, width in columns:
        parts.append(svg_text(x + width / 2, 37.2, title, "head", "middle"))
    row_height = 14.2
    y0 = 41.0
    for index, room in enumerate(report["rooms"]):
        y = y0 + index * row_height
        if index % 2:
            parts.append(f'<rect class="rowalt" x="8" y="{y:.2f}" width="386" height="{row_height:.2f}"/>')
        parts.append(f'<line class="grid" x1="8" y1="{y:.2f}" x2="394" y2="{y:.2f}"/>')
        control = room["lighting_control_groups_min"]
        power = room["general_power_groups_min"]
        equipment = room["modelled_equipment_power_demands"]
        manual = room["manual_dedicated_power_min"]
        network = room["network_data_groups_min"]
        cabinet = room["cabinet_light_feed_zones_min"]
        developer = room["developer_power_references"] + room["developer_control_references"] + room["developer_weak_references"]
        round1 = room["round1_bedside_lights"] + room["round1_new_socket_candidates"] + room["round1_cabinet_power_zones"] + room["round1_socket_rechecks"]
        values = [
            room["room_reference"], room["room_name"], control, power, equipment, manual, network, cabinet,
            room["current_light_instances"], room["current_socket_instances"], developer, round1,
        ]
        for (title, x, width), value in zip(columns, values):
            css = "cell" if title in {"Ref", "房间"} else "num"
            parts.append(svg_text(x + (1.2 if css == "cell" else width / 2), y + 8.8, value, css, "start" if css == "cell" else "middle"))
    table_bottom = y0 + len(report["rooms"]) * row_height
    parts.append(f'<line class="grid" x1="8" y1="{table_bottom:.2f}" x2="394" y2="{table_bottom:.2f}"/>')
    for _title, x, _width in columns:
        parts.append(f'<line class="grid" x1="{x:.2f}" y1="30" x2="{x:.2f}" y2="{table_bottom:.2f}"/>')
    parts.extend([
        '<rect class="panel" x="402" y="7" width="93" height="386"/>',
        svg_text(407, 16, "程序汇总", "title"),
        svg_text(407, 27, f'房间/Space：{report["summary"]["room_count"]}', "cell"),
        svg_text(407, 34, f'控制最小功能组：{report["summary"]["lighting_control_groups_min"]}', "cell"),
        svg_text(407, 41, f'普通用电最小功能组：{report["summary"]["general_power_groups_min"]}', "cell"),
        svg_text(407, 48, f'弱电/数据最小功能组：{report["summary"]["network_data_groups_min"]}', "cell"),
        svg_text(407, 55, f'模型设备供电需求：{report["summary"]["modelled_equipment_power_demands"]}', "cell"),
        svg_text(407, 62, f'手动补充专用需求：{report["summary"]["manual_dedicated_power_min"]}', "cell"),
        svg_text(407, 69, f'柜体灯供电协调区：{report["summary"]["cabinet_light_feed_zones_min"]}', "cell"),
        svg_text(407, 84, "当前正式 IFC：", "small"),
        svg_text(407, 91, "开关 0｜网络 0｜回路/端口 0", "warn"),
        svg_text(407, 105, "设备中心仅表示供电需求", "warn"),
        svg_text(407, 112, "不等于插座、出线口或回路位置", "warn"),
        svg_text(407, 126, "已确认设计规则：", "small"),
        svg_text(407, 133, "只用实体有线开关", "cell"),
        svg_text(407, 140, "卧室2个AP｜客厅/书房2个路由器", "cell"),
        svg_text(407, 147, "客厅/书房/餐厅：入户↔主卧门双控", "cell"),
        svg_text(407, 154, "面板底边：插座300｜电视600｜开关1300", "cell"),
        svg_text(407, 168, "下一道门：设备功率、柜体立面", "small"),
        svg_text(407, 175, "网络点位、端口拓扑与检修条件", "cell"),
        svg_text(407, 382, f'IFC SHA {report["source_ifc_sha256"][:12]}…', "small"),
        '</svg>',
    ])
    return "".join(parts)


def main() -> int:
    args = parse_args()
    ifc_hash = sha256(args.ifc)
    if ifc_hash != EXPECTED_IFC_SHA256:
        raise RuntimeError(f"formal IFC hash changed: {ifc_hash}")
    program = read_program(args.program)
    design_rules = read_design_rules(args.design_rules)
    design_rule_summary = compile_design_rule_summary(design_rules)
    delta = read_json(args.delta)
    elec = read_json(args.elec_existing)
    positioning = read_json(args.elec_positioning)
    round1 = read_json(args.round1)
    int1 = read_json(args.int1)
    route = read_json(args.route_readiness)
    hashes = {
        delta["source_ifc_sha256"], elec["source"]["sha256"], positioning["source_ifc_sha256"],
        round1["source_ifc_sha256"], int1["source"]["ifc_sha256"], route["source"]["ifc_sha256"],
    }
    if hashes != {ifc_hash}:
        raise RuntimeError(f"stale room-program inputs: {hashes}")
    equipment = typed_equipment_demands(elec)
    equipment.extend(proxy_equipment_demands(elec, positioning))
    equipment.extend(supplemental_equipment_demands(int1, route))
    equipment.sort(key=lambda row: (row["placement_room_reference"], row["equipment_role"], row["source_id"]))
    for index, row in enumerate(equipment, 1):
        row["demand_id"] = f"EP-{index:03d}"
    rooms = compile_rooms(program, delta, equipment, round1)
    summary = {
        "room_count": len(rooms),
        "lighting_control_groups_min": sum(row["lighting_control_groups_min"] for row in rooms),
        "general_power_groups_min": sum(row["general_power_groups_min"] for row in rooms),
        "network_data_groups_min": sum(row["network_data_groups_min"] for row in rooms),
        "manual_dedicated_power_min": sum(row["manual_dedicated_power_min"] for row in rooms),
        "cabinet_light_feed_zones_min": sum(row["cabinet_light_feed_zones_min"] for row in rooms),
        "modelled_equipment_power_demands": len(equipment),
        "formal_switch_instances": 0,
        "formal_network_instances": 0,
    }
    report = {
        "mode": "read_only_room_electrical_program_candidate",
        "source_ifc_sha256": ifc_hash,
        "program_path": str(args.program.resolve()),
        "design_rules_path": str(args.design_rules.resolve()),
        "summary": summary,
        "design_rule_summary": design_rule_summary,
        "design_rules": design_rules,
        "rooms": rooms,
        "equipment_power_demands": equipment,
        "gates": {
            "all_22_spaces_programmed": len(rooms) == 22,
            "all_equipment_demands_have_evidence": len(equipment) == 16 and all(row["basis"] for row in equipment),
            "developer_references_are_not_final_design": True,
            "confirmed_design_rules_are_complete": design_rule_summary == {
                "physical_wired_switch_only": True,
                "panel_bottom_AFF_mm": {
                    "ordinary_socket": 300,
                    "television_point": 600,
                    "physical_switch": 1300,
                },
                "bedroom_AP_count": 2,
                "living_study_router_count": 1,
                "paired_two_way_control_groups": 3,
                "entrance_master_lighting_switches": 1,
                "confirmed_island_dishwashers": 2,
                "demolition_walls_hidden_in_general_reviews": True,
            },
            "equipment_centres_are_not_socket_positions": all("only" in row["position_status"] or "pending" in row["position_status"] for row in equipment),
            "whole_home_switch_positioning_complete": False,
            "whole_home_network_positioning_complete": False,
            "automatic_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    args.output_svg.write_text(render_svg(report), encoding="utf-8")
    print(json.dumps({"summary": summary, "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
