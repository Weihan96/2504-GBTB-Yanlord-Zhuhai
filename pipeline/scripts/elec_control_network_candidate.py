#!/usr/bin/env python3
"""Generate read-only doorway control and room-level network coordination zones."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.placement

from svg_audit_underlay import validate_wall_plan_source
from sync_owner_inputs import DECISION_HEADERS, normalized_inputs, read_csv as read_owner_csv, validate_decisions

EXPECTED_IFC_SHA256 = "6c2fd8da9e9ad7ddbc2b63415a27f1c979e8995b880d8fce210a2dda2ef2aab6"
KITCHEN_FIRE_SENSOR_GLOBAL_ID = "2fwceKahvBqQXqal2ZcIUF"
SCALE_DENOMINATOR = 50.0
SVG_WORLD_OFFSET_MM = 10000.0


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument("--rules", type=Path, default=root / "pipeline/decisions/elec-design-rules.csv")
    parser.add_argument("--ceiling-devices", type=Path, default=root / "pipeline/decisions/a106-ceiling-device-review.csv")
    parser.add_argument("--doors", type=Path, default=root / "pipeline/decisions/a104-door-window-review.csv")
    parser.add_argument("--owner-decisions", type=Path, default=root / "pipeline/decisions/owner-input-register.csv")
    parser.add_argument("--spaces", type=Path, default=root / "pipeline/decisions/space-reference-review.csv")
    parser.add_argument("--router-evidence", type=Path, default=root / "build/elec/e304-router-cad-evidence.json")
    parser.add_argument("--source-svg", type=Path, default=root / "drawings/Wall Plan.svg")
    parser.add_argument("--output", type=Path, default=root / "build/elec/elec-control-network-candidate.json")
    parser.add_argument("--output-svg", type=Path, default=root / "drawings/E302-E304-control-network-candidate.svg")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def world_to_svg(position_mm: list[float]) -> tuple[float, float]:
    return (
        (position_mm[0] + SVG_WORLD_OFFSET_MM) / SCALE_DENOMINATOR,
        (SVG_WORLD_OFFSET_MM - position_mm[1]) / SCALE_DENOMINATOR,
    )


def device_position(row: dict[str, str]) -> list[float]:
    return [float(row[axis]) for axis in ("x_mm", "y_mm", "z_mm")]


def kitchen_fire_sensor(ifc_path: Path) -> dict[str, Any]:
    model = ifcopenshell.open(ifc_path)
    sensor = model.by_guid(KITCHEN_FIRE_SENSOR_GLOBAL_ID)
    if (
        sensor is None
        or not sensor.is_a("IfcSensor")
        or sensor.Name != "A106-FIRE-R04"
        or sensor.PredefinedType != "FIRESENSOR"
        or sensor.Representation is not None
    ):
        raise RuntimeError("formal kitchen fire sensor semantic gate failed")
    placement = ifcopenshell.util.placement.get_local_placement(sensor.ObjectPlacement)
    references = [
        relation.RelatingStructure.GlobalId
        for relation in model.by_type("IfcRelReferencedInSpatialStructure")
        if sensor in (relation.RelatedElements or ())
    ]
    if references != ["2fhEbDfK1EkhJwlPikNm$b"]:
        raise RuntimeError(f"formal kitchen fire sensor R04 reference drift: {references}")
    return {
        "global_id": sensor.GlobalId,
        "position_mm": [float(placement[axis][3]) for axis in range(3)],
        "predefined_type": sensor.PredefinedType,
        "representation": None,
        "space_global_id": references[0],
    }


def control_zones(doors: list[dict[str, str]]) -> list[dict[str, Any]]:
    door_by_id = {row["global_id"]: row for row in doors if row["ifc_class"] == "IfcDoor"}
    definitions = [
        ("CTRL-ENTRY", "2dxcvMre5F6fk90OTXSFJH", "R01", "入户门口", ["客厅", "书房", "餐厅"], True),
        ("CTRL-MASTER", "1TW6$_GfnABRZusYvx0zZG", "R10", "主卧门口", ["客厅", "书房", "餐厅"], False),
    ]
    rows = []
    for candidate_id, global_id, reference, location_name, groups, master in definitions:
        door = door_by_id[global_id]
        rows.append({
            "candidate_id": candidate_id,
            "kind": "doorway_control_coordination_zone",
            "room_reference": reference,
            "location_name": location_name,
            "source_global_ids": [global_id],
            "position_mm": [float(door["center_x_mm"]), float(door["center_y_mm"]), 1300.0],
            "vertical_datum": "panel_bottom_AFF_mm",
            "controlled_groups": groups,
            "entrance_master_lighting_switch": master,
            "control_method": "physical_wired_only",
            "position_basis": "confirmed doorway identity and door threshold centre; marker identifies the doorway coordination zone only",
            "coordinate_status": "doorway_zone_only_exact_jamb_side_pending",
            "confidence": 1.0,
            "review_required": True,
            "automatic_ifc_write_allowed": False,
        })
    return rows


def effective_owner_decisions(path: Path) -> dict[str, dict[str, str]]:
    rows = read_owner_csv(path, DECISION_HEADERS)
    errors = validate_decisions(rows, rows)
    if errors:
        raise RuntimeError("invalid owner decision register: " + "; ".join(errors))
    normalized = normalized_inputs(rows, [])["decisions"]
    return {row["input_id"]: row for row in normalized}


def control_wall_side_options(
    ifc_path: Path,
    controls: list[dict[str, Any]],
    owner_decisions: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    model = ifcopenshell.open(ifc_path)
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, False)
    result: list[dict[str, Any]] = []
    decision_ids = {
        "CTRL-ENTRY": "E302-ENTRY-SIDE",
        "CTRL-MASTER": "E302-MASTER-SIDE",
    }
    for control in controls:
        door = model.by_guid(control["source_global_ids"][0])
        shape = ifcopenshell.geom.create_shape(settings, door)
        vertices = list(shape.geometry.verts)
        local_x = vertices[0::3]
        local_width_mm = (max(local_x) - min(local_x)) * 1000.0
        if not 700.0 <= local_width_mm <= 1200.0:
            raise RuntimeError(
                f"unexpected door opening width for {control['candidate_id']}: {local_width_mm:.3f} mm"
            )
        placement = ifcopenshell.util.placement.get_local_placement(door.ObjectPlacement)
        axis_x, axis_y = float(placement[0][0]), float(placement[1][0])
        axis_length = (axis_x * axis_x + axis_y * axis_y) ** 0.5
        axis_x, axis_y = axis_x / axis_length, axis_y / axis_length
        offset_mm = local_width_mm / 2.0 + 150.0
        decision = owner_decisions[decision_ids[control["candidate_id"]]]
        candidate_suffix = "A" if " A" in decision["candidate_value"] else None
        effective_suffix = "A" if " A" in decision["effective_value"] else "B" if " B" in decision["effective_value"] else None
        for suffix, direction in (("A", -1.0), ("B", 1.0)):
            preferred = suffix == candidate_suffix
            if control["candidate_id"] == "CTRL-MASTER":
                review_status = "preferred_pending_a104" if preferred else "alternate_pending_a104"
            else:
                review_status = "candidate_preferred_unconfirmed" if preferred else "alternate_unconfirmed"
            result.append({
                "candidate_id": f"{control['candidate_id']}-{suffix}",
                "source_candidate_id": control["candidate_id"],
                "kind": "doorway_control_wall_side_option",
                "room_reference": control["room_reference"],
                "position_mm": [
                    round(control["position_mm"][0] + axis_x * offset_mm * direction, 6),
                    round(control["position_mm"][1] + axis_y * offset_mm * direction, 6),
                    control["position_mm"][2],
                ],
                "jamb_clearance_mm": 150.0,
                "review_status": review_status,
                "effective_selection": suffix == effective_suffix,
                "coordinate_status": "wall_side_option_not_final" if effective_suffix is None else "owner_selected_pending_geometry_evidence",
                "automatic_ifc_write_allowed": False,
            })
    return result


def network_zones(
    spaces: list[dict[str, str]], ceiling_devices: list[dict[str, str]], router_evidence: dict[str, Any]
) -> list[dict[str, Any]]:
    by_reference = {row["candidate_reference"]: row for row in spaces}
    by_candidate_id = {row["candidate_id"]: row for row in ceiling_devices}
    rows = []
    for candidate_id, reference, room_name in (
        ("A106-AP-R09", "R09", "主卧"),
        ("A106-AP-R14", "R14", "次卧"),
    ):
        candidate = by_candidate_id[candidate_id]
        space = by_reference[reference]
        rows.append({
            "candidate_id": candidate_id,
            "kind": "network_device_mechanical_position_candidate",
            "room_reference": reference,
            "room_name": room_name,
            "source_global_ids": [space["space_global_id"]],
            "position_mm": device_position(candidate),
            "network_role": "wireless_access_point",
            "wired_backhaul": True,
            "position_basis": candidate["position_basis"],
            "coordinate_status": "mechanical_position_candidate_product_power_data_review_pending",
            "confidence": float(candidate["confidence"]),
            "review_required": True,
            "automatic_ifc_write_allowed": False,
        })
    entry = by_reference["R01"]
    living = by_reference["R20"]
    study = by_reference["R22"]
    router_decision = router_evidence["router_decision"]
    weak_current_box = router_evidence["weak_current_box"]
    if router_decision["candidate_id"] != "NET-ROUTER-LIVING-STUDY" or router_decision["installation_z_mm"] is not None:
        raise RuntimeError("router CAD evidence must retain a confirmed plan position with installation Z pending")
    plan_position = [float(value) for value in router_decision["plan_position_mm"]]
    rows.append({
        "candidate_id": "NET-ROUTER-LIVING-STUDY",
        "kind": "network_device_entry_cabinet_plan_position",
        "room_reference": "R01",
        "served_room_references": router_decision["served_room_references"],
        "room_name": "玄关／过道高柜弱电箱柜位（服务客厅＋书房开放空间）",
        "source_global_ids": [entry["space_global_id"], living["space_global_id"], study["space_global_id"]],
        "plan_position_mm": plan_position,
        "installation_z_mm": None,
        "review_overlay_position_mm": [plan_position[0], plan_position[1], 3150.0],
        "cabinet_bbox_ifc_mm": weak_current_box["cabinet_bbox_ifc_mm"],
        "weak_current_box_bottom_aff_mm": weak_current_box["weak_current_box_bottom_aff_mm"],
        "network_role": "router_no_AP",
        "wired_backhaul": True,
        "position_basis": "user-confirmed entry-cabinet location; official E-2 weak-current plan LEADER #224271 maps through VIEWPORT #224238 to the right high-cabinet bay",
        "coordinate_status": "entry_cabinet_plan_xy_confirmed_router_z_product_power_data_thermal_service_pending",
        "source_evidence": {
            "official_plan_dwg_sha256": router_evidence["source"]["official_plan_dwg"]["sha256"],
            "verified_conversion_dxf_sha256": router_evidence["source"]["verified_conversion_dxf"]["sha256"],
            "text_handles": [row["handle"] for row in router_evidence["text_evidence"]],
            "leader_handle": router_evidence["leader_evidence"]["handle"],
            "viewport_handle": router_evidence["source"]["coordinate_transform"]["viewport_handle"],
        },
        "confidence": 1.0,
        "review_required": True,
        "automatic_ifc_write_allowed": False,
    })
    return rows


def safety_device_zones(
    spaces: list[dict[str, str]],
    ceiling_devices: list[dict[str, str]],
    fire_sensor: dict[str, Any],
) -> list[dict[str, Any]]:
    by_reference = {row["candidate_reference"]: row for row in spaces}
    by_candidate_id = {row["candidate_id"]: row for row in ceiling_devices}
    rows = []
    for candidate_id, reference, room_name in (
        ("A106-SMOKE-R09", "R09", "主卧"),
        ("A106-SMOKE-R14", "R14", "次卧"),
        ("A106-SMOKE-R20", "R20", "客厅"),
    ):
        space = by_reference[reference]
        candidate = by_candidate_id[candidate_id]
        rows.append({
            "candidate_id": candidate_id,
            "kind": "life_safety_device_mechanical_position_candidate",
            "room_reference": reference,
            "room_name": room_name,
            "source_global_ids": [space["space_global_id"]],
            "position_mm": device_position(candidate),
            "device_role": "smoke_alarm",
            "product_candidate": "Xiaomi smoke alarm with recessed mount",
            "recess_candidate": "opening_d105_counterbore_d140x2_embed_d30_mm",
            "approval_status": "fire_safety_device_approval_pending",
            "position_basis": candidate["position_basis"],
            "coordinate_status": "mechanical_position_candidate_product_supply_air_review_pending",
            "confidence": float(candidate["confidence"]),
            "review_required": True,
            "automatic_ifc_write_allowed": False,
            "mechanical_positioning_constraints_mm": {
                "minimum_wall_or_beam_clearance": 500,
                "minimum_unobstructed_radius": 500,
                "minimum_supply_air_edge_clearance": 1500,
                "minimum_perforated_supply_clearance": 500,
                "minimum_light_clearance": 200,
            },
            "constraint_basis": "GB 50116-2013 and DB4403/T 137-2021 10.4.3.3; final product and fire-safety review still govern",
        })
    kitchen_space = by_reference["R04"]
    kitchen_candidate = by_candidate_id["A106-FIRE-R04"]
    if device_position(kitchen_candidate) != fire_sensor["position_mm"]:
        raise RuntimeError("A-106 kitchen fire position differs from formal IFC")
    rows.append({
        "candidate_id": "A106-FIRE-R04",
        "kind": "formal_ifc_life_safety_position",
        "room_reference": "R04",
        "room_name": "中厨",
        "source_global_ids": [fire_sensor["global_id"], fire_sensor["space_global_id"]],
        "position_mm": fire_sensor["position_mm"],
        "device_role": "kitchen_fire_sensor",
        "product_candidate": None,
        "recess_candidate": None,
        "approval_status": "position_confirmed_type_product_pending",
        "position_basis": kitchen_candidate["position_basis"],
        "coordinate_status": "confirmed_ifc_position_type_product_power_communication_pending",
        "confidence": 1.0,
        "review_required": False,
        "automatic_ifc_write_allowed": False,
        "mechanical_positioning_constraints_mm": {
            "minimum_wall_or_beam_clearance": 500,
            "minimum_unobstructed_radius": 500,
            "minimum_light_clearance": 500,
        },
        "constraint_basis": "confirmed four-light-array centre and existing A-106 known-geometry checks; final type, product and manufacturer conditions pending",
    })
    rows.append({
        "candidate_id": "SAFE-GAS-KITCHEN",
        "kind": "life_safety_device_room_coordination_zone",
        "room_reference": "R04",
        "room_name": "中厨",
        "source_global_ids": [kitchen_space["space_global_id"]],
        "position_mm": [
            float(kitchen_space["centre_x_mm"]),
            float(kitchen_space["centre_y_mm"]),
            float(kitchen_space["centre_z_mm"]),
        ],
        "device_role": "combustible_gas_alarm",
        "product_candidate": "Xiaomi gas alarm with recessed mount",
        "recess_candidate": "opening_d94_counterbore_d120x2_embed_d43_mm",
        "approval_status": "gas_authority_model_approval_pending",
        "position_basis": "confirmed Space bbox centre used only to identify the room-level coordination zone",
        "coordinate_status": "room_zone_only_final_xy_z_pending",
        "confidence": 1.0,
        "review_required": True,
        "automatic_ifc_write_allowed": False,
        "mechanical_positioning_constraints_mm": None,
        "constraint_basis": "gas type, authority-approved device and manufacturer installation instructions pending",
    })
    return rows


def render_svg(source: str, report: dict[str, Any]) -> str:
    source = re.sub(r'width="400(?:\.0+)?mm"', 'width="500mm"', source, count=1)
    source = re.sub(r'viewBox="0 0 400(?:\.0+)? 400(?:\.0+)?"', 'viewBox="0 0 500 400"', source, count=1)
    if "</svg>" not in source or 'viewBox="0 0 500 400"' not in source:
        raise RuntimeError("could not prepare 500x400 control/network SVG")
    markup = []
    for index, row in enumerate(report["control_wall_side_options"]):
        x, y = world_to_svg(row["position_mm"])
        candidate_id = html.escape(row["candidate_id"])
        status = "PREFERRED" if "preferred" in row["review_status"] else "ALT"
        markup.append(
            f'<g data-candidate-id="{candidate_id}"><rect class="cn-control" x="{x-2.2:.3f}" y="{y-2.2:.3f}" width="4.4" height="4.4"/>'
            f'<text class="cn-label" x="{x+3.2:.3f}" y="{y+(-3 if index % 2 else 5):.3f}">{candidate_id} {status}</text></g>'
        )
    for index, row in enumerate(report["network_coordination_zones"]):
        position = row["plan_position_mm"] if "plan_position_mm" in row else row["position_mm"]
        x, y = world_to_svg(position)
        css = "cn-ap" if row["network_role"] == "wireless_access_point" else "cn-router"
        candidate_id = html.escape(row["candidate_id"])
        markup.append(
            f'<g data-candidate-id="{candidate_id}"><circle class="{css}" cx="{x:.3f}" cy="{y:.3f}" r="2.2"/>'
            f'<text class="cn-label" x="{x+3.4:.3f}" y="{y+(5 if index % 2 else -3):.3f}">{candidate_id}</text></g>'
        )
    for index, row in enumerate(report["safety_device_coordination_zones"]):
        x, y = world_to_svg(row["position_mm"])
        css = (
            "cn-gas"
            if row["device_role"] == "combustible_gas_alarm"
            else "cn-fire"
            if row["device_role"] == "kitchen_fire_sensor"
            else "cn-smoke"
        )
        candidate_id = html.escape(row["candidate_id"])
        markup.append(
            f'<g data-candidate-id="{candidate_id}"><circle class="{css}" cx="{x:.3f}" cy="{y:.3f}" r="2.2"/>'
            f'<text class="cn-label" x="{x+3.4:.3f}" y="{y+(5 if index % 2 else -3):.3f}">{candidate_id}</text></g>'
        )
    markup.extend([
        '<g><rect class="cn-panel" x="402" y="7" width="93" height="386"/>',
        '<text class="cn-title" x="407" y="16">控制、网络与安全设备协调区</text>',
        '<text class="cn-note" x="407" y="24">Bonsai 同批材质底图｜A-106 点位联动｜非施工发布</text>',
        '<text class="cn-text" x="407" y="39">洋红方块：4 个门口墙侧 A/B 候选</text>',
        '<text class="cn-text" x="407" y="47">蓝点：主卧/次卧 AP 机械候选点</text>',
        '<text class="cn-text" x="407" y="55">紫点：玄关高柜路由器平面柜位</text>',
        '<text class="cn-text" x="407" y="63">橙点：3 个烟感机械候选点</text>',
        '<text class="cn-text" x="407" y="71">深红点：厨房火灾正式 IFC 定位点</text>',
        '<text class="cn-text" x="407" y="79">红点：厨房燃气报警器房间区</text>',
        '<text class="cn-text" x="407" y="93">双控：客厅＋书房＋餐厅｜仅实体有线</text>',
        '<text class="cn-text" x="407" y="101">开关面板底边：1300 mm AFF</text>',
        '<text class="cn-warn" x="407" y="118">Master A 为待 A-104 取证优选；Entry A/B 未关闭</text>',
        '<text class="cn-warn" x="407" y="126">路由器 Z/产品/散热检修待设备深化</text>',
        '<text class="cn-warn" x="407" y="134">燃气型号待燃气公司确认，不锁定开孔</text>',
        f'<text class="cn-note" x="407" y="382">IFC SHA {report["source_ifc_sha256"][:12]}…</text></g>',
    ])
    style = """
@page{size:500mm 400mm;margin:0}.cn-control{fill:#d63384;stroke:#6b1742;stroke-width:.7}.cn-ap{fill:#228be6;stroke:#0b477d;stroke-width:.7}.cn-router{fill:#7048e8;stroke:#35206f;stroke-width:.7}.cn-smoke{fill:#f59f00;stroke:#7a4d00;stroke-width:.7}.cn-fire{fill:#8f1020;stroke:#4a0710;stroke-width:.7}.cn-gas{fill:#e03131;stroke:#751414;stroke-width:.7}.cn-label,.cn-title,.cn-note,.cn-text,.cn-warn{font-family:Arial,'Noto Sans CJK SC',sans-serif;fill:#102f43}.cn-label{font-size:2.2px;font-weight:700;paint-order:stroke;stroke:#fff;stroke-width:.8px}.cn-panel{fill:#fbfcfd;stroke:#102f43;stroke-width:.5}.cn-title{font-size:3.7px;font-weight:700}.cn-note{font-size:2.15px;fill:#526777}.cn-text{font-size:2.3px}.cn-warn{font-size:2.15px;fill:#c92a2a;font-weight:700}
"""
    return source.replace("</svg>", f'<style id="elec-control-network-style">{style}</style><g id="elec-control-network">{"".join(markup)}</g></svg>', 1)


def main() -> int:
    args = parse_args()
    source_hash = sha256(args.ifc)
    if source_hash != EXPECTED_IFC_SHA256:
        raise RuntimeError(f"formal IFC hash changed: {source_hash}")
    rules = read_csv(args.rules)
    controls = control_zones(read_csv(args.doors))
    owner_decisions = effective_owner_decisions(args.owner_decisions)
    control_options = control_wall_side_options(args.ifc, controls, owner_decisions)
    spaces = read_csv(args.spaces)
    ceiling_devices = read_csv(args.ceiling_devices)
    router_evidence = json.loads(args.router_evidence.read_text(encoding="utf-8"))
    if len(ceiling_devices) != 6 or len({row["candidate_id"] for row in ceiling_devices}) != 6:
        raise RuntimeError("A-106 ceiling-device register must contain six unique candidates")
    fire_sensor = kitchen_fire_sensor(args.ifc)
    networks = network_zones(spaces, ceiling_devices, router_evidence)
    safety_devices = safety_device_zones(spaces, ceiling_devices, fire_sensor)
    router_zones = [row for row in networks if row["network_role"] == "router_no_AP"]
    ap_candidates = [row for row in networks if row["network_role"] == "wireless_access_point"]
    smoke_zones = [row for row in safety_devices if row["device_role"] == "smoke_alarm"]
    fire_positions = [row for row in safety_devices if row["device_role"] == "kitchen_fire_sensor"]
    gas_zones = [row for row in safety_devices if row["device_role"] == "combustible_gas_alarm"]
    report = {
        "mode": "read_only_control_network_coordination_candidate",
        "source_ifc_sha256": source_hash,
        "rules_path": str(args.rules.resolve()),
        "owner_decisions": {
            "path": str(args.owner_decisions.resolve()),
            "sha256": sha256(args.owner_decisions),
            "effective": {
                input_id: owner_decisions[input_id]["effective_value"]
                for input_id in ("E302-ENTRY-SIDE", "E302-ENTRY-PANEL", "E302-MASTER-SIDE", "E302-MASTER-PANEL")
            },
        },
        "router_evidence_path": str(args.router_evidence.resolve()),
        "summary": {
            "control_coordination_zones": len(controls),
            "control_wall_side_options": len(control_options),
            "paired_two_way_control_groups": 3,
            "entrance_master_lighting_switches": 1,
            "bedroom_AP_candidates": len(ap_candidates),
            "living_study_router_zones": len(router_zones),
            "smoke_alarm_candidates": len(smoke_zones),
            "kitchen_fire_sensor_positions": len(fire_positions),
            "kitchen_gas_alarm_room_zones": len(gas_zones),
            "confirmed_rules": sum(row["status"] == "confirmed" for row in rules),
        },
        "control_coordination_zones": controls,
        "control_wall_side_options": control_options,
        "network_coordination_zones": networks,
        "safety_device_coordination_zones": safety_devices,
        "gates": {
            "two_doorway_zones_present": len(controls) == 2,
            "four_wall_side_options_present": len(control_options) == 4,
            "master_a_preferred_pending_a104": next(
                row for row in control_options if row["candidate_id"] == "CTRL-MASTER-A"
            )["review_status"] == "preferred_pending_a104",
            "entry_wall_side_not_auto_closed": all(
                not row["effective_selection"]
                for row in control_options if row["source_candidate_id"] == "CTRL-ENTRY"
            ),
            "three_two_way_groups_present": all(row["controlled_groups"] == ["客厅", "书房", "餐厅"] for row in controls),
            "entrance_master_switch_present": sum(row["entrance_master_lighting_switch"] for row in controls) == 1,
            "two_bedroom_AP_candidates_present": len(ap_candidates) == 2
            and {row["candidate_id"] for row in ap_candidates} == {"A106-AP-R09", "A106-AP-R14"},
            "one_shared_router_no_AP_zone_present": len(router_zones) == 1
            and router_zones[0]["room_reference"] == "R01"
            and router_zones[0]["served_room_references"] == ["R20", "R22"],
            "router_entry_cabinet_plan_position_verified": len(router_zones) == 1
            and router_zones[0]["plan_position_mm"] == router_evidence["router_decision"]["plan_position_mm"]
            and router_zones[0]["installation_z_mm"] is None
            and router_zones[0]["weak_current_box_bottom_aff_mm"] == 350.0,
            "three_smoke_candidates_present": len(smoke_zones) == 3 and {row["candidate_id"] for row in smoke_zones} == {"A106-SMOKE-R09", "A106-SMOKE-R14", "A106-SMOKE-R20"},
            "smoke_positioning_constraints_present": all(row["mechanical_positioning_constraints_mm"] == {
                "minimum_wall_or_beam_clearance": 500,
                "minimum_unobstructed_radius": 500,
                "minimum_supply_air_edge_clearance": 1500,
                "minimum_perforated_supply_clearance": 500,
                "minimum_light_clearance": 200,
            } for row in smoke_zones),
            "one_formal_kitchen_fire_position_present": len(fire_positions) == 1
            and fire_positions[0]["source_global_ids"][0] == KITCHEN_FIRE_SENSOR_GLOBAL_ID
            and fire_positions[0]["position_mm"] == [1800.0, -4576.0, 2400.0],
            "one_kitchen_gas_alarm_zone_present": len(gas_zones) == 1 and gas_zones[0]["room_reference"] == "R04",
            "gas_alarm_model_authority_pending": len(gas_zones) == 1 and gas_zones[0]["approval_status"] == "gas_authority_model_approval_pending",
            "a106_exact_positions_imported": all(
                "candidate" in row["coordinate_status"] or "confirmed_ifc_position" in row["coordinate_status"]
                for row in ap_candidates + smoke_zones + fire_positions
            ),
            "unconfirmed_positions_remain_coordination_zones": all(
                "zone_only" in row["coordinate_status"] for row in controls + gas_zones
            ) and all("pending" in row["coordinate_status"] for row in router_zones),
            "automatic_ifc_write_allowed": False,
            "construction_release_ready": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    source = args.source_svg.read_text(encoding="utf-8")
    validate_wall_plan_source(source, args.source_svg, args.ifc)
    args.output_svg.write_text(render_svg(source, report), encoding="utf-8")
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
