#!/usr/bin/env python3
"""Compile the first read-only renovation electrical demand candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
from pathlib import Path
from typing import Any, Callable

from svg_audit_underlay import validate_wall_plan_source
from sync_owner_inputs import (
    APPLIANCE_HEADERS,
    DECISION_HEADERS,
    normalized_inputs,
    read_csv as read_owner_csv,
    validate_appliances,
    validate_decisions,
)
from equipment_ssot import appliance_projection_rows, validate as validate_equipment_ssot

SCALE_DENOMINATOR = 50.0
SVG_WORLD_OFFSET_MM = 10000.0
LOAD_SCENARIO_HEADERS = [
    "scenario_id",
    "socket_id",
    "owner_input_id",
    "scenario_label",
    "appliance_ids",
]
CIRCUIT_DESIGN_HEADERS = [
    "circuit_design_id", "socket_id", "design_basis", "design_load_w",
    "design_current_a", "connection_positions", "total_branch_circuits",
    "dedicated_branch_circuits", "auxiliary_branch_circuits",
    "breaker_curve_rating", "conductor_cross_section_mm2", "rcbo_required",
    "wall_socket_strategy", "product_power_status", "status", "source_id", "notes",
]
E303_CIRCUIT_INPUT_BY_SOCKET = {
    "NS-01": "E303-NS01-CIRCUIT",
    "NS-02": "E303-NS02-CIRCUIT",
}
USE_CONFIRMED_STATUSES = {"采用候选", "自定义确认"}
PLANNING_VOLTAGE_V = 220.0
SINGLE_SOCKET_CLASS_CURRENT_A = 16.0
SINGLE_SOCKET_CLASS_CAPACITY_W = PLANNING_VOLTAGE_V * SINGLE_SOCKET_CLASS_CURRENT_A


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument("--expected-ifc-sha256", help="Optional caller-frozen formal IFC hash")
    parser.add_argument(
        "--requirements",
        type=Path,
        default=root / "pipeline/decisions/elec-renovation-requirements.csv",
    )
    parser.add_argument(
        "--elec-existing",
        type=Path,
        default=root / "build/elec/elec-existing-candidate.json",
    )
    parser.add_argument(
        "--elec-positioning",
        type=Path,
        default=root / "build/elec/elec-positioning-candidate.json",
    )
    parser.add_argument(
        "--int1",
        type=Path,
        default=root / "build/int1/int1-existing-report.json",
    )
    parser.add_argument("--appliances", type=Path)
    parser.add_argument(
        "--owner-decisions",
        type=Path,
        default=root / "pipeline/decisions/owner-input-register.csv",
    )
    parser.add_argument(
        "--load-scenarios",
        type=Path,
        default=root / "pipeline/decisions/e303-load-scenarios.csv",
    )
    parser.add_argument(
        "--circuit-design",
        type=Path,
        default=root / "pipeline/decisions/e303-circuit-design.csv",
    )
    parser.add_argument("--source-svg", type=Path, default=root / "drawings/Wall Plan.svg")
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "build/elec/elec-renovation-round1-candidate.json",
    )
    parser.add_argument(
        "--output-svg",
        type=Path,
        default=root / "drawings/E301-E303-renovation-round1-candidate.svg",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_requirements(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        records = list(csv.DictReader(handle))
    ids = [row["requirement_id"] for row in records]
    if len(ids) != len(set(ids)):
        raise RuntimeError("renovation electrical requirement IDs are not unique")
    return records


def record_by_global_id(records: list[dict[str, Any]], global_id: str) -> dict[str, Any]:
    matches = [row for row in records if row.get("global_id") == global_id]
    if len(matches) != 1:
        raise RuntimeError(f"expected one INT1 record for {global_id}, got {len(matches)}")
    return matches[0]


def bbox_centre(record: dict[str, Any]) -> list[float]:
    return [
        (float(record["bbox_min_mm"][axis]) + float(record["bbox_max_mm"][axis])) / 2.0
        for axis in range(3)
    ]


def rounded(values: list[float]) -> list[float]:
    return [round(float(value), 6) for value in values]


def bedside_candidates(int1_records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    main_bed = record_by_global_id(int1_records, "3IQBEqO5vDI8Z9k1Ltge_N")
    main_tables = [
        record_by_global_id(int1_records, "3eic1dzkn5heTIn4PhF37v"),
        record_by_global_id(int1_records, "3hA0vKpcn44u4Tsx4tqiUz"),
    ]
    guest_bed = record_by_global_id(int1_records, "1i_pqgLv9A7uuV7MjaArBW")
    guest_table = record_by_global_id(int1_records, "1HZoxe$df4cBb4UXXH5J2S")

    main_head_x = round(float(main_bed["bbox_min_mm"][0]) / 100.0) * 100.0
    guest_head_x = round(float(guest_bed["bbox_min_mm"][0]) / 100.0) * 100.0
    guest_centre_y = bbox_centre(guest_bed)[1]
    guest_table_y = bbox_centre(guest_table)[1]
    guest_mirror_y = 2.0 * guest_centre_y - guest_table_y

    rows = []
    for index, table in enumerate(main_tables, 1):
        rows.append(
            {
                "candidate_id": f"BL-{index:02d}",
                "kind": "bedside_wall_light",
                "room_reference": "R09",
                "room_name": "主卧",
                "position_mm": rounded([main_head_x, bbox_centre(table)[1], 1200.0]),
                "position_basis": "head-wall face rounded from BED01 plus existing BST01/BST02 centreline",
                "source_global_ids": [main_bed["global_id"], table["global_id"]],
                "confidence": 0.85,
                "review_required": True,
                "height_status": "1200 mm AFF review candidate; final fixture and reading ergonomics pending",
                "automatic_ifc_write_allowed": False,
            }
        )
    rows.extend(
        [
            {
                "candidate_id": "BL-03",
                "kind": "bedside_wall_light",
                "room_reference": "R14",
                "room_name": "次卧",
                "position_mm": rounded([guest_head_x, guest_table_y, 1200.0]),
                "position_basis": "head-wall face rounded from BED02 plus existing BST03 centreline",
                "source_global_ids": [guest_bed["global_id"], guest_table["global_id"]],
                "confidence": 0.80,
                "review_required": True,
                "height_status": "1200 mm AFF review candidate; final fixture and reading ergonomics pending",
                "automatic_ifc_write_allowed": False,
            },
            {
                "candidate_id": "BL-04",
                "kind": "bedside_wall_light",
                "room_reference": "R14",
                "room_name": "次卧",
                "position_mm": rounded([guest_head_x, guest_mirror_y, 1200.0]),
                "position_basis": "mirror of existing BST03 centreline about the BED02 plan centreline",
                "source_global_ids": [guest_bed["global_id"], guest_table["global_id"]],
                "confidence": 0.65,
                "review_required": True,
                "height_status": "1200 mm AFF review candidate; missing second bedside table and wall segment require later reverse check",
                "automatic_ifc_write_allowed": False,
            },
        ]
    )
    return rows


def e303_circuit_decisions(path: Path) -> dict[str, dict[str, str]]:
    rows = read_owner_csv(path, DECISION_HEADERS)
    errors = validate_decisions(rows, rows)
    if errors:
        raise RuntimeError("invalid owner decision register: " + "; ".join(errors))
    by_id = {row["input_id"]: row for row in rows}
    missing = sorted(set(E303_CIRCUIT_INPUT_BY_SOCKET.values()) - set(by_id))
    if missing:
        raise RuntimeError(f"owner decision register missing E-303 circuit inputs: {missing}")
    return {
        input_id: by_id[input_id]
        for input_id in E303_CIRCUIT_INPUT_BY_SOCKET.values()
    }


def read_load_scenarios(path: Path) -> list[dict[str, Any]]:
    rows = read_owner_csv(path, LOAD_SCENARIO_HEADERS)
    scenario_ids = [row["scenario_id"] for row in rows]
    if not rows or any(not value for value in scenario_ids):
        raise RuntimeError("E-303 load scenarios require non-blank scenario IDs")
    if len(scenario_ids) != len(set(scenario_ids)):
        raise RuntimeError("E-303 load scenario IDs are not unique")
    scenarios: list[dict[str, Any]] = []
    for row in rows:
        appliance_ids = [value.strip() for value in row["appliance_ids"].split(";") if value.strip()]
        if not row["scenario_label"] or not appliance_ids:
            raise RuntimeError(f"{row['scenario_id']}: scenario label and appliance IDs are required")
        if len(appliance_ids) != len(set(appliance_ids)):
            raise RuntimeError(f"{row['scenario_id']}: duplicate appliance IDs")
        expected_input_id = E303_CIRCUIT_INPUT_BY_SOCKET.get(row["socket_id"])
        if expected_input_id is None or row["owner_input_id"] != expected_input_id:
            raise RuntimeError(
                f"{row['scenario_id']}: socket and owner-input semantic identity mismatch"
            )
        scenarios.append({**row, "appliance_ids": appliance_ids})
    if {row["socket_id"] for row in scenarios} != set(E303_CIRCUIT_INPUT_BY_SOCKET):
        raise RuntimeError("E-303 load scenarios must cover NS-01 and NS-02")
    return scenarios


def read_circuit_design(path: Path) -> dict[str, dict[str, str]]:
    rows = read_owner_csv(path, CIRCUIT_DESIGN_HEADERS)
    by_socket = {row["socket_id"]: row for row in rows}
    if set(by_socket) != {"NS-01", "NS-02"} or len(rows) != 2:
        raise RuntimeError("E-303 circuit design must contain exactly NS-01 and NS-02")
    for row in rows:
        if row["status"] != "research_conclusion":
            raise RuntimeError(f"{row['circuit_design_id']}: circuit design is not adopted")
        if row["rcbo_required"] != "yes" or row["breaker_curve_rating"] != "C16":
            raise RuntimeError(f"{row['circuit_design_id']}: C16 RCBO invariant failed")
    return by_socket


def appliance_socket_context(
    appliance_input: Path | list[dict[str, str]],
    load_scenario_path: Path,
    circuit_decisions: dict[str, dict[str, str]],
    circuit_design_path: Path,
) -> dict[str, dict[str, Any]]:
    rows = read_owner_csv(appliance_input, APPLIANCE_HEADERS) if isinstance(appliance_input, Path) else appliance_input
    errors = validate_appliances(rows, rows)
    if errors:
        raise RuntimeError("invalid appliance input register: " + "; ".join(errors))
    normalized_by_id = {
        row["appliance_id"]: row
        for row in normalized_inputs([], rows)["appliances"]
    }
    appliance_by_id = {row["appliance_id"]: row for row in rows}
    scenarios = read_load_scenarios(load_scenario_path)
    circuit_design = read_circuit_design(circuit_design_path)
    contexts: dict[str, dict[str, Any]] = {}
    for socket_id in ("NS-01", "NS-02"):
        selected = [
            row
            for row in rows
            if socket_id
            in normalized_by_id[row["appliance_id"]]["effective"]["use_location_confirmed"]
        ]
        known_load_w = 0.0
        items = []
        for row in selected:
            effective = normalized_by_id[row["appliance_id"]]["effective"]
            if effective["rated_power_w"]:
                known_load_w += float(effective["rated_power_w"]) * int(effective["quantity"])
            group = row["simultaneous_group"] or "UNASSIGNED"
            items.append({
                "appliance_id": row["appliance_id"],
                "appliance_name": row["appliance_name"],
                "status": row["status"],
                "effective_use_location": effective["use_location_confirmed"],
                "effective_rated_power_w": effective["rated_power_w"],
                "candidate_power_range_w": None,
                "simultaneous_group": group,
            })
        selected_ids = {row["appliance_id"] for row in selected}
        circuit_input_id = E303_CIRCUIT_INPUT_BY_SOCKET[socket_id]
        circuit_decision = circuit_decisions[circuit_input_id]
        use_confirmed = circuit_decision["status"] in USE_CONFIRMED_STATUSES
        design = circuit_design[socket_id]
        socket_scenarios = [row for row in scenarios if row["socket_id"] == socket_id]
        scenario_loads = []
        for scenario in socket_scenarios:
            unknown_ids = sorted(set(scenario["appliance_ids"]) - selected_ids)
            if unknown_ids:
                raise RuntimeError(
                    f"{scenario['scenario_id']}: invalid scenario appliances; "
                    f"not_at_socket={unknown_ids}"
                )
            scenario_loads.append({
                **scenario,
                "appliance_names": [
                    appliance_by_id[appliance_id]["appliance_name"]
                    for appliance_id in scenario["appliance_ids"]
                ],
                "load_status": "product_nameplate_values_may_remain_unknown",
            })
        if set().union(*(set(row["appliance_ids"]) for row in scenario_loads)) != selected_ids:
            raise RuntimeError(f"{socket_id}: explicit scenarios do not cover its appliance use list")
        semantic_gates = {
            "use_confirmed": use_confirmed,
            "planning_envelope_compiled": True,
            "product_and_circuit_fixed": False,
        }
        contexts[socket_id] = {
            "items": items,
            "circuit_owner_input": {
                "input_id": circuit_input_id,
                "status": circuit_decision["status"],
            },
            "explicit_load_scenarios": scenario_loads,
            "known_connected_load_w": known_load_w,
            "candidate_load_ranges_by_simultaneous_group": {},
            "planning_envelope": {
                "design_load_w": float(design["design_load_w"]),
                "design_current_a": float(design["design_current_a"]),
                "design_basis": design["design_basis"],
                "planning_voltage_v": PLANNING_VOLTAGE_V,
                "total_branch_circuits": int(design["total_branch_circuits"]),
                "dedicated_branch_circuits": int(design["dedicated_branch_circuits"]),
                "auxiliary_branch_circuits": int(design["auxiliary_branch_circuits"]),
                "minimum_independent_circuit_count_candidate": int(design["total_branch_circuits"]),
                "minimum_connection_positions_candidate": int(design["connection_positions"]),
                "breaker_curve_rating": design["breaker_curve_rating"],
                "conductor_cross_section_mm2": float(design["conductor_cross_section_mm2"]),
                "rcbo_required": design["rcbo_required"] == "yes",
                "wall_socket_strategy": design["wall_socket_strategy"],
                "product_power_status": design["product_power_status"],
                "source_id": design["source_id"],
                "calculation_status": "research_conclusion_not_product_nameplate_load",
                "final_conductor_protection_rcd_pending": False,
            },
            "semantic_gates": semantic_gates,
            "candidate_ranges_are_not_confirmed_loads": True,
            "socket_form_and_circuit_sizing_ready": True,
        }
    return contexts


def new_socket_candidates(
    int1_records: list[dict[str, Any]], appliance_context: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    island_end = next(
        row
        for row in int1_records
        if row.get("type_name") == "SB02" and row.get("object_name") == "Island"
    )
    island_y = bbox_centre(island_end)[1]
    return [
        {
            "candidate_id": "NS-01",
            "kind": "new_socket",
            "candidate_role": "island_end_panel_socket",
            "room_reference": "R03",
            "room_name": "西厨",
            "position_mm": rounded([float(island_end["bbox_min_mm"][0]), island_y, 650.0]),
            "position_basis": "outer face centreline of the existing SB02 island end panel",
            "source_global_ids": [island_end["global_id"]],
            "confidence": 0.70,
            "review_required": True,
            "height_status": "650 mm AFF review candidate; socket type, splash protection and final elevation pending",
            "socket_form_candidate": "three concealed covered socket positions, one on each island face below the countertop/table overlap; final product and panel cut-out pending",
            "circuit_strategy_candidate": (
                "two independent C16 RCBO socket circuits; separate the hot-pot load from the second simultaneous mobile appliance"
                if appliance_context["NS-01"]["semantic_gates"]["use_confirmed"]
                else "no circuit-count candidate; APP-001 through APP-003 exact rated powers and owner use scenario remain pending"
            ),
            "automatic_ifc_write_allowed": False,
            "appliance_context": appliance_context["NS-01"],
        },
        {
            "candidate_id": "NS-02",
            "kind": "new_socket",
            "candidate_role": "dining_bay_socket",
            "room_reference": "R06",
            "room_name": "餐厅飘窗",
            "position_mm": [-1250.0, -4400.0, 300.0],
            "vertical_datum": "panel_bottom_AFF_mm",
            "position_basis": "R06 plan centreline on the room-side face of the existing bay wall band",
            "source_global_ids": ["0vMU_9TZXEV8GPr46AYxhu", "34MraFd1nFsep2QMO2kY48"],
            "confidence": 0.75,
            "review_required": True,
            "height_status": "panel bottom 300 mm AFF confirmed general datum; intended appliance and wall-side position remain pending",
            "socket_form_candidate": "minimal linear track socket or custom concealed flip-up assembly; product, module count and cabinet detail pending",
            "circuit_strategy_candidate": (
                "one dedicated C16 RCBO coffee circuit plus one auxiliary small-appliance circuit; final coffee faceplate follows actual plug"
                if appliance_context["NS-02"]["semantic_gates"]["use_confirmed"]
                else "no circuit-count candidate; final coffee machine, grinder and kettle rated powers remain pending"
            ),
            "automatic_ifc_write_allowed": False,
            "appliance_context": appliance_context["NS-02"],
        },
    ]


def union_bbox(records: list[dict[str, Any]]) -> tuple[list[float], list[float]]:
    if not records:
        raise RuntimeError("cabinet power zone has no source objects")
    minimum = [min(float(row["bbox_min_mm"][axis]) for row in records) for axis in range(3)]
    maximum = [max(float(row["bbox_max_mm"][axis]) for row in records) for axis in range(3)]
    return minimum, maximum


def cabinet_power_zones(int1_records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    definitions: list[tuple[str, str, str, Callable[[dict[str, Any]], bool]]] = [
        ("CP-01", "R04", "中厨柜体", lambda row: row.get("container") == "KITCHEN" and row.get("installation_role") == "fixed_furniture"),
        ("CP-02", "R03", "西厨/岛台柜体", lambda row: row.get("container") == "VVD" and row.get("installation_role") == "fixed_furniture"),
        ("CP-03", "R09", "主卧衣柜", lambda row: row.get("type_name") in {"WD01", "WD01L"}),
        ("CP-04", "R14", "次卧衣柜", lambda row: row.get("type_name") == "WD03"),
        ("CP-05", "R20", "客厅书架", lambda row: row.get("type_name") == "SHE01"),
        ("CP-06", "R20", "客厅高柜", lambda row: row.get("type_name") == "WD02"),
        ("CP-07", "R12", "主卫固定柜", lambda row: row.get("type_name") == "HIMA01"),
    ]
    rows = []
    for candidate_id, reference, name, predicate in definitions:
        sources = [row for row in int1_records if row.get("bbox_min_mm") and predicate(row)]
        minimum, maximum = union_bbox(sources)
        marker = [(minimum[0] + maximum[0]) / 2.0, (minimum[1] + maximum[1]) / 2.0, 1200.0]
        rows.append(
            {
                "candidate_id": candidate_id,
                "kind": "cabinet_lighting_power_zone",
                "room_reference": reference,
                "room_name": name,
                "marker_position_mm": rounded(marker),
                "bbox_min_mm": rounded(minimum),
                "bbox_max_mm": rounded(maximum),
                "source_global_ids": sorted(row["global_id"] for row in sources),
                "position_basis": "union world bbox of the existing fixed cabinet assembly; marker denotes coordination scope, not a final cable outlet",
                "confidence": 0.90,
                "review_required": True,
                "coordinate_status": "assembly_zone_only",
                "automatic_ifc_write_allowed": False,
            }
        )
    return rows


def kitchen_socket_rechecks(
    existing: dict[str, Any], positioning: dict[str, Any]
) -> list[dict[str, Any]]:
    roles = {row["candidate_id"]: row for row in positioning["socket_candidates"]}
    room_references = {"餐厅": "R07", "西厨": "R03", "中厨": "R04", "中厨飘窗": "R05"}
    rows = []
    for socket in existing["sheets"]["E-303"]["sockets"]:
        role = roles[socket["candidate_id"]]
        rows.append(
            {
                "candidate_id": socket["candidate_id"],
                "kind": "existing_kitchen_socket_recheck",
                "global_id": socket["global_id"],
                "room_reference": room_references[socket["candidate_space"]["long_name"]],
                "room_name": socket["candidate_space"]["long_name"],
                "position_mm": rounded(socket["bbox"]["centre_mm"]),
                "candidate_role": role["candidate_role"],
                "position_basis": "current IFC socket centre retained only as a recheck reference against final appliance and cabinet elevations",
                "confidence": 1.0,
                "review_required": True,
                "coordinate_status": "current_position_not_frozen",
                "automatic_ifc_write_allowed": False,
            }
        )
    return rows


def world_to_svg(position_mm: list[float]) -> tuple[float, float]:
    return (
        (position_mm[0] + SVG_WORLD_OFFSET_MM) / SCALE_DENOMINATOR,
        (SVG_WORLD_OFFSET_MM - position_mm[1]) / SCALE_DENOMINATOR,
    )


def inject_svg(source: str, generated: str, style: str) -> str:
    source = re.sub(r'width="400(?:\.0+)?mm"', 'width="500mm"', source, count=1)
    source = re.sub(r'viewBox="0 0 400(?:\.0+)? 400(?:\.0+)?"', 'viewBox="0 0 500 400"', source, count=1)
    if "</svg>" not in source or 'viewBox="0 0 500 400"' not in source:
        raise RuntimeError("could not prepare the 500x400 renovation electrical SVG")
    return source.replace(
        "</svg>",
        f'<style id="elec-renovation-round1-style">{style}</style><g id="elec-renovation-round1">{generated}</g></svg>',
        1,
    )


def label_layout(
    items: list[tuple[dict[str, Any], str, str]]
) -> tuple[dict[str, tuple[float, float]], int]:
    occupied: list[tuple[float, float, float, float]] = []
    positions: dict[str, tuple[float, float]] = {}
    collision_count = 0
    offsets = [(3, -2), (3, 5), (-11, -2), (-11, 5), (0, -7), (0, 9), (9, -7), (-16, -7)]
    for row, _css, position_key in items:
        x, y = world_to_svg(row[position_key])
        width = max(6.0, len(row["candidate_id"]) * 1.45)
        chosen = None
        for dx, dy in offsets:
            tx, ty = x + dx, y + dy
            box = (tx - 0.5, ty - 2.5, tx + width, ty + 0.5)
            if any(
                not (box[2] + 0.5 <= other[0] or other[2] + 0.5 <= box[0] or box[3] + 0.5 <= other[1] or other[3] + 0.5 <= box[1])
                for other in occupied
            ):
                continue
            chosen = (tx, ty, box)
            break
        if chosen is None:
            collision_count += 1
            tx, ty = x + offsets[-1][0], y + offsets[-1][1]
            chosen = (tx, ty, (tx - 0.5, ty - 2.5, tx + width, ty + 0.5))
        positions[row["candidate_id"]] = (chosen[0], chosen[1])
        occupied.append(chosen[2])
    return positions, collision_count


def svg_marker(
    row: dict[str, Any], css: str, labels: dict[str, tuple[float, float]], position_key: str = "position_mm"
) -> str:
    x, y = world_to_svg(row[position_key])
    candidate_id = html.escape(row["candidate_id"])
    label_x, label_y = labels[row["candidate_id"]]
    return (
        f'<g data-candidate-id="{candidate_id}" data-elec-kind="{html.escape(row["kind"])}">'
        f'<circle class="{css}" cx="{x:.3f}" cy="{y:.3f}" r="2.0"/>'
        f'<text class="r1-label" x="{label_x:.3f}" y="{label_y:.3f}">{candidate_id}</text>'
        f'<title>{candidate_id} | {html.escape(row.get("room_name", ""))}</title></g>'
    )


def render_svg(source: str, report: dict[str, Any]) -> str:
    items = render_items(report)
    labels, _collision_count = label_layout(items)
    markup = []
    markup.extend(svg_marker(row, css, labels, position_key) for row, css, position_key in items)
    summary = report["summary"]
    markup.append('<g><rect class="r1-panel" x="402" y="7" width="93" height="386"/>')
    markup.append('<text class="r1-title" x="407" y="16">E-301/E-303 装修用电深化第一轮</text>')
    markup.append('<text class="r1-note" x="407" y="24">Bonsai 同批材质底图｜只读候选｜不写 IFC</text>')
    markup.append(f'<text class="r1-text" x="407" y="37">青色床头灯候选：{summary["bedside_light_candidates"]}</text>')
    markup.append(f'<text class="r1-text" x="407" y="44">绿色新增插座候选：{summary["new_socket_candidates"]}</text>')
    markup.append(f'<text class="r1-text" x="407" y="51">琥珀色柜体供电协调区：{summary["cabinet_power_zones"]}</text>')
    markup.append(f'<text class="r1-text" x="407" y="58">红色厨房现有插座复核：{summary["kitchen_socket_rechecks"]}</text>')
    ns01 = next(row for row in report["new_socket_candidates"] if row["candidate_id"] == "NS-01")
    ns02 = next(row for row in report["new_socket_candidates"] if row["candidate_id"] == "NS-02")
    ns01_plan = ns01["appliance_context"]["planning_envelope"]
    ns02_plan = ns02["appliance_context"]["planning_envelope"]
    circuit_candidates_authorized = report["gates"]["use_confirmed"]
    heading = "E-303 容量研究结论｜220V / C16 RCBO"
    ns01_circuit = f'{ns01_plan["total_branch_circuits"]} 回路已闭合'
    ns02_circuit = f'{ns02_plan["total_branch_circuits"]} 回路（1咖啡专用+1辅助）'
    markup.append(f'<text class="r1-heading" x="407" y="72">{heading}</text>')
    markup.append(f'<text class="r1-text" x="407" y="81">NS-01 两个 10A 端口同时：{ns01_plan["design_current_a"]:.0f}A / {ns01_plan["design_load_w"] / 1000:.1f}kW</text>')
    markup.append(f'<text class="r1-text" x="407" y="88">三面各 1 个隐藏盖板位｜{ns01_circuit}</text>')
    markup.append(f'<text class="r1-text" x="407" y="99">NS-02 咖啡候选上限：{ns02_plan["design_current_a"]:.2f}A / {ns02_plan["design_load_w"] / 1000:.1f}kW</text>')
    markup.append(f'<text class="r1-text" x="407" y="106">线性轨道或自制翻盖｜{ns02_circuit}</text>')
    markup.append('<text class="r1-warn" x="407" y="120">设计容量 ≠ 未购设备铭牌功率</text>')
    markup.append('<text class="r1-warn" x="407" y="127">NS-02 最终面板按实购插头匹配</text>')
    markup.append('<text class="r1-note" x="407" y="141">柜体标记是供电范围，不是最终出线口</text>')
    markup.append('<text class="r1-note" x="407" y="148">开发商红色旧点位已降级为参考</text>')
    markup.append(f'<text class="r1-note" x="407" y="382">IFC SHA {report["source_ifc_sha256"][:12]}…</text></g>')
    style = """
@page{size:500mm 400mm;margin:0}html,body{margin:0;width:500mm;height:400mm;overflow:hidden}
.r1-bedside{fill:#00d5e8;stroke:#006672;stroke-width:.7}.r1-new-socket{fill:#37b24d;stroke:#145c23;stroke-width:.7}
.r1-cabinet{fill:#ffb000;stroke:#7a4a00;stroke-width:.7}.r1-recheck{fill:#fa2b2b;stroke:#7b0000;stroke-width:.7}
.r1-label{font-family:Arial,'Noto Sans CJK SC',sans-serif;font-size:2.15px;font-weight:700;fill:#102f43;paint-order:stroke;stroke:#fff;stroke-width:.8px}
.r1-panel{fill:#fbfcfd;stroke:#102f43;stroke-width:.5}.r1-title,.r1-note,.r1-text,.r1-warn{font-family:Arial,'Noto Sans CJK SC',sans-serif;fill:#102f43}
.r1-title{font-size:3.7px;font-weight:700}.r1-heading{font-family:Arial,'Noto Sans CJK SC',sans-serif;font-size:2.45px;font-weight:700;fill:#0f4c81}.r1-note{font-size:2.15px;fill:#526777}.r1-text{font-size:2.3px}.r1-warn{font-size:2.15px;fill:#c92a2a;font-weight:700}
"""
    return inject_svg(source, "".join(markup), style)


def render_items(report: dict[str, Any]) -> list[tuple[dict[str, Any], str, str]]:
    return [
        *((row, "r1-bedside", "position_mm") for row in report["bedside_light_candidates"]),
        *((row, "r1-new-socket", "position_mm") for row in report["new_socket_candidates"]),
        *((row, "r1-cabinet", "marker_position_mm") for row in report["cabinet_power_zones"]),
        *((row, "r1-recheck", "position_mm") for row in report["kitchen_socket_rechecks"]),
    ]


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    ifc_hash = sha256(args.ifc)
    if args.expected_ifc_sha256 and ifc_hash != args.expected_ifc_sha256:
        raise RuntimeError(
            f"formal IFC hash changed: expected {args.expected_ifc_sha256}, got {ifc_hash}"
        )
    requirements = read_requirements(args.requirements)
    existing = read_json(args.elec_existing)
    positioning = read_json(args.elec_positioning)
    int1 = read_json(args.int1)
    source_hashes = {
        existing["source"]["sha256"],
        positioning["source_ifc_sha256"],
        int1["source"]["ifc_sha256"],
    }
    if source_hashes != {ifc_hash}:
        raise RuntimeError(f"stale ELEC/INT1 inputs: {source_hashes}")

    bedside = bedside_candidates(int1["records"])
    circuit_decisions = e303_circuit_decisions(args.owner_decisions)
    if args.appliances:
        appliance_rows = read_owner_csv(args.appliances, APPLIANCE_HEADERS)
        appliance_source = {"path": str(args.appliances.resolve()), "sha256": sha256(args.appliances)}
        appliance_dependency_paths = [args.appliances]
    else:
        validate_equipment_ssot(root)
        appliance_rows = appliance_projection_rows(root)
        canonical_paths = [root / "pipeline/decisions/equipment-register.csv", root / "pipeline/decisions/equipment-installation-requirements.csv", root / "pipeline/decisions/source-evidence-register.csv"]
        appliance_source = {"canonical_tables": {str(path.relative_to(root)): sha256(path) for path in canonical_paths}}
        appliance_dependency_paths = canonical_paths
    appliance_context = appliance_socket_context(
        appliance_rows,
        args.load_scenarios,
        circuit_decisions,
        args.circuit_design,
    )
    new_sockets = new_socket_candidates(int1["records"], appliance_context)
    cabinet_zones = cabinet_power_zones(int1["records"])
    socket_rechecks = kitchen_socket_rechecks(existing, positioning)
    requirement_status = {row["requirement_id"]: row["status"] for row in requirements}
    report = {
        "mode": "read_only_renovation_electrical_round1_candidate",
        "source_ifc_sha256": ifc_hash,
        "requirements_path": str(args.requirements.resolve()),
        "source_dependencies": [
            {"path": str(path.resolve()), "sha256": sha256(path)}
            for path in [
                args.requirements,
                args.elec_existing,
                args.elec_positioning,
                args.int1,
                args.owner_decisions,
                args.load_scenarios,
                args.circuit_design,
                *appliance_dependency_paths,
            ]
        ],
        "appliance_inputs": appliance_source,
        "owner_inputs": {
            "path": str(args.owner_decisions.resolve()),
            "sha256": sha256(args.owner_decisions),
            "e303_circuit_decisions": {
                input_id: {"status": row["status"]}
                for input_id, row in circuit_decisions.items()
            },
        },
        "load_scenarios": {
            "path": str(args.load_scenarios.resolve()),
            "sha256": sha256(args.load_scenarios),
        },
        "circuit_design": {
            "path": str(args.circuit_design.resolve()),
            "sha256": sha256(args.circuit_design),
        },
        "summary": {
            "confirmed_requirements": sum(value.startswith("confirmed") for value in requirement_status.values()),
            "bedside_light_candidates": len(bedside),
            "new_socket_candidates": len(new_sockets),
            "new_socket_known_connected_load_w": {
                row["candidate_id"]: row["appliance_context"]["known_connected_load_w"]
                for row in new_sockets
            },
            "new_socket_planning_envelopes": {
                row["candidate_id"]: row["appliance_context"]["planning_envelope"]
                for row in new_sockets
            },
            "cabinet_power_zones": len(cabinet_zones),
            "kitchen_socket_rechecks": len(socket_rechecks),
            "developer_red_points_are_reference_only": requirement_status.get("ELEC-R1-001") == "confirmed",
        },
        "requirements": requirements,
        "bedside_light_candidates": bedside,
        "new_socket_candidates": new_sockets,
        "e303_circuit_semantic_gates": {
            row["candidate_id"]: row["appliance_context"]["semantic_gates"]
            for row in new_sockets
        },
        "cabinet_power_zones": cabinet_zones,
        "kitchen_socket_rechecks": socket_rechecks,
        "gates": {
            "four_bedside_lights_present": len(bedside) == 4,
            "island_and_dining_bay_socket_present": {row["candidate_role"] for row in new_sockets}
            == {"island_end_panel_socket", "dining_bay_socket"},
            "socket_use_lists_compiled_without_fabricated_load": {
                row["candidate_id"]: len(row["appliance_context"]["items"])
                for row in new_sockets
            } == {"NS-01": 3, "NS-02": 3}
            and {
                row["candidate_id"]: row["appliance_context"]["known_connected_load_w"]
                for row in new_sockets
            } == {"NS-01": 0.0, "NS-02": 1200.0}
            and all(
                row["appliance_context"]["socket_form_and_circuit_sizing_ready"]
                for row in new_sockets
            ),
            "circuit_planning_candidates_match_confirmed_use": {
                row["candidate_id"]: row["appliance_context"]["planning_envelope"]["minimum_independent_circuit_count_candidate"]
                for row in new_sockets
            } == {"NS-01": 2, "NS-02": 2}
            and all(
                row["appliance_context"]["planning_envelope"]["minimum_connection_positions_candidate"] == 3
                and row["appliance_context"]["planning_envelope"]["calculation_status"]
                == "research_conclusion_not_product_nameplate_load"
                for row in new_sockets
            ),
            "use_confirmed": all(
                row["appliance_context"]["semantic_gates"]["use_confirmed"]
                for row in new_sockets
            ),
            "planning_envelope_compiled": all(
                row["appliance_context"]["semantic_gates"]["planning_envelope_compiled"]
                for row in new_sockets
            ),
            "product_and_circuit_fixed": False,
            "illuminated_cabinet_power_is_grouped_not_fabricated": len(cabinet_zones) == 7
            and all(row["coordinate_status"] == "assembly_zone_only" for row in cabinet_zones),
            "all_current_kitchen_sockets_reopened_for_review": len(socket_rechecks) == 11,
            "reverse_requirement_audit_planned": requirement_status.get("ELEC-R1-007") == "planned",
            "whole_home_electrical_positioning_complete": False,
            "automatic_ifc_write_allowed": False,
        },
    }
    _labels, collision_count = label_layout(render_items(report))
    report["summary"]["label_collision_count"] = collision_count
    report["gates"]["label_collision_free"] = collision_count == 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    source = args.source_svg.read_text(encoding="utf-8")
    validate_wall_plan_source(source, args.source_svg, args.ifc)
    args.output_svg.write_text(
        render_svg(source, report),
        encoding="utf-8",
    )
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
