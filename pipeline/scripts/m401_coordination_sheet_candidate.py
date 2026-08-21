#!/usr/bin/env python3
"""Build an evidence-bounded M-401 HVAC and safety coordination sheet candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell


EXPECTED_ROLE_COUNTS = {
    "assigned_ac_equipment_instance": 5,
    "placement_only_ac_equipment_instance": 1,
    "high_level_ceiling_or_led_coordination_context": 15,
    "legacy_base_condensate_geometry": 1,
    "legacy_base_refrigerant_gas_geometry": 2,
    "legacy_base_refrigerant_liquid_geometry": 2,
    "named_embedded_ac_diffuser_proxy": 2,
    "named_flue_check_valve_proxy": 1,
}
EXPECTED_ROUTE_IDS = {
    "RCP1-AIRSIDE-A02",
    "RCP1-AIRSIDE-A06",
    "RCP1-SERVICE-A02",
    "RCP1-SERVICE-A03",
    "RCP1-CONDENSATE-ENDPOINT",
    "RCP1-OUTDOOR-ENDPOINT",
}
EXPECTED_WAYPOINTS = {
    "RCP1-SERVICE-A02": ["A02", "H03"],
    "RCP1-SERVICE-A03": ["A03", "H04", "H02"],
    "RCP1-CONDENSATE-ENDPOINT": ["H07"],
    "RCP1-OUTDOOR-ENDPOINT": ["H01"],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def find_chrome(explicit: Path | None) -> Path:
    candidates = [
        explicit,
        Path(os.environ["CHROME_BIN"]) if os.environ.get("CHROME_BIN") else None,
        Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"),
        Path("/Applications/Chromium.app/Contents/MacOS/Chromium"),
    ]
    for name in ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser"):
        binary = shutil.which(name)
        if binary:
            candidates.append(Path(binary))
    for candidate in candidates:
        if candidate and candidate.is_file():
            return candidate
    raise RuntimeError("Chrome/Chromium not found; pass --chrome or set CHROME_BIN")


def validate_inventory(
    report: dict[str, Any],
    review_rows: list[dict[str, str]],
    ifc_hash: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if report.get("source_ifc_sha256") != ifc_hash:
        raise RuntimeError("M-401 inventory report is stale against the formal IFC")
    summary = report.get("summary", {})
    if summary.get("actual_instances") != 29:
        raise RuntimeError("M-401 actual-instance count drifted from 29")
    if summary.get("instance_role_counts") != EXPECTED_ROLE_COUNTS:
        raise RuntimeError("M-401 observable role counts drifted")
    if summary.get("type_definitions") != 3 or summary.get("missing_input_blocks") != 5:
        raise RuntimeError("M-401 type/BLOCK boundary drifted")
    if summary.get("ifc_sensor_instances") != 1:
        raise RuntimeError("M-401 inventory no longer reports exactly one IfcSensor")
    if len(review_rows) != 37 or any(row.get("source_ifc_sha256") != ifc_hash for row in review_rows):
        raise RuntimeError("M-401 review CSV is incomplete or stale")
    if Counter(row["record_kind"] for row in review_rows) != {
        "actual_instance": 29,
        "type_definition": 3,
        "missing_input": 5,
    }:
        raise RuntimeError("M-401 review CSV record-kind counts drifted")
    if sum(row["review_status"] == "BLOCK" for row in review_rows) != 5:
        raise RuntimeError("M-401 review CSV must retain five BLOCK records")

    records = report.get("records", [])
    types = [record for record in records if record.get("record_kind") == "type_definition"]
    blockers = [record for record in records if record.get("review_status") == "BLOCK"]
    type_counts = {record.get("type_name"): record.get("type_occurrence_count") for record in types}
    if type_counts != {"AC1180": 1, "AC700": 2, "AC700F": 2}:
        raise RuntimeError(f"M-401 AC type assignment drifted: {type_counts}")
    if len(blockers) != 5:
        raise RuntimeError("M-401 report must retain five BLOCK records")
    return types, blockers


def validate_routes(
    route_rows: list[dict[str, str]],
    waypoint_rows: list[dict[str, str]],
) -> list[dict[str, Any]]:
    if {row.get("route_id") for row in route_rows} != EXPECTED_ROUTE_IDS:
        raise RuntimeError("RCP1 route register IDs drifted")
    for row in route_rows:
        if row.get("status") != "user_confirmed" or row.get("formal_ifc_write_allowed") != "no":
            raise RuntimeError(f"RCP1 route boundary is not confirmed/read-only: {row.get('route_id')}")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in waypoint_rows:
        if row.get("status") != "user_confirmed":
            raise RuntimeError(f"RCP1 waypoint is not user-confirmed: {row.get('route_id')}")
        grouped[row["route_id"]].append(row)
    actual = {
        route_id: [row["anchor_id"] for row in sorted(rows, key=lambda item: int(item["sequence"]))]
        for route_id, rows in grouped.items()
    }
    if actual != EXPECTED_WAYPOINTS:
        raise RuntimeError(f"RCP1 route waypoint order drifted: {actual}")
    result = []
    for row in route_rows:
        result.append(
            {
                "route_id": row["route_id"],
                "route_kind": row["route_kind"],
                "equipment_id": row["equipment_id"],
                "served_space_reference": row["served_space_reference"],
                "served_space_long_name": row["served_space_long_name"],
                "terminal_anchor_id": row["terminal_anchor_id"],
                "waypoint_order": actual.get(row["route_id"], []),
                "status": row["status"],
                "basis": row["basis"],
                "formal_ifc_write_allowed": False,
                "final_remodel_route": False,
            }
        )
    return result


def validate_mapping_coverage(
    report: dict[str, Any], ifc_hash: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if report.get("source", {}).get("ifc_sha256") != ifc_hash:
        raise RuntimeError("RCP1 route-readiness mapping is stale against the formal IFC")
    mappings = report.get("equipment_opening_mapping", [])
    interfaces = report.get("interface_coverage", [])
    if [row.get("equipment_id") for row in mappings] != [f"A0{index}" for index in range(1, 7)]:
        raise RuntimeError("RCP1 equipment-opening mapping must cover A01–A06 exactly once")
    if [row.get("opening_id") for row in interfaces] != [f"H0{index}" for index in range(1, 8)]:
        raise RuntimeError("RCP1 interface coverage must cover H01–H07 exactly once")
    if Counter(row.get("status") for row in mappings) != {
        "user_confirmed": 2,
        "geometry_candidate": 4,
    }:
        raise RuntimeError("RCP1 mapping confirmation/candidate boundary drifted")
    if any(row.get("formal_ifc_write_allowed") is not False for row in mappings + interfaces):
        raise RuntimeError("RCP1 mapping coverage unexpectedly grants IFC write authority")
    required_gates = (
        "six_equipment_opening_mappings_registered",
        "seven_interfaces_have_controlled_roles",
        "confirmed_and_candidate_mapping_split_preserved",
    )
    if any(report.get("gates", {}).get(name) is not True for name in required_gates):
        raise RuntimeError("RCP1 mapping coverage gates are incomplete")
    return mappings, interfaces


def validate_safety(
    ifc_path: Path,
    a106_report: dict[str, Any],
    ifc_hash: str,
) -> dict[str, Any]:
    if a106_report.get("source_ifc_sha256") != ifc_hash:
        raise RuntimeError("A-106 ceiling-device evidence is stale against the formal IFC")
    fire_candidates = [
        item for item in a106_report.get("candidates", [])
        if item.get("candidate_id") == "A106-FIRE-R04"
    ]
    if len(fire_candidates) != 1:
        raise RuntimeError("A-106 must contain exactly one A106-FIRE-R04 record")
    fire = fire_candidates[0]
    if fire.get("review_status") != "position_confirmed_type_pending":
        raise RuntimeError("A106-FIRE-R04 must remain position-confirmed/type-pending")

    model = ifcopenshell.open(str(ifc_path))
    sensors = model.by_type("IfcSensor")
    alarms = model.by_type("IfcAlarm")
    if len(sensors) != 1:
        raise RuntimeError(f"formal IFC must contain exactly one IfcSensor, found {len(sensors)}")
    sensor = sensors[0]
    if sensor.Name != "A106-FIRE-R04" or str(sensor.PredefinedType) != "FIRESENSOR":
        raise RuntimeError("the one formal IfcSensor is not the confirmed A106-FIRE-R04 FIRESENSOR")
    return {
        "ifc_sensor_instances": len(sensors),
        "ifc_alarm_instances": len(alarms),
        "sensor_global_id": sensor.GlobalId,
        "sensor_name": sensor.Name,
        "sensor_predefined_type": str(sensor.PredefinedType),
        "confirmed_position_mm": fire.get("position_mm"),
        "space_global_id": fire.get("space_global_id"),
        "review_status": fire.get("review_status"),
        "final_type_product_power_communication_pending": True,
        "final_release_pass": bool(fire.get("final_release_pass")),
        "release_blocker": fire.get("release_blocker"),
    }


def validate_l_supply_design_direction(
    owner_rows: list[dict[str, str]],
    requirement_rows: list[dict[str, str]],
) -> dict[str, Any]:
    owner = {row["input_id"]: row for row in owner_rows}
    row = owner.get("LEGACY-M401-L-SUPPLY-AIR")
    if row is None:
        raise RuntimeError("M-401 L-shaped supply-air project direction is missing")
    if (
        row["status"] != "需证据"
        or "L 形转角风口" not in row["candidate_value"]
        or "侧出风＋下回风" not in row["candidate_value"]
        or "项目先计算并画图" not in row["user_value"]
    ):
        raise RuntimeError("M-401 L-shaped supply-air project direction changed")
    requirements = {row["requirement_id"]: row for row in requirement_rows}
    insulation = requirements.get("REQ-HVACINS-021")
    if (
        insulation is None
        or insulation["status"] != "pending"
        or insulation["blocks_release"] != "yes"
    ):
        raise RuntimeError("M-401 final insulation schedule boundary changed")
    return {
        "input_id": row["input_id"],
        "design_direction": row["candidate_value"],
        "owner_direction": row["user_value"],
        "project_calculation_method": "A_eff >= Q / v_design",
        "project_drawing_scope": (
            "M-401 plan + section + grille schedule: airflow, effective free area, "
            "face velocity, external static pressure, corner loss, noise, return-air access"
        ),
        "external_review_scope": (
            "verify actual grille free-area ratio, corner-loss data, selected unit duty, "
            "noise and maintenance access against the issued project drawing"
        ),
        "unknown_until_external_data": [
            "final room airflow",
            "available external static pressure",
            "selected grille free-area ratio",
            "corner pressure-loss coefficient",
            "manufacturer noise data",
            "final segment insulation schedule",
        ],
    }


def esc(value: object) -> str:
    return html.escape(str(value))


def render_svg(
    ifc_hash: str,
    inventory_summary: dict[str, Any],
    types: list[dict[str, Any]],
    routes: list[dict[str, Any]],
    blockers: list[dict[str, Any]],
    safety: dict[str, Any],
    mappings: list[dict[str, Any]],
    l_supply: dict[str, Any],
) -> str:
    type_cards = []
    for index, item in enumerate(sorted(types, key=lambda row: row["type_name"])):
        x = 70 + index * 300
        type_cards.append(
            f'<rect class="card" x="{x}" y="268" width="270" height="96" rx="14"/>'
            f'<text class="type" x="{x + 22}" y="304">{esc(item["type_name"])}</text>'
            f'<text class="count" x="{x + 230}" y="304">×{esc(item["type_occurrence_count"])}</text>'
            f'<text class="small" x="{x + 22}" y="330">{esc(item["type_description"])}</text>'
            f'<text class="warn" x="{x + 22}" y="352">型号可查 ≠ 最终设备选型</text>'
        )

    route_labels = {
        "RCP1-AIRSIDE-A02": "A02 → R07 餐厅（服务对象）",
        "RCP1-AIRSIDE-A06": "A06 → R07 餐厅（正式 placement-only 身份）",
        "RCP1-SERVICE-A02": "A02 → H03",
        "RCP1-SERVICE-A03": "A03 → H04 → H02",
        "RCP1-CONDENSATE-ENDPOINT": "H07 = 冷凝水排放接口位置",
        "RCP1-OUTDOOR-ENDPOINT": "H01 = 室外机接口位置",
    }
    route_cards = []
    for index, route in enumerate(routes):
        y = 425 + index * 67
        route_cards.append(
            f'<rect class="route" x="70" y="{y}" width="860" height="52" rx="9"/>'
            f'<text class="route-id" x="92" y="{y + 22}">{esc(route["route_id"])}</text>'
            f'<text class="route-label" x="360" y="{y + 22}">{esc(route_labels[route["route_id"]])}</text>'
            f'<text class="route-note" x="92" y="{y + 42}">用户确认约束 · 非最终管线 · 不写 IFC</text>'
        )

    blocker_names = {
        "M401-MISS-001": "正式风口 / 系统 / 端口连接",
        "M401-MISS-002": "厨房可燃气体报警器",
        "M401-MISS-003": "卧室 / 客厅火警与厨房最终产品",
        "M401-MISS-004": "卫生间暖风 / 通风设备",
        "M401-MISS-005": "风量、管径、坡度、电气、控制与检修",
    }
    blocker_cards = []
    for index, blocker in enumerate(blockers):
        y = 425 + index * 82
        blocker_cards.append(
            f'<rect class="block" x="995" y="{y}" width="530" height="70" rx="10"/>'
            f'<text class="block-id" x="1018" y="{y + 24}">{esc(blocker["queue_id"])} · BLOCK</text>'
            f'<text class="block-name" x="1018" y="{y + 50}">{esc(blocker_names[blocker["queue_id"]])}</text>'
        )

    mapping_cards = []
    for index, mapping in enumerate(mappings):
        column = index % 3
        row = index // 3
        x = 70 + column * 500
        y = 905 + row * 63
        confirmed = mapping["status"] == "user_confirmed"
        status_label = "用户确认" if confirmed else "几何候选"
        card_class = "map-confirmed" if confirmed else "map-candidate"
        mapping_cards.append(
            f'<rect class="{card_class}" x="{x}" y="{y}" width="470" height="50" rx="9"/>'
            f'<text class="map-id" x="{x + 18}" y="{y + 22}">{esc(mapping["equipment_id"])} → {esc(mapping["opening_id"])}</text>'
            f'<text class="map-status" x="{x + 160}" y="{y + 22}">{status_label}</text>'
            f'<text class="map-note" x="{x + 18}" y="{y + 41}">最小包围盒净距候选 {float(mapping["minimum_clearance_candidate_mm"]):.1f} mm · 不写 IFC</text>'
        )

    roles = inventory_summary["instance_role_counts"]
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="1400" viewBox="0 0 1600 1400">
<style>
.bg{{fill:#f3f6fa}}.panel{{fill:#fff;stroke:#cbd5e1;stroke-width:2}}.card{{fill:#eef6ff;stroke:#93c5fd;stroke-width:2}}
.title{{font:700 38px -apple-system,"PingFang SC",sans-serif;fill:#0f172a}}.subtitle{{font:18px -apple-system,"PingFang SC",sans-serif;fill:#475569}}
.kpi{{font:700 31px ui-monospace,SFMono-Regular,monospace;fill:#0f4c81}}.kpi-label{{font:15px -apple-system,"PingFang SC",sans-serif;fill:#475569}}
.section{{font:700 24px -apple-system,"PingFang SC",sans-serif;fill:#1e293b}}.type{{font:700 26px ui-monospace,monospace;fill:#1e3a8a}}
.count{{font:700 24px ui-monospace,monospace;fill:#0369a1}}.small{{font:14px -apple-system,"PingFang SC",sans-serif;fill:#334155}}.warn{{font:700 15px -apple-system,"PingFang SC",sans-serif;fill:#a16207}}
.route{{fill:#f8fafc;stroke:#cbd5e1}}.route-id{{font:700 15px ui-monospace,monospace;fill:#334155}}.route-label{{font:700 17px -apple-system,"PingFang SC",sans-serif;fill:#0f4c81}}
.route-note{{font:13px -apple-system,"PingFang SC",sans-serif;fill:#a16207}}.block{{fill:#fff7ed;stroke:#fdba74;stroke-width:2}}
.block-id{{font:700 15px ui-monospace,monospace;fill:#9a3412}}.block-name{{font:16px -apple-system,"PingFang SC",sans-serif;fill:#7c2d12}}
.map-confirmed{{fill:#eff6ff;stroke:#60a5fa;stroke-width:2}}.map-candidate{{fill:#f8fafc;stroke:#94a3b8;stroke-width:2}}
.map-id{{font:700 17px ui-monospace,monospace;fill:#0f4c81}}.map-status{{font:700 15px -apple-system,"PingFang SC",sans-serif;fill:#475569}}
.map-note{{font:13px -apple-system,"PingFang SC",sans-serif;fill:#64748b}}
.safety{{fill:#ecfdf5;stroke:#6ee7b7;stroke-width:2}}.safety-title{{font:700 19px -apple-system,"PingFang SC",sans-serif;fill:#065f46}}
.safety-text{{font:15px -apple-system,"PingFang SC",sans-serif;fill:#064e3b}}.gate{{font:700 16px ui-monospace,monospace;fill:#991b1b}}
.footer{{font:14px ui-monospace,monospace;fill:#64748b}}
</style>
<rect class="bg" width="1600" height="1240"/>
<text class="title" x="55" y="62">M-401 空调、通风及安全设备协调候选</text>
<text class="subtitle" x="55" y="100">证据化索引 + 路线约束图 · 只读生成 · 旧紫色管线不是装修后最终路线 · 非施工发布</text>
<rect class="panel" x="55" y="125" width="1490" height="90" rx="14"/>
<text class="kpi" x="90" y="168">29</text><text class="kpi-label" x="90" y="194">可追踪实例</text>
<text class="kpi" x="250" y="168">5 + 1</text><text class="kpi-label" x="250" y="194">实体 AC / 点位 AC</text>
<text class="kpi" x="500" y="168">5</text><text class="kpi-label" x="500" y="194">旧冷媒/冷凝几何</text>
<text class="kpi" x="675" y="168">2 + 1</text><text class="kpi-label" x="675" y="194">风口 Proxy / 止回阀 Proxy</text>
<text class="kpi" x="970" y="168">6 / 7</text><text class="kpi-label" x="970" y="194">确认约束 / 顺序锚点</text>
<text class="kpi" x="1240" y="168">5 BLOCK</text><text class="kpi-label" x="1240" y="194">施工发布未关闭</text>

<text class="section" x="70" y="250">当前 IFC 空调类型身份（不是最终选型）</text>
{''.join(type_cards)}

<rect class="panel" x="55" y="390" width="900" height="475" rx="16"/>
<text class="section" x="78" y="418">用户确认的服务 / 端点 / 经过顺序约束</text>
{''.join(route_cards)}
<text class="warn" x="78" y="838">只表达对象、端点和顺序；不表达介质、管径、保温、坡度、端口连接或施工走向。</text>

<rect class="panel" x="975" y="390" width="570" height="475" rx="16"/>
<text class="section" x="998" y="418">不得冒充完成的发布阻塞</text>
{''.join(blocker_cards)}

<text class="section" x="70" y="892">6 台设备—7 个既有接口覆盖矩阵（确认关系与几何候选分层）</text>
{''.join(mapping_cards)}

<rect class="safety" x="55" y="1040" width="1490" height="125" rx="16"/>
<text class="safety-title" x="82" y="1075">安全设备现状：正式 IFC 仅 1 个 IfcSensor / FIRESENSOR</text>
<text class="safety-text" x="82" y="1105">A106-FIRE-R04 · 中厨 R04 · 位置 {esc(safety["confirmed_position_mm"])} mm 已确认；最终感温/感烟/复合类型、产品、供电通信和厂家安装条件仍待确认。</text>
<text class="gate" x="82" y="1137">automatic_ifc_write_allowed=false · construction_release_ready=false · 未建模内容不得由近接、名称或旧几何推断</text>

<rect class="panel" x="55" y="1185" width="1490" height="145" rx="16"/>
<text class="section" x="82" y="1220">客厅 L 形转角送风口｜项目设计方向已定，性能数据待对图复核</text>
<text class="route-label" x="82" y="1250">侧出风＋下回风；可拆回风口兼检修口；L 形风口用于弱化窗帘盒阳角</text>
<text class="small" x="82" y="1278">项目计算：有效开口 Aeff ≥ Q / 设计面风速；M-401 须列风量、有效开口、面风速、机外静压、转角压损、噪声与检修路线。</text>
<text class="small" x="82" y="1304">厂家／安装方只按项目图核对真实风口开孔率、转角损失、机组余压和噪声；没有数据时保持 unknown，不按效果图猜尺寸。</text>
<text class="warn" x="82" y="1325">当前不可施工放样：最终房间风量、可用余压、风口开孔率、转角损失系数及逐段保温表仍待关闭。</text>
<text class="footer" x="55" y="1373">IFC SHA-256 {ifc_hash} · 6/6 equipment mappings · 7/7 controlled interface roles · 5 BLOCK</text>
</svg>
'''


def render_png(svg_path: Path, png_path: Path, chrome: Path) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(chrome),
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            "--force-device-scale-factor=1",
            f"--screenshot={png_path}",
            "--window-size=1600,1400",
            svg_path.resolve().as_uri(),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if not png_path.is_file() or png_path.stat().st_size == 0:
        raise RuntimeError("M-401 proof PNG was not created")
    if png_path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError("M-401 proof output is not a PNG")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--input-ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument("--expected-ifc-sha256")
    parser.add_argument("--m401-report", type=Path, default=Path("build/rcp1/m401-existing-report.json"))
    parser.add_argument("--m401-review", type=Path, default=Path("pipeline/decisions/m401-existing-review.csv"))
    parser.add_argument("--route-register", type=Path, default=Path("pipeline/decisions/rcp1-hvac-route-register.csv"))
    parser.add_argument("--route-waypoints", type=Path, default=Path("pipeline/decisions/rcp1-hvac-route-waypoints.csv"))
    parser.add_argument("--route-readiness", type=Path, default=Path("build/rcp1/route-readiness-candidate.json"))
    parser.add_argument("--a106-report", type=Path, default=Path("build/elec/a106-ceiling-device-candidate.json"))
    parser.add_argument("--owner-inputs", type=Path, default=Path("pipeline/decisions/owner-input-register.csv"))
    parser.add_argument(
        "--installation-requirements",
        type=Path,
        default=Path("pipeline/decisions/equipment-installation-requirements.csv"),
    )
    parser.add_argument("--output-svg", type=Path, default=Path("drawings/M-401-hvac-safety-coordination-candidate.svg"))
    parser.add_argument("--proof-png", type=Path, default=Path("build/rcp1/M-401-hvac-safety-coordination-candidate.png"))
    parser.add_argument("--report", type=Path, default=Path("build/rcp1/m401-coordination-sheet-candidate.json"))
    parser.add_argument("--chrome", type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    resolve = lambda path: path if path.is_absolute() else root / path
    paths = {
        "ifc": resolve(args.input_ifc),
        "m401_report": resolve(args.m401_report),
        "m401_review": resolve(args.m401_review),
        "route_register": resolve(args.route_register),
        "route_waypoints": resolve(args.route_waypoints),
        "route_readiness": resolve(args.route_readiness),
        "a106_report": resolve(args.a106_report),
        "owner_inputs": resolve(args.owner_inputs),
        "installation_requirements": resolve(args.installation_requirements),
    }
    for label, path in paths.items():
        if not path.is_file():
            raise RuntimeError(f"missing M-401 coordination input {label}: {path}")

    ifc_hash = sha256(paths["ifc"])
    if args.expected_ifc_sha256 and args.expected_ifc_sha256 != ifc_hash:
        raise RuntimeError(
            f"formal IFC SHA mismatch: expected {args.expected_ifc_sha256}, got {ifc_hash}"
        )
    m401_report = read_json(paths["m401_report"])
    types, blockers = validate_inventory(
        m401_report,
        read_csv(paths["m401_review"]),
        ifc_hash,
    )
    routes = validate_routes(
        read_csv(paths["route_register"]),
        read_csv(paths["route_waypoints"]),
    )
    mappings, interfaces = validate_mapping_coverage(
        read_json(paths["route_readiness"]),
        ifc_hash,
    )
    safety = validate_safety(paths["ifc"], read_json(paths["a106_report"]), ifc_hash)
    l_supply = validate_l_supply_design_direction(
        read_csv(paths["owner_inputs"]),
        read_csv(paths["installation_requirements"]),
    )

    output_svg = resolve(args.output_svg)
    proof_png = resolve(args.proof_png)
    report_path = resolve(args.report)
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    output_svg.write_text(
        render_svg(
            ifc_hash,
            m401_report["summary"],
            types,
            routes,
            blockers,
            safety,
            mappings,
            l_supply,
        ),
        encoding="utf-8",
    )
    render_png(output_svg, proof_png, find_chrome(args.chrome))

    input_hashes = {
        label: {"path": str(path), "sha256": sha256(path)} for label, path in paths.items()
    }
    report: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_m401_hvac_safety_coordination_sheet_candidate",
        "source_ifc_sha256": ifc_hash,
        "source": {"ifc_sha256": ifc_hash, "inputs": input_hashes},
        "summary": {
            "actual_instance_count": m401_report["summary"]["actual_instances"],
            "instance_role_counts": m401_report["summary"]["instance_role_counts"],
            "ac_type_count": len(types),
            "assigned_ac_instance_count": sum(int(item["type_occurrence_count"]) for item in types),
            "placement_only_ac_instance_count": m401_report["summary"]["instance_role_counts"]["placement_only_ac_equipment_instance"],
            "missing_input_block_count": len(blockers),
            "confirmed_route_constraint_count": len(routes),
            "confirmed_waypoint_count": sum(len(items) for items in EXPECTED_WAYPOINTS.values()),
            "equipment_opening_mapping_count": len(mappings),
            "controlled_interface_role_count": len(interfaces),
            "user_confirmed_equipment_mapping_count": sum(
                item["status"] == "user_confirmed" for item in mappings
            ),
            "geometry_candidate_equipment_mapping_count": sum(
                item["status"] == "geometry_candidate" for item in mappings
            ),
            "ifc_sensor_instances": safety["ifc_sensor_instances"],
            "ifc_alarm_instances": safety["ifc_alarm_instances"],
            "legacy_routes_declared_final_count": 0,
        },
        "equipment_types": types,
        "route_constraints": routes,
        "equipment_opening_mapping": mappings,
        "interface_coverage": interfaces,
        "release_blockers": blockers,
        "safety_context": safety,
        "l_shaped_supply_air_design": l_supply,
        "outputs": {
            "svg": {"path": str(output_svg), "sha256": sha256(output_svg)},
            "proof_png": {"path": str(proof_png), "sha256": sha256(proof_png)},
        },
        "gates": {
            "caller_frozen_ifc_hash_match": (
                None
                if args.expected_ifc_sha256 is None
                else args.expected_ifc_sha256 == ifc_hash
            ),
            "source_hashes_current": True,
            "inventory_counts_match": True,
            "route_constraints_confirmed": True,
            "six_equipment_opening_mappings_registered": len(mappings) == 6,
            "seven_interfaces_have_controlled_roles": (
                len(interfaces) == 7 and all(item["role"] for item in interfaces)
            ),
            "mapping_confirmation_boundary_preserved": (
                sum(item["status"] == "user_confirmed" for item in mappings) == 2
                and sum(item["status"] == "geometry_candidate" for item in mappings) == 4
            ),
            "legacy_routes_marked_nonfinal": all(not item["final_remodel_route"] for item in routes),
            "one_ifc_sensor_preserved": safety["ifc_sensor_instances"] == 1,
            "kitchen_fire_position_only_closed": safety["review_status"] == "position_confirmed_type_pending",
            "five_release_blockers_explicit": len(blockers) == 5,
            "automatic_ifc_write_allowed": False,
            "construction_release_ready": False,
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {"report": str(report_path), "summary": report["summary"], "gates": report["gates"]},
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
