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


EXPECTED_IFC_SHA256 = "7521c09991f3d0c7b7d91ca2324fd55ad961d8e32e9e3e9a9777a4cc19b06e81"
SCALE_DENOMINATOR = 50.0
SVG_WORLD_OFFSET_MM = 10000.0


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument("--rules", type=Path, default=root / "pipeline/decisions/elec-design-rules.csv")
    parser.add_argument("--doors", type=Path, default=root / "pipeline/decisions/a104-door-window-review.csv")
    parser.add_argument("--spaces", type=Path, default=root / "pipeline/decisions/space-reference-review.csv")
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


def network_zones(spaces: list[dict[str, str]]) -> list[dict[str, Any]]:
    by_reference = {row["candidate_reference"]: row for row in spaces}
    definitions = [
        ("NET-AP-MASTER", ["R09"], "主卧", "wireless_access_point"),
        ("NET-AP-GUEST", ["R14"], "次卧", "wireless_access_point"),
        ("NET-ROUTER-LIVING-STUDY", ["R20", "R22"], "客厅＋书房开放空间", "router_no_AP"),
    ]
    rows = []
    for candidate_id, references, room_name, role in definitions:
        source_spaces = [by_reference[reference] for reference in references]
        centres = [
            [float(space["centre_x_mm"]), float(space["centre_y_mm"]), float(space["centre_z_mm"])]
            for space in source_spaces
        ]
        rows.append({
            "candidate_id": candidate_id,
            "kind": "network_device_room_coordination_zone",
            "room_reference": ";".join(references),
            "room_name": room_name,
            "source_global_ids": [space["space_global_id"] for space in source_spaces],
            "position_mm": [sum(values) / len(values) for values in zip(*centres)],
            "network_role": role,
            "wired_backhaul": True,
            "position_basis": "confirmed Space bbox centre average used only to identify the shared open-space coordination zone" if len(references) > 1 else "confirmed Space bbox centre used only to identify the room-level coordination zone",
            "coordinate_status": "shared_open_space_zone_only_final_xy_z_pending" if len(references) > 1 else "room_zone_only_final_xy_z_pending",
            "confidence": 1.0,
            "review_required": True,
            "automatic_ifc_write_allowed": False,
        })
    return rows


def render_svg(source: str, report: dict[str, Any]) -> str:
    source = re.sub(r'width="400(?:\.0+)?mm"', 'width="500mm"', source, count=1)
    source = re.sub(r'viewBox="0 0 400(?:\.0+)? 400(?:\.0+)?"', 'viewBox="0 0 500 400"', source, count=1)
    if "</svg>" not in source or 'viewBox="0 0 500 400"' not in source:
        raise RuntimeError("could not prepare 500x400 control/network SVG")
    markup = []
    for index, row in enumerate(report["control_coordination_zones"]):
        x, y = world_to_svg(row["position_mm"])
        candidate_id = html.escape(row["candidate_id"])
        markup.append(
            f'<g data-candidate-id="{candidate_id}"><rect class="cn-control" x="{x-2.2:.3f}" y="{y-2.2:.3f}" width="4.4" height="4.4"/>'
            f'<text class="cn-label" x="{x+3.2:.3f}" y="{y+(-3 if index else 5):.3f}">{candidate_id}</text></g>'
        )
    for index, row in enumerate(report["network_coordination_zones"]):
        x, y = world_to_svg(row["position_mm"])
        css = "cn-ap" if row["network_role"] == "wireless_access_point" else "cn-router"
        candidate_id = html.escape(row["candidate_id"])
        markup.append(
            f'<g data-candidate-id="{candidate_id}"><circle class="{css}" cx="{x:.3f}" cy="{y:.3f}" r="2.2"/>'
            f'<text class="cn-label" x="{x+3.4:.3f}" y="{y+(5 if index % 2 else -3):.3f}">{candidate_id}</text></g>'
        )
    markup.extend([
        '<g><rect class="cn-panel" x="402" y="7" width="93" height="386"/>',
        '<text class="cn-title" x="407" y="16">E-302/E-304 控制与网络协调区</text>',
        '<text class="cn-note" x="407" y="24">只读候选｜非最终安装点｜不写 IFC</text>',
        '<text class="cn-text" x="407" y="39">洋红方块：2 个门口控制面板区</text>',
        '<text class="cn-text" x="407" y="47">蓝点：主卧/次卧 AP 房间区</text>',
        '<text class="cn-text" x="407" y="55">紫点：客厅＋书房共享路由器区</text>',
        '<text class="cn-text" x="407" y="69">双控：客厅＋书房＋餐厅</text>',
        '<text class="cn-text" x="407" y="77">控制方式：仅实体有线</text>',
        '<text class="cn-text" x="407" y="85">开关面板底边：1300 mm AFF</text>',
        '<text class="cn-warn" x="407" y="102">门中心只标识门口，不是最终墙上点</text>',
        '<text class="cn-warn" x="407" y="110">房间中心只标识房间，不是设备点</text>',
        f'<text class="cn-note" x="407" y="382">IFC SHA {report["source_ifc_sha256"][:12]}…</text></g>',
    ])
    style = """
@page{size:500mm 400mm;margin:0}.cn-control{fill:#d63384;stroke:#6b1742;stroke-width:.7}.cn-ap{fill:#228be6;stroke:#0b477d;stroke-width:.7}.cn-router{fill:#7048e8;stroke:#35206f;stroke-width:.7}.cn-label,.cn-title,.cn-note,.cn-text,.cn-warn{font-family:Arial,'Noto Sans CJK SC',sans-serif;fill:#102f43}.cn-label{font-size:2.2px;font-weight:700;paint-order:stroke;stroke:#fff;stroke-width:.8px}.cn-panel{fill:#fbfcfd;stroke:#102f43;stroke-width:.5}.cn-title{font-size:3.7px;font-weight:700}.cn-note{font-size:2.15px;fill:#526777}.cn-text{font-size:2.3px}.cn-warn{font-size:2.15px;fill:#c92a2a;font-weight:700}
"""
    return source.replace("</svg>", f'<style id="elec-control-network-style">{style}</style><g id="elec-control-network">{"".join(markup)}</g></svg>', 1)


def main() -> int:
    args = parse_args()
    source_hash = sha256(args.ifc)
    if source_hash != EXPECTED_IFC_SHA256:
        raise RuntimeError(f"formal IFC hash changed: {source_hash}")
    rules = read_csv(args.rules)
    controls = control_zones(read_csv(args.doors))
    networks = network_zones(read_csv(args.spaces))
    router_zones = [row for row in networks if row["network_role"] == "router_no_AP"]
    report = {
        "mode": "read_only_control_network_coordination_candidate",
        "source_ifc_sha256": source_hash,
        "rules_path": str(args.rules.resolve()),
        "summary": {
            "control_coordination_zones": len(controls),
            "paired_two_way_control_groups": 3,
            "entrance_master_lighting_switches": 1,
            "bedroom_AP_zones": sum(row["network_role"] == "wireless_access_point" for row in networks),
            "living_study_router_zones": len(router_zones),
            "confirmed_rules": sum(row["status"] == "confirmed" for row in rules),
        },
        "control_coordination_zones": controls,
        "network_coordination_zones": networks,
        "gates": {
            "two_doorway_zones_present": len(controls) == 2,
            "three_two_way_groups_present": all(row["controlled_groups"] == ["客厅", "书房", "餐厅"] for row in controls),
            "entrance_master_switch_present": sum(row["entrance_master_lighting_switch"] for row in controls) == 1,
            "two_bedroom_AP_zones_present": sum(row["network_role"] == "wireless_access_point" for row in networks) == 2,
            "one_shared_router_no_AP_zone_present": len(router_zones) == 1 and router_zones[0]["room_reference"] == "R20;R22",
            "all_positions_are_coordination_zones": all("zone_only" in row["coordinate_status"] for row in controls + networks),
            "automatic_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    args.output_svg.write_text(render_svg(args.source_svg.read_text(encoding="utf-8"), report), encoding="utf-8")
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
