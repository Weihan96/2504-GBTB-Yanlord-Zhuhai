#!/usr/bin/env python3
"""Extract the entry-cabinet router location evidence from the official handover CAD."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
from typing import Any


EXPECTED_PLAN_DWG_SHA256 = "ba355f6a90732ad07f843d59e8bab5e1da9daffe7aed74889d1d265b2ce22d7e"
EXPECTED_PLAN_DXF_SHA256 = "706483de83d90e526d7d7cf4f0b097902a4f97868c16aba5f6fac9d9e942f0f0"
EXPECTED_ELEVATION_DWG_SHA256 = "80c0e01a3e713941f729c3b01750181c18ea4d0ff4a5ae8d65800686f19646dc"
SOURCE_TO_IFC_X_OFFSET_MM = 129253.2415
SOURCE_TO_IFC_Y_OFFSET_MM = -50000.611869

GENERAL_NOTE_HANDLE = "2598C0"
WEAK_BOX_DATUM_HANDLE = "224270"
WEAK_BOX_SOCKET_HANDLE = "224295"
WEAK_BOX_LEADER_HANDLE = "224271"
WEAK_PLAN_VIEWPORT_HANDLE = "224238"
CABINET_HANDLES = ("22017E", "220182")
EXPECTED_TEXT = {
    GENERAL_NOTE_HANDLE: "4. 强电箱、弱电箱配置于过道高柜内,见机电综合点位图",
    WEAK_BOX_DATUM_HANDLE: "弱电配电箱(底边距地H+350mm)",
    WEAK_BOX_SOCKET_HANDLE: "弱电箱箱内预留五孔插座",
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    source_dir = root.parent / "图纸" / "矩阵纵横"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dwg", type=Path, default=source_dir / "01 成品房(D1户型-115)平面系统图.dwg")
    parser.add_argument("--plan-dxf", type=Path, default=root / "tmp/dwg/d1-handover-plan.dxf")
    parser.add_argument("--elevation-dwg", type=Path, default=source_dir / "02 成品房(D1户型-115)立面图.dwg")
    parser.add_argument("--evidence-register", type=Path, default=root / "pipeline/decisions/elec-source-evidence.csv")
    parser.add_argument("--output", type=Path, default=root / "build/elec/e304-router-cad-evidence.json")
    parser.add_argument("--output-svg", type=Path, default=root / "drawings/E304-router-entry-cabinet-evidence.svg")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_evidence_register(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["sheet_id"] == "E-304"]
    required_ids = {"E304-CAD-001", "E304-CAD-002", "E304-USER-001"}
    evidence_ids = [row["evidence_id"] for row in rows]
    if not required_ids.issubset(evidence_ids) or len(evidence_ids) != len(set(evidence_ids)):
        raise RuntimeError("E-304 evidence register is missing required evidence or has duplicate IDs")
    required = ("source_kind", "source_document", "source_sha256", "source_locator", "evidence", "proves", "does_not_prove", "status")
    if any(not row[field] for row in rows for field in required):
        raise RuntimeError("E-304 evidence register contains an incomplete evidence row")
    for row in rows:
        try:
            confidence = float(row["confidence"])
        except (TypeError, ValueError):
            raise RuntimeError(f"{row['evidence_id']}: confidence is not numeric")
        if not 0 <= confidence <= 1:
            raise RuntimeError(f"{row['evidence_id']}: confidence is outside 0..1")
        if row["review_required"] not in {"yes", "no"}:
            raise RuntimeError(f"{row['evidence_id']}: review_required is not yes/no")
        if row["formal_ifc_write_allowed"] != "no":
            raise RuntimeError(f"{row['evidence_id']}: evidence must not authorize IFC writes")
    return rows


def parse_dxf(path: Path) -> list[dict[str, Any]]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) % 2:
        raise RuntimeError(f"DXF has an odd group-code line count: {path}")
    pairs = [(lines[index].strip(), lines[index + 1].strip()) for index in range(0, len(lines), 2)]
    section: str | None = None
    current: dict[str, Any] | None = None
    block_name: str | None = None
    entities: list[dict[str, Any]] = []
    for code, value in pairs:
        if code == "0" and value == "SECTION":
            section = "WAITING_FOR_NAME"
            continue
        if section == "WAITING_FOR_NAME" and code == "2":
            section = value
            continue
        if code == "0" and value == "ENDSEC":
            section = None
            current = None
            block_name = None
            continue
        if section not in {"BLOCKS", "ENTITIES"}:
            continue
        if code == "0":
            if section == "BLOCKS" and value == "BLOCK":
                current = {"type": value, "section": section, "block_name": None, "attributes": {}}
                block_name = "__pending__"
                continue
            if section == "BLOCKS" and value == "ENDBLK":
                current = None
                block_name = None
                continue
            current = {"type": value, "section": section, "block_name": None if block_name == "__pending__" else block_name, "attributes": {}}
            entities.append(current)
        elif current is not None:
            current["attributes"].setdefault(code, []).append(value)
            if section == "BLOCKS" and current["type"] == "BLOCK" and code == "2":
                block_name = value
    return entities


def first(entity: dict[str, Any], code: str, default: str = "") -> str:
    return entity["attributes"].get(code, [default])[0]


def text_value(entity: dict[str, Any]) -> str:
    attributes = entity["attributes"]
    return "".join(attributes.get("3", []) + attributes.get("1", []))


def entity_by_handle(entities: list[dict[str, Any]], handle: str) -> dict[str, Any]:
    matches = [entity for entity in entities if first(entity, "5") == handle]
    if len(matches) != 1:
        raise RuntimeError(f"expected one DXF entity for handle {handle}, found {len(matches)}")
    return matches[0]


def point(entity: dict[str, Any], index: int = 0) -> list[float]:
    attributes = entity["attributes"]
    return [float(attributes["10"][index]), float(attributes["20"][index])]


def source_to_ifc(source_position_mm: list[float]) -> list[float]:
    return [
        source_position_mm[0] - SOURCE_TO_IFC_X_OFFSET_MM,
        source_position_mm[1] - SOURCE_TO_IFC_Y_OFFSET_MM,
    ]


def paper_to_model(paper_position: list[float], viewport: dict[str, Any]) -> list[float]:
    attributes = viewport["attributes"]
    paper_centre = [float(first(viewport, "10")), float(first(viewport, "20"))]
    view_centre = [float(first(viewport, "12")), float(first(viewport, "22"))]
    view_target = [float(first(viewport, "17")), float(first(viewport, "27"))]
    paper_height = float(first(viewport, "41"))
    view_height = float(first(viewport, "45"))
    twist = float(first(viewport, "51", "0"))
    if abs(twist) > 1e-9:
        raise RuntimeError(f"rotated viewport is not supported: {twist}")
    scale = view_height / paper_height
    return [
        view_target[0] + view_centre[0] + (paper_position[0] - paper_centre[0]) * scale,
        view_target[1] + view_centre[1] + (paper_position[1] - paper_centre[1]) * scale,
    ]


def polyline_points(entity: dict[str, Any]) -> list[list[float]]:
    attributes = entity["attributes"]
    return [[float(x), float(y)] for x, y in zip(attributes.get("10", []), attributes.get("20", []))]


def evidence_record(entity: dict[str, Any]) -> dict[str, Any]:
    return {
        "handle": first(entity, "5"),
        "entity_type": entity["type"],
        "section": entity["section"],
        "block_name": entity.get("block_name"),
        "layer": first(entity, "8"),
        "owner_handle": first(entity, "330"),
        "paper_space": first(entity, "67", "0") == "1",
        "insertion_point": point(entity) if first(entity, "10") and first(entity, "20") else None,
        "text": text_value(entity) or None,
    }


def render_svg(report: dict[str, Any], entities: list[dict[str, Any]]) -> str:
    width, height = 1400, 800
    plan_left, plan_top, plan_width, plan_height = 60, 130, 790, 560
    bounds = {"min_x": 2500.0, "max_x": 6500.0, "min_y": -1750.0, "max_y": 850.0}
    sx = plan_width / (bounds["max_x"] - bounds["min_x"])
    sy = plan_height / (bounds["max_y"] - bounds["min_y"])

    def screen(position: list[float]) -> tuple[float, float]:
        x, y = source_to_ifc(position)
        return plan_left + (x - bounds["min_x"]) * sx, plan_top + (bounds["max_y"] - y) * sy

    geometry: list[str] = []
    relevant_layers = ("FF1-固定家具", "FF-DOOR", "EM-ELEC", "FF1-平面中英文字")
    for entity in entities:
        if entity["section"] != "ENTITIES" or first(entity, "67", "0") != "0":
            continue
        layer = first(entity, "8")
        if not any(layer.startswith(prefix) for prefix in relevant_layers):
            continue
        points = polyline_points(entity)
        if entity["type"] == "LINE" and first(entity, "11") and first(entity, "21"):
            points.append([float(first(entity, "11")), float(first(entity, "21"))])
        if not points:
            continue
        if not any(bounds["min_x"] <= source_to_ifc(item)[0] <= bounds["max_x"] and bounds["min_y"] <= source_to_ifc(item)[1] <= bounds["max_y"] for item in points):
            continue
        screen_points = [screen(item) for item in points]
        colour = "#14213d" if layer.startswith("FF1-固定家具") else "#5c677d" if layer.startswith("FF-DOOR") else "#d97706" if layer == "EM-ELEC" else "#334155"
        if entity["type"] in {"LINE", "LWPOLYLINE"}:
            close = entity["type"] == "LWPOLYLINE" and int(first(entity, "70", "0")) & 1
            command = "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in screen_points) + (" Z" if close else "")
            geometry.append(f'<path d="{command}" fill="none" stroke="{colour}" stroke-width="2"/>')
        elif entity["type"] in {"TEXT", "MTEXT"} and text_value(entity):
            x, y = screen_points[0]
            geometry.append(f'<text x="{x:.2f}" y="{y:.2f}" class="plan-text">{html.escape(text_value(entity))}</text>')

    cabinet_markup: list[str] = []
    for cabinet in report["cabinet_bays"]:
        points = [screen(item) for item in cabinet["source_polygon_mm"]]
        command = "M " + " L ".join(f"{x:.2f} {y:.2f}" for x, y in points) + " Z"
        cabinet_markup.append(f'<path d="{command}" class="cabinet" data-handle="{cabinet["handle"]}"/>')
    marker_source = report["weak_current_box"]["source_plan_position_mm"]
    marker_x, marker_y = screen(marker_source)
    source = report["source"]
    exact = report["weak_current_box"]
    hashes = {
        "plan": source["official_plan_dwg"]["sha256"][:16],
        "dxf": source["verified_conversion_dxf"]["sha256"][:16],
        "elevation": source["user_supplied_elevation_dwg"]["sha256"][:16],
    }
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="800" viewBox="0 0 {width} {height}">
<rect width="1400" height="800" fill="#f8fafc"/>
<text x="60" y="60" class="title">E-304｜玄关柜路由器定位取证</text>
<text x="60" y="96" class="subtitle">官方 CAD 机械提取｜只确认平面柜位，不虚构产品与安装高度</text>
<rect x="40" y="115" width="830" height="610" rx="16" class="panel"/>
<text x="65" y="155" class="section">玄关 / 过道高柜平面局部</text>
{''.join(geometry)}
{''.join(cabinet_markup)}
<circle cx="{marker_x:.2f}" cy="{marker_y:.2f}" r="14" class="marker"/>
<path d="M {marker_x + 16:.2f} {marker_y - 8:.2f} L {marker_x + 125:.2f} {marker_y - 70:.2f}" class="leader"/>
<text x="{marker_x + 133:.2f}" y="{marker_y - 78:.2f}" class="marker-label">弱电箱 / 路由器柜位</text>
<text x="{marker_x + 133:.2f}" y="{marker_y - 50:.2f}" class="marker-note">IFC XY ≈ ({exact["ifc_plan_position_mm"][0]:.3f}, {exact["ifc_plan_position_mm"][1]:.3f}) mm</text>
<text x="65" y="700" class="footnote">蓝色轮廓：玄关高柜两柜格｜紫色：弱电标注引线落点｜橙色：原 CAD 电气图层</text>
<rect x="900" y="115" width="460" height="610" rx="16" class="panel"/>
<text x="930" y="160" class="section">取证结论</text>
<text x="930" y="205" class="body">1. “强电箱、弱电箱配置于过道高柜内”</text>
<text x="950" y="234" class="meta">TEXT #{GENERAL_NOTE_HANDLE}｜DIML【标注引线】</text>
<text x="930" y="285" class="body">2. “弱电配电箱(底边距地H+350mm)”</text>
<text x="950" y="314" class="meta">TEXT #{WEAK_BOX_DATUM_HANDLE}｜EM-TEXT</text>
<text x="930" y="365" class="body">3. “弱电箱箱内预留五孔插座”</text>
<text x="950" y="394" class="meta">TEXT #{WEAK_BOX_SOCKET_HANDLE}｜EM-TEXT</text>
<text x="930" y="448" class="decision">已确认：路由器放玄关/过道高柜弱电箱柜位</text>
<text x="930" y="482" class="warning">未确认：路由器 Z、层板、产品、散热及检修净空</text>
<text x="930" y="516" class="warning">H+350 是弱电箱底边，不是路由器安装高度</text>
<text x="930" y="570" class="meta">01 平面系统图 DWG  {hashes["plan"]}…</text>
<text x="930" y="599" class="meta">验证转换 DXF       {hashes["dxf"]}…</text>
<text x="930" y="628" class="meta">02 立面图 DWG       {hashes["elevation"]}…</text>
<text x="930" y="669" class="meta">02 文件已登记哈希；本平面定位证据取自 01 的 E-2 视口</text>
<style>
.title{{font:700 32px Arial,'Noto Sans CJK SC',sans-serif;fill:#0f172a}}.subtitle{{font:18px Arial,'Noto Sans CJK SC',sans-serif;fill:#475569}}
.panel{{fill:#fff;stroke:#cbd5e1;stroke-width:2}}.section{{font:700 22px Arial,'Noto Sans CJK SC',sans-serif;fill:#0f172a}}
.plan-text{{font:16px Arial,'Noto Sans CJK SC',sans-serif;fill:#334155}}.cabinet{{fill:#dbeafe;fill-opacity:.55;stroke:#2563eb;stroke-width:5}}
.marker{{fill:#7c3aed;stroke:#4c1d95;stroke-width:5}}.leader{{fill:none;stroke:#7c3aed;stroke-width:4}}.marker-label{{font:700 19px Arial,'Noto Sans CJK SC',sans-serif;fill:#4c1d95}}
.marker-note,.footnote{{font:16px Arial,'Noto Sans CJK SC',sans-serif;fill:#475569}}.body{{font:18px Arial,'Noto Sans CJK SC',sans-serif;fill:#1e293b}}
.meta{{font:15px Menlo,'Noto Sans CJK SC',monospace;fill:#64748b}}.decision{{font:700 18px Arial,'Noto Sans CJK SC',sans-serif;fill:#166534}}
.warning{{font:700 17px Arial,'Noto Sans CJK SC',sans-serif;fill:#b91c1c}}
</style></svg>'''


def main() -> int:
    args = parse_args()
    source_hashes = {
        "plan_dwg": sha256(args.plan_dwg),
        "plan_dxf": sha256(args.plan_dxf),
        "elevation_dwg": sha256(args.elevation_dwg),
    }
    expected = {
        "plan_dwg": EXPECTED_PLAN_DWG_SHA256,
        "plan_dxf": EXPECTED_PLAN_DXF_SHA256,
        "elevation_dwg": EXPECTED_ELEVATION_DWG_SHA256,
    }
    if source_hashes != expected:
        raise RuntimeError(f"official CAD evidence hashes changed: {source_hashes}")
    evidence_register = read_evidence_register(args.evidence_register)

    entities = parse_dxf(args.plan_dxf)
    selected = {handle: entity_by_handle(entities, handle) for handle in (*EXPECTED_TEXT, WEAK_BOX_LEADER_HANDLE, WEAK_PLAN_VIEWPORT_HANDLE, *CABINET_HANDLES)}
    for handle, expected_text in EXPECTED_TEXT.items():
        if text_value(selected[handle]) != expected_text:
            raise RuntimeError(f"DXF evidence text changed at handle {handle}")
    leader = selected[WEAK_BOX_LEADER_HANDLE]
    if leader["type"] != "LEADER" or first(leader, "67") != "1":
        raise RuntimeError("weak-current-box leader is no longer a paper-space LEADER")
    viewport = selected[WEAK_PLAN_VIEWPORT_HANDLE]
    if viewport["type"] != "VIEWPORT" or first(viewport, "69") != "4":
        raise RuntimeError("weak-current plan viewport identity changed")
    source_position = paper_to_model(point(leader), viewport)
    ifc_position = source_to_ifc(source_position)
    expected_ifc_position = [4600.016493, -735.368749]
    if max(abs(actual - target) for actual, target in zip(ifc_position, expected_ifc_position)) > 0.001:
        raise RuntimeError(f"weak-current-box coordinate changed: {ifc_position}")

    cabinet_bays = []
    for handle in CABINET_HANDLES:
        entity = selected[handle]
        source_polygon = polyline_points(entity)
        if entity["type"] != "LWPOLYLINE" or len(source_polygon) != 4:
            raise RuntimeError(f"entry cabinet geometry changed at handle {handle}")
        ifc_polygon = [source_to_ifc(item) for item in source_polygon]
        cabinet_bays.append({
            "handle": handle,
            "layer": first(entity, "8"),
            "source_polygon_mm": source_polygon,
            "ifc_polygon_mm": ifc_polygon,
            "ifc_bbox_mm": [
                min(item[0] for item in ifc_polygon), min(item[1] for item in ifc_polygon),
                max(item[0] for item in ifc_polygon), max(item[1] for item in ifc_polygon),
            ],
        })
    router_bay = cabinet_bays[1]
    report = {
        "mode": "read_only_official_cad_evidence",
        "source": {
            "official_plan_dwg": {"path": str(args.plan_dwg.resolve()), "sha256": source_hashes["plan_dwg"]},
            "verified_conversion_dxf": {"path": str(args.plan_dxf.resolve()), "sha256": source_hashes["plan_dxf"]},
            "user_supplied_elevation_dwg": {
                "path": str(args.elevation_dwg.resolve()),
                "sha256": source_hashes["elevation_dwg"],
                "evidence_role": "file identity recorded; not used as the plan-coordinate source",
            },
            "coordinate_transform": {
                "ifc_x_mm": "source_x_mm - 129253.2415",
                "ifc_y_mm": "source_y_mm + 50000.611869",
                "viewport_handle": WEAK_PLAN_VIEWPORT_HANDLE,
                "viewport_number": 4,
                "viewport_scale": 50.0,
            },
        },
        "evidence_register": {
            "path": str(args.evidence_register.resolve()),
            "rows": evidence_register,
        },
        "text_evidence": [evidence_record(selected[handle]) for handle in EXPECTED_TEXT],
        "leader_evidence": evidence_record(leader),
        "cabinet_bays": cabinet_bays,
        "weak_current_box": {
            "source_plan_position_mm": source_position,
            "ifc_plan_position_mm": ifc_position,
            "weak_current_box_bottom_aff_mm": 350.0,
            "cabinet_bbox_ifc_mm": router_bay["ifc_bbox_mm"],
        },
        "router_decision": {
            "candidate_id": "NET-ROUTER-LIVING-STUDY",
            "location": "entry/corridor high cabinet at the official weak-current-box bay",
            "plan_position_mm": ifc_position,
            "installation_z_mm": None,
            "served_room_references": ["R20", "R22"],
            "wired_backhaul": True,
            "status": "plan_position_confirmed_product_z_power_data_thermal_service_pending",
            "automatic_ifc_write_allowed": False,
        },
        "gates": {
            "official_source_hashes_match": True,
            "three_text_entities_match": True,
            "paper_viewport_transform_verified": True,
            "entry_cabinet_geometry_verified": True,
            "router_z_not_inferred_from_weak_box_datum": True,
            "automatic_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(render_svg(report, entities), encoding="utf-8")
    print(json.dumps({"weak_current_box": report["weak_current_box"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
