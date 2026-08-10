#!/usr/bin/env python3
"""Compile the first read-only renovation electrical demand candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
from pathlib import Path
from typing import Any, Callable

from svg_audit_underlay import strip_wall_plan_raster_underlay

EXPECTED_IFC_SHA256 = "6c2fd8da9e9ad7ddbc2b63415a27f1c979e8995b880d8fce210a2dda2ef2aab6"
SCALE_DENOMINATOR = 50.0
SVG_WORLD_OFFSET_MM = 10000.0


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
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


def new_socket_candidates(int1_records: list[dict[str, Any]]) -> list[dict[str, Any]]:
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
            "automatic_ifc_write_allowed": False,
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
            "automatic_ifc_write_allowed": False,
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
    source = strip_wall_plan_raster_underlay(source, "E-301/E-303 renovation round 1")
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
    markup.append('<text class="r1-note" x="407" y="24">纯矢量审核底图｜只读候选｜不写 IFC</text>')
    markup.append(f'<text class="r1-text" x="407" y="37">青色床头灯候选：{summary["bedside_light_candidates"]}</text>')
    markup.append(f'<text class="r1-text" x="407" y="44">绿色新增插座候选：{summary["new_socket_candidates"]}</text>')
    markup.append(f'<text class="r1-text" x="407" y="51">琥珀色柜体供电协调区：{summary["cabinet_power_zones"]}</text>')
    markup.append(f'<text class="r1-text" x="407" y="58">红色厨房现有插座复核：{summary["kitchen_socket_rechecks"]}</text>')
    markup.append('<text class="r1-warn" x="407" y="72">柜体标记是供电范围，不是最终出线口</text>')
    markup.append('<text class="r1-warn" x="407" y="79">开发商红色旧点位已降级为参考</text>')
    markup.append('<text class="r1-note" x="407" y="92">高度、回路、防水、设备功率和柜体开孔未冻结</text>')
    markup.append(f'<text class="r1-note" x="407" y="382">IFC SHA {report["source_ifc_sha256"][:12]}…</text></g>')
    style = """
@page{size:500mm 400mm;margin:0}html,body{margin:0;width:500mm;height:400mm;overflow:hidden}
.r1-bedside{fill:#00d5e8;stroke:#006672;stroke-width:.7}.r1-new-socket{fill:#37b24d;stroke:#145c23;stroke-width:.7}
.r1-cabinet{fill:#ffb000;stroke:#7a4a00;stroke-width:.7}.r1-recheck{fill:#fa2b2b;stroke:#7b0000;stroke-width:.7}
.r1-label{font-family:Arial,'Noto Sans CJK SC',sans-serif;font-size:2.15px;font-weight:700;fill:#102f43;paint-order:stroke;stroke:#fff;stroke-width:.8px}
.r1-panel{fill:#fbfcfd;stroke:#102f43;stroke-width:.5}.r1-title,.r1-note,.r1-text,.r1-warn{font-family:Arial,'Noto Sans CJK SC',sans-serif;fill:#102f43}
.r1-title{font-size:3.7px;font-weight:700}.r1-note{font-size:2.15px;fill:#526777}.r1-text{font-size:2.3px}.r1-warn{font-size:2.15px;fill:#c92a2a;font-weight:700}
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
    ifc_hash = sha256(args.ifc)
    if ifc_hash != EXPECTED_IFC_SHA256:
        raise RuntimeError(f"formal IFC hash changed: {ifc_hash}")
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
    new_sockets = new_socket_candidates(int1["records"])
    cabinet_zones = cabinet_power_zones(int1["records"])
    socket_rechecks = kitchen_socket_rechecks(existing, positioning)
    requirement_status = {row["requirement_id"]: row["status"] for row in requirements}
    report = {
        "mode": "read_only_renovation_electrical_round1_candidate",
        "source_ifc_sha256": ifc_hash,
        "requirements_path": str(args.requirements.resolve()),
        "summary": {
            "confirmed_requirements": sum(value.startswith("confirmed") for value in requirement_status.values()),
            "bedside_light_candidates": len(bedside),
            "new_socket_candidates": len(new_sockets),
            "cabinet_power_zones": len(cabinet_zones),
            "kitchen_socket_rechecks": len(socket_rechecks),
            "developer_red_points_are_reference_only": requirement_status.get("ELEC-R1-001") == "confirmed",
        },
        "requirements": requirements,
        "bedside_light_candidates": bedside,
        "new_socket_candidates": new_sockets,
        "cabinet_power_zones": cabinet_zones,
        "kitchen_socket_rechecks": socket_rechecks,
        "gates": {
            "four_bedside_lights_present": len(bedside) == 4,
            "island_and_dining_bay_socket_present": {row["candidate_role"] for row in new_sockets}
            == {"island_end_panel_socket", "dining_bay_socket"},
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
    args.output_svg.write_text(
        render_svg(args.source_svg.read_text(encoding="utf-8"), report),
        encoding="utf-8",
    )
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
