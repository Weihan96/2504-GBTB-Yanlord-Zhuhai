#!/usr/bin/env python3
"""Extract the developer handover MEP points into the current IFC coordinate system."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import numpy as np
from ifcopenshell.util.element import get_psets
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

from svg_audit_underlay import validate_electrical_coordination_source


EXPECTED_SOURCE_DXF_SHA256 = "51dc592477cb772da36e643d507877f0cc28c253797c7bae0c913dcaeaabf94f"
EXPECTED_DWG_CONVERSION_SHA256 = "706483de83d90e526d7d7cf4f0b097902a4f97868c16aba5f6fac9d9e942f0f0"

# The source plan repeats the project Grid in millimetres.  Grid 07 is x=-100
# in IFC and x=129153.2415 in the source; y=-5300 is y=-55300.611869.
SOURCE_TO_IFC_X_OFFSET_MM = 129253.2415
SOURCE_TO_IFC_Y_OFFSET_MM = -50000.611869
SOURCE_BOUNDS = (120000.0, -57000.0, 138000.0, -44000.0)
SOURCE_LAYERS = {
    "EM-ELEC",
    "CW-SWITCH",
    "DG-DRAINAGE",
    "FC-FLOOR DRAIN【地漏】",
}
EXPECTED_COUNTS = {
    "developer_electrical_point": 66,
    "developer_switch_or_control": 11,
    "developer_water_or_wc_drain": 20,
    "developer_floor_drain": 6,
}
SCALE_DENOMINATOR = 50.0
SVG_WORLD_OFFSET_MM = 10000.0


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    source_root = root.parent / "图纸" / "矩阵纵横"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument(
        "--source-dxf",
        type=Path,
        default=source_root / "(D1户型-115)建筑资料平面图.dxf",
    )
    parser.add_argument("--converted-dxf", type=Path, default=root / "tmp/dwg/d1-handover-plan.dxf")
    parser.add_argument(
        "--source-svg", type=Path, default=root / "drawings/Electrical Coordination Plan.svg"
    )
    parser.add_argument(
        "--elec-svg",
        type=Path,
        default=root / "drawings/E302-E304-developer-handover-reference.svg",
    )
    parser.add_argument(
        "--plum-svg",
        type=Path,
        default=root / "drawings/P201-developer-handover-reference.svg",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "build/mep-positioning/developer-handover-mep-candidate.json",
    )
    parser.add_argument("--expected-ifc-sha256", help="Optional caller-frozen source hash")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_ascii_dxf(path: Path) -> list[dict[str, Any]]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) % 2:
        raise RuntimeError(f"DXF has an odd group-code line count: {path}")
    pairs = [(lines[index].strip(), lines[index + 1].strip()) for index in range(0, len(lines), 2)]
    section: str | None = None
    current: dict[str, Any] | None = None
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
            continue
        if section != "ENTITIES":
            continue
        if code == "0":
            current = {"type": value, "attributes": defaultdict(list)}
            entities.append(current)
        elif current is not None:
            current["attributes"][code].append(value)
    return entities


def insert_records(path: Path) -> list[dict[str, Any]]:
    minimum_x, minimum_y, maximum_x, maximum_y = SOURCE_BOUNDS
    records = []
    for entity in parse_ascii_dxf(path):
        attributes = entity["attributes"]
        if entity["type"] != "INSERT" or attributes.get("67", ["0"])[0] != "0":
            continue
        layer = attributes.get("8", [""])[0]
        if layer not in SOURCE_LAYERS:
            continue
        x = float(attributes.get("10", ["0"])[0])
        y = float(attributes.get("20", ["0"])[0])
        if not (minimum_x < x < maximum_x and minimum_y < y < maximum_y):
            continue
        records.append(
            {
                "source_layer": layer,
                "source_block": attributes.get("2", [""])[0],
                "source_position_mm": [x, y],
                "source_rotation_deg": float(attributes.get("50", ["0"])[0]),
            }
        )
    return records


def source_text_records(path: Path) -> list[dict[str, Any]]:
    records = []
    for entity in parse_ascii_dxf(path):
        attributes = entity["attributes"]
        if entity["type"] not in {"TEXT", "MTEXT"} or attributes.get("67", ["0"])[0] != "0":
            continue
        if attributes.get("8", [""])[0] != "EM-ELEC":
            continue
        x = float(attributes.get("10", ["0"])[0])
        y = float(attributes.get("20", ["0"])[0])
        if not (SOURCE_BOUNDS[0] < x < SOURCE_BOUNDS[2] and SOURCE_BOUNDS[1] < y < SOURCE_BOUNDS[3]):
            continue
        value = "".join(attributes.get("3", []) + attributes.get("1", []))
        if re.fullmatch(r"\+[0-9]+", value):
            records.append({"position_mm": [x, y], "height_mm": float(value[1:]), "text": value})
    return records


def source_to_ifc(source_position_mm: list[float]) -> list[float]:
    return [
        source_position_mm[0] - SOURCE_TO_IFC_X_OFFSET_MM,
        source_position_mm[1] - SOURCE_TO_IFC_Y_OFFSET_MM,
    ]


def classify(record: dict[str, Any]) -> tuple[str, str, float]:
    layer = record["source_layer"]
    block = record["source_block"]
    if layer == "DG-DRAINAGE":
        return {
            "冷水": ("cold_water_rough_in_reference", "developer symbol explicitly names cold water", 1.0),
            "热水": ("hot_water_rough_in_reference", "developer symbol explicitly names hot water", 1.0),
            "马桶落水": ("wc_drain_reference", "developer symbol explicitly names WC drainage", 1.0),
        }[block]
    if layer == "FC-FLOOR DRAIN【地漏】":
        return "floor_drain_reference", "developer layer explicitly identifies floor drains", 1.0
    if layer == "CW-SWITCH":
        if block in {"A$C596C36D7", "A$C0A3239C7"}:
            return "ac_control_panel_reference", "block contains the literal AC", 0.95
        if block == "A$C16FC1894":
            return "bathroom_warm_air_panel_reference", "block text identifies an integrated warm-air panel", 0.95
        return "switch_or_control_panel_reference", "developer switch layer and symbol position", 0.65
    if block == "FL11DJGVE11电力平面图$0$E-BS201":
        return "general_power_point_reference", "repeated developer E-BS201 power symbol", 0.85
    if block in {"A$C349922F2", "A$C77E06134"}:
        return "grouped_power_point_reference", "wrapper contains the same E-BS201 power symbol", 0.80
    if block == "无线AP墙插面板":
        return "wireless_ap_wall_point_reference", "block name explicitly identifies a wireless AP wall panel", 1.0
    if block == "A$C3B9C162C":
        return "tv_point_reference", "block contains the literal TV", 0.95
    if block == "$equip$00003047":
        return "tv_telephone_box_reference", "block attribute definition identifies TV and telephone", 0.95
    if block == "A$C34D428B3":
        return "equipotential_terminal_reference", "block tag is LEB", 0.90
    if block == "$equip_U$00000028":
        return "button_reference", "block attribute definition identifies a button", 0.80
    return "unresolved_developer_electrical_symbol", "source layer and position are reliable but block identity is unresolved", 0.35


def point_kind(layer: str) -> str:
    if layer == "EM-ELEC":
        return "developer_electrical_point"
    if layer == "CW-SWITCH":
        return "developer_switch_or_control"
    if layer == "DG-DRAINAGE":
        return "developer_water_or_wc_drain"
    return "developer_floor_drain"


def space_footprints(model: ifcopenshell.file) -> list[tuple[Any, Any, str]]:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    records = []
    for space in model.by_type("IfcSpace"):
        shape = ifcopenshell.geom.create_shape(settings, space)
        vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3))
        faces = np.asarray(shape.geometry.faces, dtype=int).reshape((-1, 3))
        polygons = [Polygon(vertices[face, :2]) for face in faces]
        footprint = unary_union([polygon for polygon in polygons if polygon.area > 1e-10]).buffer(1e-8)
        reference = str(get_psets(space).get("Pset_SpaceCommon", {}).get("Reference", ""))
        records.append((space, footprint, reference))
    return records


def candidate_space(position_mm: list[float], footprints: list[tuple[Any, Any, str]]) -> dict[str, Any]:
    point = Point(position_mm[0] / 1000.0, position_mm[1] / 1000.0)
    covering = [(space, reference) for space, footprint, reference in footprints if footprint.covers(point)]
    distances = sorted(
        ((footprint.distance(point) * 1000.0, space, reference) for space, footprint, reference in footprints),
        key=lambda row: (row[0], row[2], row[1].GlobalId),
    )
    nearest_distance, nearest, nearest_reference = distances[0]
    return {
        "candidate_reference": covering[0][1] if len(covering) == 1 else nearest_reference,
        "candidate_global_id": covering[0][0].GlobalId if len(covering) == 1 else nearest.GlobalId,
        "candidate_long_name": str((covering[0][0] if len(covering) == 1 else nearest).LongName or ""),
        "relationship": "contained" if len(covering) == 1 else ("boundary_ambiguous" if covering else "nearest_only"),
        "covering_space_count": len(covering),
        "nearest_distance_mm": round(nearest_distance, 6),
        "review_required": True,
    }


def nearest_height(record: dict[str, Any], texts: list[dict[str, Any]]) -> dict[str, Any]:
    point = record["source_position_mm"]
    candidates = sorted(
        (
            (math.dist(point, row["position_mm"]), row)
            for row in texts
        ),
        key=lambda item: item[0],
    )
    if not candidates or candidates[0][0] > 260.0:
        return {"installation_height_mm": None, "height_basis": "no unique nearby height text", "height_confidence": 0.0}
    distance, row = candidates[0]
    return {
        "installation_height_mm": row["height_mm"],
        "height_basis": f"nearest developer +height text at {distance:.3f} mm in source coordinates",
        "height_confidence": 0.70,
    }


def assign_ids(records: list[dict[str, Any]]) -> None:
    prefixes = {
        "developer_electrical_point": "DEV-E",
        "developer_switch_or_control": "DEV-S",
        "developer_water_or_wc_drain": "DEV-W",
        "developer_floor_drain": "DEV-D",
    }
    grouped: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[record["kind"]].append(record)
    for kind, rows in grouped.items():
        rows.sort(key=lambda row: (-row["ifc_position_mm"][1], row["ifc_position_mm"][0], row["source_block"]))
        for index, row in enumerate(rows, 1):
            row["candidate_id"] = f"{prefixes[kind]}{index:03d}"


def world_to_svg(position_mm: list[float]) -> tuple[float, float]:
    return (
        (position_mm[0] + SVG_WORLD_OFFSET_MM) / SCALE_DENOMINATOR,
        (SVG_WORLD_OFFSET_MM - position_mm[1]) / SCALE_DENOMINATOR,
    )


def inject_svg(source: str, generated: str, style: str, element_id: str) -> str:
    source = re.sub(r'width="400(?:\.0+)?mm"', 'width="500mm"', source, count=1)
    source = re.sub(r'viewBox="0 0 400(?:\.0+)? 400(?:\.0+)?"', 'viewBox="0 0 500 400"', source, count=1)
    if "</svg>" not in source or 'viewBox="0 0 500 400"' not in source:
        raise RuntimeError("could not prepare the 500x400 MEP reference SVG")
    return source.replace(
        "</svg>",
        f'<style id="{element_id}-style">{style}</style><g id="{element_id}">{generated}</g></svg>',
        1,
    )


def render_reference_svg(source: str, records: list[dict[str, Any]], discipline: str, ifc_hash: str) -> str:
    if discipline == "ELEC":
        rows = [row for row in records if row["kind"] in {"developer_electrical_point", "developer_switch_or_control"}]
        title = "E-302/E-304 开发商交付水电参考层"
    else:
        rows = [row for row in records if row["kind"] in {"developer_water_or_wc_drain", "developer_floor_drain"}]
        title = "P-201 开发商交付给排水参考层"
    markup = []
    for row in rows:
        x, y = world_to_svg(row["ifc_position_mm"])
        css = "source-" + row["kind"].replace("developer_", "").replace("_", "-")
        markup.append(
            f'<g data-candidate-id="{html.escape(row["candidate_id"])}" data-source-block="{html.escape(row["source_block"])}">'
            f'<circle class="{css}" cx="{x:.3f}" cy="{y:.3f}" r="1.8"/>'
            f'<title>{html.escape(row["candidate_id"] + " | " + row["candidate_role"] + " | " + row["candidate_space"]["candidate_reference"])}</title>'
            "</g>"
        )
    counts = Counter(row["kind"] for row in rows)
    markup.append('<g><rect class="source-panel" x="402" y="7" width="93" height="386"/>')
    markup.append(f'<text class="source-title" x="407" y="16">{html.escape(title)}</text>')
    markup.append('<text class="source-note" x="407" y="24">改造前参考｜非施工发布｜不写 IFC</text>')
    y = 36.0
    for kind, count in sorted(counts.items()):
        markup.append(f'<text class="source-text" x="407" y="{y:.1f}">{html.escape(kind)}：{count}</text>')
        y += 6.0
    markup.append('<text class="source-warn" x="407" y="82">只证明交付位置，不等于装修后保留位置</text>')
    markup.append('<text class="source-warn" x="407" y="89">未解析符号、回路、管径和连接关系保持待审</text>')
    markup.append(f'<text class="source-note" x="407" y="382">IFC SHA {ifc_hash[:12]}…</text></g>')
    style = """
@page{size:500mm 400mm;margin:0}html,body{margin:0;width:500mm;height:400mm;overflow:hidden}
.source-electrical-point{fill:#845ef7;stroke:#3b1c8c;stroke-width:.55}.source-switch-or-control{fill:#ff922b;stroke:#a94800;stroke-width:.55}
.source-water-or-wc-drain{fill:#228be6;stroke:#0b4f8a;stroke-width:.55}.source-floor-drain{fill:#22b8cf;stroke:#075b68;stroke-width:.55}
.source-panel{fill:#fbfcfd;stroke:#102f43;stroke-width:.5}.source-title,.source-note,.source-text,.source-warn{font-family:Arial,'Noto Sans CJK SC',sans-serif;fill:#102f43}
.source-title{font-size:3.7px;font-weight:700}.source-note{font-size:2.15px;fill:#526777}.source-text{font-size:2.3px}.source-warn{font-size:2.15px;fill:#c92a2a;font-weight:700}
"""
    return inject_svg(source, "".join(markup), style, f"developer-{discipline.lower()}-reference")


def main() -> int:
    args = parse_args()
    ifc_hash = sha256(args.ifc)
    source_hash = sha256(args.source_dxf)
    if args.expected_ifc_sha256 and ifc_hash != args.expected_ifc_sha256:
        raise RuntimeError(
            f"formal IFC hash changed: expected {args.expected_ifc_sha256}, got {ifc_hash}"
        )
    if source_hash != EXPECTED_SOURCE_DXF_SHA256:
        raise RuntimeError(f"developer DXF hash changed: {source_hash}")

    source_inserts = insert_records(args.source_dxf)
    converted_equivalence = None
    if args.converted_dxf.exists():
        converted_hash = sha256(args.converted_dxf)
        converted_records = insert_records(args.converted_dxf)
        source_keys = [
            (row["source_layer"], row["source_block"], *[round(value, 6) for value in row["source_position_mm"]], round(row["source_rotation_deg"], 6))
            for row in source_inserts
        ]
        converted_keys = [
            (row["source_layer"], row["source_block"], *[round(value, 6) for value in row["source_position_mm"]], round(row["source_rotation_deg"], 6))
            for row in converted_records
        ]
        converted_equivalence = {
            "sha256": converted_hash,
            "expected_sha256": EXPECTED_DWG_CONVERSION_SHA256,
            "in_scope_records_equal": Counter(source_keys) == Counter(converted_keys),
            "in_scope_record_count": len(converted_keys),
        }

    model = ifcopenshell.open(args.ifc)
    footprints = space_footprints(model)
    height_texts = source_text_records(args.source_dxf)
    records = []
    for source_record in source_inserts:
        role, basis, confidence = classify(source_record)
        ifc_position = source_to_ifc(source_record["source_position_mm"])
        kind = point_kind(source_record["source_layer"])
        record = {
            **source_record,
            "kind": kind,
            "ifc_position_mm": [round(value, 6) for value in ifc_position],
            "candidate_role": role,
            "role_basis": basis,
            "role_confidence": confidence,
            "candidate_space": candidate_space(ifc_position, footprints),
            "source_status": "developer_handover_existing_reference",
            "review_required": True,
            "automatic_ifc_write_allowed": False,
        }
        if kind in {"developer_electrical_point", "developer_switch_or_control"}:
            record.update(nearest_height(source_record, height_texts))
        records.append(record)
    assign_ids(records)

    counts = Counter(row["kind"] for row in records)
    room_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    for row in records:
        room_counts[row["candidate_space"]["candidate_reference"]][row["kind"]] += 1
    report = {
        "source_ifc_sha256": ifc_hash,
        "source": {
            "dxf_path": str(args.source_dxf.resolve()),
            "dxf_sha256": source_hash,
            "dwg_conversion_evidence": converted_equivalence,
            "coordinate_transform": {
                "ifc_x_mm": "source_x_mm - 129253.2415",
                "ifc_y_mm": "source_y_mm + 50000.611869",
                "basis": "source Grid 07 x=129153.2415 maps to IFC x=-100; source Grid J y=-55300.611869 maps to IFC y=-5300",
            },
        },
        "summary": {
            "point_counts": dict(sorted(counts.items())),
            "unresolved_electrical_symbols": sum(row["candidate_role"] == "unresolved_developer_electrical_symbol" for row in records),
            "points_with_source_height": sum(row.get("installation_height_mm") is not None for row in records),
            "points_with_unique_current_space": sum(row["candidate_space"]["relationship"] == "contained" for row in records),
            "room_counts": {reference: dict(sorted(value.items())) for reference, value in sorted(room_counts.items())},
        },
        "records": records,
        "gates": {
            "expected_point_counts_pass": dict(counts) == EXPECTED_COUNTS,
            "source_matches_dwg_conversion": bool(converted_equivalence and converted_equivalence["in_scope_records_equal"]),
            "coordinate_transform_documented": True,
            "all_points_have_space_candidate": all(row["candidate_space"]["candidate_reference"] for row in records),
            "developer_reference_is_renovation_design": False,
            "automatic_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    source_svg = args.source_svg.read_text(encoding="utf-8")
    validate_electrical_coordination_source(source_svg, args.source_svg, args.ifc)
    args.elec_svg.write_text(render_reference_svg(source_svg, records, "ELEC", ifc_hash), encoding="utf-8")
    args.plum_svg.write_text(render_reference_svg(source_svg, records, "PLUM", ifc_hash), encoding="utf-8")
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
