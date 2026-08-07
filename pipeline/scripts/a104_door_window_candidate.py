#!/usr/bin/env python3
"""Generate the read-only A-104 door/window plan and review register.

The generator inventories every IfcDoor and IfcWindow, validates hosted
opening dimensions against IFC nominal sizes, assigns deterministic candidate
numbers, and records every inference with its evidence and review boundary.
It never writes the IFC.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import ifcopenshell
import ifcopenshell.util.element

from a103_wall_plan_candidate import (
    Box,
    format_mm,
    geometry_settings,
    place_label,
    sha256,
    world_bbox_mm,
    world_to_svg,
)


MASTER_BEDROOM_DOOR_CANDIDATE = "3xKBbA2CT9mfby2$MODnzM"
CONFIRMED_GUEST_BATHROOM_DOOR = "1To6gRjrf2nBU5GQuuD8FV"
UNHOSTED_DOORS = {
    "1TW6$_GfnABRZusYvx0zZG",
    "2D5BPoo2XFSvhTdfPenCh7",
    "0zjVS5FBbBewgUkk0fdfiv",
}
BATHROOM_UNHOSTED_PAIR = {
    "2D5BPoo2XFSvhTdfPenCh7",
    "0zjVS5FBbBewgUkk0fdfiv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--source-svg", required=True, type=Path)
    parser.add_argument("--output-svg", required=True, type=Path)
    parser.add_argument("--review-register", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def svg_escape(value: Any) -> str:
    return html.escape(str(value), quote=True)


def read_bbim_data(product: ifcopenshell.entity_instance) -> dict[str, Any]:
    psets = ifcopenshell.util.element.get_psets(product)
    key = "BBIM_Door" if product.is_a("IfcDoor") else "BBIM_Window"
    raw = psets.get(key, {}).get("Data")
    if not raw:
        return {}
    try:
        value = json.loads(str(raw))
    except json.JSONDecodeError as error:
        raise RuntimeError(f"invalid {key}.Data JSON on {product.GlobalId}") from error
    return value if isinstance(value, dict) else {}


def product_container(product: ifcopenshell.entity_instance) -> dict[str, str]:
    relations = list(getattr(product, "ContainedInStructure", ()) or ())
    if len(relations) != 1:
        return {"global_id": "", "name": "", "class": ""}
    container = relations[0].RelatingStructure
    return {
        "global_id": str(container.GlobalId),
        "name": str(container.Name or ""),
        "class": str(container.is_a()),
    }


def fill_chain(product: ifcopenshell.entity_instance) -> dict[str, Any]:
    fills = list(getattr(product, "FillsVoids", ()) or ())
    opening = fills[0].RelatingOpeningElement if len(fills) == 1 else None
    voids = list(getattr(opening, "VoidsElements", ()) or ()) if opening else []
    host = voids[0].RelatingBuildingElement if len(voids) == 1 else None
    return {
        "fill_relationship_count": len(fills),
        "void_relationship_count": len(voids),
        "opening": opening,
        "host": host,
    }


def width_axis(bbox: dict[str, list[float]]) -> str:
    dx, dy = bbox["dimensions_mm"][:2]
    return "X" if dx >= dy else "Y"


def interval_for_bbox(bbox: dict[str, list[float]], axis: str) -> tuple[float, float]:
    index = 0 if axis == "X" else 1
    return float(bbox["min_mm"][index]), float(bbox["max_mm"][index])


def interval_relationship(first: tuple[float, float], second: tuple[float, float]) -> dict[str, float | str]:
    overlap = max(0.0, min(first[1], second[1]) - max(first[0], second[0]))
    gap = max(0.0, max(first[0], second[0]) - min(first[1], second[1]))
    smaller = min(first[1] - first[0], second[1] - second[0])
    return {
        "kind": "OVERLAP" if overlap > 0.0 else "GAP_OR_TOUCH",
        "overlap_mm": overlap,
        "gap_mm": gap,
        "overlap_ratio_of_smaller": overlap / smaller if smaller > 0 else 0.0,
    }


def point_in_bbox(x: float, y: float, bbox: dict[str, list[float]], tolerance: float = 0.1) -> bool:
    return (
        bbox["min_mm"][0] - tolerance <= x <= bbox["max_mm"][0] + tolerance
        and bbox["min_mm"][1] - tolerance <= y <= bbox["max_mm"][1] + tolerance
    )


def space_candidates(
    center: Sequence[float],
    axis: str,
    spaces: Sequence[dict[str, Any]],
) -> list[str]:
    """Return evidence candidates only; never assert a room relationship."""

    x, y = float(center[0]), float(center[1])
    points = [(x, y)]
    normal_axis = 1 if axis == "X" else 0
    for offset in (150.0, 300.0):
        first = [x, y]
        second = [x, y]
        first[normal_axis] -= offset
        second[normal_axis] += offset
        points.extend((tuple(first), tuple(second)))
    names: set[str] = set()
    for px, py in points:
        for space in spaces:
            if point_in_bbox(px, py, space["bbox"]):
                names.add(space["long_name"])
    return sorted(names)


def format_operation(value: str) -> str:
    mapping = {
        "SINGLE_SWING_LEFT": "单扇左开",
        "SINGLE_SWING_RIGHT": "单扇右开",
        "SLIDING_TO_LEFT": "左向推拉",
        "SINGLE_PANEL": "单扇窗",
        "DOUBLE_PANEL_HORIZONTAL": "双扇横分",
        "DOUBLE_PANEL_VERTICAL": "双扇竖分",
        "TRIPLE_PANEL_LEFT": "三扇左分格",
        "TRIPLE_PANEL_RIGHT": "三扇右分格",
    }
    return mapping.get(value, value or "待确认")


def review_metadata(record: dict[str, Any]) -> dict[str, Any]:
    global_id = record["global_id"]
    if global_id == MASTER_BEDROOM_DOOR_CANDIDATE:
        return {
            "review_group": "A104-R01",
            "review_required": "yes",
            "review_question": "该 900×2000 左开门位于次卧—客卫干区之间；确认它的真实房间归属，不能直接写为主卧门。",
            "confidence": 1.0,
        }
    if global_id == "1TW6$_GfnABRZusYvx0zZG":
        return {
            "review_group": "A104-R01",
            "review_required": "yes",
            "review_question": "该无宿主三角开启门位于主卧入口范围；确认是否为主卧门，并核定洞口与名义尺寸。",
            "confidence": 0.75,
        }
    if global_id in BATHROOM_UNHOSTED_PAIR:
        return {
            "review_group": "A104-R02",
            "review_required": "yes",
            "review_question": "两对象共用放置：一块 1000×38×2400 门板和一条 1200×61×7 顶部构件；确认是否属于同一樘门。",
            "confidence": 0.9,
        }
    return {
        "review_group": "",
        "review_required": "no",
        "review_question": "",
        "confidence": 1.0,
    }


def plan_sort_key(record: dict[str, Any]) -> tuple[float, float, str]:
    center = record["bbox"]["centre_mm"]
    return (-round(float(center[1]), 3), round(float(center[0]), 3), record["global_id"])


def build_inventory(
    model: ifcopenshell.file,
    tolerance_mm: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    settings = geometry_settings()
    spaces: list[dict[str, Any]] = []
    for space in model.by_type("IfcSpace"):
        spaces.append(
            {
                "global_id": space.GlobalId,
                "long_name": str(space.LongName or space.Name or space.GlobalId),
                "bbox": world_bbox_mm(settings, space),
            }
        )

    products: list[dict[str, Any]] = []
    for product in list(model.by_type("IfcDoor")) + list(model.by_type("IfcWindow")):
        bbox = world_bbox_mm(settings, product)
        chain = fill_chain(product)
        opening = chain["opening"]
        host = chain["host"]
        opening_bbox = world_bbox_mm(settings, opening) if opening else None
        bbim = read_bbim_data(product)
        axis = width_axis(bbox)
        nominal_width = float(product.OverallWidth) if product.OverallWidth is not None else None
        nominal_height = float(product.OverallHeight) if product.OverallHeight is not None else None
        opening_width = None
        opening_height = None
        width_delta = None
        height_delta = None
        if opening_bbox and nominal_width is not None and nominal_height is not None:
            axis_index = 0 if axis == "X" else 1
            opening_width = float(opening_bbox["dimensions_mm"][axis_index])
            opening_height = float(opening_bbox["dimensions_mm"][2])
            width_delta = abs(opening_width - nominal_width)
            height_delta = abs(opening_height - nominal_height)
        container = product_container(product)
        operation = str(bbim.get("door_type") or bbim.get("window_type") or "")
        record = {
            "ifc_class": product.is_a(),
            "global_id": product.GlobalId,
            "name": str(product.Name or ""),
            "current_tag": str(product.Tag or ""),
            "current_predefined_type": str(product.PredefinedType or ""),
            "candidate_predefined_type": "DOOR" if product.is_a("IfcDoor") else "WINDOW",
            "container_global_id": container["global_id"],
            "container_name": container["name"],
            "bbox": bbox,
            "width_axis": axis,
            "nominal_width_mm": nominal_width,
            "nominal_height_mm": nominal_height,
            "operation_type": operation,
            "operation_label": format_operation(operation),
            "fill_relationship_count": chain["fill_relationship_count"],
            "void_relationship_count": chain["void_relationship_count"],
            "opening_global_id": opening.GlobalId if opening else "",
            "host_wall_global_id": host.GlobalId if host else "",
            "host_wall_name": str(host.Name or "") if host else "",
            "opening_bbox": opening_bbox,
            "opening_width_mm": opening_width,
            "opening_height_mm": opening_height,
            "opening_width_delta_mm": width_delta,
            "opening_height_delta_mm": height_delta,
            "geometry_relation_pass": bool(
                opening
                and host
                and width_delta is not None
                and height_delta is not None
                and width_delta <= tolerance_mm
                and height_delta <= tolerance_mm
            ),
            "space_candidates": space_candidates(bbox["centre_mm"], axis, spaces),
            "basis": (
                "IfcRelFillsElement + IfcRelVoidsElement + nominal size/opening geometry"
                if opening and host
                else "product world geometry + spatial container only; no host inferred"
            ),
        }
        record.update(review_metadata(record))
        products.append(record)

    doors = sorted((record for record in products if record["ifc_class"] == "IfcDoor"), key=plan_sort_key)
    windows = sorted((record for record in products if record["ifc_class"] == "IfcWindow"), key=plan_sort_key)
    for index, record in enumerate(doors, start=1):
        record["candidate_id"] = f"M{index:02d}"
    for index, record in enumerate(windows, start=1):
        record["candidate_id"] = f"W{index:02d}"

    pair_relations: list[dict[str, Any]] = []
    by_host: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in products:
        if record["host_wall_global_id"] and record["opening_bbox"]:
            by_host[record["host_wall_global_id"]].append(record)
    for host_id, records in sorted(by_host.items()):
        for index, first in enumerate(records):
            for second in records[index + 1 :]:
                if first["width_axis"] != second["width_axis"]:
                    continue
                first_interval = interval_for_bbox(first["opening_bbox"], first["width_axis"])
                second_interval = interval_for_bbox(second["opening_bbox"], second["width_axis"])
                relation = interval_relationship(first_interval, second_interval)
                if relation["overlap_mm"] <= tolerance_mm and relation["gap_mm"] > tolerance_mm:
                    continue
                pair_relations.append(
                    {
                        "host_wall_global_id": host_id,
                        "first_candidate_id": first["candidate_id"],
                        "first_global_id": first["global_id"],
                        "second_candidate_id": second["candidate_id"],
                        "second_global_id": second["global_id"],
                        "width_axis": first["width_axis"],
                        **relation,
                        "review_required": relation["overlap_mm"] > tolerance_mm,
                    }
                )
                if relation["overlap_mm"] > tolerance_mm:
                    for member in (first, second):
                        member["review_group"] = member["review_group"] or f"A104-R{2 + len([p for p in pair_relations if p['review_required']]):02d}"
                        member["review_required"] = "yes"
                        member["review_question"] = (
                            f"与 {second['candidate_id'] if member is first else first['candidate_id']} 的宿主洞口沿宽度方向重叠 "
                            f"{relation['overlap_mm']:.3f} mm；确认是同一窗组的重叠表达还是应保留两樘编号。"
                        )
    return doors, windows, pair_relations


def register_rows(records: Iterable[dict[str, Any]], source_sha: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        bbox = record["bbox"]
        rows.append(
            {
                "candidate_id": record["candidate_id"],
                "candidate_tag": record["candidate_id"],
                "ifc_class": record["ifc_class"],
                "global_id": record["global_id"],
                "name": record["name"],
                "current_tag": record["current_tag"],
                "current_predefined_type": record["current_predefined_type"],
                "candidate_predefined_type": record["candidate_predefined_type"],
                "container_name": record["container_name"],
                "host_relation": "HOSTED" if record["host_wall_global_id"] else "UNHOSTED",
                "opening_global_id": record["opening_global_id"],
                "host_wall_global_id": record["host_wall_global_id"],
                "nominal_width_mm": "" if record["nominal_width_mm"] is None else format_mm(record["nominal_width_mm"]),
                "nominal_height_mm": "" if record["nominal_height_mm"] is None else format_mm(record["nominal_height_mm"]),
                "operation_type": record["operation_type"],
                "width_axis": record["width_axis"],
                "opening_width_mm": "" if record["opening_width_mm"] is None else f"{record['opening_width_mm']:.6f}",
                "opening_height_mm": "" if record["opening_height_mm"] is None else f"{record['opening_height_mm']:.6f}",
                "opening_width_delta_mm": "" if record["opening_width_delta_mm"] is None else f"{record['opening_width_delta_mm']:.6f}",
                "opening_height_delta_mm": "" if record["opening_height_delta_mm"] is None else f"{record['opening_height_delta_mm']:.6f}",
                "geometry_relation_pass": str(record["geometry_relation_pass"]).lower(),
                "center_x_mm": f"{bbox['centre_mm'][0]:.6f}",
                "center_y_mm": f"{bbox['centre_mm'][1]:.6f}",
                "sill_or_threshold_z_mm": f"{bbox['min_mm'][2]:.6f}",
                "space_candidates": "; ".join(record["space_candidates"]),
                "basis": record["basis"],
                "confidence": f"{record['confidence']:.2f}",
                "review_group": record["review_group"],
                "review_required": record["review_required"],
                "review_question": record["review_question"],
                "formal_ifc_write_allowed": "no",
                "source_ifc_sha256": source_sha,
            }
        )
    return rows


def write_register(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError("cannot write an empty A-104 register")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def render_plan_markers(records: Sequence[dict[str, Any]], occupied: list[Box]) -> tuple[str, list[dict[str, Any]]]:
    pieces = ['<g id="a104-product-markers">']
    labels: list[dict[str, Any]] = []
    offsets = ((0, -7), (0, 8), (11, 0), (-11, 0), (11, -7), (-11, -7), (11, 8), (-11, 8))
    for record in records:
        bbox = record["opening_bbox"] or record["bbox"]
        cx, cy = world_to_svg(bbox["centre_mm"][0], bbox["centre_mm"][1])
        axis = record["width_axis"]
        if record["opening_bbox"]:
            lower, upper = interval_for_bbox(bbox, axis)
            if axis == "X":
                x1, _ = world_to_svg(lower, 0)
                x2, _ = world_to_svg(upper, 0)
                y1 = y2 = cy
            else:
                _, y1 = world_to_svg(0, upper)
                _, y2 = world_to_svg(0, lower)
                x1 = x2 = cx
            pieces.append(
                f'<line class="a104-location a104-{record["ifc_class"].lower()}" '
                f'data-ifc-guid="{svg_escape(record["global_id"])}" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}"/>'
            )
        else:
            x, y = world_to_svg(bbox["min_mm"][0], bbox["max_mm"][1])
            w = max(1.2, bbox["dimensions_mm"][0] / 50.0)
            h = max(1.2, bbox["dimensions_mm"][1] / 50.0)
            pieces.append(
                f'<rect class="a104-unhosted" data-ifc-guid="{svg_escape(record["global_id"])}" '
                f'x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}"/>'
            )
        size = (
            f"{format_mm(record['nominal_width_mm'])}×{format_mm(record['nominal_height_mm'])}"
            if record["nominal_width_mm"] is not None and record["nominal_height_mm"] is not None
            else "尺寸待定"
        )
        label = f"{record['candidate_id']} {size}"
        tx, ty, box = place_label(cx, cy, label, occupied, offsets=offsets, font_size=2.15)
        issue = " a104-review" if record["review_required"] == "yes" else ""
        pieces.append(
            f'<line class="a104-leader{issue}" x1="{cx:.3f}" y1="{cy:.3f}" x2="{tx:.3f}" y2="{ty - 0.8:.3f}"/>'
            f'<rect class="a104-label-bg" x="{box.x0 - 0.7:.3f}" y="{box.y0 - 0.4:.3f}" '
            f'width="{box.x1 - box.x0 + 1.4:.3f}" height="{box.y1 - box.y0 + 0.8:.3f}"/>'
            f'<text class="a104-label{issue}" x="{tx:.3f}" y="{ty:.3f}">{svg_escape(label)}</text>'
        )
        labels.append({"candidate_id": record["candidate_id"], "global_id": record["global_id"], "text": label, "box": [box.x0, box.y0, box.x1, box.y1]})
    pieces.append("</g>")
    return "".join(pieces), labels


def render_side_panel(
    source_sha: str,
    doors: Sequence[dict[str, Any]],
    windows: Sequence[dict[str, Any]],
    pair_relations: Sequence[dict[str, Any]],
    collision_count: int,
) -> str:
    lines = ['<g id="a104-side-panel"><rect class="a104-panel" x="402" y="7" width="93" height="386"/>']
    y = 16.0

    def line(value: str, css: str = "a104-panel-text", gap: float = 4.7) -> None:
        nonlocal y
        lines.append(f'<text class="{css}" x="407" y="{y:.1f}">{svg_escape(value)}</text>')
        y += gap

    line("A-104 门窗定位图及门窗表", "a104-panel-title", 7.5)
    line("候选编号｜未写 IFC｜1:50", "a104-panel-note", 7.0)
    line("门表", "a104-panel-heading", 5.5)
    for record in doors:
        size = (
            f"{format_mm(record['nominal_width_mm'])}×{format_mm(record['nominal_height_mm'])}"
            if record["nominal_width_mm"] is not None else "尺寸待定"
        )
        marker = "※" if record["review_required"] == "yes" else ""
        line(f"{record['candidate_id']}{marker}  {size}  {record['operation_label']}", gap=4.4)
    y += 2.0
    line("窗表", "a104-panel-heading", 5.5)
    for record in windows:
        marker = "※" if record["review_required"] == "yes" else ""
        line(
            f"{record['candidate_id']}{marker}  {format_mm(record['nominal_width_mm'])}×{format_mm(record['nominal_height_mm'])}  {record['operation_label']}",
            gap=4.4,
        )
    y += 2.0
    line("机械 QA", "a104-panel-heading", 5.5)
    hosted = [record for record in [*doors, *windows] if record["host_wall_global_id"]]
    line(f"门 8／窗 11；有宿主 {len(hosted)}", gap=4.3)
    line(f"洞口名义宽高 0.1mm：{sum(r['geometry_relation_pass'] for r in hosted)}/{len(hosted)}", gap=4.3)
    line(f"候选编号重复：0；标注碰撞：{collision_count}", gap=4.3)
    line(f"重叠洞口组：{sum(bool(r['review_required']) for r in pair_relations)}", gap=4.3)
    line(f"IFC SHA {source_sha[:12]}…", "a104-panel-note", 5.0)
    line("※ 仅表示必须在 Blender 确认", "a104-panel-note", 4.5)
    lines.append("</g>")
    return "".join(lines)


def inject_candidate_svg(source: str, generated: str) -> str:
    if 'data-scale="1:50"' not in source:
        raise RuntimeError("Wall Plan SVG is not at expected 1:50 scale")
    source = re.sub(r'width="400(?:\.0+)?mm"', 'width="500mm"', source, count=1)
    source = re.sub(r'viewBox="0 0 400(?:\.0+)? 400(?:\.0+)?"', 'viewBox="0 0 500 400"', source, count=1)
    style = """
<style id="a104-candidate-style">
  @page { size:500mm 400mm; margin:0; }
  html,body { margin:0; width:500mm; height:400mm; overflow:hidden; }
  .a104-location { fill:none; stroke-width:1.2; stroke-linecap:round; }
  .a104-ifcdoor { stroke:#008b8b; }
  .a104-ifcwindow { stroke:#2563a6; }
  .a104-unhosted { fill:#ff6b6b55; stroke:#c92a2a; stroke-width:0.8; }
  .a104-leader { stroke:#526777; stroke-width:0.3; }
  .a104-leader.a104-review { stroke:#c92a2a; stroke-width:0.55; }
  .a104-label-bg { fill:#fff; stroke:#d4dde4; stroke-width:0.2; }
  .a104-label { font-family:Arial,'Noto Sans CJK SC',sans-serif; fill:#102f43; text-anchor:middle; font-size:2.15px; font-weight:700; }
  .a104-label.a104-review { fill:#c92a2a; }
  .a104-panel { fill:#fbfcfd; stroke:#102f43; stroke-width:0.5; }
  .a104-panel-title,.a104-panel-heading,.a104-panel-text,.a104-panel-note { font-family:Arial,'Noto Sans CJK SC',sans-serif; fill:#102f43; }
  .a104-panel-title { font-size:4px; font-weight:700; }
  .a104-panel-heading { font-size:3px; font-weight:700; }
  .a104-panel-text { font-size:2.25px; }
  .a104-panel-note { font-size:2.15px; fill:#526777; }
</style>
"""
    if "</svg>" not in source:
        raise RuntimeError("invalid SVG: closing tag missing")
    return source.replace("</svg>", f'{style}<g id="a104-generated">{generated}</g></svg>', 1)


def main() -> None:
    args = parse_args()
    model = ifcopenshell.open(args.input)
    source_svg = args.source_svg.read_text(encoding="utf-8")
    source_sha = sha256(args.input)
    doors, windows, pair_relations = build_inventory(model, args.tolerance_mm)
    if len(doors) != 8 or len(windows) != 11:
        raise RuntimeError(f"expected 8 doors and 11 windows, found {len(doors)} and {len(windows)}")
    if {record["global_id"] for record in doors if not record["host_wall_global_id"]} != UNHOSTED_DOORS:
        raise RuntimeError("A-104 unhosted door set drifted from the C003 handoff")
    if len({record["candidate_id"] for record in [*doors, *windows]}) != 19:
        raise RuntimeError("duplicate A-104 candidate identifier")

    rows = register_rows([*doors, *windows], source_sha)
    write_register(args.review_register, rows)
    occupied: list[Box] = []
    marker_markup, labels = render_plan_markers([*doors, *windows], occupied)
    collisions = [
        [first.label, second.label]
        for index, first in enumerate(occupied)
        for second in occupied[index + 1 :]
        if first.intersects(second)
    ]
    if collisions:
        raise RuntimeError(f"generated label collisions remain: {collisions}")
    panel = render_side_panel(source_sha, doors, windows, pair_relations, len(collisions))
    output_svg = inject_candidate_svg(source_svg, marker_markup + panel)
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(output_svg, encoding="utf-8")

    hosted = [record for record in [*doors, *windows] if record["host_wall_global_id"]]
    review_queue = [
        {
            key: record[key]
            for key in (
                "candidate_id",
                "ifc_class",
                "global_id",
                "review_group",
                "review_question",
                "space_candidates",
            )
        }
        for record in [*doors, *windows]
        if record["review_required"] == "yes"
    ]
    gates = {
        "ifc_schema_is_ifc4": model.schema == "IFC4",
        "door_count": len(doors),
        "window_count": len(windows),
        "hosted_product_count": len(hosted),
        "unhosted_door_count": sum(not record["host_wall_global_id"] for record in doors),
        "hosted_nominal_opening_checks_passed": sum(record["geometry_relation_pass"] for record in hosted),
        "candidate_identifier_count": len({record["candidate_id"] for record in [*doors, *windows]}),
        "candidate_identifier_duplicate_count": 19 - len({record["candidate_id"] for record in [*doors, *windows]}),
        "generated_label_count": len(labels),
        "generated_label_collision_count": len(collisions),
        "review_queue_count": len(review_queue),
        "automatic_ifc_write_allowed": False,
    }
    gates["mechanical_pass"] = (
        gates["ifc_schema_is_ifc4"]
        and gates["door_count"] == 8
        and gates["window_count"] == 11
        and gates["hosted_product_count"] == 16
        and gates["unhosted_door_count"] == 3
        and gates["hosted_nominal_opening_checks_passed"] == 16
        and gates["candidate_identifier_duplicate_count"] == 0
        and gates["generated_label_collision_count"] == 0
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a104-door-window-candidate",
        "source": {
            "ifc": str(args.input.resolve()),
            "ifc_sha256": source_sha,
            "schema": model.schema,
            "wall_plan_svg": str(args.source_svg.resolve()),
            "wall_plan_svg_sha256": sha256(args.source_svg),
        },
        "tolerance_mm": args.tolerance_mm,
        "numbering_rule": "M=door and W=window; plan scan north-to-south, then west-to-east; candidate only until review; D prefix is reserved by A-102 demolition walls",
        "doors": doors,
        "windows": windows,
        "shared_host_opening_relations": pair_relations,
        "review_queue": review_queue,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_svg": str(args.output_svg), "review_register": str(args.review_register), "report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["mechanical_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
