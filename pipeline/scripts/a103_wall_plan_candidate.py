#!/usr/bin/env python3
"""Generate the read-only A-103 wall positioning candidate.

The formal IFC is never modified.  The script combines the current Bonsai Wall
Plan SVG with mechanically derived grid chains, confirmed-new-wall dimensions,
hosted-door opening marks, a wall status register, and machine-readable QA.

Wall construction status comes from the confirmed IFC semantics:

* ``Status=NEW`` is shown green;
* ``Status=EXISTING`` with ``LoadBearing=false`` is shown orange;
* ``Status=EXISTING`` with ``LoadBearing=true`` is shown blue-gray.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.placement


SCALE_DENOMINATOR = 50.0
SVG_WORLD_OFFSET_MM = 10000.0
SVG_PLAN_SIZE_MM = 400.0
SVG_SHEET_WIDTH_MM = 500.0
CONFIRMED_NEW_GUIDS = {
    "04rs0EDjn2EvxytEQSxWRB",
    "3jha5L04zBjh0$pMl_tHLy",
    "2Ca2tGerPBj8kvUTjlUtDl",
    "3pAfMJYxPBlwR30CstZ4kK",
}


@dataclass(frozen=True)
class Box:
    x0: float
    y0: float
    x1: float
    y1: float
    label: str

    def intersects(self, other: "Box", padding: float = 0.4) -> bool:
        return not (
            self.x1 + padding <= other.x0
            or other.x1 + padding <= self.x0
            or self.y1 + padding <= other.y0
            or other.y1 + padding <= self.y0
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def world_to_svg(x_mm: float, y_mm: float) -> tuple[float, float]:
    return (
        (x_mm + SVG_WORLD_OFFSET_MM) / SCALE_DENOMINATOR,
        (SVG_WORLD_OFFSET_MM - y_mm) / SCALE_DENOMINATOR,
    )


def svg_escape(value: Any) -> str:
    return html.escape(str(value), quote=True)


def format_mm(value: float, decimals: int = 3) -> str:
    if abs(value - round(value)) <= 1e-6:
        return str(int(round(value)))
    return f"{value:.{decimals}f}".rstrip("0").rstrip(".")


def geometry_settings() -> ifcopenshell.geom.settings:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    return settings


def world_bbox_mm(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> dict[str, list[float]]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    raw = list(shape.geometry.verts)
    points = [
        tuple(float(raw[index + offset]) * 1000.0 for offset in range(3))
        for index in range(0, len(raw), 3)
    ]
    if not points:
        raise RuntimeError(f"{product.GlobalId} has no shape vertices")
    minimum = [min(point[axis] for point in points) for axis in range(3)]
    maximum = [max(point[axis] for point in points) for axis in range(3)]
    return {
        "min_mm": minimum,
        "max_mm": maximum,
        "dimensions_mm": [maximum[axis] - minimum[axis] for axis in range(3)],
        "centre_mm": [(minimum[axis] + maximum[axis]) / 2.0 for axis in range(3)],
    }


def material_names(product: ifcopenshell.entity_instance) -> list[str]:
    material = ifcopenshell.util.element.get_material(product)
    if material is None:
        return []
    names: set[str] = set()

    def visit(entity: Any) -> None:
        if entity is None or not hasattr(entity, "is_a"):
            return
        if entity.is_a("IfcMaterial") and entity.Name:
            names.add(str(entity.Name))
        for attribute in (
            "ForLayerSet",
            "MaterialLayers",
            "Material",
            "ForProfileSet",
            "MaterialProfiles",
        ):
            child = getattr(entity, attribute, None)
            if isinstance(child, tuple):
                for item in child:
                    visit(item)
            else:
                visit(child)

    visit(material)
    if not names and getattr(material, "Name", None):
        names.add(str(material.Name))
    return sorted(names)


def wall_status_candidate(
    current_status: str | None,
    current_load_bearing: bool | None,
    materials: Sequence[str],
) -> tuple[str, str, float, str]:
    normalized = (current_status or "").upper()
    if normalized == "NEW":
        return (
            "CONFIRMED_NEW",
            "IFC Pset_WallCommon.Status=NEW",
            1.0,
            "no",
        )
    if normalized == "EXISTING" and current_load_bearing is False and "Aircrete" in materials:
        return (
            "CONFIRMED_EXISTING_NON_LOAD_BEARING",
            "IFC Pset_WallCommon.Status=EXISTING and LoadBearing=false; user-confirmed orange group",
            1.0,
            "no",
        )
    if normalized == "EXISTING" and current_load_bearing is True and "Concrete" in materials:
        return (
            "CONFIRMED_EXISTING_LOAD_BEARING",
            "IFC Pset_WallCommon.Status=EXISTING and LoadBearing=true; user-confirmed blue-gray group",
            1.0,
            "no",
        )
    return (
        "INVALID_SEMANTICS",
        "Current IFC status, LoadBearing and material do not match the confirmed A-103 classification",
        0.0,
        "yes",
    )


def current_wall_status(wall: ifcopenshell.entity_instance) -> str | None:
    psets = ifcopenshell.util.element.get_psets(wall)
    value = psets.get("Pset_WallCommon", {}).get("Status")
    return str(value) if value else None


def current_wall_load_bearing(wall: ifcopenshell.entity_instance) -> bool | None:
    psets = ifcopenshell.util.element.get_psets(wall)
    value = psets.get("Pset_WallCommon", {}).get("LoadBearing")
    return bool(value) if value is not None else None


def grid_axis_coordinate(axis: ifcopenshell.entity_instance) -> float:
    curve = axis.AxisCurve
    if not curve or not curve.is_a("IfcPolyline") or len(curve.Points) < 2:
        raise RuntimeError(f"unsupported GridAxis curve for {axis.AxisTag}")
    first = curve.Points[0].Coordinates
    second = curve.Points[-1].Coordinates
    if abs(float(first[0]) - float(second[0])) <= abs(float(first[1]) - float(second[1])):
        return (float(first[0]) + float(second[0])) / 2.0
    return (float(first[1]) + float(second[1])) / 2.0


def grid_chain(axes: Iterable[ifcopenshell.entity_instance]) -> list[dict[str, Any]]:
    records = sorted(
        ({"tag": str(axis.AxisTag), "coordinate_mm": grid_axis_coordinate(axis)} for axis in axes),
        key=lambda record: record["coordinate_mm"],
    )
    for index, record in enumerate(records):
        if index == 0:
            record["segment_from_previous_mm"] = None
        else:
            record["segment_from_previous_mm"] = (
                record["coordinate_mm"] - records[index - 1]["coordinate_mm"]
            )
    return records


def chain_closure(chain: Sequence[dict[str, Any]]) -> dict[str, Any]:
    segments = [
        float(record["segment_from_previous_mm"])
        for record in chain
        if record["segment_from_previous_mm"] is not None
    ]
    overall = float(chain[-1]["coordinate_mm"] - chain[0]["coordinate_mm"])
    residual = sum(segments) - overall
    return {
        "segment_sum_mm": sum(segments),
        "overall_mm": overall,
        "residual_mm": residual,
        "pass": abs(residual) <= 1e-6,
    }


def visible_wall_guids(svg: str, wall_guids: Iterable[str]) -> set[str]:
    return {global_id for global_id in wall_guids if global_id in svg}


def add_wall_status_classes(svg: str, status_by_guid: dict[str, str]) -> str:
    pattern = re.compile(r'<g\b[^>]*\bclass="[^"]*\bIfcWall\b[^"]*"[^>]*>')

    def replace_tag(match: re.Match[str]) -> str:
        tag = match.group(0)
        statuses = {
            status
            for global_id, status in status_by_guid.items()
            if global_id in tag
        }
        if not statuses:
            return tag
        if len(statuses) != 1:
            raise RuntimeError(f"one SVG wall group mixes status candidates: {sorted(statuses)}")
        css_class = f"a103-{next(iter(statuses)).lower().replace('_', '-')}"
        return re.sub(
            r'class="([^"]*)"',
            lambda class_match: f'class="{class_match.group(1)} {css_class}"',
            tag,
            count=1,
        )

    return pattern.sub(replace_tag, svg)


def make_label_box(x: float, y: float, text: str, font_size: float = 2.6) -> Box:
    width = max(font_size * 0.65 * len(text), font_size)
    return Box(x - width / 2.0, y - font_size, x + width / 2.0, y + 0.4, text)


def place_label(
    x: float,
    y: float,
    text: str,
    occupied: list[Box],
    offsets: Sequence[tuple[float, float]] | None = None,
    font_size: float = 2.6,
) -> tuple[float, float, Box]:
    offsets = offsets or ((0, 0), (0, -5), (0, 5), (7, 0), (-7, 0), (7, -5), (-7, -5))
    for dx, dy in offsets:
        candidate = make_label_box(x + dx, y + dy, text, font_size)
        if not any(candidate.intersects(existing) for existing in occupied):
            occupied.append(candidate)
            return x + dx, y + dy, candidate
    raise RuntimeError(f"unable to place generated label without collision: {text}")


def dimension_line(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    label: str,
    text_x: float,
    text_y: float,
    rotate: bool = False,
) -> str:
    transform = f' transform="rotate(-90 {text_x:.3f} {text_y:.3f})"' if rotate else ""
    return (
        f'<line class="a103-dimension PredefinedType-DIMENSION" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}"/>'
        f'<line class="a103-tick" x1="{x1 - 1.2:.3f}" y1="{y1 - 1.2:.3f}" x2="{x1 + 1.2:.3f}" y2="{y1 + 1.2:.3f}"/>'
        f'<line class="a103-tick" x1="{x2 - 1.2:.3f}" y1="{y2 - 1.2:.3f}" x2="{x2 + 1.2:.3f}" y2="{y2 + 1.2:.3f}"/>'
        f'<text class="a103-dimension-text" x="{text_x:.3f}" y="{text_y:.3f}"{transform}>{svg_escape(label)}</text>'
    )


def render_grid_dimensions(
    v_chain: Sequence[dict[str, Any]],
    u_chain: Sequence[dict[str, Any]],
    occupied: list[Box],
) -> tuple[str, list[dict[str, Any]]]:
    pieces = ['<g id="a103-grid-dimensions">']
    records: list[dict[str, Any]] = []

    x_positions = [world_to_svg(record["coordinate_mm"], 0.0)[0] for record in v_chain]
    for index in range(1, len(v_chain)):
        x1, x2 = x_positions[index - 1], x_positions[index]
        value = v_chain[index]["segment_from_previous_mm"]
        label = format_mm(value)
        tx, ty, _ = place_label((x1 + x2) / 2.0, 343.0, label, occupied, font_size=2.4)
        pieces.append(dimension_line(x1, 339.0, x2, 339.0, label, tx, ty))
        records.append({"axis_from": v_chain[index - 1]["tag"], "axis_to": v_chain[index]["tag"], "value_mm": value})
    overall = v_chain[-1]["coordinate_mm"] - v_chain[0]["coordinate_mm"]
    tx, ty, _ = place_label((x_positions[0] + x_positions[-1]) / 2.0, 351.0, format_mm(overall), occupied, font_size=2.4)
    pieces.append(dimension_line(x_positions[0], 347.0, x_positions[-1], 347.0, format_mm(overall), tx, ty))

    y_positions = [world_to_svg(0.0, record["coordinate_mm"])[1] for record in u_chain]
    for index in range(1, len(u_chain)):
        y1, y2 = y_positions[index - 1], y_positions[index]
        value = u_chain[index]["segment_from_previous_mm"]
        label = format_mm(value)
        tx, ty, _ = place_label(34.0, (y1 + y2) / 2.0, label, occupied, font_size=2.4)
        pieces.append(dimension_line(39.0, y1, 39.0, y2, label, tx, ty, rotate=True))
        records.append({"axis_from": u_chain[index - 1]["tag"], "axis_to": u_chain[index]["tag"], "value_mm": value})
    overall_y = u_chain[-1]["coordinate_mm"] - u_chain[0]["coordinate_mm"]
    tx, ty, _ = place_label(25.0, (y_positions[0] + y_positions[-1]) / 2.0, format_mm(overall_y), occupied, font_size=2.4)
    pieces.append(dimension_line(30.0, y_positions[0], 30.0, y_positions[-1], format_mm(overall_y), tx, ty, rotate=True))
    pieces.append("</g>")
    return "".join(pieces), records


def render_new_wall_annotations(
    records: Sequence[dict[str, Any]], occupied: list[Box]
) -> str:
    pieces = ['<g id="a103-confirmed-new-wall-dimensions">']
    for index, record in enumerate(records, start=1):
        bbox = record["bbox"]
        cx, cy = world_to_svg(*bbox["centre_mm"][:2])
        label = f"N{index:02d}"
        tx, ty, _ = place_label(cx, cy, label, occupied, font_size=2.8)
        pieces.append(f'<circle class="a103-new-tag" cx="{tx:.3f}" cy="{ty - 0.9:.3f}" r="3.2"/>')
        pieces.append(f'<text class="a103-tag-text" x="{tx:.3f}" y="{ty:.3f}">{label}</text>')

        dx, dy, _ = bbox["dimensions_mm"]
        long_axis = "x" if dx >= dy else "y"
        long_value = max(dx, dy)
        if long_axis == "x":
            x1 = world_to_svg(bbox["min_mm"][0], 0)[0]
            x2 = world_to_svg(bbox["max_mm"][0], 0)[0]
            y = world_to_svg(0, bbox["min_mm"][1])[1] + 5.0
            dim_label = format_mm(long_value)
            tdx, tdy, _ = place_label((x1 + x2) / 2.0, y - 1.2, dim_label, occupied, font_size=2.2)
            pieces.append(dimension_line(x1, y, x2, y, dim_label, tdx, tdy))
        else:
            y1 = world_to_svg(0, bbox["max_mm"][1])[1]
            y2 = world_to_svg(0, bbox["min_mm"][1])[1]
            x = world_to_svg(bbox["max_mm"][0], 0)[0] + 5.0
            dim_label = format_mm(long_value)
            tdx, tdy, _ = place_label(x + 1.3, (y1 + y2) / 2.0, dim_label, occupied, font_size=2.2)
            pieces.append(dimension_line(x, y1, x, y2, dim_label, tdx, tdy, rotate=True))
    pieces.append("</g>")
    return "".join(pieces)


def render_door_annotations(
    records: Sequence[dict[str, Any]], occupied: list[Box]
) -> str:
    pieces = ['<g id="a103-hosted-door-openings">']
    for index, record in enumerate(records, start=1):
        bbox = record["opening_bbox"]
        cx, cy = world_to_svg(*bbox["centre_mm"][:2])
        dx, dy, _ = bbox["dimensions_mm"]
        if dx >= dy:
            x1 = world_to_svg(bbox["min_mm"][0], 0)[0]
            x2 = world_to_svg(bbox["max_mm"][0], 0)[0]
            y1 = y2 = cy
        else:
            y1 = world_to_svg(0, bbox["max_mm"][1])[1]
            y2 = world_to_svg(0, bbox["min_mm"][1])[1]
            x1 = x2 = cx
        pieces.append(f'<line class="a103-door-opening PredefinedType-DIMENSION" x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}"/>')
        label = f"D{index:02d} {format_mm(record['overall_width_mm'])}×{format_mm(record['overall_height_mm'])}"
        tx, ty, _ = place_label(cx, cy - 5.0, label, occupied, font_size=2.2)
        pieces.append(f'<text class="a103-door-text" x="{tx:.3f}" y="{ty:.3f}">{svg_escape(label)}</text>')
    pieces.append("</g>")
    return "".join(pieces)


def render_side_panel(
    source_sha: str,
    wall_counts: Counter[str],
    confirmed_new: Sequence[dict[str, Any]],
    hosted_doors: Sequence[dict[str, Any]],
    absent_wall_count: int,
    qa_pass: bool,
) -> str:
    lines: list[str] = []
    y = 17.0

    def text_line(value: str, css_class: str = "a103-panel-text", gap: float = 6.0) -> None:
        nonlocal y
        lines.append(f'<text class="{css_class}" x="408" y="{y:.1f}">{svg_escape(value)}</text>')
        y += gap

    lines.append('<g id="a103-side-panel">')
    lines.append('<rect class="a103-panel" x="402" y="7" width="93" height="386"/>')
    text_line("A-103 新建墙体及施工定位图", "a103-panel-title", 8.0)
    text_line("候选版｜比例 1:50｜不写 IFC", "a103-panel-subtitle", 10.0)
    text_line("墙体状态审核图例", "a103-panel-heading", 6.5)
    legend = [
        ("#19a974", f"绿色 新建墙 {wall_counts['CONFIRMED_NEW']}"),
        ("#e68619", f"橙色 现状非承重墙 {wall_counts['CONFIRMED_EXISTING_NON_LOAD_BEARING']}"),
        ("#6d7f96", f"蓝灰 现状承重墙 {wall_counts['CONFIRMED_EXISTING_LOAD_BEARING']}"),
    ]
    for color, label in legend:
        lines.append(f'<rect x="408" y="{y - 3.6:.1f}" width="4" height="4" fill="{color}"/>')
        lines.append(f'<text class="a103-panel-text" x="415" y="{y:.1f}">{svg_escape(label)}</text>')
        y += 5.5
    y += 4.0
    text_line("已确认新建墙", "a103-panel-heading", 6.5)
    for index, record in enumerate(confirmed_new, start=1):
        dims = record["bbox"]["dimensions_mm"]
        dimension_text = "×".join(format_mm(value) for value in dims)
        text_line(f"N{index:02d} {record['type_name'] or '无类型'} {dimension_text}", gap=5.0)
    y += 4.0
    text_line("有宿主门洞定位", "a103-panel-heading", 6.5)
    for index, record in enumerate(hosted_doors, start=1):
        text_line(f"D{index:02d} {format_mm(record['overall_width_mm'])}×{format_mm(record['overall_height_mm'])}", gap=5.0)
    text_line("另 3 樘无宿主门转 A-104", "a103-panel-note", 7.0)
    text_line("机械 QA", "a103-panel-heading", 6.5)
    text_line(f"Grid：20 轴，含 08；闭合 {'PASS' if qa_pass else 'FAIL'}", gap=5.0)
    text_line(f"Wall：88；本切面显示 {88 - absent_wall_count}", gap=5.0)
    text_line(f"本切面外/不相交墙：{absent_wall_count}", gap=5.0)
    text_line("生成标注碰撞：0", gap=5.0)
    text_line(f"IFC SHA {source_sha[:12]}…", "a103-panel-note", 5.0)
    lines.append("</g>")
    return "".join(lines)


def inject_candidate_svg(svg: str, status_by_guid: dict[str, str], generated: str) -> str:
    if 'data-scale="1:50"' not in svg:
        raise RuntimeError("Wall Plan SVG is not at expected 1:50 scale")
    svg = re.sub(r'width="400(?:\.0+)?mm"', 'width="500mm"', svg, count=1)
    svg = re.sub(r'viewBox="0 0 400(?:\.0+)? 400(?:\.0+)?"', 'viewBox="0 0 500 400"', svg, count=1)
    svg = add_wall_status_classes(svg, status_by_guid)
    style = """
<style id="a103-candidate-style">
  @page { size: 500mm 400mm; margin: 0; }
  html, body { margin:0; width:500mm; height:400mm; overflow:hidden; }
  g.a103-confirmed-new path { fill:#19a974 !important; stroke:#0b6b49 !important; }
  g.a103-confirmed-existing-non-load-bearing path { fill:#e68619 !important; stroke:#9c5008 !important; }
  g.a103-confirmed-existing-load-bearing path { fill:#6d7f96 !important; stroke:#34465d !important; }
  g.a103-invalid-semantics path { fill:#c43a84 !important; stroke:#7d1750 !important; }
  g.a103-excluded-demolish { display:none !important; }
  .a103-dimension,.a103-tick { stroke:#143a52; stroke-width:0.35; fill:none; }
  .a103-dimension-text,.a103-tag-text,.a103-door-text { font-family:Arial,'Noto Sans CJK SC',sans-serif; fill:#102f43; text-anchor:middle; font-size:2.4px; }
  .a103-new-tag { fill:#ffffff; stroke:#087f5b; stroke-width:0.6; }
  .a103-tag-text { font-size:2.5px; font-weight:700; }
  .a103-door-opening { stroke:#00a6a6 !important; stroke-width:1.1 !important; stroke-linecap:round; marker-start:none !important; marker-end:none !important; }
  .a103-door-text { fill:#006d77; font-size:2.2px; font-weight:700; paint-order:stroke; stroke:#fff; stroke-width:1px; }
  .a103-panel { fill:#fbfcfd; stroke:#102f43; stroke-width:0.5; }
  .a103-panel-title,.a103-panel-heading,.a103-panel-subtitle,.a103-panel-text,.a103-panel-note { font-family:Arial,'Noto Sans CJK SC',sans-serif; fill:#102f43; }
  .a103-panel-title { font-size:4px; font-weight:700; }
  .a103-panel-heading { font-size:3px; font-weight:700; }
  .a103-panel-subtitle,.a103-panel-text { font-size:2.5px; }
  .a103-panel-note { font-size:2.2px; fill:#526777; }
</style>
"""
    insertion = f"{style}<g id=\"a103-generated\">{generated}</g>"
    if "</svg>" not in svg:
        raise RuntimeError("invalid SVG: closing tag missing")
    candidate = svg.replace("</svg>", f"{insertion}</svg>", 1)
    return re.sub(r"[ \t]+\r?\n", "\n", candidate)


def build_inventory(
    model: ifcopenshell.file,
    source_svg: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    settings = geometry_settings()
    walls: list[dict[str, Any]] = []
    wall_guids = [wall.GlobalId for wall in model.by_type("IfcWall")]
    visible = visible_wall_guids(source_svg, wall_guids)
    demolition_guids: list[str] = []
    for wall in sorted(model.by_type("IfcWall"), key=lambda entity: entity.GlobalId):
        current_status = current_wall_status(wall)
        if (current_status or "").upper() == "DEMOLISH":
            demolition_guids.append(wall.GlobalId)
            continue
        wall_type = ifcopenshell.util.element.get_type(wall)
        materials = material_names(wall)
        current_load_bearing = current_wall_load_bearing(wall)
        candidate, basis, confidence, review_required = wall_status_candidate(
            current_status, current_load_bearing, materials
        )
        record = {
            "global_id": wall.GlobalId,
            "name": wall.Name or "",
            "object_type": wall.ObjectType or "",
            "occurrence_predefined_type": wall.PredefinedType or "",
            "type_name": wall_type.Name if wall_type else "",
            "type_predefined_type": wall_type.PredefinedType if wall_type else "",
            "materials": materials,
            "current_status": current_status or "",
            "current_load_bearing": current_load_bearing,
            "candidate_status": candidate,
            "basis": basis,
            "confidence": confidence,
            "review_required": review_required,
            "visible_in_wall_plan": wall.GlobalId in visible,
            "bbox": world_bbox_mm(settings, wall),
        }
        walls.append(record)

    confirmed_new = [record for record in walls if record["candidate_status"] == "CONFIRMED_NEW"]
    if {record["global_id"] for record in confirmed_new} != CONFIRMED_NEW_GUIDS:
        raise RuntimeError("confirmed NEW wall set drifted from the reviewed A-103 baseline")

    hosted_doors: list[dict[str, Any]] = []
    unhosted_doors: list[dict[str, Any]] = []
    for door in sorted(model.by_type("IfcDoor"), key=lambda entity: entity.GlobalId):
        fills = list(door.FillsVoids or ())
        opening = fills[0].RelatingOpeningElement if len(fills) == 1 else None
        voids = list(opening.VoidsElements or ()) if opening else []
        host = voids[0].RelatingBuildingElement if len(voids) == 1 else None
        record = {
            "door_global_id": door.GlobalId,
            "door_name": door.Name or "",
            "overall_width_mm": float(door.OverallWidth) if door.OverallWidth is not None else None,
            "overall_height_mm": float(door.OverallHeight) if door.OverallHeight is not None else None,
            "opening_global_id": opening.GlobalId if opening else "",
            "host_wall_global_id": host.GlobalId if host else "",
            "fill_relationship_count": len(fills),
            "void_relationship_count": len(voids),
        }
        if (
            opening
            and host
            and door.OverallWidth is not None
            and door.OverallHeight is not None
        ):
            record["opening_bbox"] = world_bbox_mm(settings, opening)
            hosted_doors.append(record)
        else:
            unhosted_doors.append(record)
    return walls, hosted_doors, unhosted_doors, demolition_guids


def write_wall_register(path: Path, walls: Sequence[dict[str, Any]], source_sha: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "global_id",
        "name",
        "object_type",
        "type_name",
        "type_predefined_type",
        "occurrence_predefined_type",
        "materials",
        "current_status",
        "current_load_bearing",
        "candidate_status",
        "basis",
        "confidence",
        "review_required",
        "visible_in_wall_plan",
        "min_x_mm",
        "min_y_mm",
        "min_z_mm",
        "max_x_mm",
        "max_y_mm",
        "max_z_mm",
        "size_x_mm",
        "size_y_mm",
        "size_z_mm",
        "source_ifc_sha256",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for record in walls:
            bbox = record["bbox"]
            writer.writerow(
                {
                    **{key: record.get(key, "") for key in fieldnames},
                    "materials": "; ".join(record["materials"]),
                    "min_x_mm": f"{bbox['min_mm'][0]:.6f}",
                    "min_y_mm": f"{bbox['min_mm'][1]:.6f}",
                    "min_z_mm": f"{bbox['min_mm'][2]:.6f}",
                    "max_x_mm": f"{bbox['max_mm'][0]:.6f}",
                    "max_y_mm": f"{bbox['max_mm'][1]:.6f}",
                    "max_z_mm": f"{bbox['max_mm'][2]:.6f}",
                    "size_x_mm": f"{bbox['dimensions_mm'][0]:.6f}",
                    "size_y_mm": f"{bbox['dimensions_mm'][1]:.6f}",
                    "size_z_mm": f"{bbox['dimensions_mm'][2]:.6f}",
                    "source_ifc_sha256": source_sha,
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--source-svg", type=Path, required=True)
    parser.add_argument("--output-svg", type=Path, required=True)
    parser.add_argument("--wall-register", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument(
        "--expected-ifc-sha256",
        help="Optional caller-frozen formal IFC hash",
    )
    args = parser.parse_args()

    source_sha = sha256(args.input)
    if args.expected_ifc_sha256 and source_sha != args.expected_ifc_sha256:
        raise RuntimeError(
            "formal IFC hash differs from caller-frozen hash: "
            f"expected {args.expected_ifc_sha256}, found {source_sha}"
        )
    model = ifcopenshell.open(args.input)
    source_svg = args.source_svg.read_text(encoding="utf-8")
    source_svg_sha = sha256(args.source_svg)
    walls, hosted_doors, unhosted_doors, demolition_guids = build_inventory(model, source_svg)
    status_counts = Counter(record["candidate_status"] for record in walls)
    if status_counts != Counter(
        {
            "CONFIRMED_NEW": 4,
            "CONFIRMED_EXISTING_NON_LOAD_BEARING": 64,
            "CONFIRMED_EXISTING_LOAD_BEARING": 20,
        }
    ):
        raise RuntimeError(f"A-103 wall grouping drifted: {dict(status_counts)}")

    grid = model.by_type("IfcGrid")
    if len(grid) != 1:
        raise RuntimeError(f"expected one IfcGrid, found {len(grid)}")
    u_chain = grid_chain(grid[0].UAxes)
    v_chain = grid_chain(grid[0].VAxes)
    u_closure = chain_closure(u_chain)
    v_closure = chain_closure(v_chain)
    if len(u_chain) + len(v_chain) != 20 or not u_closure["pass"] or not v_closure["pass"]:
        raise RuntimeError("Grid inventory or dimension-chain closure failed")
    if "08" not in {record["tag"] for record in v_chain}:
        raise RuntimeError("Grid 08 is missing from A-103 dimension chain")

    confirmed_new = [record for record in walls if record["candidate_status"] == "CONFIRMED_NEW"]
    occupied: list[Box] = []
    grid_markup, grid_dimensions = render_grid_dimensions(v_chain, u_chain, occupied)
    wall_markup = render_new_wall_annotations(confirmed_new, occupied)
    door_markup = render_door_annotations(hosted_doors, occupied)
    collisions = [
        [first.label, second.label]
        for index, first in enumerate(occupied)
        for second in occupied[index + 1 :]
        if first.intersects(second)
    ]
    if collisions:
        raise RuntimeError(f"generated annotation collisions remain: {collisions}")

    visible_count = sum(bool(record["visible_in_wall_plan"]) for record in walls)
    qa_pass = u_closure["pass"] and v_closure["pass"] and not collisions
    panel_markup = render_side_panel(
        source_sha,
        status_counts,
        confirmed_new,
        hosted_doors,
        len(walls) - visible_count,
        qa_pass,
    )
    generated = grid_markup + wall_markup + door_markup + panel_markup
    status_by_guid = {record["global_id"]: record["candidate_status"] for record in walls}
    status_by_guid.update({global_id: "EXCLUDED_DEMOLISH" for global_id in demolition_guids})
    candidate_svg = inject_candidate_svg(source_svg, status_by_guid, generated)
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(candidate_svg, encoding="utf-8")
    write_wall_register(args.wall_register, walls, source_sha)

    new_wall_dimension_checks = []
    for record in confirmed_new:
        dimensions = record["bbox"]["dimensions_mm"]
        nearest = [round(value) for value in dimensions]
        deltas = [abs(value - rounded) for value, rounded in zip(dimensions, nearest)]
        new_wall_dimension_checks.append(
            {
                "global_id": record["global_id"],
                "dimensions_mm": dimensions,
                "nearest_integer_mm": nearest,
                "maximum_integer_delta_mm": max(deltas),
                "within_0_1_mm": max(deltas) <= 0.1,
            }
        )

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "read_only": True,
        "source": {
            "ifc": str(args.input),
            "ifc_sha256": source_sha,
            "svg": str(args.source_svg),
            "svg_sha256": source_svg_sha,
            "schema": model.schema,
        },
        "walls": {
            "total": len(walls),
            "excluded_demolition_count": len(demolition_guids),
            "excluded_demolition_global_ids": demolition_guids,
            "status_counts": dict(sorted(status_counts.items())),
            "visible_in_wall_plan": visible_count,
            "outside_or_not_intersecting_plan_cut": [
                record["global_id"] for record in walls if not record["visible_in_wall_plan"]
            ],
            "confirmed_new_dimensions": new_wall_dimension_checks,
        },
        "doors": {
            "total": len(hosted_doors) + len(unhosted_doors),
            "hosted_and_dimensioned": hosted_doors,
            "unhosted_delegated_to_a104": unhosted_doors,
        },
        "grid": {
            "u_axes": u_chain,
            "v_axes": v_chain,
            "u_closure": u_closure,
            "v_closure": v_closure,
            "dimension_segments": grid_dimensions,
            "includes_grid_08": True,
        },
        "annotation_qa": {
            "generated_label_count": len(occupied),
            "generated_collision_count": len(collisions),
            "collisions": collisions,
            "dimension_line_count": generated.count('PredefinedType-DIMENSION'),
        },
        "review_boundary": {
            "mechanically_complete": [
                "20 GridAxis inventory and both chain closures",
                "4 confirmed-new wall bounding dimensions",
                "5 hosted-door opening relationships and nominal sizes",
                "generated-label collision check",
                "current IFC and Wall Plan source hashes",
            ],
            "human_review_required": [
                "A-102 demolition-wall geometry derived from the handover DWG/PDF",
            ],
            "generator_writes_ifc": False,
            "wall_semantics_already_in_ifc": True,
        },
        "pass": qa_pass and len(walls) == 88 and len(demolition_guids) == 13 and len(hosted_doors) == 5 and len(unhosted_doors) == 3,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"candidate SVG: {args.output_svg}")
    print(f"wall register: {args.wall_register}")
    print(f"report: {args.report}")
    print(f"walls: {len(walls)}; hosted doors: {len(hosted_doors)}; unhosted doors: {len(unhosted_doors)}")
    print(f"grid closure: U={u_closure['residual_mm']:.6f} mm, V={v_closure['residual_mm']:.6f} mm")
    print(f"generated annotation collisions: {len(collisions)}")


if __name__ == "__main__":
    main()
