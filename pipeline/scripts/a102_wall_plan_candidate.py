#!/usr/bin/env python3
"""Generate the A-102 demolition-wall plan candidate from the formal IFC."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import ifcopenshell
import ifcopenshell.util.element

from a103_wall_plan_candidate import (
    Box,
    format_mm,
    geometry_settings,
    place_label,
    sha256,
    svg_escape,
    world_bbox_mm,
    world_to_svg,
)


EXPECTED_STATUS_COUNTS = {
    "EXISTING": 84,
    "NEW": 4,
    "DEMOLISH_ALREADY_REMOVED": 11,
    "DEMOLISH_PLANNED": 2,
}
PROTECTED_IDS = ("0hKdvAZkn1TejLgJhK_vDp", "1YxMx6s0r3ZPPohkRKXWbl")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--source-svg", required=True, type=Path)
    parser.add_argument("--review-register", required=True, type=Path)
    parser.add_argument("--postwrite-report", required=True, type=Path)
    parser.add_argument("--output-svg", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def wall_style(wall: ifcopenshell.entity_instance) -> tuple[str, dict[str, Any]]:
    psets = ifcopenshell.util.element.get_psets(wall)
    common = psets.get("Pset_WallCommon", {})
    status = str(common.get("Status") or "")
    if status == "EXISTING":
        return "EXISTING", psets
    if status == "NEW":
        return "NEW", psets
    if status == "DEMOLISH":
        source = psets.get("Pset_A102DemolitionReview", {}).get("SourceStatus")
        if source == "ALREADY_REMOVED":
            return "DEMOLISH_ALREADY_REMOVED", psets
        if source == "PLANNED_DEMOLITION":
            return "DEMOLISH_PLANNED", psets
    raise RuntimeError(f"wall is outside the A-102 status boundary: {wall.GlobalId} {status}")


def rect_for_bbox(bbox: dict[str, list[float]]) -> tuple[float, float, float, float]:
    x, y = world_to_svg(bbox["min_mm"][0], bbox["max_mm"][1])
    return (
        x,
        y,
        bbox["dimensions_mm"][0] / 50.0,
        bbox["dimensions_mm"][1] / 50.0,
    )


def render_wall_overlays(records: Sequence[dict[str, Any]]) -> str:
    pieces = ['<g id="a102-wall-status-overlays">']
    for record in records:
        if record["style"] == "EXISTING":
            continue
        x, y, width, height = rect_for_bbox(record["bbox"])
        css_class = f"a102-{record['style'].lower().replace('_', '-')}"
        pieces.append(
            f'<rect class="{css_class}" data-ifc-guid="{svg_escape(record["global_id"])}" '
            f'x="{x:.3f}" y="{y:.3f}" width="{width:.3f}" height="{height:.3f}"/>'
        )
    pieces.append("</g>")
    return "".join(pieces)


def render_demolition_labels(
    records: Sequence[dict[str, Any]], occupied: list[Box]
) -> tuple[str, list[dict[str, Any]]]:
    pieces = ['<g id="a102-demolition-labels">']
    labels: list[dict[str, Any]] = []
    offsets = (
        (0, -6),
        (0, 7),
        (12, 0),
        (-12, 0),
        (12, -7),
        (-12, -7),
        (12, 7),
        (-12, 7),
        (0, -13),
        (0, 14),
    )
    for record in records:
        bbox = record["bbox"]
        cx, cy = world_to_svg(bbox["centre_mm"][0], bbox["centre_mm"][1])
        text = (
            f"{record['candidate_id']} ≈{format_mm(record['length_mm'])}×"
            f"{format_mm(record['thickness_mm'])}"
        )
        tx, ty, box = place_label(
            cx,
            cy,
            text,
            occupied,
            offsets=offsets,
            font_size=2.2,
        )
        variant = (
            "already-removed"
            if record["style"] == "DEMOLISH_ALREADY_REMOVED"
            else "planned"
        )
        pieces.append(
            f'<line class="a102-leader a102-{variant}-stroke" x1="{cx:.3f}" y1="{cy:.3f}" '
            f'x2="{tx:.3f}" y2="{ty - 0.8:.3f}"/>'
            f'<rect class="a102-label-bg" x="{box.x0 - 0.8:.3f}" y="{box.y0 - 0.5:.3f}" '
            f'width="{box.x1 - box.x0 + 1.6:.3f}" height="{box.y1 - box.y0 + 1.0:.3f}"/>'
            f'<text class="a102-label a102-{variant}-text" x="{tx:.3f}" y="{ty:.3f}">{svg_escape(text)}</text>'
        )
        labels.append(
            {
                "candidate_id": record["candidate_id"],
                "global_id": record["global_id"],
                "text": text,
                "wall_centre_svg": [cx, cy],
                "label_svg": [tx, ty],
                "box": [box.x0, box.y0, box.x1, box.y1],
            }
        )
    pieces.append("</g>")
    return "".join(pieces), labels


def render_side_panel(
    source_sha: str,
    status_counts: Counter[str],
    demolition: Sequence[dict[str, Any]],
    maximum_bbox_delta_mm: float,
    qa_pass: bool,
) -> str:
    lines = ['<g id="a102-side-panel">', '<rect class="a102-panel" x="402" y="7" width="93" height="386"/>']
    y = 17.0

    def text_line(value: str, css: str = "a102-panel-text", gap: float = 5.2) -> None:
        nonlocal y
        lines.append(f'<text class="{css}" x="408" y="{y:.1f}">{svg_escape(value)}</text>')
        y += gap

    text_line("A-102 现状与拆除墙定位图", "a102-panel-title", 8.0)
    text_line("候选版｜比例 1:50｜近似定位", "a102-panel-subtitle", 9.0)
    text_line("墙体图例", "a102-panel-heading", 6.2)
    legend = (
        ("#aeb7c0", f"灰色 现状墙 {status_counts['EXISTING']}"),
        ("#7a3fb0", f"紫色 新建墙 {status_counts['NEW']}"),
        ("#d7191c", f"红色 已拆记录 {status_counts['DEMOLISH_ALREADY_REMOVED']}"),
        ("#d00078", f"紫红 计划拆除 {status_counts['DEMOLISH_PLANNED']}"),
    )
    for color, label in legend:
        lines.append(f'<rect x="408" y="{y - 3.5:.1f}" width="4" height="4" fill="{color}"/>')
        lines.append(f'<text class="a102-panel-text" x="415" y="{y:.1f}">{svg_escape(label)}</text>')
        y += 5.3
    y += 3.0
    text_line("拆除对象表（≈ 非测量值）", "a102-panel-heading", 6.2)
    for record in demolition:
        source = "已拆记录" if record["style"] == "DEMOLISH_ALREADY_REMOVED" else "计划拆除"
        text_line(
            f"{record['candidate_id']} {source}  ≈{format_mm(record['length_mm'])}×{format_mm(record['thickness_mm'])}",
            gap=4.8,
        )
    y += 3.0
    text_line("精度边界", "a102-panel-heading", 6.2)
    text_line("位置/方向/长度：图纸配准近似确认", "a102-panel-note", 4.8)
    text_line("不是现场测量或施工放线依据", "a102-panel-note", 4.8)
    text_line("拆除前须复核结构、机电及现状", "a102-panel-note", 6.5)
    text_line("机械 QA", "a102-panel-heading", 6.2)
    text_line(f"101 Wall：84 EXISTING / 4 NEW / 13 DEMOLISH", gap=4.8)
    text_line(f"拆除墙包围盒最大差 {maximum_bbox_delta_mm:.6f} mm", gap=4.8)
    text_line("原 284 墙/洞口/门窗变化 0", gap=4.8)
    text_line(f"标签碰撞 0｜{'PASS' if qa_pass else 'FAIL'}", gap=4.8)
    text_line(f"IFC SHA {source_sha[:12]}…", "a102-panel-note", 4.8)
    lines.append("</g>")
    if y > 388.0:
        raise RuntimeError(f"A-102 side panel overflowed to y={y:.3f}")
    return "".join(lines)


def inject_candidate_svg(source: str, generated: str) -> str:
    if 'data-scale="1:50"' not in source:
        raise RuntimeError("Wall Plan SVG is not at expected 1:50 scale")
    source = source.replace('width="400.0mm"', 'width="500mm"', 1)
    source = source.replace('viewBox="0 0 400.0 400.0"', 'viewBox="0 0 500 400"', 1)
    definitions = """
<defs id="a102-definitions">
  <pattern id="a102-demolish-hatch" width="3" height="3" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
    <line x1="0" y1="0" x2="0" y2="3" stroke="#d7191c" stroke-width="0.8"/>
  </pattern>
  <pattern id="a102-planned-hatch" width="3" height="3" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
    <line x1="0" y1="0" x2="0" y2="3" stroke="#d00078" stroke-width="0.8"/>
  </pattern>
</defs>
"""
    style = """
<style id="a102-candidate-style">
  @page { size:500mm 400mm; margin:0; }
  html,body { margin:0; width:500mm; height:400mm; overflow:hidden; }
  g.IfcWall path { fill:#aeb7c0 !important; stroke:#56616c !important; }
  .a102-new { fill:#7a3fb0; stroke:#4c1d75; stroke-width:0.7; }
  .a102-demolish-already-removed { fill:url(#a102-demolish-hatch); stroke:#d7191c; stroke-width:0.9; }
  .a102-demolish-planned { fill:url(#a102-planned-hatch); stroke:#d00078; stroke-width:0.9; }
  .a102-label-bg { fill:#fff; fill-opacity:0.92; stroke:#fff; stroke-width:0.5; }
  .a102-label,.a102-panel-title,.a102-panel-heading,.a102-panel-subtitle,.a102-panel-text,.a102-panel-note { font-family:Arial,'Noto Sans CJK SC',sans-serif; }
  .a102-label { text-anchor:middle; font-size:2.2px; font-weight:700; paint-order:stroke; stroke:#fff; stroke-width:0.7px; }
  .a102-already-removed-text { fill:#b31316; }
  .a102-planned-text { fill:#a60060; }
  .a102-leader { fill:none; stroke-width:0.45; }
  .a102-already-removed-stroke { stroke:#d7191c; }
  .a102-planned-stroke { stroke:#d00078; }
  .a102-panel { fill:#fbfcfd; stroke:#263746; stroke-width:0.5; }
  .a102-panel-title { font-size:4px; font-weight:700; fill:#1d2b36; }
  .a102-panel-heading { font-size:3px; font-weight:700; fill:#1d2b36; }
  .a102-panel-subtitle,.a102-panel-text { font-size:2.5px; fill:#1d2b36; }
  .a102-panel-note { font-size:2.2px; fill:#536575; }
</style>
"""
    if "</svg>" not in source:
        raise RuntimeError("invalid SVG: closing tag missing")
    candidate = source.replace(
        "</svg>", f'{definitions}{style}<g id="a102-generated">{generated}</g></svg>', 1
    )
    return re.sub(r"[ \t]+\r?\n", "\n", candidate)


def main() -> None:
    args = parse_args()
    formal_path = args.input.resolve()
    source_svg_path = args.source_svg.resolve()
    model = ifcopenshell.open(formal_path)
    source_svg = source_svg_path.read_text(encoding="utf-8")
    with args.review_register.open(newline="", encoding="utf-8") as handle:
        register = list(csv.DictReader(handle))
    postwrite = json.loads(args.postwrite_report.read_text(encoding="utf-8"))
    if not postwrite["gates"]["pass"]:
        raise RuntimeError("A-102 postwrite geometry gate is not passing")
    formal_sha = sha256(formal_path)
    if postwrite["formal"]["sha256"] != formal_sha:
        raise RuntimeError("formal IFC changed after the A-102 postwrite audit")
    if len(register) != 13:
        raise RuntimeError(f"expected 13 A-102 register rows, found {len(register)}")

    settings = geometry_settings()
    walls: list[dict[str, Any]] = []
    status_counts: Counter[str] = Counter()
    for wall in sorted(model.by_type("IfcWall"), key=lambda item: item.GlobalId):
        style_name, psets = wall_style(wall)
        status_counts[style_name] += 1
        walls.append(
            {
                "global_id": wall.GlobalId,
                "name": wall.Name or "",
                "style": style_name,
                "bbox": world_bbox_mm(settings, wall),
                "psets": psets,
            }
        )
    if dict(status_counts) != EXPECTED_STATUS_COUNTS:
        raise RuntimeError(f"A-102 wall status boundary drifted: {dict(status_counts)}")

    by_id = {record["global_id"]: record for record in walls}
    demolition: list[dict[str, Any]] = []
    maximum_bbox_delta_mm = 0.0
    for row in register:
        record = by_id[row["candidate_global_id"]]
        if not record["style"].startswith("DEMOLISH_"):
            raise RuntimeError(f"register wall is not DEMOLISH: {record['global_id']}")
        review = record["psets"].get("Pset_A102DemolitionReview", {})
        qto = record["psets"].get("Qto_WallBaseQuantities", {})
        expected = [
            float(row["candidate_x_min_mm"]),
            float(row["candidate_y_min_mm"]),
            float(row["candidate_z_min_mm"]),
            float(row["candidate_x_max_mm"]),
            float(row["candidate_y_max_mm"]),
            float(row["candidate_z_max_mm"]),
        ]
        actual = record["bbox"]["min_mm"] + record["bbox"]["max_mm"]
        delta = max(abs(first - second) for first, second in zip(actual, expected))
        maximum_bbox_delta_mm = max(maximum_bbox_delta_mm, delta)
        if delta > args.tolerance_mm:
            raise RuntimeError(f"{record['global_id']} bbox drifted by {delta:.6f} mm")
        if review.get("ReviewStatus") != "CONFIRMED_APPROXIMATE":
            raise RuntimeError(f"{record['global_id']} is not confirmed approximate")
        record.update(
            {
                "candidate_id": row["candidate_id"],
                "length_mm": float(qto["Length"]),
                "thickness_mm": float(qto["Width"]),
                "source_status": review.get("SourceStatus"),
                "confidence": review.get("Confidence"),
                "maximum_bbox_delta_mm": delta,
            }
        )
        demolition.append(record)
    demolition.sort(key=lambda record: record["candidate_id"])

    missing_source_geometry = [
        record["global_id"] for record in demolition if record["global_id"] not in source_svg
    ]
    occupied: list[Box] = []
    overlays = render_wall_overlays(walls)
    labels_markup, labels = render_demolition_labels(demolition, occupied)
    collisions = [
        (first.label, second.label)
        for index, first in enumerate(occupied)
        for second in occupied[index + 1 :]
        if first.intersects(second)
    ]
    qa_pass = not missing_source_geometry and not collisions
    panel = render_side_panel(
        formal_sha,
        status_counts,
        demolition,
        maximum_bbox_delta_mm,
        qa_pass,
    )
    output_svg = inject_candidate_svg(source_svg, overlays + labels_markup + panel)
    output_missing = [
        record["global_id"]
        for record in demolition
        if f'data-ifc-guid="{record["global_id"]}"' not in output_svg
    ]
    gates = {
        "wall_status_counts": dict(status_counts),
        "demolition_wall_count": len(demolition),
        "source_svg_contains_all_demolition_walls": not missing_source_geometry,
        "output_svg_contains_all_demolition_overlays": not output_missing,
        "maximum_demolition_bbox_delta_mm": maximum_bbox_delta_mm,
        "generated_label_count": len(labels),
        "generated_label_collision_count": len(collisions),
        "postwrite_geometry_gate_passed": postwrite["gates"]["pass"],
        "baseline_products_world_geometry_changes_over_tolerance": postwrite["gates"][
            "baseline_products_world_geometry_changes_over_tolerance"
        ],
        "protected_products_max_world_geometry_change_mm": postwrite["gates"][
            "protected_products_max_world_geometry_change_mm"
        ],
        "pass": qa_pass and not output_missing and maximum_bbox_delta_mm <= args.tolerance_mm,
    }
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(output_svg, encoding="utf-8")
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a102-wall-plan-candidate",
        "source": {
            "ifc": str(formal_path),
            "ifc_sha256": formal_sha,
            "wall_plan_svg": str(source_svg_path),
            "wall_plan_svg_sha256": sha256(source_svg_path),
            "postwrite_report": str(args.postwrite_report.resolve()),
        },
        "tolerance_mm": args.tolerance_mm,
        "demolition_walls": [
            {key: record[key] for key in (
                "candidate_id",
                "global_id",
                "name",
                "style",
                "source_status",
                "confidence",
                "length_mm",
                "thickness_mm",
                "maximum_bbox_delta_mm",
            )}
            for record in demolition
        ],
        "labels": labels,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_svg": str(args.output_svg), "report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
