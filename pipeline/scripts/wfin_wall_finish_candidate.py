#!/usr/bin/env python3
"""Generate the read-only WFIN wall-finish register and plan candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element
import numpy as np


TADELAKT_SPACE_NAMES = {
    "次卧",
    "次卧飘窗",
    "主卫干区",
    "主卫湿区",
    "主卫湿区飘窗",
    "客卫干区",
    "客卫",
    "客卫飘窗",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--source-svg", required=True, type=Path)
    parser.add_argument("--output-svg", required=True, type=Path)
    parser.add_argument("--register", required=True, type=Path)
    parser.add_argument("--segment-register", required=True, type=Path)
    parser.add_argument("--boundary-decisions", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--containment-tolerance-mm", type=float, default=0.5)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def shape_bbox(settings: ifcopenshell.geom.settings, product: ifcopenshell.entity_instance) -> dict[str, list[float]]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3)) * 1000.0
    if len(vertices) == 0:
        raise RuntimeError(f"{product.GlobalId} has no Body geometry")
    minimum = vertices.min(axis=0)
    maximum = vertices.max(axis=0)
    return {
        "min_mm": minimum.tolist(),
        "max_mm": maximum.tolist(),
        "centre_mm": ((minimum + maximum) / 2.0).tolist(),
        "dimensions_mm": (maximum - minimum).tolist(),
    }


def effective_predefined_type(product: ifcopenshell.entity_instance) -> str:
    value = str(product.PredefinedType or "")
    if value and value != "NOTDEFINED":
        return value
    product_type = ifcopenshell.util.element.get_type(product)
    return str(product_type.PredefinedType or "") if product_type is not None else ""


def material_name(product: ifcopenshell.entity_instance) -> str:
    material = ifcopenshell.util.element.get_material(product, should_skip_usage=True)
    if material is None:
        return ""
    if material.is_a("IfcMaterial"):
        return str(material.Name or "")
    if material.is_a("IfcMaterialLayerSet"):
        return str(material.LayerSetName or "")
    return str(getattr(material, "Name", "") or material)


def centre_in_space_bbox(centre: list[float], bbox: dict[str, list[float]], tolerance: float) -> bool:
    return all(
        bbox["min_mm"][axis] - tolerance <= centre[axis] <= bbox["max_mm"][axis] + tolerance
        for axis in range(3)
    )


def finish_segments(
    covering_bbox: dict[str, list[float]],
    spaces: list[dict[str, Any]],
    tolerance: float,
) -> tuple[str, list[dict[str, Any]]]:
    """Partition a thin vertical CLADDING along aligned Space boundaries."""
    minimum = covering_bbox["min_mm"]
    maximum = covering_bbox["max_mm"]
    dimensions = covering_bbox["dimensions_mm"]
    axis = 0 if dimensions[0] >= dimensions[1] else 1
    perpendicular = 1 - axis
    centre_perpendicular = covering_bbox["centre_mm"][perpendicular]
    centre_z = covering_bbox["centre_mm"][2]
    raw = []
    for space in spaces:
        bbox = space["bbox"]
        if not (
            bbox["min_mm"][perpendicular] - tolerance
            <= centre_perpendicular
            <= bbox["max_mm"][perpendicular] + tolerance
        ):
            continue
        if not bbox["min_mm"][2] - tolerance <= centre_z <= bbox["max_mm"][2] + tolerance:
            continue
        start = max(minimum[axis], bbox["min_mm"][axis])
        end = min(maximum[axis], bbox["max_mm"][axis])
        if end - start <= tolerance:
            continue
        raw.append(
            {
                "axis": "X" if axis == 0 else "Y",
                "constant_mm": centre_perpendicular,
                "start_mm": start,
                "end_mm": end,
                "space_global_ids": [space["global_id"]],
                "space_long_names": [space["long_name"]],
                "candidate_finish_code": "TADELAKT" if space["long_name"] in TADELAKT_SPACE_NAMES else "WHITE_WALL",
                "candidate_finish": "Tadelakt" if space["long_name"] in TADELAKT_SPACE_NAMES else "大白墙",
            }
        )
    raw.sort(key=lambda item: (item["start_mm"], item["end_mm"], item["space_long_names"]))
    merged: list[dict[str, Any]] = []
    for item in raw:
        if merged and item["start_mm"] < merged[-1]["end_mm"] - tolerance:
            raise RuntimeError("overlapping Space intervals make the CLADDING finish boundary ambiguous")
        if (
            merged
            and item["candidate_finish_code"] == merged[-1]["candidate_finish_code"]
            and abs(item["start_mm"] - merged[-1]["end_mm"]) <= tolerance
        ):
            merged[-1]["end_mm"] = max(merged[-1]["end_mm"], item["end_mm"])
            merged[-1]["space_global_ids"].extend(item["space_global_ids"])
            merged[-1]["space_long_names"].extend(item["space_long_names"])
        else:
            merged.append(item)
    covered = sum(item["end_mm"] - item["start_mm"] for item in merged)
    if not merged or dimensions[axis] - covered > tolerance * 2:
        raise RuntimeError(
            f"CLADDING segment coverage drift: length={dimensions[axis]:.3f} mm, covered={covered:.3f} mm"
        )
    return ("X" if axis == 0 else "Y"), merged


def read_boundary_decisions(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    decisions = {row["covering_global_id"]: row for row in rows}
    expected = {
        "0e0XOb$L18ZBVYJiQJrQ1p": "FULL_OBJECT_WHITE_WALL",
        "3bMoS7bIT8wBmjqRvdPunu": "FULL_OBJECT_WHITE_WALL",
        "06wFwLoDD6ie5iCTnc_yad": "DEFER_TO_FINAL_MATERIAL_REVIEW",
    }
    if {global_id: row["decision"] for global_id, row in decisions.items()} != expected:
        raise RuntimeError("WFIN material-boundary decision boundary drift")
    return decisions


def apply_boundary_decision(
    covering_bbox: dict[str, list[float]],
    segments: list[dict[str, Any]],
    decision: dict[str, str] | None,
) -> list[dict[str, Any]]:
    if not decision or decision["decision"] != "FULL_OBJECT_WHITE_WALL":
        return segments
    axis = segments[0]["axis"]
    axis_index = 0 if axis == "X" else 1
    return [
        {
            "axis": axis,
            "constant_mm": segments[0]["constant_mm"],
            "start_mm": covering_bbox["min_mm"][axis_index],
            "end_mm": covering_bbox["max_mm"][axis_index],
            "space_global_ids": [space_id for segment in segments for space_id in segment["space_global_ids"]],
            "space_long_names": [name for segment in segments for name in segment["space_long_names"]],
            "candidate_finish_code": "WHITE_WALL",
            "candidate_finish": "大白墙",
        }
    ]


def add_class(svg: str, global_id: str, class_name: str) -> tuple[str, int]:
    pattern = re.compile(
        rf'(<g\b[^>]*\bclass=")([^"]*)("[^>]*\bifc:guid="{re.escape(global_id)}"[^>]*>)'
    )

    def replace(match: re.Match[str]) -> str:
        classes = match.group(2).split()
        if class_name not in classes:
            classes.append(class_name)
        return f'{match.group(1)}{" ".join(classes)}{match.group(3)}'

    return pattern.subn(replace, svg)


def side_panel(records: list[dict[str, Any]], segments: list[dict[str, Any]], source_sha: str) -> str:
    object_counts = Counter(record["candidate_finish_code"] for record in records)
    segment_counts = Counter(segment["candidate_finish_code"] for segment in segments)
    tadelakt_spaces = sorted({name for segment in segments if segment["candidate_finish_code"] == "TADELAKT" for name in segment["space_long_names"]})
    rows = [
        ("IFC 饰面对象", str(len(records))),
        ("跨材料边界对象", str(object_counts["MULTI_FINISH_SPLIT_REQUIRED"])),
        ("Tadelakt 分段", str(segment_counts["TADELAKT"])),
        ("大白墙分段", str(segment_counts["WHITE_WALL"])),
    ]
    y = 22.0
    parts = [
        '<g id="wfin-side-panel">',
        '<rect x="402" y="2" width="96" height="396" rx="1.5" fill="#fbfbf8" stroke="#2b2b2b" stroke-width="0.35"/>',
        '<text x="406" y="10" font-size="5.2" font-weight="700">WFIN 墙面材料候选</text>',
        '<text x="406" y="15" font-size="2.8" fill="#555">只读候选｜不写正式 IFC</text>',
    ]
    for label, value in rows:
        parts.append(f'<text x="406" y="{y:.1f}" font-size="3.0">{html.escape(label)}：{html.escape(value)}</text>')
        y += 5.0
    y += 3.0
    parts.extend(
        [
            f'<rect x="406" y="{y:.1f}" width="8" height="4" fill="#dce3e8" stroke="#68727a" stroke-width="0.25"/>',
            f'<text x="417" y="{y + 3.2:.1f}" font-size="3.0">大白墙（默认）</text>',
            f'<rect x="406" y="{y + 7:.1f}" width="8" height="4" fill="#65b9a9" stroke="#156b60" stroke-width="0.25"/>',
            f'<text x="417" y="{y + 10.2:.1f}" font-size="3.0">Tadelakt</text>',
        ]
    )
    y += 18.0
    parts.append(f'<text x="406" y="{y:.1f}" font-size="3.2" font-weight="700">Tadelakt 房间范围</text>')
    y += 5.0
    for name in tadelakt_spaces:
        parts.append(f'<text x="408" y="{y:.1f}" font-size="2.8">• {html.escape(name)}</text>')
        y += 4.2
    y += 3.0
    notes = [
        "规则：次卧及其飘窗、两卫干/湿区及飘窗",
        "当前画面 2 件已确认大白墙；1 件延后判断",
        "待定：品牌系统、颜色、厚度、基层、防水节点",
        "参考链接只作材质方向证据，不作施工参数",
        f"IFC SHA {source_sha[:12]}…",
    ]
    parts.append(f'<text x="406" y="{y:.1f}" font-size="3.2" font-weight="700">边界与停止条件</text>')
    y += 5.0
    for note in notes:
        parts.append(f'<text x="408" y="{y:.1f}" font-size="2.6">{html.escape(note)}</text>')
        y += 4.0
    parts.append('</g>')
    return "".join(parts)


def fallback_plan_rect(record: dict[str, Any], class_name: str) -> str:
    """Draw a top-view bbox only when the Bonsai source SVG is older than the IFC."""
    bbox = record["bbox"]
    x = (float(bbox["min_mm"][0]) + 10000.0) / 50.0
    y = (10000.0 - float(bbox["max_mm"][1])) / 50.0
    width = float(bbox["dimensions_mm"][0]) / 50.0
    height = float(bbox["dimensions_mm"][1]) / 50.0
    if width < 0.4:
        x -= (0.4 - width) / 2.0
        width = 0.4
    if height < 0.4:
        y -= (0.4 - height) / 2.0
        height = 0.4
    return (
        f'<g class="IfcCovering projection {class_name} wfin-fallback" '
        f'ifc:guid="{html.escape(record["covering_global_id"], quote=True)}">'
        f'<rect x="{x:.6f}" y="{y:.6f}" width="{width:.6f}" height="{height:.6f}"/>'
        "</g>"
    )


def plan_indicator(segment: dict[str, Any], class_name: str) -> str:
    if segment["axis"] == "X":
        x1, y1 = (float(segment["start_mm"]) + 10000.0) / 50.0, (10000.0 - float(segment["constant_mm"])) / 50.0
        x2, y2 = (float(segment["end_mm"]) + 10000.0) / 50.0, y1
    else:
        x1, y1 = (float(segment["constant_mm"]) + 10000.0) / 50.0, (10000.0 - float(segment["start_mm"])) / 50.0
        x2, y2 = x1, (10000.0 - float(segment["end_mm"])) / 50.0
    return (
        f'<line class="wfin-indicator {class_name}" data-guid="{html.escape(segment["covering_global_id"], quote=True)}" '
        f'x1="{x1:.6f}" y1="{y1:.6f}" x2="{x2:.6f}" y2="{y2:.6f}"/>'
    )


def main() -> None:
    args = parse_args()
    source_sha = sha256(args.input)
    model = ifcopenshell.open(args.input)
    boundary_decisions = read_boundary_decisions(args.boundary_decisions)
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)

    spaces = []
    for space in model.by_type("IfcSpace"):
        spaces.append(
            {
                "global_id": space.GlobalId,
                "long_name": str(space.LongName or space.Name or ""),
                "bbox": shape_bbox(settings, space),
            }
        )

    records = []
    for covering in model.by_type("IfcCovering"):
        if effective_predefined_type(covering) != "CLADDING":
            continue
        bbox = shape_bbox(settings, covering)
        centre_matches = [
            space
            for space in spaces
            if centre_in_space_bbox(bbox["centre_mm"], space["bbox"], args.containment_tolerance_mm)
        ]
        if len(centre_matches) != 1:
            raise RuntimeError(f"{covering.GlobalId} centre matches {len(centre_matches)} Space bboxes; expected exactly one")
        _, segments = finish_segments(bbox, spaces, args.containment_tolerance_mm)
        boundary_decision = boundary_decisions.get(covering.GlobalId)
        segments = apply_boundary_decision(bbox, segments, boundary_decision)
        finish_codes = {segment["candidate_finish_code"] for segment in segments}
        mixed = len(finish_codes) > 1
        deferred = bool(boundary_decision and boundary_decision["decision"] == "DEFER_TO_FINAL_MATERIAL_REVIEW")
        space_ids = [space_id for segment in segments for space_id in segment["space_global_ids"]]
        space_names = [name for segment in segments for name in segment["space_long_names"]]
        for segment in segments:
            segment["covering_global_id"] = covering.GlobalId
        records.append(
            {
                "covering_global_id": covering.GlobalId,
                "covering_name": str(covering.Name or ""),
                "effective_predefined_type": "CLADDING",
                "current_material": material_name(covering),
                "space_global_id": "; ".join(space_ids),
                "space_long_name": "; ".join(space_names),
                "space_match_count": len(space_ids),
                "candidate_finish_code": "MULTI_FINISH_SPLIT_REQUIRED" if mixed else next(iter(finish_codes)),
                "candidate_finish": "按 Space 边界拆分" if mixed else segments[0]["candidate_finish"],
                "basis": "用户确认房间级材料规则；薄型垂直饰面沿 Grid Space 边界机械分段",
                "confidence": 1.0,
                "review_required": "yes" if deferred else "no",
                "status": "deferred_material_review" if deferred else "confirmed_candidate",
                "formal_ifc_write_allowed": "no",
                "bbox": bbox,
                "segments": segments,
            }
        )
    records.sort(key=lambda record: (record["space_long_name"], record["covering_global_id"]))
    if len(records) != 51:
        raise RuntimeError(f"expected 51 effective CLADDING objects, got {len(records)}")
    counts = Counter(record["candidate_finish_code"] for record in records)
    all_segments = [segment for record in records for segment in record["segments"]]
    segment_counts = Counter(segment["candidate_finish_code"] for segment in all_segments)
    if counts != Counter({"WHITE_WALL": 27, "TADELAKT": 23, "MULTI_FINISH_SPLIT_REQUIRED": 1}):
        raise RuntimeError(f"unexpected object finish counts: {dict(counts)}")
    if segment_counts != Counter({"WHITE_WALL": 28, "TADELAKT": 24}) or len(all_segments) != 52:
        raise RuntimeError(f"unexpected finish segments: count={len(all_segments)}, finishes={dict(segment_counts)}")

    args.register.parent.mkdir(parents=True, exist_ok=True)
    with args.register.open("w", encoding="utf-8-sig", newline="") as handle:
        fieldnames = [
            "covering_global_id",
            "covering_name",
            "effective_predefined_type",
            "current_material",
            "space_global_id",
            "space_long_name",
            "space_match_count",
            "candidate_finish_code",
            "candidate_finish",
            "basis",
            "confidence",
            "review_required",
            "status",
            "formal_ifc_write_allowed",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({key: record[key] for key in fieldnames})

    args.segment_register.parent.mkdir(parents=True, exist_ok=True)
    segment_fields = [
        "segment_id", "covering_global_id", "axis", "constant_mm", "start_mm", "end_mm",
        "space_global_ids", "space_long_names", "candidate_finish_code", "candidate_finish",
        "review_required", "status", "formal_ifc_write_allowed",
    ]
    with args.segment_register.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=segment_fields, lineterminator="\n")
        writer.writeheader()
        for index, segment in enumerate(all_segments, 1):
            row = dict(segment)
            row.update(
                {
                    "segment_id": f"WFIN-S{index:03d}",
                    "space_global_ids": "; ".join(segment["space_global_ids"]),
                    "space_long_names": "; ".join(segment["space_long_names"]),
                }
            )
            source_record = next(record for record in records if record["covering_global_id"] == segment["covering_global_id"])
            row["review_required"] = source_record["review_required"]
            row["status"] = "deferred_material_review" if source_record["review_required"] == "yes" else "confirmed_candidate"
            row["formal_ifc_write_allowed"] = "no"
            writer.writerow({key: row[key] for key in segment_fields})

    svg = args.source_svg.read_text(encoding="utf-8")
    svg = re.sub(r'width="400(?:\.00003814697266)?mm"', 'width="500mm"', svg, count=1)
    svg = re.sub(r'height="400(?:\.00003814697266)?mm"', 'height="400mm"', svg, count=1)
    svg = re.sub(r'viewBox="0 0 400(?:\.00003814697266)? 400(?:\.00003814697266)?"', 'viewBox="0 0 500 400"', svg, count=1)
    css = """
/* WFIN read-only candidate overlay */
@page { size: 500mm 400mm; margin: 0; }
html, body { margin: 0; padding: 0; width: 500mm; height: 400mm; }
svg { display: block; }
* { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
.wfin-white-wall path { fill: #dce3e8 !important; stroke: #68727a !important; stroke-width: 0.35 !important; }
.wfin-tadelakt path { fill: #65b9a9 !important; stroke: #156b60 !important; stroke-width: 0.45 !important; }
.wfin-split-required path { fill: #f1d4a9 !important; stroke: #c65f00 !important; stroke-width: 0.55 !important; }
.wfin-white-wall rect { fill: #dce3e8 !important; stroke: #68727a !important; stroke-width: 0.35 !important; }
.wfin-tadelakt rect { fill: #65b9a9 !important; stroke: #156b60 !important; stroke-width: 0.45 !important; }
.wfin-split-required rect { fill: #f1d4a9 !important; stroke: #c65f00 !important; stroke-width: 0.55 !important; }
.wfin-fallback rect { stroke-dasharray: 1.2 0.7; }
.wfin-indicator { fill: none !important; stroke-width: 1.3 !important; stroke-linecap: butt; opacity: 0.92; }
.wfin-indicator.wfin-white-wall { stroke: #8e9ba5 !important; }
.wfin-indicator.wfin-tadelakt { stroke: #00a88f !important; }
"""
    svg = svg.replace("</style>", css + "</style>", 1)
    fragment_count = 0
    fragment_object_count = 0
    fallback_fragments = []
    indicators = []
    for record in records:
        class_name = {
            "TADELAKT": "wfin-tadelakt",
            "WHITE_WALL": "wfin-white-wall",
            "MULTI_FINISH_SPLIT_REQUIRED": "wfin-split-required",
        }[record["candidate_finish_code"]]
        svg, count = add_class(svg, record["covering_global_id"], class_name)
        if count == 0:
            fallback_fragments.append(fallback_plan_rect(record, class_name))
        else:
            fragment_object_count += 1
        fragment_count += count
    for segment in all_segments:
        class_name = "wfin-tadelakt" if segment["candidate_finish_code"] == "TADELAKT" else "wfin-white-wall"
        indicators.append(plan_indicator(segment, class_name))
    svg = svg.replace(
        "</svg>",
        '<g id="wfin-ifc-fallbacks">' + "".join(fallback_fragments) + "</g>"
        + '<g id="wfin-plan-indicators">' + "".join(indicators) + "</g>"
        + side_panel(records, all_segments, source_sha)
        + "</svg>",
        1,
    )
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(svg, encoding="utf-8")

    main_bath_dry_segments = [segment for segment in all_segments if "主卫干区" in segment["space_long_names"]]
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_ifc_sha256": source_sha,
        "source": {"path": str(args.input), "sha256": source_sha, "schema": model.schema},
        "source_svg": {"path": str(args.source_svg), "sha256": sha256(args.source_svg)},
        "automatic_ifc_write_allowed": False,
        "cladding_count": len(records),
        "single_finish_object_count": counts["TADELAKT"] + counts["WHITE_WALL"],
        "mixed_finish_object_count": counts["MULTI_FINISH_SPLIT_REQUIRED"],
        "finish_segment_count": len(all_segments),
        "tadelakt_segment_count": segment_counts["TADELAKT"],
        "white_wall_segment_count": segment_counts["WHITE_WALL"],
        "source_svg_fragment_count": fragment_count,
        "source_svg_object_count": fragment_object_count,
        "fallback_object_count": len(fallback_fragments),
        "main_bath_dry_segment_count": len(main_bath_dry_segments),
        "mixed_finish_object_ids": [record["covering_global_id"] for record in records if record["candidate_finish_code"] == "MULTI_FINISH_SPLIT_REQUIRED"],
        "controlled_gap": "当前画面两块橙色饰面已确认整件为大白墙；剩余 1 个跨界对象延后至最终统一选材",
        "unresolved_parameters": ["品牌/产品系统", "颜色/样板", "完成面总厚度", "基层", "湿区防水与收口节点"],
        "references": [
            "https://mp.weixin.qq.com/s/roqN3h93WZ0AER71lJzaUQ",
            "https://www.modamuri.com/article195.html",
            "https://www.modamuri.com/ysq.html",
            "https://us.tadelakt.com/tadelakt/",
        ],
        "qa": {
            "effective_cladding_count_51": len(records) == 51,
            "finish_segment_partition_complete": segment_counts["TADELAKT"] + segment_counts["WHITE_WALL"] == len(all_segments),
            "mixed_finish_objects_blocked_from_whole_object_assignment": counts["MULTI_FINISH_SPLIT_REQUIRED"] == 1,
            "candidate_visual_contains_every_cladding": fragment_object_count + len(fallback_fragments) == len(records),
            "formal_ifc_unchanged": True,
        },
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        f"WFIN candidate: {len(records)} CLADDING, {len(all_segments)} finish segments "
        f"({segment_counts['TADELAKT']} Tadelakt, {segment_counts['WHITE_WALL']} white wall), "
        f"1 deferred split-review object, IFC unchanged"
    )


if __name__ == "__main__":
    main()
