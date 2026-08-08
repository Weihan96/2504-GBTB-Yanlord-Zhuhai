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


def side_panel(records: list[dict[str, Any]], source_sha: str) -> str:
    counts = Counter(record["candidate_finish"] for record in records)
    tadelakt_spaces = sorted({record["space_long_name"] for record in records if record["candidate_finish_code"] == "TADELAKT"})
    rows = [
        ("IFC 饰面对象", str(len(records))),
        ("Tadelakt 候选", str(counts["Tadelakt"])),
        ("大白墙候选", str(counts["大白墙"])),
        ("唯一房间归属", f'{sum(record["space_match_count"] == 1 for record in records)}/{len(records)}'),
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
        "主卫干区当前无可归属的 CLADDING 对象",
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


def plan_indicator(record: dict[str, Any], class_name: str) -> str:
    bbox = record["bbox"]
    minimum = bbox["min_mm"]
    maximum = bbox["max_mm"]
    centre = bbox["centre_mm"]
    if bbox["dimensions_mm"][0] >= bbox["dimensions_mm"][1]:
        x1, y1 = (float(minimum[0]) + 10000.0) / 50.0, (10000.0 - float(centre[1])) / 50.0
        x2, y2 = (float(maximum[0]) + 10000.0) / 50.0, y1
    else:
        x1, y1 = (float(centre[0]) + 10000.0) / 50.0, (10000.0 - float(minimum[1])) / 50.0
        x2, y2 = x1, (10000.0 - float(maximum[1])) / 50.0
    return (
        f'<line class="wfin-indicator {class_name}" data-guid="{html.escape(record["covering_global_id"], quote=True)}" '
        f'x1="{x1:.6f}" y1="{y1:.6f}" x2="{x2:.6f}" y2="{y2:.6f}"/>'
    )


def main() -> None:
    args = parse_args()
    source_sha = sha256(args.input)
    model = ifcopenshell.open(args.input)
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
        matches = [
            space
            for space in spaces
            if centre_in_space_bbox(bbox["centre_mm"], space["bbox"], args.containment_tolerance_mm)
        ]
        if len(matches) != 1:
            raise RuntimeError(f"{covering.GlobalId} matches {len(matches)} Space bboxes; expected exactly one")
        space = matches[0]
        is_tadelakt = space["long_name"] in TADELAKT_SPACE_NAMES
        records.append(
            {
                "covering_global_id": covering.GlobalId,
                "covering_name": str(covering.Name or ""),
                "effective_predefined_type": "CLADDING",
                "current_material": material_name(covering),
                "space_global_id": space["global_id"],
                "space_long_name": space["long_name"],
                "space_match_count": len(matches),
                "candidate_finish_code": "TADELAKT" if is_tadelakt else "WHITE_WALL",
                "candidate_finish": "Tadelakt" if is_tadelakt else "大白墙",
                "basis": "用户确认房间级材料规则；饰面世界包围盒中心唯一落入当前 Space 世界包围盒",
                "confidence": 1.0,
                "review_required": "no",
                "status": "confirmed_candidate",
                "formal_ifc_write_allowed": "no",
                "bbox": bbox,
            }
        )
    records.sort(key=lambda record: (record["space_long_name"], record["covering_global_id"]))
    if len(records) != 51:
        raise RuntimeError(f"expected 51 effective CLADDING objects, got {len(records)}")
    counts = Counter(record["candidate_finish_code"] for record in records)
    if counts != Counter({"WHITE_WALL": 27, "TADELAKT": 24}):
        raise RuntimeError(f"unexpected finish counts: {dict(counts)}")

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
.wfin-white-wall rect { fill: #dce3e8 !important; stroke: #68727a !important; stroke-width: 0.35 !important; }
.wfin-tadelakt rect { fill: #65b9a9 !important; stroke: #156b60 !important; stroke-width: 0.45 !important; }
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
        class_name = "wfin-tadelakt" if record["candidate_finish_code"] == "TADELAKT" else "wfin-white-wall"
        indicators.append(plan_indicator(record, class_name))
        svg, count = add_class(svg, record["covering_global_id"], class_name)
        if count == 0:
            fallback_fragments.append(fallback_plan_rect(record, class_name))
        else:
            fragment_object_count += 1
        fragment_count += count
    svg = svg.replace(
        "</svg>",
        '<g id="wfin-ifc-fallbacks">' + "".join(fallback_fragments) + "</g>"
        + '<g id="wfin-plan-indicators">' + "".join(indicators) + "</g>"
        + side_panel(records, source_sha)
        + "</svg>",
        1,
    )
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(svg, encoding="utf-8")

    main_bath_dry_objects = [record for record in records if record["space_long_name"] == "主卫干区"]
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {"path": str(args.input), "sha256": source_sha, "schema": model.schema},
        "source_svg": {"path": str(args.source_svg), "sha256": sha256(args.source_svg)},
        "automatic_ifc_write_allowed": False,
        "cladding_count": len(records),
        "unique_space_assignment_count": sum(record["space_match_count"] == 1 for record in records),
        "tadelakt_count": counts["TADELAKT"],
        "white_wall_count": counts["WHITE_WALL"],
        "source_svg_fragment_count": fragment_count,
        "source_svg_object_count": fragment_object_count,
        "fallback_object_count": len(fallback_fragments),
        "main_bath_dry_cladding_count": len(main_bath_dry_objects),
        "controlled_gap": "主卫干区当前没有可按中心唯一归属的 CLADDING 对象；不自动新增饰面几何",
        "unresolved_parameters": ["品牌/产品系统", "颜色/样板", "完成面总厚度", "基层", "湿区防水与收口节点"],
        "references": [
            "https://mp.weixin.qq.com/s/roqN3h93WZ0AER71lJzaUQ",
            "https://www.modamuri.com/article195.html",
            "https://www.modamuri.com/ysq.html",
            "https://us.tadelakt.com/tadelakt/",
        ],
        "qa": {
            "effective_cladding_count_51": len(records) == 51,
            "every_cladding_has_one_space": all(record["space_match_count"] == 1 for record in records),
            "finish_partition_complete": counts["TADELAKT"] + counts["WHITE_WALL"] == len(records),
            "candidate_visual_contains_every_cladding": fragment_object_count + len(fallback_fragments) == len(records),
            "formal_ifc_unchanged": True,
        },
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        f"WFIN candidate: {len(records)} CLADDING, {counts['TADELAKT']} Tadelakt, "
        f"{counts['WHITE_WALL']} white wall, IFC unchanged"
    )


if __name__ == "__main__":
    main()
