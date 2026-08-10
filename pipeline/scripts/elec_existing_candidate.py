#!/usr/bin/env python3
"""Generate a read-only ELEC existing-condition candidate from the formal IFC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import numpy as np
from ifcopenshell.util.element import get_container, get_psets, get_type
from ifcopenshell.util.placement import get_local_placement
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union


EXPECTED_SOURCE_SHA256 = "6c2fd8da9e9ad7ddbc2b63415a27f1c979e8995b880d8fce210a2dda2ef2aab6"
LIGHT_TYPE_GLOBAL_ID = "26DmeZ15D7ZeHYavba3FdK"
EXPECTED_PROXY_IDS = {
    "1faflkXXH6M9cnYPE9Liir",
    "1WMWifzpjFOxatyFb9stFb",
    "3ufPf53o1CbQlrIxF$Oogk",
    "2qvK0LPkPAr8uEI8Oa_U1h",
    "3Irt7GGfb5MP7Q6g$zFVaX",
    "2xvG7uDof7OwbTH_w0edo0",
    "2SdKxc_q96SBGrAUSYsWM9",
    "1SNXhKeZb7r9Y0MmzcQnc3",
    "2d2Vw3ZSn0exH1seMBMiVf",
}
CONTROLLED_SOCKET_IDS = {
    "0laejMoxn8Lu_X3FZaCXmi",
    "27MTenki57DQsfMryX_1U0",
    "2OOjqQDMHDjRcXQCniWXnp",
    "3KXtmVvejA78j_iAS$kydj",
}
NETWORK_TYPE_CATEGORIES = {"LAN", "LANFLUSH", "LANSOCKET"}
CONTROL_TYPE_CATEGORIES = {"SWITCHPANEL", "CONTROLPANEL"}
SCALE_DENOMINATOR = 50.0
SVG_WORLD_OFFSET_MM = 10000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument(
        "--review",
        type=Path,
        default=Path("pipeline/decisions/elec-existing-review.csv"),
    )
    parser.add_argument("--output", type=Path, default=Path("build/elec/elec-existing-candidate.json"))
    parser.add_argument("--source-svg", type=Path)
    parser.add_argument("--e301-svg", type=Path)
    parser.add_argument("--e303-svg", type=Path)
    parser.add_argument("--expected-sha256", default=EXPECTED_SOURCE_SHA256)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


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


def shape_data(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> tuple[np.ndarray, np.ndarray]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3))
    faces = np.asarray(shape.geometry.faces, dtype=int).reshape((-1, 3))
    if not len(vertices):
        raise RuntimeError(f"{product.GlobalId} has no renderable world geometry")
    return vertices, faces


def bbox_mm(vertices_m: np.ndarray) -> dict[str, list[float]]:
    vertices = vertices_m * 1000.0
    minimum = vertices.min(axis=0)
    maximum = vertices.max(axis=0)
    return {
        "min_mm": minimum.tolist(),
        "max_mm": maximum.tolist(),
        "dimensions_mm": (maximum - minimum).tolist(),
        "centre_mm": ((minimum + maximum) / 2.0).tolist(),
    }


def placement_record(product: ifcopenshell.entity_instance) -> dict[str, Any]:
    if product.ObjectPlacement is None:
        raise RuntimeError(f"{product.GlobalId} has no ObjectPlacement")
    origin = np.asarray(get_local_placement(product.ObjectPlacement)[:3, 3], dtype=float)
    residual = np.abs(origin - np.round(origin))
    return {
        "origin_mm": origin.tolist(),
        "maximum_integer_residual_mm": float(residual.max()),
    }


def read_review(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        records = list(csv.DictReader(handle))
    ids = [row["review_id"] for row in records]
    if len(ids) != len(set(ids)):
        raise RuntimeError("ELEC review IDs are not unique")
    return records


def review_index(records: list[dict[str, str]], scope: str) -> dict[str, dict[str, str]]:
    return {row["global_id"]: row for row in records if row["scope"] == scope}


def space_footprints(
    model: ifcopenshell.file,
    settings: ifcopenshell.geom.settings,
) -> list[tuple[ifcopenshell.entity_instance, Any]]:
    records: list[tuple[ifcopenshell.entity_instance, Any]] = []
    for space in model.by_type("IfcSpace"):
        vertices, faces = shape_data(settings, space)
        polygons = []
        for face in faces:
            polygon = Polygon(vertices[face, :2])
            if polygon.area > 1e-10:
                polygons.append(polygon)
        if not polygons:
            raise RuntimeError(f"Space {space.GlobalId} has no plan footprint")
        records.append((space, unary_union(polygons).buffer(1e-8)))
    return records


def candidate_space(
    centre_mm: list[float],
    footprints: list[tuple[ifcopenshell.entity_instance, Any]],
) -> dict[str, Any]:
    point = Point(float(centre_mm[0]) / 1000.0, float(centre_mm[1]) / 1000.0)
    matches = [space for space, footprint in footprints if footprint.covers(point)]
    if len(matches) != 1:
        raise RuntimeError(
            f"electrical point at {centre_mm[:2]} has {len(matches)} candidate Spaces"
        )
    space = matches[0]
    return {
        "global_id": space.GlobalId,
        "name": str(space.Name or ""),
        "long_name": str(space.LongName or ""),
        "basis": "world Body bbox centre contained by current IfcSpace plan footprint",
        "confidence": 0.95,
        "review_required": "yes",
    }


def instance_record(
    product: ifcopenshell.entity_instance,
    settings: ifcopenshell.geom.settings,
    footprints: list[tuple[ifcopenshell.entity_instance, Any]],
) -> dict[str, Any]:
    vertices, _ = shape_data(settings, product)
    bbox = bbox_mm(vertices)
    assigned_type = get_type(product)
    return {
        "global_id": product.GlobalId,
        "ifc_class": product.is_a(),
        "name": str(product.Name or ""),
        "object_type": str(product.ObjectType or "") if hasattr(product, "ObjectType") else "",
        "container": str(getattr(get_container(product), "Name", "") or ""),
        "assigned_type": {
            "global_id": str(getattr(assigned_type, "GlobalId", "") or ""),
            "name": str(getattr(assigned_type, "Name", "") or ""),
            "element_type": str(getattr(assigned_type, "ElementType", "") or ""),
            "predefined_type": str(getattr(assigned_type, "PredefinedType", "") or ""),
        },
        "bbox": bbox,
        "placement": placement_record(product),
        "candidate_space": candidate_space(bbox["centre_mm"], footprints),
        "has_psets": bool(get_psets(product)),
        "formal_ifc_write_allowed": "no",
    }


def plan_sort_key(record: dict[str, Any]) -> tuple[float, float, str]:
    centre = record["bbox"]["centre_mm"]
    return (-round(float(centre[1]), 3), round(float(centre[0]), 3), record["global_id"])


def assign_candidate_ids(records: list[dict[str, Any]], prefix: str) -> None:
    records.sort(key=plan_sort_key)
    for index, record in enumerate(records, 1):
        record["candidate_id"] = f"{prefix}{index:03d}"


def duplicate_centres(records: list[dict[str, Any]], precision_mm: int = 3) -> list[dict[str, Any]]:
    groups: defaultdict[tuple[float, float, float], list[str]] = defaultdict(list)
    for record in records:
        centre = record["bbox"]["centre_mm"]
        key = tuple(round(float(value), precision_mm) for value in centre)
        groups[key].append(record["global_id"])
    return [
        {"centre_mm": list(centre), "global_ids": ids}
        for centre, ids in sorted(groups.items())
        if len(ids) > 1
    ]


def type_record(type_object: ifcopenshell.entity_instance, instance_count: int) -> dict[str, Any]:
    return {
        "global_id": type_object.GlobalId,
        "name": str(type_object.Name or ""),
        "element_type": str(type_object.ElementType or ""),
        "predefined_type": str(type_object.PredefinedType or ""),
        "instance_count": instance_count,
        "is_instance_point": False,
    }


def opening_record(opening: ifcopenshell.entity_instance) -> dict[str, Any]:
    hosts = [relation.RelatingBuildingElement for relation in opening.VoidsElements]
    return {
        "global_id": opening.GlobalId,
        "name": str(opening.Name or ""),
        "hosts": [
            {"ifc_class": host.is_a(), "global_id": host.GlobalId, "name": str(host.Name or "")}
            for host in hosts
        ],
        "formal_ifc_write_allowed": "no",
    }


def place_plan_labels(
    records: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], int]:
    occupied: list[tuple[float, float, float, float]] = []
    labels: list[dict[str, Any]] = []
    offsets = sorted(
        (
            (dx, dy)
            for dx in range(-27, 28, 3)
            for dy in range(-27, 28, 3)
            if abs(dx) + abs(dy) >= 4
        ),
        key=lambda item: (abs(item[0]) + abs(item[1]), abs(item[1]), abs(item[0])),
    )
    collision_count = 0
    for record in records:
        cx, cy = world_to_svg(*record["bbox"]["centre_mm"][:2])
        text = record["candidate_id"]
        width = max(5.8, len(text) * 1.45)
        chosen = None
        for dx, dy in offsets:
            tx = cx + dx
            ty = cy + dy
            box = (tx - width / 2, ty - 2.6, tx + width / 2, ty + 0.8)
            if box[0] < 1 or box[1] < 1 or box[2] > 399 or box[3] > 399:
                continue
            if any(
                not (
                    box[2] + 0.35 <= other[0]
                    or other[2] + 0.35 <= box[0]
                    or box[3] + 0.35 <= other[1]
                    or other[3] + 0.35 <= box[1]
                )
                for other in occupied
            ):
                continue
            chosen = (tx, ty, box)
            break
        if chosen is None:
            collision_count += 1
            chosen = (cx, cy, (cx - width / 2, cy - 2.6, cx + width / 2, cy + 0.8))
        tx, ty, box = chosen
        occupied.append(box)
        labels.append(
            {
                "global_id": record["global_id"],
                "candidate_id": text,
                "marker": [cx, cy],
                "text": [tx, ty],
                "box": list(box),
            }
        )
    return labels, collision_count


def inject_elec_svg(source: str, generated: str, sheet: str, style: str) -> str:
    if 'data-scale="1:50"' not in source:
        raise RuntimeError("Wall Plan SVG is not at expected 1:50 scale")
    source = re.sub(r'width="400(?:\.0+)?mm"', 'width="500mm"', source, count=1)
    source = re.sub(
        r'viewBox="0 0 400(?:\.0+)? 400(?:\.0+)?"',
        'viewBox="0 0 500 400"',
        source,
        count=1,
    )
    if 'width="500mm"' not in source or 'viewBox="0 0 500 400"' not in source:
        raise RuntimeError("failed to expand ELEC candidate sheet to 500x400 mm")
    if "</svg>" not in source:
        raise RuntimeError("invalid SVG: closing tag missing")
    return source.replace(
        "</svg>",
        f'<style id="{sheet.lower()}-candidate-style">{style}</style>'
        f'<g id="{sheet.lower()}-generated">{generated}</g></svg>',
        1,
    )


def panel_line(markup: list[str], text: str, y: float, css: str = "elec-text") -> None:
    markup.append(f'<text class="{css}" x="407" y="{y:.1f}">{svg_escape(text)}</text>')


def render_e301(
    source: str,
    source_sha: str,
    lights: list[dict[str, Any]],
    room_counts: Counter[str],
) -> tuple[str, dict[str, Any]]:
    labels, collisions = place_plan_labels(lights)
    by_id = {record["global_id"]: record for record in lights}
    markup = ['<g id="e301-light-markers">']
    for label in labels:
        record = by_id[label["global_id"]]
        cx, cy = label["marker"]
        tx, ty = label["text"]
        box = label["box"]
        markup.append(
            f'<g data-elec-kind="light" data-ifc-guid="{svg_escape(record["global_id"])}" '
            f'data-candidate-id="{record["candidate_id"]}">'
            f'<circle class="e301-light-ring" cx="{cx:.3f}" cy="{cy:.3f}" r="2.0"/>'
            f'<line class="e301-light-cross" x1="{cx-1.3:.3f}" y1="{cy:.3f}" x2="{cx+1.3:.3f}" y2="{cy:.3f}"/>'
            f'<line class="e301-light-cross" x1="{cx:.3f}" y1="{cy-1.3:.3f}" x2="{cx:.3f}" y2="{cy+1.3:.3f}"/>'
            f'<line class="elec-leader" x1="{cx:.3f}" y1="{cy:.3f}" x2="{tx:.3f}" y2="{ty-1.0:.3f}"/>'
            f'<rect class="elec-label-bg" x="{box[0]:.3f}" y="{box[1]:.3f}" width="{box[2]-box[0]:.3f}" height="{box[3]-box[1]:.3f}"/>'
            f'<text class="elec-label" x="{tx:.3f}" y="{ty:.3f}">{record["candidate_id"]}</text>'
            '</g>'
        )
    markup.append('</g><g id="e301-side-panel"><rect class="elec-panel" x="402" y="7" width="93" height="386"/>')
    panel_line(markup, "E-301 灯具定位候选", 16, "elec-title")
    panel_line(markup, "当前正式 IFC 派生｜1:50｜非施工发布", 23, "elec-note")
    panel_line(markup, "仅表达既有位置与 Space 候选", 30)
    panel_line(markup, "不含灯组、回路、功率或控制推断", 36, "elec-warn")
    panel_line(markup, "灯具统计", 46, "elec-heading")
    panel_line(markup, "RA.LP / DIRECTIONSOURCE：79", 52)
    panel_line(markup, "安装中心 0.1 mm：79/79", 58)
    panel_line(markup, "重复中心：0｜标注碰撞：0", 64)
    panel_line(markup, "Space 几何候选", 75, "elec-heading")
    y = 81.0
    for room, count in sorted(room_counts.items(), key=lambda item: (-item[1], item[0])):
        panel_line(markup, f"{room}：{count}", y)
        y += 4.5
    y += 3
    panel_line(markup, "施工发布停止条件", y, "elec-heading")
    panel_line(markup, "产品、安装方式、灯组与控制未确认", y + 6, "elec-warn")
    panel_line(markup, f"IFC SHA {source_sha[:12]}…", 382, "elec-note")
    markup.append("</g>")
    style = """
@page{size:500mm 400mm;margin:0}html,body{margin:0;width:500mm;height:400mm;overflow:hidden}
.e301-light-ring{fill:#b2f2f2;fill-opacity:.7;stroke:#087f8c;stroke-width:.75}
.e301-light-cross{stroke:#087f8c;stroke-width:.55;stroke-linecap:round}
.elec-leader{stroke:#526777;stroke-width:.25}.elec-label-bg{fill:#fff;fill-opacity:.9;stroke:#8fa3b1;stroke-width:.18}
.elec-label,.elec-title,.elec-heading,.elec-text,.elec-note,.elec-warn{font-family:Arial,'Noto Sans CJK SC',sans-serif;fill:#102f43}
.elec-label{font-size:2.05px;font-weight:700;text-anchor:middle}.elec-panel{fill:#fbfcfd;stroke:#102f43;stroke-width:.5}
.elec-title{font-size:4px;font-weight:700}.elec-heading{font-size:3px;font-weight:700}.elec-text{font-size:2.3px}
.elec-note{font-size:2.15px;fill:#526777}.elec-warn{font-size:2.25px;fill:#c92a2a;font-weight:700}
"""
    gates = {
        "marker_count": len(lights),
        "label_count": len(labels),
        "label_collisions": collisions,
        "candidate_ids_unique": len({record["candidate_id"] for record in lights}) == len(lights),
        "markers_in_plan_bounds": all(
            0 <= label["marker"][0] <= 400 and 0 <= label["marker"][1] <= 400
            for label in labels
        ),
    }
    gates["pass"] = all(
        [
            gates["marker_count"] == 79,
            gates["label_count"] == 79,
            gates["label_collisions"] == 0,
            gates["candidate_ids_unique"],
            gates["markers_in_plan_bounds"],
        ]
    )
    return inject_elec_svg(source, "".join(markup), "E301", style), gates


def render_e303(
    source: str,
    source_sha: str,
    sockets: list[dict[str, Any]],
    equipment: list[dict[str, Any]],
    proxies: list[dict[str, Any]],
) -> tuple[str, dict[str, Any]]:
    records = [*sockets, *equipment, *proxies]
    labels, collisions = place_plan_labels(records)
    by_id = {record["global_id"]: record for record in records}
    socket_ids = {record["global_id"] for record in sockets}
    equipment_ids = {record["global_id"] for record in equipment}
    markup = ['<g id="e303-location-markers">']
    for label in labels:
        record = by_id[label["global_id"]]
        cx, cy = label["marker"]
        tx, ty = label["text"]
        box = label["box"]
        if record["global_id"] in socket_ids:
            kind = "socket-exception" if record["controlled_installation_exception"] else "socket"
            symbol = (
                f'<rect class="e303-{kind}" x="{cx-1.7:.3f}" y="{cy-1.7:.3f}" width="3.4" height="3.4" '
                f'transform="rotate(45 {cx:.3f} {cy:.3f})"/>'
            )
        else:
            kind = "equipment" if record["global_id"] in equipment_ids else "proxy"
            bbox = record["bbox"]
            x, y = world_to_svg(bbox["min_mm"][0], bbox["max_mm"][1])
            width = max(1.2, bbox["dimensions_mm"][0] / SCALE_DENOMINATOR)
            height = max(1.2, bbox["dimensions_mm"][1] / SCALE_DENOMINATOR)
            symbol = (
                f'<rect class="e303-{kind}" x="{x:.3f}" y="{y:.3f}" width="{width:.3f}" height="{height:.3f}"/>'
            )
        markup.append(
            f'<g data-elec-kind="{kind}" data-ifc-guid="{svg_escape(record["global_id"])}" '
            f'data-candidate-id="{record["candidate_id"]}">{symbol}'
            f'<line class="elec-leader" x1="{cx:.3f}" y1="{cy:.3f}" x2="{tx:.3f}" y2="{ty-1.0:.3f}"/>'
            f'<rect class="elec-label-bg" x="{box[0]:.3f}" y="{box[1]:.3f}" width="{box[2]-box[0]:.3f}" height="{box[3]-box[1]:.3f}"/>'
            f'<text class="elec-label" x="{tx:.3f}" y="{ty:.3f}">{record["candidate_id"]}</text>'
            '</g>'
        )
    markup.append('</g><g id="e303-side-panel"><rect class="elec-panel" x="402" y="7" width="93" height="386"/>')
    panel_line(markup, "E-303 插座与设备定位候选", 16, "elec-title")
    panel_line(markup, "当前正式 IFC 派生｜1:50｜非施工发布", 23, "elec-note")
    panel_line(markup, "仅表达既有点位、轮廓与待确认项", 30)
    panel_line(markup, "不含回路、功率、防水或接口推断", 36, "elec-warn")
    panel_line(markup, "既有实例", 47, "elec-heading")
    panel_line(markup, "插座 11：SOC01×2 / SOC04×6 / SOCF04×3", 53)
    panel_line(markup, "类型化设备 8：空调 5 / 厨房设备 3", 59)
    panel_line(markup, "ELEC Proxy 9：身份或锚点待确认", 65)
    panel_line(markup, "定位 QA", 76, "elec-heading")
    panel_line(markup, "插座 0.1 mm：7/11｜受控例外 4", 82)
    panel_line(markup, "设备整数原点：8/8", 88)
    panel_line(markup, "Proxy 超 0.1 mm：9/9（不移动）", 94)
    panel_line(markup, "重复中心：0｜标注碰撞：0", 100)
    controlled_ids = [record["candidate_id"] for record in sockets if record["controlled_installation_exception"]]
    panel_line(markup, "受控插座例外", 111, "elec-heading")
    panel_line(markup, " / ".join(controlled_ids), 117, "elec-warn")
    panel_line(markup, "设备类型", 128, "elec-heading")
    panel_line(markup, "AC1180×1 / AC700×2 / AC700F×2", 134)
    panel_line(markup, "HD01 / WD01 / OV01 各 1", 140)
    panel_line(markup, "施工发布停止条件", 151, "elec-heading")
    panel_line(markup, "设备参数、专用回路与 Proxy 身份未确认", 157, "elec-warn")
    panel_line(markup, f"IFC SHA {source_sha[:12]}…", 382, "elec-note")
    markup.append("</g>")
    style = """
@page{size:500mm 400mm;margin:0}html,body{margin:0;width:500mm;height:400mm;overflow:hidden}
.e303-socket{fill:#d0ebff;stroke:#1864ab;stroke-width:.65}.e303-socket-exception{fill:#ffc9c9;stroke:#c92a2a;stroke-width:.8}
.e303-equipment{fill:#d0bfff44;stroke:#7048e8;stroke-width:.55}.e303-proxy{fill:#ffd8a844;stroke:#e67700;stroke-width:.65;stroke-dasharray:2 1}
.elec-leader{stroke:#526777;stroke-width:.25}.elec-label-bg{fill:#fff;fill-opacity:.9;stroke:#8fa3b1;stroke-width:.18}
.elec-label,.elec-title,.elec-heading,.elec-text,.elec-note,.elec-warn{font-family:Arial,'Noto Sans CJK SC',sans-serif;fill:#102f43}
.elec-label{font-size:2.05px;font-weight:700;text-anchor:middle}.elec-panel{fill:#fbfcfd;stroke:#102f43;stroke-width:.5}
.elec-title{font-size:4px;font-weight:700}.elec-heading{font-size:3px;font-weight:700}.elec-text{font-size:2.3px}
.elec-note{font-size:2.15px;fill:#526777}.elec-warn{font-size:2.25px;fill:#c92a2a;font-weight:700}
"""
    gates = {
        "socket_markers": len(sockets),
        "typed_equipment_markers": len(equipment),
        "proxy_markers": len(proxies),
        "label_count": len(labels),
        "label_collisions": collisions,
        "candidate_ids_unique": len({record["candidate_id"] for record in records}) == len(records),
        "markers_in_plan_bounds": all(
            0 <= label["marker"][0] <= 400 and 0 <= label["marker"][1] <= 400
            for label in labels
        ),
    }
    gates["pass"] = all(
        [
            gates["socket_markers"] == 11,
            gates["typed_equipment_markers"] == 8,
            gates["proxy_markers"] == 9,
            gates["label_count"] == 28,
            gates["label_collisions"] == 0,
            gates["candidate_ids_unique"],
            gates["markers_in_plan_bounds"],
        ]
    )
    return inject_elec_svg(source, "".join(markup), "E303", style), gates


def main() -> int:
    args = parse_args()
    drawing_paths = [args.source_svg, args.e301_svg, args.e303_svg]
    if any(drawing_paths) and not all(drawing_paths):
        raise RuntimeError("--source-svg, --e301-svg, and --e303-svg must be supplied together")
    source_sha = sha256(args.input)
    if source_sha != args.expected_sha256:
        raise RuntimeError(
            f"formal IFC SHA drift: {source_sha} != {args.expected_sha256}"
        )

    model = ifcopenshell.open(args.input)
    if model.schema != "IFC4":
        raise RuntimeError(f"unexpected schema {model.schema}")
    reviews = read_review(args.review)
    proxy_reviews = review_index(reviews, "proxy_handoff")
    socket_exception_reviews = review_index(reviews, "socket_exception")
    if set(proxy_reviews) != EXPECTED_PROXY_IDS:
        raise RuntimeError("proxy handoff review set drift")
    if set(socket_exception_reviews) != CONTROLLED_SOCKET_IDS:
        raise RuntimeError("controlled socket review set drift")

    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    footprints = space_footprints(model, settings)

    lights = [instance_record(product, settings, footprints) for product in model.by_type("IfcLightFixture")]
    for record in lights:
        record["review_basis"] = "existing RA.LP world installation centre; room is a geometric candidate"
        record["review_required"] = "yes"
    assign_candidate_ids(lights, "L")

    appliances = model.by_type("IfcElectricAppliance")
    sockets = [
        instance_record(product, settings, footprints)
        for product in appliances
        if str(getattr(get_type(product), "ElementType", "") or "") == "SOCKET"
    ]
    for record in sockets:
        record["controlled_installation_exception"] = record["global_id"] in CONTROLLED_SOCKET_IDS
        record["review_id"] = socket_exception_reviews.get(record["global_id"], {}).get("review_id", "")
    assign_candidate_ids(sockets, "P")

    equipment = [
        instance_record(product, settings, footprints)
        for product in appliances
        if str(getattr(get_type(product), "ElementType", "") or "") != "SOCKET"
    ]
    for record in equipment:
        record["power_and_interface_status"] = "missing"
        record["review_required"] = "yes"
    assign_candidate_ids(equipment, "EQ")

    proxies = []
    for global_id in EXPECTED_PROXY_IDS:
        product = model.by_guid(global_id)
        if product is None or not product.is_a("IfcBuildingElementProxy"):
            raise RuntimeError(f"ELEC proxy handoff missing or class drift: {global_id}")
        record = instance_record(product, settings, footprints)
        review = proxy_reviews[global_id]
        record.update(
            {
                "candidate_role": review["candidate_role"],
                "review_id": review["review_id"],
                "review_required": review["review_required"],
                "review_status": review["status"],
                "protected_action": review["protected_action"],
                "required_confirmation": review["required_confirmation"],
                "related_openings": [
                    relation.RelatedOpeningElement.GlobalId for relation in product.HasOpenings
                ],
            }
        )
        proxies.append(record)
    assign_candidate_ids(proxies, "PX")

    openings = []
    for product in proxies:
        for global_id in product["related_openings"]:
            openings.append(opening_record(model.by_guid(global_id)))
    openings.sort(key=lambda record: record["global_id"])

    type_use_counts = Counter(
        get_type(product).GlobalId for product in appliances if get_type(product) is not None
    )
    appliance_types = model.by_type("IfcElectricApplianceType")
    used_types = [type_record(item, type_use_counts[item.GlobalId]) for item in appliance_types if type_use_counts[item.GlobalId]]
    unused_types = [type_record(item, 0) for item in appliance_types if not type_use_counts[item.GlobalId]]
    used_types.sort(key=lambda record: (record["element_type"], record["name"], record["global_id"]))
    unused_types.sort(key=lambda record: (record["element_type"], record["name"], record["global_id"]))
    control_types = [record for record in unused_types if record["element_type"] in CONTROL_TYPE_CATEGORIES]
    network_types = [record for record in unused_types if record["element_type"] in NETWORK_TYPE_CATEGORIES]

    root_ids = [root.GlobalId for root in model.by_type("IfcRoot")]
    light_type_ids = {record["assigned_type"]["global_id"] for record in lights}
    socket_ids = {record["global_id"] for record in sockets}
    proxy_ids = {record["global_id"] for record in proxies}
    light_room_counts = Counter(record["candidate_space"]["long_name"] for record in lights)
    light_residuals = [record["placement"]["maximum_integer_residual_mm"] for record in lights]
    socket_residuals = [record["placement"]["maximum_integer_residual_mm"] for record in sockets]
    proxy_residuals = [record["placement"]["maximum_integer_residual_mm"] for record in proxies]

    topology = {
        "ifc_systems": len(model.by_type("IfcSystem")),
        "ifc_distribution_systems": len(model.by_type("IfcDistributionSystem")),
        "distribution_ports": len(model.by_type("IfcDistributionPort")),
        "port_connections": len(model.by_type("IfcRelConnectsPorts")),
        "port_to_element_connections": len(model.by_type("IfcRelConnectsPortToElement")),
        "control_assignments": len(model.by_type("IfcRelAssignsToControl")),
    }
    e302_instances = len(model.by_type("IfcSwitchingDevice"))
    e304_instances = len(model.by_type("IfcCommunicationsAppliance"))

    gates = {
        "root_global_ids_unique": len(root_ids) == len(set(root_ids)),
        "light_count": len(lights),
        "light_type_ids": sorted(light_type_ids),
        "light_geometry_and_type_complete": all(record["assigned_type"]["global_id"] for record in lights),
        "light_space_candidates_single": all(record["candidate_space"]["global_id"] for record in lights),
        "light_duplicate_centres": duplicate_centres(lights),
        "light_origins_within_0_1_mm": sum(value <= args.tolerance_mm + 1e-9 for value in light_residuals),
        "light_maximum_integer_residual_mm": max(light_residuals),
        "socket_count": len(sockets),
        "socket_duplicate_centres": duplicate_centres(sockets),
        "socket_origins_within_0_1_mm": sum(value <= args.tolerance_mm + 1e-9 for value in socket_residuals),
        "controlled_socket_ids": sorted(record["global_id"] for record in sockets if record["controlled_installation_exception"]),
        "typed_equipment_count": len(equipment),
        "typed_equipment_integer_origins": sum(
            record["placement"]["maximum_integer_residual_mm"] <= 1e-9 for record in equipment
        ),
        "proxy_handoff_count": len(proxies),
        "proxy_handoff_ids": sorted(proxy_ids),
        "proxy_handoff_over_0_1_mm": sum(value > args.tolerance_mm + 1e-9 for value in proxy_residuals),
        "proxy_handoff_maximum_integer_residual_mm": max(proxy_residuals),
        "related_opening_count": len(openings),
        "related_openings_have_one_host": all(len(record["hosts"]) == 1 for record in openings),
        "appliance_type_count": len(appliance_types),
        "used_appliance_type_count": len(used_types),
        "unused_appliance_type_count": len(unused_types),
        "type_library_separated_from_instances": not ({record["global_id"] for record in used_types} & {record["global_id"] for record in unused_types}),
        "e302_switch_instances": e302_instances,
        "e302_unused_control_type_definitions": len(control_types),
        "e304_network_instances": e304_instances,
        "e304_unused_network_type_definitions": len(network_types),
        "topology": topology,
        "automatic_ifc_write_allowed": False,
    }
    gates["candidate_pass"] = all(
        [
            gates["root_global_ids_unique"],
            gates["light_count"] == 79,
            gates["light_type_ids"] == [LIGHT_TYPE_GLOBAL_ID],
            gates["light_geometry_and_type_complete"],
            gates["light_space_candidates_single"],
            gates["light_duplicate_centres"] == [],
            gates["light_origins_within_0_1_mm"] == 79,
            gates["socket_count"] == 11,
            gates["socket_duplicate_centres"] == [],
            gates["socket_origins_within_0_1_mm"] == 7,
            set(gates["controlled_socket_ids"]) == CONTROLLED_SOCKET_IDS,
            gates["typed_equipment_count"] == 8,
            gates["typed_equipment_integer_origins"] == 8,
            gates["proxy_handoff_count"] == 9,
            set(gates["proxy_handoff_ids"]) == EXPECTED_PROXY_IDS,
            gates["proxy_handoff_over_0_1_mm"] == 9,
            gates["related_opening_count"] == 4,
            gates["related_openings_have_one_host"],
            gates["appliance_type_count"] == 42,
            gates["used_appliance_type_count"] == 9,
            gates["unused_appliance_type_count"] == 33,
            gates["type_library_separated_from_instances"],
            gates["e302_switch_instances"] == 0,
            gates["e302_unused_control_type_definitions"] == 13,
            gates["e304_network_instances"] == 0,
            gates["e304_unused_network_type_definitions"] == 7,
            all(value == 0 for value in topology.values()),
        ]
    )
    gates["construction_release_ready"] = False
    gates["release_blocks"] = [
        "E-301 fixture product, mounting method, circuit, and control group are not confirmed",
        "E-302 has no switch instances or control relations",
        "E-303 has no equipment power, circuit, waterproofing, or interface properties",
        "E-304 has no network instances, ports, systems, or topology",
        "nine ELEC proxy handoff objects have no confirmed installation anchors",
    ]

    drawing_source = None
    drawing_gates: dict[str, Any] = {}
    if args.source_svg is not None:
        source_svg = args.source_svg.read_text(encoding="utf-8")
        e301_svg, e301_gates = render_e301(source_svg, source_sha, lights, light_room_counts)
        e303_svg, e303_gates = render_e303(source_svg, source_sha, sockets, equipment, proxies)
        drawing_gates = {"E-301": e301_gates, "E-303": e303_gates}
        gates["drawing_candidate_pass"] = e301_gates["pass"] and e303_gates["pass"]
        gates["candidate_pass"] = gates["candidate_pass"] and gates["drawing_candidate_pass"]
        if not gates["candidate_pass"]:
            raise RuntimeError("ELEC drawing candidate gates failed")
        args.e301_svg.parent.mkdir(parents=True, exist_ok=True)
        args.e303_svg.parent.mkdir(parents=True, exist_ok=True)
        args.e301_svg.write_text(e301_svg, encoding="utf-8")
        args.e303_svg.write_text(e303_svg, encoding="utf-8")
        drawing_source = {
            "path": str(args.source_svg.resolve()),
            "sha256": sha256(args.source_svg),
            "scale": "1:50",
        }

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_existing_condition_candidate",
        "source": {"path": str(args.input.resolve()), "sha256": source_sha, "schema": model.schema},
        "tolerance_mm": args.tolerance_mm,
        "sheets": {
            "E-301": {
                "status": "existing_points_candidate",
                "lights": lights,
                "room_counts": dict(sorted(light_room_counts.items())),
                "missing": ["fixture product", "mounting method", "circuit", "control group"],
            },
            "E-302": {
                "status": "pending_layout_no_instances",
                "instances": [],
                "available_uninstantiated_type_definitions": control_types,
                "confirmed_requirements_only": [
                    "entry lighting master control",
                    "entry-to-master-bedroom corridor/living dual control",
                    "text or icon label for every key",
                ],
            },
            "E-303": {
                "status": "existing_points_and_footprints_candidate",
                "sockets": sockets,
                "typed_equipment": equipment,
                "proxy_handoffs": proxies,
                "related_openings": openings,
                "missing": ["power", "voltage", "circuit", "waterproofing", "interface"],
            },
            "E-304": {
                "status": "pending_layout_no_instances_or_topology",
                "instances": [],
                "available_uninstantiated_type_definitions": network_types,
                "topology": topology,
            },
        },
        "type_library": {"used": used_types, "unused": unused_types},
        "review_register": {"path": str(args.review), "records": len(reviews)},
        "drawing_source": drawing_source,
        "drawing_gates": drawing_gates,
        "gates": gates,
    }
    if not gates["candidate_pass"]:
        raise RuntimeError("ELEC existing-condition candidate gates failed")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "gates": gates}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
