#!/usr/bin/env python3
"""Generate the read-only A-105 floor finish, slope, and threshold candidate."""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element
import numpy as np

from a103_wall_plan_candidate import Box, format_mm, sha256, world_to_svg


A105_SLAB_GLOBAL_ID = "3ARl_CqPrCWQrA6_$W073W"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--source-svg", required=True, type=Path)
    parser.add_argument("--surface-audit", required=True, type=Path)
    parser.add_argument("--a104-register", required=True, type=Path)
    parser.add_argument("--output-svg", required=True, type=Path)
    parser.add_argument("--floor-register", required=True, type=Path)
    parser.add_argument("--threshold-register", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--plane-residual-mm", type=float, default=0.001)
    return parser.parse_args()


def svg_escape(value: Any) -> str:
    return html.escape(str(value), quote=True)


def shape_data(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> tuple[np.ndarray, np.ndarray]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3)) * 1000.0
    faces = np.asarray(shape.geometry.faces, dtype=int).reshape((-1, 3))
    if len(vertices) == 0:
        raise RuntimeError(f"{product.GlobalId} has no world vertices")
    return vertices, faces


def bbox_record(vertices: np.ndarray) -> dict[str, list[float]]:
    minimum = vertices.min(axis=0)
    maximum = vertices.max(axis=0)
    return {
        "min_mm": minimum.tolist(),
        "max_mm": maximum.tolist(),
        "dimensions_mm": (maximum - minimum).tolist(),
        "centre_mm": ((minimum + maximum) / 2.0).tolist(),
    }


def top_plane(vertices: np.ndarray, faces: np.ndarray) -> dict[str, Any] | None:
    indices: list[int] = []
    for face in faces:
        points = vertices[face]
        normal = np.cross(points[1] - points[0], points[2] - points[0])
        length = float(np.linalg.norm(normal))
        if length and float(normal[2] / length) > 0.7:
            indices.extend(int(index) for index in face)
    if not indices:
        return None
    points = np.unique(np.round(vertices[indices], 9), axis=0)
    matrix = np.column_stack((points[:, 0], points[:, 1], np.ones(len(points))))
    coefficients, *_ = np.linalg.lstsq(matrix, points[:, 2], rcond=None)
    fitted = matrix @ coefficients
    residual = float(np.max(np.abs(fitted - points[:, 2])))
    a, b, c = (float(value) for value in coefficients)
    slope_percent = math.hypot(a, b) * 100.0
    downhill = np.array([-a, -b], dtype=float)
    direction = "LEVEL"
    if float(np.linalg.norm(downhill)) > 1e-12:
        azimuth = math.degrees(math.atan2(float(downhill[0]), float(downhill[1]))) % 360.0
        labels = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
        direction = labels[int((azimuth + 22.5) // 45.0) % 8]
    return {
        "point_count": len(points),
        "a_dz_dx": a,
        "b_dz_dy": b,
        "c_mm": c,
        "slope_percent": slope_percent,
        "downhill_direction": direction,
        "top_min_z_mm": float(points[:, 2].min()),
        "top_max_z_mm": float(points[:, 2].max()),
        "maximum_fit_residual_mm": residual,
    }


def material_name(product: ifcopenshell.entity_instance) -> str:
    material = ifcopenshell.util.element.get_material(product)
    if material is None:
        return ""
    if material.is_a("IfcMaterial"):
        return str(material.Name or "")
    return str(material)


def plan_sort_key(record: dict[str, Any]) -> tuple[float, float, str]:
    centre = record["bbox"]["centre_mm"]
    return (-round(float(centre[1]), 3), round(float(centre[0]), 3), record["global_id"])


def read_surface_audit(path: Path, source_sha: str) -> dict[str, dict[str, Any]]:
    audit = json.loads(path.read_text(encoding="utf-8"))
    if audit.get("source", {}).get("sha256") != source_sha:
        raise RuntimeError("construction-surface audit does not match the formal IFC")
    return {record["global_id"]: record for record in audit.get("records", [])}


def floor_inventory(
    model: ifcopenshell.file,
    surface_records: dict[str, dict[str, Any]],
    plane_residual_mm: float,
) -> list[dict[str, Any]]:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    records: list[dict[str, Any]] = []
    for product in model.by_type("IfcCovering"):
        if str(product.PredefinedType or "") != "FLOORING":
            continue
        vertices, faces = shape_data(settings, product)
        bbox = bbox_record(vertices)
        material = material_name(product)
        plane = top_plane(vertices, faces)
        audit = surface_records.get(product.GlobalId, {})
        if material == "TerrazzoMosaicTile":
            kind = "sloped_wet_tile"
            review_group = "A105-R02"
            review_question = "确认客卫/主卫 1.047% 找坡方向、地漏关系和铺贴分区是否符合设计与安装要求。"
            if plane is None or plane["maximum_fit_residual_mm"] > plane_residual_mm:
                raise RuntimeError(f"wet tile {product.GlobalId} has no reliable planar top face")
        elif not material and bbox["dimensions_mm"][2] <= plane_residual_mm:
            kind = "finish_reference_plane_material_pending"
            review_group = "A105-R01"
            review_question = "确定该干区地坪的材料、完成面构造总厚度及房间分界；当前仅为零厚度 FFL 参考面。"
        else:
            kind = "unclassified_flooring"
            review_group = "A105-R05"
            review_question = "现有 IFC 证据不足以确认该地坪材料和构造。"
        records.append(
            {
                "global_id": product.GlobalId,
                "name": str(product.Name or ""),
                "current_tag": str(product.Tag or ""),
                "material": material,
                "kind": kind,
                "container": str(audit.get("container") or ""),
                "primary_space": str(audit.get("primary_space_long_name") or ""),
                "primary_space_overlap_ratio": audit.get("primary_space_overlap_ratio"),
                "space_overlaps": audit.get("space_overlaps", []),
                "bbox": bbox,
                "top_plane": plane,
                "review_group": review_group,
                "review_required": "yes",
                "review_question": review_question,
                "basis": "formal IFC Body geometry + material association + current Space overlap audit",
                "confidence": 1.0 if kind != "unclassified_flooring" else 0.5,
                "formal_ifc_write_allowed": "no",
            }
        )
    records.sort(key=plan_sort_key)
    for index, record in enumerate(records, 1):
        record["candidate_id"] = f"F{index:02d}"
    return records


def point_in_bbox(x: float, y: float, bbox: dict[str, list[float]], tolerance: float = 1.0) -> bool:
    return (
        bbox["min_mm"][0] - tolerance <= x <= bbox["max_mm"][0] + tolerance
        and bbox["min_mm"][1] - tolerance <= y <= bbox["max_mm"][1] + tolerance
    )


def plane_z_at(record: dict[str, Any], x: float, y: float) -> float | None:
    plane = record.get("top_plane")
    if plane is None:
        return None
    return plane["a_dz_dx"] * x + plane["b_dz_dy"] * y + plane["c_mm"]


def threshold_inventory(
    path: Path,
    floors: Sequence[dict[str, Any]],
    source_sha: str,
) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["ifc_class"] == "IfcDoor"]
    if len(rows) != 8 or any(row["source_ifc_sha256"] != source_sha for row in rows):
        raise RuntimeError("A-104 door register does not match the formal IFC")
    records = []
    for row in rows:
        x = float(row["center_x_mm"])
        y = float(row["center_y_mm"])
        door_z = float(row["sill_or_threshold_z_mm"])
        candidates = []
        for floor in floors:
            if point_in_bbox(x, y, floor["bbox"]):
                floor_z = plane_z_at(floor, x, y)
                candidates.append(
                    {
                        "candidate_id": floor["candidate_id"],
                        "global_id": floor["global_id"],
                        "material": floor["material"],
                        "floor_z_mm": floor_z,
                        "door_to_floor_delta_mm": door_z - floor_z if floor_z is not None else None,
                    }
                )
        records.append(
            {
                "door_candidate_id": row["candidate_id"],
                "door_global_id": row["global_id"],
                "host_relation": row["host_relation"],
                "space_candidates": row["space_candidates"],
                "center_x_mm": x,
                "center_y_mm": y,
                "door_geometry_min_z_mm": door_z,
                "floor_candidates": candidates,
                "a104_review_group": row["review_group"],
                "review_group": "A105-R03",
                "review_required": "yes",
                "review_question": "确定门槛两侧材料、完成面高差、防水收口和门下净空；无宿主门须先完成 A-104 身份确认。",
                "basis": "A-104 door world geometry + A-105 floor top planes at the door plan centre",
                "confidence": 1.0 if row["host_relation"] == "HOSTED" else 0.5,
                "formal_ifc_write_allowed": "no",
            }
        )
    return records


def slab_inventory(
    model: ifcopenshell.file,
    settings: ifcopenshell.geom.settings,
) -> dict[str, Any]:
    slab = model.by_guid(A105_SLAB_GLOBAL_ID)
    if slab is None or not slab.is_a("IfcSlab"):
        raise RuntimeError(f"missing delegated A-105 slab {A105_SLAB_GLOBAL_ID}")
    vertices, _ = shape_data(settings, slab)
    return {
        "global_id": slab.GlobalId,
        "name": str(slab.Name or ""),
        "predefined_type": str(slab.PredefinedType or ""),
        "material": material_name(slab),
        "bbox": bbox_record(vertices),
        "review_group": "A105-R04",
        "review_required": "yes",
        "review_question": "确认该 1800×2200×550 mm Aircrete 板是否为湿区降板/构造层，以及应保留的施工语义。",
        "formal_ifc_write_allowed": "no",
    }


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in rows:
            row = dict(record)
            for key, value in list(row.items()):
                if isinstance(value, (list, dict)):
                    row[key] = json.dumps(value, ensure_ascii=False, separators=(",", ":"))
            writer.writerow({field: row.get(field, "") for field in fields})


def render_overlay(floors: Sequence[dict[str, Any]], source_sha: str) -> tuple[str, int]:
    pieces = ['<g id="a105-floor-overlay">']
    occupied: list[Box] = []
    collision_count = 0
    for record in floors:
        bbox = record["bbox"]
        x, y = world_to_svg(bbox["min_mm"][0], bbox["max_mm"][1])
        width = max(0.5, bbox["dimensions_mm"][0] / 50.0)
        height = max(0.5, bbox["dimensions_mm"][1] / 50.0)
        css = "a105-wet" if record["kind"] == "sloped_wet_tile" else "a105-material-pending"
        pieces.append(
            f'<rect class="a105-floor {css}" data-ifc-guid="{svg_escape(record["global_id"])}" '
            f'x="{x:.3f}" y="{y:.3f}" width="{width:.3f}" height="{height:.3f}"/>'
        )
        cx, cy = world_to_svg(bbox["centre_mm"][0], bbox["centre_mm"][1])
        label = record["candidate_id"]
        box = Box(cx - 3.2, cy - 2.2, cx + 3.2, cy + 1.2, label)
        if any(box.intersects(other) for other in occupied):
            collision_count += 1
            for offset in (4.5, -4.5, 9.0, -9.0):
                candidate = Box(box.x0, box.y0 + offset, box.x1, box.y1 + offset, label)
                if not any(candidate.intersects(other) for other in occupied):
                    box = candidate
                    collision_count -= 1
                    break
        occupied.append(box)
        pieces.append(
            f'<rect class="a105-label-bg" x="{box.x0:.3f}" y="{box.y0:.3f}" width="{box.x1-box.x0:.3f}" height="{box.y1-box.y0:.3f}"/>'
            f'<text class="a105-label" x="{(box.x0+box.x1)/2:.3f}" y="{box.y1-0.7:.3f}">{label}</text>'
        )
    pieces.append("</g>")
    panel = [
        '<g id="a105-side-panel"><rect class="a105-panel" x="402" y="7" width="93" height="386"/>',
        '<text class="a105-title" x="407" y="16">A-105 地坪完成面候选</text>',
        '<text class="a105-note" x="407" y="23">只读｜未写 IFC｜1:50</text>',
        f'<text class="a105-text" x="407" y="32">地坪 {len(floors)}：湿区砖 18／材料待定 3</text>',
        '<text class="a105-text" x="407" y="38">湿区坡度：18/18 = 1.047%</text>',
        '<text class="a105-text" x="407" y="44">地砖分缝 2.0 mm</text>',
        '<text class="a105-text" x="407" y="50">线性地漏留缝 14.6 mm</text>',
        '<text class="a105-heading" x="407" y="60">集中审核</text>',
        '<text class="a105-text" x="407" y="67">R01 干区材料与构造厚度</text>',
        '<text class="a105-text" x="407" y="73">R02 两个卫生间找坡与地漏</text>',
        '<text class="a105-text" x="407" y="79">R03 8 樘门门槛/高差/防水</text>',
        '<text class="a105-text" x="407" y="85">R04 Aircrete 板施工语义</text>',
        f'<text class="a105-note" x="407" y="96">IFC SHA {source_sha[:12]}…</text>',
        '</g>',
    ]
    return "".join(pieces + panel), collision_count


def inject_svg(source: str, generated: str) -> str:
    if 'data-scale="1:50"' not in source:
        raise RuntimeError("FFL PLAN SVG is not at expected 1:50 scale")
    source = re.sub(r'width="400(?:\.0+)?mm"', 'width="500mm"', source, count=1)
    source = re.sub(r'viewBox="0 0 400(?:\.0+)? 400(?:\.0+)?"', 'viewBox="0 0 500 400"', source, count=1)
    style = """
<style id="a105-candidate-style">
  @page { size:500mm 400mm; margin:0; }
  html,body { margin:0; width:500mm; height:400mm; overflow:hidden; }
  .a105-floor { stroke-width:0.45; }
  .a105-wet { fill:#2b8aaf33; stroke:#0b7285; }
  .a105-material-pending { fill:#f59f0030; stroke:#e67700; stroke-dasharray:2 1; }
  .a105-label-bg { fill:#fff; stroke:#526777; stroke-width:0.18; }
  .a105-label { font-family:Arial,'Noto Sans CJK SC',sans-serif; fill:#102f43; text-anchor:middle; font-size:2.4px; font-weight:700; }
  .a105-panel { fill:#fbfcfd; stroke:#102f43; stroke-width:0.5; }
  .a105-title,.a105-heading,.a105-text,.a105-note { font-family:Arial,'Noto Sans CJK SC',sans-serif; fill:#102f43; }
  .a105-title { font-size:4px; font-weight:700; }
  .a105-heading { font-size:3px; font-weight:700; }
  .a105-text { font-size:2.35px; }
  .a105-note { font-size:2.15px; fill:#526777; }
</style>
"""
    if "</svg>" not in source:
        raise RuntimeError("invalid SVG: closing tag missing")
    return source.replace("</svg>", f'{style}<g id="a105-generated">{generated}</g></svg>', 1)


def main() -> None:
    args = parse_args()
    source_sha = sha256(args.input)
    model = ifcopenshell.open(args.input)
    if model.schema != "IFC4":
        raise RuntimeError(f"expected IFC4, found {model.schema}")
    surface_records = read_surface_audit(args.surface_audit, source_sha)
    floors = floor_inventory(model, surface_records, args.plane_residual_mm)
    thresholds = threshold_inventory(args.a104_register, floors, source_sha)
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    slab = slab_inventory(model, settings)

    wet = [record for record in floors if record["kind"] == "sloped_wet_tile"]
    material_pending = [record for record in floors if record["kind"] == "finish_reference_plane_material_pending"]
    generated, collision_count = render_overlay(floors, source_sha)
    output = inject_svg(args.source_svg.read_text(encoding="utf-8"), generated)
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(output, encoding="utf-8")

    write_csv(
        args.floor_register,
        floors,
        [
            "candidate_id", "global_id", "name", "current_tag", "material", "kind", "container",
            "primary_space", "primary_space_overlap_ratio", "space_overlaps", "bbox", "top_plane",
            "review_group", "review_required", "review_question", "basis", "confidence",
            "formal_ifc_write_allowed",
        ],
    )
    write_csv(
        args.threshold_register,
        thresholds,
        [
            "door_candidate_id", "door_global_id", "host_relation", "space_candidates", "center_x_mm",
            "center_y_mm", "door_geometry_min_z_mm", "floor_candidates", "a104_review_group", "review_group",
            "review_required", "review_question", "basis", "confidence", "formal_ifc_write_allowed",
        ],
    )

    slopes = [record["top_plane"]["slope_percent"] for record in wet]
    maximum_residual = max(record["top_plane"]["maximum_fit_residual_mm"] for record in wet)
    gates = {
        "ifc_schema_is_ifc4": model.schema == "IFC4",
        "flooring_count": len(floors),
        "wet_tile_count": len(wet),
        "material_pending_count": len(material_pending),
        "wet_tile_planar_top_count": sum(record["top_plane"] is not None for record in wet),
        "wet_tile_maximum_plane_residual_mm": maximum_residual,
        "wet_tile_slope_min_percent": min(slopes),
        "wet_tile_slope_max_percent": max(slopes),
        "door_threshold_count": len(thresholds),
        "candidate_identifier_duplicate_count": len(floors) - len({record["candidate_id"] for record in floors}),
        "generated_label_collision_count": collision_count,
        "automatic_ifc_write_allowed": False,
    }
    gates["mechanical_pass"] = (
        gates["flooring_count"] == 21
        and gates["wet_tile_count"] == 18
        and gates["material_pending_count"] == 3
        and gates["wet_tile_planar_top_count"] == 18
        and gates["wet_tile_maximum_plane_residual_mm"] <= args.plane_residual_mm
        and gates["door_threshold_count"] == 8
        and gates["candidate_identifier_duplicate_count"] == 0
        and gates["generated_label_collision_count"] == 0
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a105-floor-finish-candidate",
        "source": {
            "ifc": str(args.input.resolve()),
            "ifc_sha256": source_sha,
            "schema": model.schema,
            "ffl_plan_svg": str(args.source_svg.resolve()),
            "ffl_plan_svg_sha256": sha256(args.source_svg),
            "surface_audit": str(args.surface_audit.resolve()),
            "a104_register": str(args.a104_register.resolve()),
        },
        "tolerance_mm": args.tolerance_mm,
        "plane_residual_mm": args.plane_residual_mm,
        "confirmed_installation_gaps": {"tile_joint_mm": 2.0, "linear_drain_gap_mm": 14.6},
        "floors": floors,
        "thresholds": thresholds,
        "delegated_slab": slab,
        "review_groups": {
            "A105-R01": "3 个零厚度干区 FFL 参考面：材料、构造厚度和分界待定",
            "A105-R02": "18 块湿区地砖：1.047% 几何找坡已机械计算，设计与安装方向待确认",
            "A105-R03": "8 樘门：门槛材料、高差、防水收口与门下净空待确认",
            "A105-R04": "Aircrete 板 3ARl_CqPrCWQrA6_$W073W 的施工语义待确认",
        },
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_svg": str(args.output_svg), "floor_register": str(args.floor_register), "threshold_register": str(args.threshold_register), "report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["mechanical_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
