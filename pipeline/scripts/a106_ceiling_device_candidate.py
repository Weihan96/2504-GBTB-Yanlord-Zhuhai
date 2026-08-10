#!/usr/bin/env python3
"""Compile read-only A-106 smoke-alarm and bedroom-AP positioning candidates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import numpy as np
from shapely.geometry import Point, Polygon, box
from shapely.ops import unary_union

from svg_audit_underlay import validate_wall_plan_source


EXPECTED_IFC_SHA256 = "6c2fd8da9e9ad7ddbc2b63415a27f1c979e8995b880d8fce210a2dda2ef2aab6"
SPACE_IDS = {
    "R04": "2fhEbDfK1EkhJwlPikNm$b",
    "R09": "3gHz6U6BfFXgV6PnRzfOf$",
    "R14": "0WyQ2Z9pX5qgwfdTOZwAkw",
    "R20": "2wgBPVUpv2DvcZCfbe6fdv",
}
MODELED_SUPPLY_OUTLET_IDS = {
    "16Ey9Flj9BK9VRun$ozzjH",
    "3Bv_Kl3jDC5RvaUMdyge1U",
}
SCALE_DENOMINATOR = 50.0
SVG_WORLD_OFFSET_MM = 10000.0


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument("--register", type=Path, default=root / "pipeline/decisions/a106-ceiling-device-review.csv")
    parser.add_argument("--demolition-register", type=Path, default=root / "pipeline/decisions/a102-demolition-review.csv")
    parser.add_argument("--elec-existing", type=Path, default=root / "build/elec/elec-existing-candidate.json")
    parser.add_argument("--rcp1-existing", type=Path, default=root / "build/rcp1/rcp1-existing-candidate.json")
    parser.add_argument("--source-svg", type=Path, default=root / "drawings/Wall Plan.svg")
    parser.add_argument("--output", type=Path, default=root / "build/elec/a106-ceiling-device-candidate.json")
    parser.add_argument("--output-svg", type=Path, default=root / "drawings/A106-ceiling-device-candidate.svg")
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ids = [row["candidate_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise RuntimeError("A-106 candidate IDs are not unique")
    return rows


def shape_data(settings: ifcopenshell.geom.settings, product: Any) -> tuple[np.ndarray, np.ndarray]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3))
    faces = np.asarray(shape.geometry.faces, dtype=int).reshape((-1, 3))
    if not len(vertices):
        raise RuntimeError(f"{product.GlobalId} has no renderable geometry")
    return vertices, faces


def plan_polygon(vertices: np.ndarray, faces: np.ndarray) -> Any:
    polygons = []
    for face in faces:
        polygon = Polygon(vertices[face, :2])
        if polygon.area > 1e-10:
            polygons.append(polygon)
    if not polygons:
        raise RuntimeError("rendered shape has no plan polygon")
    return unary_union(polygons).buffer(1e-8)


def bbox_record(vertices: np.ndarray) -> dict[str, list[float]]:
    minimum = vertices.min(axis=0) * 1000.0
    maximum = vertices.max(axis=0) * 1000.0
    return {
        "min_mm": minimum.tolist(),
        "max_mm": maximum.tolist(),
        "dimensions_mm": (maximum - minimum).tolist(),
    }


def polygon_bbox_mm(polygon: Any) -> list[float]:
    return [float(value * 1000.0) for value in polygon.bounds]


def nearest_distance_mm(point: Point, polygons: list[Any]) -> float | None:
    if not polygons:
        return None
    return float(min(point.distance(polygon) for polygon in polygons) * 1000.0)


def world_to_svg(position_mm: list[float]) -> tuple[float, float]:
    return (
        (position_mm[0] + SVG_WORLD_OFFSET_MM) / SCALE_DENOMINATOR,
        (SVG_WORLD_OFFSET_MM - position_mm[1]) / SCALE_DENOMINATOR,
    )


def compile_report(args: argparse.Namespace) -> dict[str, Any]:
    source_hash = sha256(args.ifc)
    if source_hash != EXPECTED_IFC_SHA256:
        raise RuntimeError(f"formal IFC hash changed: {source_hash}")
    existing = json.loads(args.elec_existing.read_text(encoding="utf-8"))
    if existing["source"]["sha256"] != source_hash:
        raise RuntimeError("ELEC existing report does not match the formal IFC")
    rcp1_existing = json.loads(args.rcp1_existing.read_text(encoding="utf-8"))
    if rcp1_existing["source"]["sha256"] != source_hash:
        raise RuntimeError("RCP1 existing report does not match the formal IFC")
    demolition_rows = read_csv(args.demolition_register)
    demolition_ids = [row["candidate_global_id"] for row in demolition_rows if row["ifc_status_candidate"] == "DEMOLISH"]
    if len(demolition_ids) != 13 or len(demolition_ids) != len(set(demolition_ids)):
        raise RuntimeError("A-106 requires exactly 13 unique DEMOLISH wall IDs")

    model = ifcopenshell.open(args.ifc)
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)

    spaces: dict[str, dict[str, Any]] = {}
    for reference, global_id in SPACE_IDS.items():
        product = model.by_guid(global_id)
        vertices, faces = shape_data(settings, product)
        polygon = plan_polygon(vertices, faces)
        spaces[reference] = {
            "global_id": global_id,
            "polygon": polygon,
            "bbox_mm": polygon_bbox_mm(polygon),
            "centre_mm": [
                float((polygon.bounds[0] + polygon.bounds[2]) * 500.0),
                float((polygon.bounds[1] + polygon.bounds[3]) * 500.0),
            ],
        }

    beams = []
    for product in model.by_type("IfcBeam"):
        vertices, faces = shape_data(settings, product)
        if float(vertices[:, 2].max() * 1000.0) < 2400.0:
            continue
        beams.append({
            "global_id": product.GlobalId,
            "polygon": plan_polygon(vertices, faces),
            "bbox": bbox_record(vertices),
        })

    ceilings = []
    for product in model.by_type("IfcBuildingElementProxy"):
        if "Ceiling" not in str(product.Name or ""):
            continue
        vertices, faces = shape_data(settings, product)
        ceilings.append({
            "global_id": product.GlobalId,
            "name": str(product.Name or ""),
            "polygon": plan_polygon(vertices, faces),
            "bbox": bbox_record(vertices),
        })

    outlets = []
    for global_id in sorted(MODELED_SUPPLY_OUTLET_IDS):
        product = model.by_guid(global_id)
        vertices, faces = shape_data(settings, product)
        outlets.append({
            "global_id": global_id,
            "name": str(product.Name or ""),
            "polygon": plan_polygon(vertices, faces),
            "bbox": bbox_record(vertices),
        })

    obstacle_categories = {
        "typed_high_equipment",
        "high_flow_segments",
        "other_high_proxies",
    }
    obstacle_ids = {
        row["global_id"]: category
        for category, rows in rcp1_existing["inventory"].items()
        if category in obstacle_categories
        for row in rows
    }
    obstacle_ids.update({
        row["global_id"]: "non_ceiling_dcl_proxy"
        for row in rcp1_existing["inventory"]["dcl_proxies"]
        if "Ceiling" not in str(row["name"] or "") and row["global_id"] not in MODELED_SUPPLY_OUTLET_IDS
    })
    high_level_obstacles = []
    for global_id, category in sorted(obstacle_ids.items()):
        product = model.by_guid(global_id)
        vertices, faces = shape_data(settings, product)
        high_level_obstacles.append({
            "global_id": global_id,
            "name": str(product.Name or ""),
            "category": category,
            "polygon": plan_polygon(vertices, faces),
            "bbox": bbox_record(vertices),
        })

    lights = []
    for row in existing["sheets"]["E-301"]["lights"]:
        minimum = row["bbox"]["min_mm"]
        maximum = row["bbox"]["max_mm"]
        lights.append({
            "global_id": row["global_id"],
            "candidate_id": row["candidate_id"],
            "room_name": row["candidate_space"]["long_name"],
            "centre_mm": row["bbox"]["centre_mm"],
            "polygon": box(minimum[0] / 1000.0, minimum[1] / 1000.0, maximum[0] / 1000.0, maximum[1] / 1000.0),
        })

    records = []
    for row in read_csv(args.register):
        reference = row["room_reference"]
        room = spaces[reference]
        position = [float(row["x_mm"]), float(row["y_mm"]), float(row["z_mm"])]
        point = Point(position[0] / 1000.0, position[1] / 1000.0)
        room_beams = [item for item in beams if item["polygon"].intersects(room["polygon"].buffer(0.5))]
        room_lights = [item for item in lights if item["room_name"] == row["room_name"]]
        room_outlets = [item for item in outlets if item["polygon"].intersects(room["polygon"])]
        room_obstacles = [item for item in high_level_obstacles if item["polygon"].intersects(room["polygon"])]
        matching_ceilings = [item for item in ceilings if item["polygon"].buffer(args.tolerance_mm / 1000.0).covers(point)]
        matching_ceilings.sort(key=lambda item: abs(item["bbox"]["min_mm"][2] - position[2]))
        ceiling = matching_ceilings[0] if matching_ceilings else None

        wall_clearance = float(point.distance(room["polygon"].boundary) * 1000.0) if room["polygon"].covers(point) else 0.0
        beam_clearance = nearest_distance_mm(point, [item["polygon"] for item in room_beams])
        light_clearance = nearest_distance_mm(point, [item["polygon"] for item in room_lights])
        outlet_clearance = nearest_distance_mm(point, [item["polygon"] for item in room_outlets])
        obstacle_distances = sorted(
            (float(point.distance(item["polygon"]) * 1000.0), item)
            for item in room_obstacles
        )
        obstacle_clearance = obstacle_distances[0][0] if obstacle_distances else None
        nearest_obstacle = obstacle_distances[0][1] if obstacle_distances else None
        ceiling_residual = None if ceiling is None else abs(float(ceiling["bbox"]["min_mm"][2]) - position[2])
        room_axis_offset = abs(position[0] - room["centre_mm"][0])
        required_light_clearance = 500.0 if row["device_role"] in {"smoke_alarm", "kitchen_fire_sensor"} else 200.0
        known_geometry_pass = (
            room["polygon"].covers(point)
            and wall_clearance + args.tolerance_mm >= 500.0
            and (beam_clearance is None or beam_clearance + args.tolerance_mm >= 500.0)
            and (light_clearance is None or light_clearance + args.tolerance_mm >= required_light_clearance)
            and (obstacle_clearance is None or obstacle_clearance + args.tolerance_mm >= 500.0)
            and ceiling_residual is not None
            and ceiling_residual <= args.tolerance_mm
        )
        supply_status = "known_modeled_clearance_pass" if outlet_clearance is not None and outlet_clearance + args.tolerance_mm >= 1500.0 else (
            "known_modeled_clearance_fail" if outlet_clearance is not None else "unknown_incomplete_supply_air_model"
        )
        records.append({
            **row,
            "position_mm": position,
            "space_global_id": room["global_id"],
            "space_bbox_mm": room["bbox_mm"],
            "room_axis_offset_mm": room_axis_offset,
            "room_center_offset_mm": float(np.linalg.norm(np.asarray(position[:2]) - np.asarray(room["centre_mm"]))),
            "wall_boundary_clearance_mm": wall_clearance,
            "nearest_beam_clearance_mm": beam_clearance,
            "nearest_light_edge_clearance_mm": light_clearance,
            "required_light_clearance_mm": required_light_clearance,
            "modeled_supply_air_edge_clearance_mm": outlet_clearance,
            "supply_air_clearance_status": supply_status,
            "nearest_high_level_obstacle_clearance_mm": obstacle_clearance,
            "nearest_high_level_obstacle": None if nearest_obstacle is None else {
                "global_id": nearest_obstacle["global_id"],
                "name": nearest_obstacle["name"],
                "category": nearest_obstacle["category"],
            },
            "ceiling_context_global_id": None if ceiling is None else ceiling["global_id"],
            "ceiling_context_name": None if ceiling is None else ceiling["name"],
            "ceiling_underside_mm": None if ceiling is None else ceiling["bbox"]["min_mm"][2],
            "ceiling_datum_residual_mm": ceiling_residual,
            "known_geometry_pass": known_geometry_pass,
            "final_release_pass": False,
            "release_blocker": (
                "confirmed supply/return air layout and fire-safety review pending"
                if row["device_role"] == "smoke_alarm"
                else "final heat/smoke/combined type, product, power, communication and manufacturer installation conditions pending"
                if row["device_role"] == "kitchen_fire_sensor"
                else "AP product, PoE route and signal/service review pending"
            ),
        })

    pair_checks = []
    for reference in ("R09", "R14"):
        pair = [record for record in records if record["room_reference"] == reference]
        if len(pair) != 2:
            raise RuntimeError(f"{reference} must have one smoke alarm and one AP")
        midpoint = [sum(record["position_mm"][axis] for record in pair) / 2.0 for axis in range(2)]
        centre = spaces[reference]["centre_mm"]
        pair_checks.append({
            "room_reference": reference,
            "candidate_ids": [record["candidate_id"] for record in pair],
            "separation_mm": float(np.linalg.norm(np.asarray(pair[0]["position_mm"][:2]) - np.asarray(pair[1]["position_mm"][:2]))),
            "midpoint_mm": midpoint,
            "room_centre_mm": centre,
            "symmetry_residual_mm": float(np.linalg.norm(np.asarray(midpoint) - np.asarray(centre))),
        })

    smoke = [record for record in records if record["device_role"] == "smoke_alarm"]
    aps = [record for record in records if record["device_role"] == "wireless_access_point"]
    kitchen_fire = [record for record in records if record["device_role"] == "kitchen_fire_sensor"]
    return {
        "mode": "read_only_a106_ceiling_device_position_candidate",
        "source_ifc_sha256": source_hash,
        "tolerance_mm": args.tolerance_mm,
        "summary": {
            "candidate_count": len(records),
            "smoke_alarm_candidates": len(smoke),
            "bedroom_AP_candidates": len(aps),
            "kitchen_fire_sensor_candidates": len(kitchen_fire),
            "known_geometry_pass_count": sum(bool(record["known_geometry_pass"]) for record in records),
            "final_release_pass_count": sum(bool(record["final_release_pass"]) for record in records),
        },
        "candidates": records,
        "bedroom_pair_checks": pair_checks,
        "modeled_supply_air_outlets": [
            {"global_id": item["global_id"], "name": item["name"], "bbox": item["bbox"]}
            for item in outlets
        ],
        "existing_light_markers": [
            {"global_id": item["global_id"], "candidate_id": item["candidate_id"], "centre_mm": item["centre_mm"]}
            for item in lights
        ],
        "excluded_demolition_wall_global_ids": demolition_ids,
        "gates": {
            "six_candidates_present": len(records) == 6,
            "three_smoke_candidates_present": len(smoke) == 3,
            "two_bedroom_AP_candidates_present": len(aps) == 2,
            "one_kitchen_fire_sensor_present": len(kitchen_fire) == 1,
            "known_geometry_pass": all(bool(record["known_geometry_pass"]) for record in records),
            "candidate_xy_coordinates_are_integer_mm": all(abs(value - round(value)) <= args.tolerance_mm for record in records for value in record["position_mm"][:2]),
            "bedroom_pair_separation_pass": all(check["separation_mm"] >= 700.0 for check in pair_checks),
            "thirteen_demolition_walls_excluded": len(demolition_ids) == 13,
            "all_smoke_supply_air_clearances_verified": all(record["supply_air_clearance_status"] == "known_modeled_clearance_pass" for record in smoke),
            "automatic_ifc_write_allowed": False,
        },
    }


def render_svg(source: str, report: dict[str, Any], output: Path) -> None:
    source = re.sub(r'width="400(?:\.0+)?mm"', 'width="500mm"', source, count=1)
    source = re.sub(r'viewBox="0 0 400(?:\.0+)? 400(?:\.0+)?"', 'viewBox="0 0 500 400"', source, count=1)
    demolition_ids = report["excluded_demolition_wall_global_ids"]
    hidden_group_counts = {global_id: 0 for global_id in demolition_ids}
    wall_group_pattern = re.compile(r'<g\b[^>]*\bclass="[^"]*\bIfcWall\b[^"]*"[^>]*>')

    def hide_demolition_group(match: re.Match[str]) -> str:
        tag = match.group(0)
        matched_ids = [global_id for global_id in demolition_ids if global_id in tag]
        if not matched_ids:
            return tag
        for global_id in matched_ids:
            hidden_group_counts[global_id] += 1
        return re.sub(
            r'class="([^"]*)"',
            lambda class_match: f'class="{class_match.group(1)} a106-excluded-demolish"',
            tag,
            count=1,
        )

    source = wall_group_pattern.sub(hide_demolition_group, source)
    missing_demolition_groups = [global_id for global_id, count in hidden_group_counts.items() if count == 0]
    if missing_demolition_groups:
        raise RuntimeError(f"could not exclude DEMOLISH walls from A-106 SVG: {missing_demolition_groups}")
    if "</svg>" not in source or 'viewBox="0 0 500 400"' not in source:
        raise RuntimeError("could not prepare 500x400 A-106 SVG")

    markup = []
    for record in report["candidates"]:
        x, y = world_to_svg(record["position_mm"])
        css = (
            "a106-smoke"
            if record["device_role"] == "smoke_alarm"
            else "a106-fire"
            if record["device_role"] == "kitchen_fire_sensor"
            else "a106-ap"
        )
        label = html.escape(record["candidate_id"])
        if record["device_role"] in {"smoke_alarm", "kitchen_fire_sensor"}:
            markup.append(f'<circle class="a106-known-clearance" cx="{x:.3f}" cy="{y:.3f}" r="{500 / SCALE_DENOMINATOR:.3f}"/>')
        markup.append(f'<circle class="{css}" cx="{x:.3f}" cy="{y:.3f}" r="2.4"/><text class="a106-label" x="{x+3.5:.3f}" y="{y-3.0:.3f}">{label}</text>')

    for light in report["existing_light_markers"]:
        x, y = world_to_svg(light["centre_mm"])
        markup.append(f'<circle class="a106-light" cx="{x:.3f}" cy="{y:.3f}" r="0.8"/>')

    markup.extend([
        '<g><rect class="a106-panel" x="402" y="7" width="93" height="386"/>',
        '<text class="a106-title" x="407" y="16">A-106 天花设备定位候选</text>',
        '<text class="a106-note" x="407" y="24">Bonsai 同批材质底图｜厨房火灾探测为 IFC 定位点｜非施工发布</text>',
        '<text class="a106-text" x="407" y="39">橙点：烟感候选（3）</text>',
        '<text class="a106-text" x="407" y="47">蓝点：卧室吸顶 AP 候选（2）</text>',
        '<text class="a106-text" x="407" y="55">红点：厨房火灾探测 IFC 定位点（1）</text>',
        '<text class="a106-text" x="407" y="63">灰点：正式 IFC 既有灯具（79）</text>',
        '<text class="a106-text" x="407" y="71">淡圈：500 mm 已知几何检查范围</text>',
        '<text class="a106-text" x="407" y="86">烟感：居中或中心轴最近安全点</text>',
        '<text class="a106-text" x="407" y="94">13 面 DEMOLISH 墙已排除</text>',
        '<text class="a106-warn" x="407" y="111">风口模型不完整：送风禁距门未关闭</text>',
        '<text class="a106-warn" x="407" y="119">厨房火灾点位已确认；最终类型/产品待确认</text>',
        '<text class="a106-warn" x="407" y="127">厨房燃气探测点位及型号待确认</text>',
        f'<text class="a106-note" x="407" y="382">IFC SHA {report["source_ifc_sha256"][:12]}…</text></g>',
    ])
    style = """
@page{size:500mm 400mm;margin:0}.a106-excluded-demolish{display:none!important}.a106-smoke{fill:#f59f00;stroke:#7a4d00;stroke-width:.7}.a106-fire{fill:#e03131;stroke:#7d1010;stroke-width:.7}.a106-ap{fill:#228be6;stroke:#0b477d;stroke-width:.7}.a106-light{fill:#868e96;stroke:#343a40;stroke-width:.25}.a106-known-clearance{fill:#f59f00;fill-opacity:.035;stroke:#f59f00;stroke-opacity:.32;stroke-width:.3;stroke-dasharray:1.2 1.2}.a106-label,.a106-title,.a106-note,.a106-text,.a106-warn{font-family:Arial,'Noto Sans CJK SC',sans-serif;fill:#102f43}.a106-label{font-size:2.2px;font-weight:700;paint-order:stroke;stroke:#fff;stroke-width:.8px}.a106-panel{fill:#fbfcfd;stroke:#102f43;stroke-width:.5}.a106-title{font-size:3.7px;font-weight:700}.a106-note{font-size:2.15px;fill:#526777}.a106-text{font-size:2.3px}.a106-warn{font-size:2.15px;fill:#c92a2a;font-weight:700}
"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(source.replace("</svg>", f'<style id="a106-style">{style}</style><g id="a106-ceiling-device-candidate">{"".join(markup)}</g></svg>', 1), encoding="utf-8")


def main() -> int:
    args = parse_args()
    report = compile_report(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    source = args.source_svg.read_text(encoding="utf-8")
    validate_wall_plan_source(source, args.source_svg, args.ifc)
    render_svg(source, report, args.output_svg)
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
