#!/usr/bin/env python3
"""Build a current-IFC electrical coordination underlay with doors and fixed joinery."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
from pathlib import Path
from typing import Any, Iterable

import ifcopenshell
import ifcopenshell.geom
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union

from plan_elevation_index import apply_official_elevation_index, update_source_manifest
from svg_audit_underlay import validate_wall_plan_source


EXPECTED_WALLS = 101
EXPECTED_DOORS = 8
EXPECTED_FIXED_FURNITURE = 70
SCALE_DENOMINATOR = 50.0
SVG_WORLD_OFFSET_MM = 10000.0


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument("--wall-plan", type=Path, default=root / "drawings/Wall Plan.svg")
    parser.add_argument("--int1", type=Path, default=root / "build/int1/int1-existing-report.json")
    parser.add_argument(
        "--output-svg", type=Path, default=root / "drawings/Electrical Coordination Plan.svg"
    )
    parser.add_argument(
        "--manifest", type=Path, default=root / "drawings/Electrical Coordination Plan-source.json"
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def project_relative(path: Path) -> str:
    root = Path(__file__).resolve().parents[2]
    resolved = path.resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError as exc:
        raise RuntimeError(f"coordination dependency resolves outside the project: {resolved}") from exc


def world_to_svg(x_mm: float, y_mm: float) -> tuple[float, float]:
    return (
        (x_mm + SVG_WORLD_OFFSET_MM) / SCALE_DENOMINATOR,
        (SVG_WORLD_OFFSET_MM - y_mm) / SCALE_DENOMINATOR,
    )


def projected_footprint(product: Any, settings: ifcopenshell.geom.settings):
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3)) * 1000.0
    faces = np.asarray(shape.geometry.faces, dtype=int).reshape((-1, 3))
    polygons = []
    for face in faces:
        polygon = Polygon(vertices[face, :2])
        if polygon.is_valid and polygon.area > 1e-6:
            polygons.append(polygon)
    geometry = unary_union(polygons) if polygons else Polygon()
    if geometry.is_empty:
        raise RuntimeError(f"{product.GlobalId}: plan footprint is empty")
    return geometry


def polygons(geometry: Any) -> Iterable[Any]:
    if geometry.geom_type == "Polygon":
        yield geometry
    elif geometry.geom_type == "MultiPolygon":
        yield from geometry.geoms
    elif geometry.geom_type == "GeometryCollection":
        for child in geometry.geoms:
            yield from polygons(child)
    else:
        raise RuntimeError(f"unsupported footprint geometry: {geometry.geom_type}")


def ring_path(coords: Iterable[tuple[float, float]]) -> str:
    points = [world_to_svg(float(x), float(y)) for x, y in coords]
    return "M " + " L ".join(f"{x:.4f} {y:.4f}" for x, y in points) + " Z"


def geometry_path(geometry: Any) -> str:
    commands: list[str] = []
    for polygon in polygons(geometry):
        commands.append(ring_path(polygon.exterior.coords))
        commands.extend(ring_path(interior.coords) for interior in polygon.interiors)
    return " ".join(commands)


def main() -> int:
    args = parse_args()
    ifc_hash = sha256(args.ifc)
    wall_source = args.wall_plan.read_text(encoding="utf-8")
    validate_wall_plan_source(wall_source, args.wall_plan, args.ifc)
    int1 = json.loads(args.int1.read_text(encoding="utf-8"))
    if int1.get("source", {}).get("ifc_sha256") != ifc_hash:
        raise RuntimeError("INT1 fixed-furniture report does not match the current formal IFC")

    model = ifcopenshell.open(args.ifc)
    walls = model.by_type("IfcWall")
    doors = model.by_type("IfcDoor")
    fixed_ids = {
        row["global_id"]
        for row in int1.get("records", [])
        if row.get("ifc_class") == "IfcFurniture"
        and row.get("installation_role") == "fixed_furniture"
    }
    fixed = [model.by_guid(global_id) for global_id in sorted(fixed_ids)]
    if (len(walls), len(doors), len(fixed)) != (
        EXPECTED_WALLS,
        EXPECTED_DOORS,
        EXPECTED_FIXED_FURNITURE,
    ):
        raise RuntimeError(
            "electrical coordination source counts changed: "
            f"walls={len(walls)}, doors={len(doors)}, fixed_furniture={len(fixed)}"
        )
    if any(product is None or not product.is_a("IfcFurniture") for product in fixed):
        raise RuntimeError("INT1 fixed-furniture selection contains a missing or non-furniture object")

    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    groups: list[str] = []
    footprint_parts = {"door": 0, "fixed_furniture": 0}
    for kind, products_in_scope in (("door", doors), ("fixed_furniture", fixed)):
        for product in products_in_scope:
            geometry = projected_footprint(product, settings)
            part_count = sum(1 for _ in polygons(geometry))
            footprint_parts[kind] += part_count
            groups.append(
                f'<g class="ec-{kind.replace("_", "-")}" '
                f'data-coordination-kind="{kind}" data-global-id="{html.escape(product.GlobalId)}">'
                f'<path d="{geometry_path(geometry)}" fill-rule="evenodd">'
                f'<title>{html.escape(str(product.Name or product.GlobalId))}</title></path></g>'
            )

    style = """
.ec-door path{fill:#f8fafc;fill-opacity:.78;stroke:#4b5563;stroke-width:.34}
.ec-fixed-furniture path{fill:#dbe8ef;fill-opacity:.34;stroke:#3d5968;stroke-width:.24}
"""
    generated = (
        f'<style id="electrical-coordination-underlay-style">{style}</style>'
        f'<g id="electrical-coordination-underlay" data-electrical-coordination="current-ifc">'
        f'{"".join(groups)}</g>'
    )
    if "</svg>" not in wall_source:
        raise RuntimeError("Wall Plan SVG has no closing root element")
    output = wall_source.replace("</svg>", generated + "</svg>", 1)
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(output, encoding="utf-8")
    elevation_index_report = apply_official_elevation_index(args.output_svg)

    manifest = {
        "mode": "current_ifc_electrical_coordination_underlay",
        "formal_ifc_sha256": ifc_hash,
        "wall_plan_svg_path": project_relative(args.wall_plan),
        "wall_plan_svg_sha256": sha256(args.wall_plan),
        "int1_report_path": project_relative(args.int1),
        "int1_report_sha256": sha256(args.int1),
        "coordination_svg_path": project_relative(args.output_svg),
        "coordination_svg_sha256": sha256(args.output_svg),
        "counts": {
            "walls": len(walls),
            "doors": len(doors),
            "fixed_furniture": len(fixed),
            "door_footprint_parts": footprint_parts["door"],
            "fixed_furniture_footprint_parts": footprint_parts["fixed_furniture"],
        },
        "gates": {
            "wall_plan_current": True,
            "doors_from_current_ifc": True,
            "fixed_furniture_from_current_ifc_and_int1_role_register": True,
            "stale_furniture_plan_not_used": True,
            "automatic_ifc_write_allowed": False,
        },
    }
    args.manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    update_source_manifest(args.output_svg, elevation_index_report)
    print(json.dumps({"counts": manifest["counts"], "gates": manifest["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
