"""Show confirmed HVAC constraint paths in Blender with true depth."""

from __future__ import annotations

import json
from pathlib import Path

import bmesh
import bpy
import bonsai.tool as tool
from mathutils import Vector


COLLECTION_NAME = "RCP1_ROUTE_CONSTRAINT_REVIEW"
CONSTRAINT_COLOR = (0.00, 0.82, 1.00, 1.0)
OUTDOOR_COLOR = (1.00, 0.16, 0.72, 1.0)
CONDENSATE_COLOR = (0.08, 0.42, 1.00, 1.0)
TEXT_COLOR = (1.00, 1.00, 1.00, 1.0)
PAIRING_REVIEW_HIDDEN_PREFIXES = ("RCP1_PIPE_", "RCP1_SERVICE_ARROW_")


def remove_collection() -> None:
    collection = bpy.data.collections.get(COLLECTION_NAME)
    if collection is None:
        return
    for obj in list(collection.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    bpy.data.collections.remove(collection)


def load_report() -> dict:
    root = Path(str(tool.Ifc.get_path())).resolve().parent
    return json.loads((root / "build/rcp1/route-readiness-candidate.json").read_text(encoding="utf-8"))


def point_m(values_mm: list[float]) -> Vector:
    return Vector(tuple(value / 1000.0 for value in values_mm))


def add_line(
    collection: bpy.types.Collection,
    name: str,
    start: Vector,
    end: Vector,
    color: tuple[float, float, float, float],
    bevel_depth: float = 0.025,
) -> bpy.types.Object:
    curve = bpy.data.curves.new(name, "CURVE")
    curve.dimensions = "3D"
    curve.bevel_depth = bevel_depth
    curve.bevel_resolution = 2
    spline = curve.splines.new("POLY")
    spline.points.add(1)
    spline.points[0].co = (*start, 1.0)
    spline.points[1].co = (*end, 1.0)
    obj = bpy.data.objects.new(name, curve)
    collection.objects.link(obj)
    obj.color = color
    obj.show_in_front = False
    return obj


def add_label(
    collection: bpy.types.Collection,
    name: str,
    body: str,
    location: Vector,
    color: tuple[float, float, float, float] = TEXT_COLOR,
    size: float = 0.12,
) -> bpy.types.Object:
    curve = bpy.data.curves.new(name, "FONT")
    curve.body = body
    curve.align_x = "CENTER"
    curve.align_y = "CENTER"
    curve.size = size
    curve.extrude = 0.003
    obj = bpy.data.objects.new(name, curve)
    collection.objects.link(obj)
    obj.location = location
    obj.color = color
    obj.show_in_front = False
    return obj


def add_marker(
    collection: bpy.types.Collection,
    name: str,
    location: Vector,
    color: tuple[float, float, float, float],
) -> bpy.types.Object:
    mesh = bpy.data.meshes.new(name)
    bm = bmesh.new()
    bmesh.ops.create_icosphere(bm, subdivisions=2, radius=0.08)
    bm.to_mesh(mesh)
    bm.free()
    obj = bpy.data.objects.new(name, mesh)
    collection.objects.link(obj)
    obj.location = location
    obj.color = color
    obj.display_type = "SOLID"
    obj.show_in_front = False
    return obj


def configure_viewport() -> None:
    for obj in bpy.context.scene.objects:
        obj.show_in_front = False
        if obj.name.startswith(PAIRING_REVIEW_HIDDEN_PREFIXES):
            obj.hide_set(True)
    if bpy.context.object is not None and bpy.context.object.mode != "OBJECT":
        bpy.ops.object.mode_set(mode="OBJECT")
    bpy.ops.object.select_all(action="DESELECT")
    if bpy.context.screen is None:
        return
    for area in bpy.context.screen.areas:
        if area.type != "VIEW_3D":
            continue
        space = area.spaces.active
        space.shading.type = "SOLID"
        space.shading.color_type = "OBJECT"
        space.shading.show_xray = False
        space.overlay.show_wireframes = False


def main() -> dict:
    report = load_report()
    remove_collection()
    collection = bpy.data.collections.new(COLLECTION_NAME)
    bpy.context.scene.collection.children.link(collection)
    path_objects = []
    segment_count = 0
    for route in report["confirmed_route_graph"]:
        if not route["segments"]:
            continue
        waypoints = {point["anchor_id"]: point for point in route["waypoints"]}
        for index, segment in enumerate(route["segments"], start=1):
            start = point_m(waypoints[segment["from_anchor_id"]]["centre_mm"])
            end = point_m(waypoints[segment["to_anchor_id"]]["centre_mm"])
            name = f"RCP1_CONSTRAINT_{route['route_id']}_{index:02d}"
            path_objects.append(add_line(collection, name, start, end, CONSTRAINT_COLOR))
            add_label(
                collection,
                f"{name}_LABEL",
                f"{segment['from_anchor_id']}→{segment['to_anchor_id']}",
                (start + end) / 2.0 + Vector((0.0, 0.0, 0.08)),
                CONSTRAINT_COLOR,
                0.13,
            )
            segment_count += 1

    opening_by_id = {item["opening_id"]: item for item in report["openings"]}
    for anchor_id, label, color in (
        ("H01", "H01 室外机接口位置", OUTDOOR_COLOR),
        ("H07", "H07 冷凝水排放接口", CONDENSATE_COLOR),
    ):
        location = point_m(opening_by_id[anchor_id]["centre_mm"])
        add_marker(collection, f"RCP1_ENDPOINT_{anchor_id}", location, color)
        add_label(collection, f"RCP1_ENDPOINT_{anchor_id}_LABEL", label, location + Vector((0.0, 0.0, 0.15)), color)

    add_label(
        collection,
        "RCP1_CONSTRAINT_LEGEND",
        "青线＝已确认锚点顺序骨架；不是最终风管、冷媒管或厂家接口",
        Vector((0.0, 3.85, 3.28)),
        TEXT_COLOR,
        0.18,
    )
    configure_viewport()
    for obj in path_objects:
        obj.select_set(True)
    return {
        "collection": COLLECTION_NAME,
        "route_segment_count": segment_count,
        "confirmed_paths": ["A02→H03", "A03→H04→H02"],
        "shared_openings_allowed": True,
        "multihop_routes_allowed": True,
        "outdoor_interface": "H01",
        "condensate_discharge_interface": "H07",
        "legacy_pipe_and_service_arrow_context_hidden": True,
        "show_in_front": False,
        "xray": False,
        "wireframes": False,
        "formal_ifc_write_allowed": False,
    }


RESULT = main()
print(RESULT)
