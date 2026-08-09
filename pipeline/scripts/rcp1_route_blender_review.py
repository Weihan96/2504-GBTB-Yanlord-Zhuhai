"""Show fixed-equipment/opening pairing evidence in Blender with true depth."""

from __future__ import annotations

import json
from pathlib import Path

import bmesh
import bpy
import bonsai.tool as tool
from mathutils import Vector


COLLECTION_NAME = "RCP1_ROUTE_READINESS_REVIEW"
PAIR_COLOR = (0.00, 0.82, 1.00, 1.0)
UNRESOLVED_COLOR = (1.00, 0.08, 0.03, 1.0)
CHAIN_COLOR = (1.00, 0.38, 0.00, 1.0)
TEXT_COLOR = (1.00, 1.00, 1.00, 1.0)
REVIEW_Z_M = 3.12
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


def point_m(values_mm: list[float], z_m: float = REVIEW_Z_M) -> Vector:
    return Vector((values_mm[0] / 1000.0, values_mm[1] / 1000.0, z_m))


def add_line(
    collection: bpy.types.Collection,
    name: str,
    start: Vector,
    end: Vector,
    color: tuple[float, float, float, float],
    bevel_depth: float = 0.038,
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
    bmesh.ops.create_icosphere(bm, subdivisions=2, radius=0.09)
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
    equipment = {item["equipment_id"]: item for item in report["equipment"]}
    openings = {item["opening_id"]: item for item in report["openings"]}
    assignment = report["global_equipment_opening_assignment"]
    remove_collection()
    collection = bpy.data.collections.new(COLLECTION_NAME)
    bpy.context.scene.collection.children.link(collection)
    pair_objects = []
    for pair in assignment["pairs"]:
        equipment_point = point_m(equipment[pair["equipment_id"]]["centre_mm"])
        opening_point = point_m(openings[pair["opening_id"]]["centre_mm"])
        name = f"RCP1_PAIR_{pair['equipment_id']}_{pair['opening_id']}"
        pair_objects.append(add_line(collection, name, equipment_point, opening_point, PAIR_COLOR))
        midpoint = (equipment_point + opening_point) / 2.0 + Vector((0.0, 0.0, 0.06))
        add_label(
            collection,
            f"{name}_LABEL",
            f"{pair['equipment_id']}–{pair['opening_id']} {pair['clearance_mm']:.0f} mm",
            midpoint,
            PAIR_COLOR,
            0.145,
        )

    h01 = point_m(openings["H01"]["centre_mm"])
    h06 = point_m(openings["H06"]["centre_mm"])
    add_marker(collection, "RCP1_H01_UNUSED_MARKER", h01, UNRESOLVED_COLOR)
    add_label(
        collection,
        "RCP1_H01_UNUSED_LABEL",
        "H01 未配对 / 外部接口链待确认",
        h01 + Vector((0.0, 0.0, 0.16)),
        UNRESOLVED_COLOR,
        0.12,
    )
    add_line(collection, "RCP1_CHAIN_H06_H01_CANDIDATE", h06, h01, CHAIN_COLOR, 0.016)
    add_label(
        collection,
        "RCP1_CHAIN_H06_H01_LABEL",
        "H06→H01？仅穿墙链候选",
        (h06 + h01) / 2.0 + Vector((0.0, 0.0, 0.06)),
        CHAIN_COLOR,
        0.145,
    )
    add_label(
        collection,
        "RCP1_PAIRING_ONLY_LEGEND",
        "青线＝设备—既有洞口最可信配对（不是风管/冷媒管）",
        Vector((0.0, 3.85, REVIEW_Z_M + 0.12)),
        TEXT_COLOR,
        0.19,
    )
    configure_viewport()
    for obj in pair_objects:
        obj.select_set(True)
    return {
        "collection": COLLECTION_NAME,
        "pair_count": len(pair_objects),
        "pairing_total_clearance_mm": assignment["total_clearance_mm"],
        "best_to_second_margin_mm": assignment["best_to_second_margin_mm"],
        "unused_opening": "H01",
        "h06_to_h01_chain_status": "candidate_only",
        "review_elevation_m": REVIEW_Z_M,
        "legacy_pipe_and_service_arrow_context_hidden": True,
        "show_in_front": False,
        "xray": False,
        "wireframes": False,
        "formal_ifc_write_allowed": False,
    }


RESULT = main()
print(RESULT)
