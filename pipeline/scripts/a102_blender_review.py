"""Build non-IFC Blender review helpers for A-102 demolition candidates."""

from __future__ import annotations

import csv
from pathlib import Path

import bpy
import ifcopenshell.util.element
import bonsai.tool as tool
from mathutils import Vector


REGISTER = Path(
    "/Users/jiaxinchen/Documents/Projects/仁恒-滨海湾/2504 GBTB Yanlord Zhuhai/"
    "pipeline/decisions/a102-demolition-review.csv"
)
COLLECTION_NAME = "A102_DEMOLITION_REVIEW"

COLORS = {
    "existing_non_load": (0.90, 0.35, 0.05, 1.0),
    "existing_load": (0.30, 0.42, 0.55, 1.0),
    "new": (0.55, 0.18, 0.75, 1.0),
    "already_removed": (1.00, 0.05, 0.02, 1.0),
    "planned_demolition": (1.00, 0.05, 0.55, 1.0),
    "door": (0.00, 0.75, 0.85, 1.0),
    "window": (0.05, 0.35, 1.00, 1.0),
    "furniture": (1.00, 0.75, 0.00, 1.0),
    "grid": (0.15, 0.80, 0.20, 1.0),
    "slab": (0.70, 0.72, 0.74, 1.0),
}


def add_prism(
    vertices: list[tuple[float, float, float]],
    faces: list[tuple[int, ...]],
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    z0: float,
    z1: float,
) -> None:
    offset = len(vertices)
    vertices.extend(
        [
            (x0, y0, z0),
            (x1, y0, z0),
            (x1, y1, z0),
            (x0, y1, z0),
            (x0, y0, z1),
            (x1, y0, z1),
            (x1, y1, z1),
            (x0, y1, z1),
        ]
    )
    faces.extend(
        [
            (offset + 0, offset + 1, offset + 2, offset + 3),
            (offset + 4, offset + 7, offset + 6, offset + 5),
            (offset + 0, offset + 4, offset + 5, offset + 1),
            (offset + 1, offset + 5, offset + 6, offset + 2),
            (offset + 2, offset + 6, offset + 7, offset + 3),
            (offset + 4, offset + 0, offset + 3, offset + 7),
        ]
    )


def build_cage(record: dict[str, str], collection: bpy.types.Collection) -> bpy.types.Object:
    x0 = float(record["candidate_x_min_mm"]) / 1000.0
    y0 = float(record["candidate_y_min_mm"]) / 1000.0
    x1 = float(record["candidate_x_max_mm"]) / 1000.0
    y1 = float(record["candidate_y_max_mm"]) / 1000.0
    z0 = float(record["candidate_z_min_mm"]) / 1000.0
    z1 = float(record["candidate_z_max_mm"]) / 1000.0
    edge = 0.025
    half = edge / 2.0
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, ...]] = []

    for z in (z0, z1 - edge):
        add_prism(vertices, faces, x0, x1, y0 - half, y0 + half, z, z + edge)
        add_prism(vertices, faces, x0, x1, y1 - half, y1 + half, z, z + edge)
        add_prism(vertices, faces, x0 - half, x0 + half, y0, y1, z, z + edge)
        add_prism(vertices, faces, x1 - half, x1 + half, y0, y1, z, z + edge)
    for x in (x0, x1):
        for y in (y0, y1):
            add_prism(vertices, faces, x - half, x + half, y - half, y + half, z0, z1)

    name = f"A102_{record['candidate_id']}_{record['source_status']}"
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    obj = bpy.data.objects.new(name, mesh)
    collection.objects.link(obj)
    obj.color = (
        COLORS["already_removed"]
        if record["source_status"] == "ALREADY_REMOVED"
        else COLORS["planned_demolition"]
    )
    obj.display_type = "SOLID"
    obj.show_in_front = False
    obj["a102_candidate_id"] = record["candidate_id"]
    obj["a102_source_status"] = record["source_status"]
    obj["a102_confidence"] = float(record["confidence"])
    obj["a102_review_status"] = record["review_status"]
    obj["a102_ifc_write_allowed"] = False
    return obj


def build_label(record: dict[str, str], collection: bpy.types.Collection) -> bpy.types.Object:
    x = (float(record["candidate_x_min_mm"]) + float(record["candidate_x_max_mm"])) / 2000.0
    y = (float(record["candidate_y_min_mm"]) + float(record["candidate_y_max_mm"])) / 2000.0
    text_curve = bpy.data.curves.new(f"A102_LABEL_{record['candidate_id']}", "FONT")
    text_curve.body = record["candidate_id"]
    text_curve.align_x = "CENTER"
    text_curve.align_y = "CENTER"
    text_curve.size = 0.22
    text_curve.extrude = 0.008
    label = bpy.data.objects.new(f"A102_LABEL_{record['candidate_id']}", text_curve)
    collection.objects.link(label)
    label.location = (x, y, 3.05)
    label.color = (
        COLORS["already_removed"]
        if record["source_status"] == "ALREADY_REMOVED"
        else COLORS["planned_demolition"]
    )
    label.show_in_front = False
    return label


def ifc_entity(obj: bpy.types.Object):
    try:
        return tool.Ifc.get_entity(obj)
    except Exception:
        return None


def restore_context() -> dict[str, int]:
    counts: dict[str, int] = {}
    for obj in bpy.context.scene.objects:
        obj.show_in_front = False
        entity = ifc_entity(obj)
        if entity is None:
            continue
        if entity.is_a("IfcSpace") or entity.is_a("IfcOpeningElement") or entity.is_a("IfcAnnotation"):
            obj.hide_set(True)
            continue
        if entity.is_a("IfcWall"):
            obj.hide_set(False)
            common = ifcopenshell.util.element.get_psets(entity).get("Pset_WallCommon", {})
            status = common.get("Status")
            load_bearing = common.get("LoadBearing")
            if status == "NEW":
                key = "new"
            elif load_bearing is True:
                key = "existing_load"
            else:
                key = "existing_non_load"
            obj.color = COLORS[key]
            counts[key] = counts.get(key, 0) + 1
        elif entity.is_a("IfcDoor"):
            obj.hide_set(False)
            obj.color = COLORS["door"]
            counts["door"] = counts.get("door", 0) + 1
        elif entity.is_a("IfcWindow"):
            obj.hide_set(False)
            obj.color = COLORS["window"]
            counts["window"] = counts.get("window", 0) + 1
        elif entity.is_a("IfcFurniture") or entity.is_a("IfcFurnishingElement"):
            obj.hide_set(False)
            obj.color = COLORS["furniture"]
            counts["furniture"] = counts.get("furniture", 0) + 1
        elif entity.is_a("IfcGridAxis"):
            obj.hide_set(False)
            obj.color = COLORS["grid"]
            counts["grid"] = counts.get("grid", 0) + 1
        elif entity.is_a("IfcSlab"):
            world_z = [(obj.matrix_world @ Vector(corner)).z for corner in obj.bound_box]
            is_floor = max(world_z) <= 0.1
            obj.hide_set(not is_floor)
            if is_floor:
                obj.color = COLORS["slab"]
                counts["floor_slab"] = counts.get("floor_slab", 0) + 1
        elif entity.is_a("IfcCovering"):
            world_z = [(obj.matrix_world @ Vector(corner)).z for corner in obj.bound_box]
            is_floor_finish = max(world_z) <= 0.15
            obj.hide_set(not is_floor_finish)
            if is_floor_finish:
                obj.color = COLORS["slab"]
                counts["floor_finish"] = counts.get("floor_finish", 0) + 1
        else:
            obj.hide_set(True)
    return counts


def configure_viewport() -> None:
    for area in bpy.context.screen.areas:
        if area.type != "VIEW_3D":
            continue
        space = area.spaces.active
        space.shading.type = "SOLID"
        space.shading.color_type = "OBJECT"
        space.shading.show_xray = False
        space.overlay.show_overlays = True
        space.overlay.show_outline_selected = True
        space.region_3d.view_perspective = "PERSP"


def main() -> None:
    if not REGISTER.exists():
        raise RuntimeError(f"A-102 register does not exist: {REGISTER}")
    with REGISTER.open(encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle))
    if len(records) != 13:
        raise RuntimeError(f"expected 13 A-102 review records, found {len(records)}")

    old = bpy.data.collections.get(COLLECTION_NAME)
    if old is not None:
        for obj in list(old.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        bpy.data.collections.remove(old)
    collection = bpy.data.collections.new(COLLECTION_NAME)
    bpy.context.scene.collection.children.link(collection)

    counts = restore_context()
    helpers = []
    for record in records:
        helpers.append(build_cage(record, collection))
        build_label(record, collection)

    bpy.ops.object.select_all(action="DESELECT")
    for helper in helpers:
        helper.select_set(True)
    bpy.context.view_layer.objects.active = helpers[0]
    configure_viewport()
    print(
        {
            "collection": COLLECTION_NAME,
            "candidate_helpers": len(helpers),
            "already_removed": sum(record["source_status"] == "ALREADY_REMOVED" for record in records),
            "planned_demolition": sum(record["source_status"] == "PLANNED_DEMOLITION" for record in records),
            "context": counts,
            "show_in_front": 0,
            "ifc_modified": False,
        }
    )


main()
