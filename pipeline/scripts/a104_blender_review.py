"""Prepare Blender for direct A-104 door/window review with true depth."""

from __future__ import annotations

import csv
from pathlib import Path

import bpy
import bonsai.tool as tool
import ifcopenshell.util.element
from mathutils import Vector


LABEL_COLLECTION = "A104_REVIEW_LABELS"
COLORS = {
    "existing_wall": (0.48, 0.52, 0.57, 1.0),
    "new_wall": (0.58, 0.26, 0.72, 1.0),
    "door": (0.00, 0.72, 0.78, 1.0),
    "window": (0.10, 0.36, 0.95, 1.0),
    "review": (1.00, 0.04, 0.03, 1.0),
    "furniture": (1.00, 0.72, 0.00, 1.0),
    "grid": (0.12, 0.78, 0.20, 1.0),
    "floor": (0.66, 0.68, 0.72, 1.0),
}


def entity_for(obj: bpy.types.Object):
    try:
        return tool.Ifc.get_entity(obj)
    except Exception:
        return None


def remove_labels() -> None:
    collection = bpy.data.collections.get(LABEL_COLLECTION)
    if collection is None:
        return
    for obj in list(collection.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    bpy.data.collections.remove(collection)


def load_register() -> dict[str, dict[str, str]]:
    ifc_path = Path(str(tool.Ifc.get_path())).resolve()
    register = ifc_path.parent / "pipeline/decisions/a104-door-window-review.csv"
    with register.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 19:
        raise RuntimeError(f"expected 19 A-104 register rows, found {len(rows)}")
    return {row["global_id"]: row for row in rows}


def world_bounds(obj: bpy.types.Object) -> tuple[Vector, Vector]:
    corners = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    minimum = Vector(tuple(min(point[axis] for point in corners) for axis in range(3)))
    maximum = Vector(tuple(max(point[axis] for point in corners) for axis in range(3)))
    return minimum, maximum


def switch_door_windows_to_model_body(register: dict[str, dict[str, str]]) -> int:
    switched = 0
    model = tool.Ifc.get()
    for global_id in register:
        entity = model.by_guid(global_id)
        obj = tool.Ifc.get_object(entity)
        if obj is None or entity.Representation is None:
            raise RuntimeError(f"cannot find Blender object/representation for {global_id}")
        body = next(
            (
                representation
                for representation in entity.Representation.Representations
                if representation.RepresentationIdentifier == "Body"
                and representation.ContextOfItems.ContextType == "Model"
            ),
            None,
        )
        if body is None:
            raise RuntimeError(f"no Model/Body representation for {global_id}")
        active = tool.Geometry.get_active_representation(obj)
        if active is None or active.id() != body.id():
            result = bpy.ops.bim.switch_representation(
                obj=obj.name,
                ifc_definition_id=body.id(),
                disable_opening_subtractions=False,
            )
            if "FINISHED" not in result:
                raise RuntimeError(f"failed to show Model/Body for {global_id}: {result}")
            switched += 1
    return switched


def add_labels(review_objects: list[bpy.types.Object], register: dict[str, dict[str, str]]) -> int:
    collection = bpy.data.collections.new(LABEL_COLLECTION)
    bpy.context.scene.collection.children.link(collection)
    for obj in review_objects:
        entity = entity_for(obj)
        row = register[entity.GlobalId]
        minimum, maximum = world_bounds(obj)
        curve = bpy.data.curves.new(f"A104_LABEL_{row['candidate_id']}", "FONT")
        curve.body = f"{row['candidate_id']}  {row['review_group']}"
        curve.align_x = "CENTER"
        curve.align_y = "CENTER"
        curve.size = 0.16
        curve.extrude = 0.004
        label = bpy.data.objects.new(f"A104_LABEL_{row['candidate_id']}", curve)
        collection.objects.link(label)
        label.location = (
            (minimum.x + maximum.x) / 2.0,
            (minimum.y + maximum.y) / 2.0,
            maximum.z + 0.08,
        )
        label.color = COLORS["review"]
        label.show_in_front = False
    return len(review_objects)


def configure_objects(register: dict[str, dict[str, str]]) -> tuple[list[bpy.types.Object], dict[str, int]]:
    review_objects: list[bpy.types.Object] = []
    counts: dict[str, int] = {}
    for obj in bpy.context.scene.objects:
        obj.show_in_front = False
        entity = entity_for(obj)
        if entity is None:
            continue
        key = ""
        visible = False
        if entity.is_a("IfcWall"):
            common = ifcopenshell.util.element.get_psets(entity).get("Pset_WallCommon", {})
            if str(common.get("Status") or "").upper() != "DEMOLISH":
                visible = True
                key = "new_wall" if str(common.get("Status") or "").upper() == "NEW" else "existing_wall"
        elif entity.is_a("IfcDoor") or entity.is_a("IfcWindow"):
            visible = True
            row = register.get(entity.GlobalId)
            if row is None:
                raise RuntimeError(f"door/window missing from A-104 register: {entity.GlobalId}")
            if row["review_required"] == "yes":
                key = "review"
                review_objects.append(obj)
            else:
                key = "door" if entity.is_a("IfcDoor") else "window"
        elif entity.is_a("IfcFurniture") or entity.is_a("IfcFurnishingElement"):
            visible = True
            key = "furniture"
        elif entity.is_a("IfcGridAxis"):
            visible = True
            key = "grid"
        elif entity.is_a("IfcSlab") or entity.is_a("IfcCovering"):
            _minimum, maximum = world_bounds(obj)
            visible = maximum.z <= 0.15
            key = "floor"
        obj.hide_set(not visible)
        if visible:
            obj.color = COLORS[key]
            obj.display_type = "SOLID"
            counts[key] = counts.get(key, 0) + 1
    return review_objects, counts


def configure_viewport() -> None:
    screen = bpy.context.screen
    if screen is None:
        return
    for area in screen.areas:
        if area.type != "VIEW_3D":
            continue
        space = area.spaces.active
        space.shading.type = "SOLID"
        space.shading.light = "STUDIO"
        space.shading.color_type = "OBJECT"
        space.shading.background_type = "VIEWPORT"
        space.shading.background_color = (0.055, 0.055, 0.055)
        space.shading.show_xray = False
        space.overlay.show_overlays = True
        space.overlay.show_outline_selected = True
        space.region_3d.view_perspective = "PERSP"


def main() -> None:
    model = tool.Ifc.get()
    if model is None:
        raise RuntimeError("no IFC is loaded")
    register = load_register()
    remove_labels()
    bpy.ops.object.select_all(action="DESELECT")
    switched_to_model_body = switch_door_windows_to_model_body(register)
    review_objects, counts = configure_objects(register)
    if len(review_objects) != 8:
        raise RuntimeError(f"expected 8 A-104 review objects, found {len(review_objects)}")
    label_count = add_labels(review_objects, register)
    for obj in review_objects:
        obj.select_set(True)
    bpy.context.view_layer.objects.active = review_objects[0]
    configure_viewport()
    rows = []
    for obj in sorted(review_objects, key=lambda item: register[entity_for(item).GlobalId]["candidate_id"]):
        entity = entity_for(obj)
        minimum, maximum = world_bounds(obj)
        row = register[entity.GlobalId]
        rows.append(
            {
                "candidate_id": row["candidate_id"],
                "review_group": row["review_group"],
                "global_id": entity.GlobalId,
                "ifc_class": entity.is_a(),
                "world_bbox_mm": [
                    *[round(value * 1000.0, 3) for value in minimum],
                    *[round(value * 1000.0, 3) for value in maximum],
                ],
            }
        )
    print(
        {
            "ifc_path": str(tool.Ifc.get_path()),
            "schema": model.schema,
            "review_objects": len(review_objects),
            "review_labels": label_count,
            "switched_to_model_body": switched_to_model_body,
            "rows": rows,
            "context": counts,
            "selected": len(bpy.context.selected_objects),
            "show_in_front": sum(obj.show_in_front for obj in bpy.context.scene.objects),
            "solid_view": True,
            "xray": False,
        }
    )


main()
