"""Prepare the A-102 candidate IFC for direct IfcWall review in Blender."""

from __future__ import annotations

import bpy
import ifcopenshell.util.element
import bonsai.tool as tool
from mathutils import Vector


OLD_HELPER_COLLECTION = "A102_DEMOLITION_REVIEW"
LABEL_COLLECTION = "A102_DIMENSION_LABELS"
COLORS = {
    "existing_non_load": (0.90, 0.35, 0.05, 1.0),
    "existing_load": (0.30, 0.42, 0.55, 1.0),
    "new": (0.55, 0.18, 0.75, 1.0),
    "already_removed": (1.00, 0.04, 0.02, 1.0),
    "planned_demolition": (1.00, 0.04, 0.55, 1.0),
    "door": (0.00, 0.75, 0.85, 1.0),
    "window": (0.05, 0.35, 1.00, 1.0),
    "furniture": (1.00, 0.75, 0.00, 1.0),
    "grid": (0.15, 0.80, 0.20, 1.0),
    "floor": (0.70, 0.72, 0.74, 1.0),
}


def entity_for(obj: bpy.types.Object):
    try:
        return tool.Ifc.get_entity(obj)
    except Exception:
        return None


def world_z_bounds(obj: bpy.types.Object) -> tuple[float, float]:
    values = [(obj.matrix_world @ Vector(corner)).z for corner in obj.bound_box]
    return min(values), max(values)


def remove_old_helpers() -> None:
    for name in (OLD_HELPER_COLLECTION, LABEL_COLLECTION):
        collection = bpy.data.collections.get(name)
        if collection is None:
            continue
        for obj in list(collection.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        bpy.data.collections.remove(collection)


def add_dimension_labels(demolition: list[bpy.types.Object]) -> int:
    collection = bpy.data.collections.new(LABEL_COLLECTION)
    bpy.context.scene.collection.children.link(collection)
    for obj in demolition:
        entity = entity_for(obj)
        quantities = ifcopenshell.util.element.get_psets(entity).get("Qto_WallBaseQuantities", {})
        length = float(quantities["Length"])
        width = float(quantities["Width"])
        corners = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
        centre_x = (min(point.x for point in corners) + max(point.x for point in corners)) / 2.0
        centre_y = (min(point.y for point in corners) + max(point.y for point in corners)) / 2.0
        top_z = max(point.z for point in corners)
        curve = bpy.data.curves.new(f"A102_LABEL_{entity.Tag}", "FONT")
        curve.body = f"{entity.Tag}  {length:.0f}x{width:.0f}"
        curve.align_x = "CENTER"
        curve.align_y = "CENTER"
        curve.size = 0.16
        curve.extrude = 0.004
        label = bpy.data.objects.new(f"A102_LABEL_{entity.Tag}", curve)
        collection.objects.link(label)
        label.location = (centre_x, centre_y, top_z + 0.03)
        label.color = (1.0, 1.0, 1.0, 1.0)
        label.show_in_front = False
    return len(demolition)


def prepare_objects() -> tuple[list[bpy.types.Object], dict[str, int]]:
    demolition: list[bpy.types.Object] = []
    counts: dict[str, int] = {}
    for obj in bpy.context.scene.objects:
        obj.show_in_front = False
        entity = entity_for(obj)
        if entity is None:
            continue
        if entity.is_a("IfcSpace") or entity.is_a("IfcOpeningElement") or entity.is_a("IfcAnnotation"):
            obj.hide_set(True)
            continue
        if entity.is_a("IfcWall"):
            psets = ifcopenshell.util.element.get_psets(entity)
            common = psets.get("Pset_WallCommon", {})
            status = common.get("Status")
            if status == "DEMOLISH":
                review = psets.get("Pset_A102DemolitionReview", {})
                source_status = review.get("SourceStatus")
                key = "already_removed" if source_status == "ALREADY_REMOVED" else "planned_demolition"
                demolition.append(obj)
            elif status == "NEW":
                key = "new"
            elif common.get("LoadBearing") is True:
                key = "existing_load"
            else:
                key = "existing_non_load"
            obj.hide_set(False)
            obj.color = COLORS[key]
            obj.display_type = "SOLID"
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
            obj.display_type = "SOLID"
            counts["grid"] = counts.get("grid", 0) + 1
        elif entity.is_a("IfcSlab"):
            _, z_max = world_z_bounds(obj)
            visible = z_max <= 0.1
            obj.hide_set(not visible)
            if visible:
                obj.color = COLORS["floor"]
                counts["floor_slab"] = counts.get("floor_slab", 0) + 1
        elif entity.is_a("IfcCovering"):
            _, z_max = world_z_bounds(obj)
            visible = z_max <= 0.15
            obj.hide_set(not visible)
            if visible:
                obj.color = COLORS["floor"]
                counts["floor_finish"] = counts.get("floor_finish", 0) + 1
        else:
            obj.hide_set(True)
    return demolition, counts


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
        space.shading.background_color = (0.05, 0.05, 0.05)
        space.shading.show_xray = False
        space.overlay.show_overlays = True
        space.overlay.show_outline_selected = True
        space.region_3d.view_perspective = "PERSP"


def main() -> None:
    model = tool.Ifc.get()
    if model is None:
        raise RuntimeError("no IFC is loaded")
    remove_old_helpers()
    demolition, counts = prepare_objects()
    if len(demolition) != 13:
        raise RuntimeError(f"expected 13 A-102 IfcWall candidates, found {len(demolition)}")
    label_count = add_dimension_labels(demolition)
    rows = []
    for obj in sorted(demolition, key=lambda item: entity_for(item).Tag or ""):
        entity = entity_for(obj)
        psets = ifcopenshell.util.element.get_psets(entity)
        review = psets.get("Pset_A102DemolitionReview", {})
        dimensions_mm = [round(float(value) * 1000.0, 3) for value in obj.dimensions]
        rows.append(
            {
                "candidate_id": entity.Tag,
                "global_id": entity.GlobalId,
                "ifc_class": entity.is_a(),
                "name": entity.Name,
                "source_status": review.get("SourceStatus"),
                "dimensions_mm": dimensions_mm,
            }
        )
    bpy.ops.object.select_all(action="DESELECT")
    configure_viewport()
    print(
        {
            "ifc_path": str(tool.Ifc.get_path()),
            "schema": model.schema,
            "ifc_walls": len(model.by_type("IfcWall")),
            "demolition_ifc_walls": len(demolition),
            "dimension_labels": label_count,
            "rows": rows,
            "context": counts,
            "selected": len(bpy.context.selected_objects),
            "show_in_front": sum(obj.show_in_front for obj in bpy.context.scene.objects),
            "solid_view": True,
            "xray": False,
        }
    )


main()
