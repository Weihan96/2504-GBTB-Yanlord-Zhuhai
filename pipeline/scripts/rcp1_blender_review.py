"""Prepare the RCP1 Blender view without treating demolition walls as built walls."""

from __future__ import annotations

import bpy
import bonsai.tool as tool
import ifcopenshell.util.element


DEMOLITION_COLLECTION = "A102_DEMOLITION_REFERENCE"
OBSOLETE_ROUTE_PREFIXES = ("RCP1_OPTION_A", "RCP1_OPTION_B")
DEMOLITION_STATUSES = {"DEMOLISH", "DEMOLISHED"}
DEMOLITION_COLOR = (0.82, 0.04, 0.03, 0.22)


def wall_status(entity) -> str:
    common = ifcopenshell.util.element.get_psets(entity).get("Pset_WallCommon", {})
    return str(common.get("Status") or "").strip().upper()


def remove_collection(name: str) -> None:
    collection = bpy.data.collections.get(name)
    if collection is None:
        return
    for obj in list(collection.objects):
        bpy.data.objects.remove(obj, do_unlink=True)
    bpy.data.collections.remove(collection)


def remove_obsolete_route_options() -> int:
    removed = 0
    for obj in list(bpy.data.objects):
        if obj.name.startswith(OBSOLETE_ROUTE_PREFIXES):
            bpy.data.objects.remove(obj, do_unlink=True)
            removed += 1
    return removed


def label_fixed_equipment_positions() -> int:
    replacements = {
        "A06 旧版东侧候选": ("A06 固定机位", "A06 固定机位 → 送回风待深化"),
        "RCP1_LEGACY_AC_LABEL": (None, "A06 固定机位（正式 placement-only 身份）"),
        "RCP1_LIVING_GAP_LABEL": (None, "R20 客厅：送回风待深化"),
    }
    updated = 0
    for object_name, (new_name, new_body) in replacements.items():
        obj = bpy.data.objects.get(object_name)
        if obj is None or obj.type != "FONT":
            continue
        if new_name:
            obj.name = new_name
        obj.data.body = new_body
        updated += 1
    legacy_mesh = bpy.data.objects.get("RCP1_A06_旧版东侧候选")
    if legacy_mesh is not None:
        legacy_mesh.name = "RCP1_A06_固定机位"
        updated += 1
    return updated


def demolition_material() -> bpy.types.Material:
    material = bpy.data.materials.get("A102_DEMOLITION_TRANSPARENT")
    if material is None:
        material = bpy.data.materials.new("A102_DEMOLITION_TRANSPARENT")
    material.diffuse_color = DEMOLITION_COLOR
    material.use_nodes = True
    principled = material.node_tree.nodes.get("Principled BSDF")
    if principled is not None:
        principled.inputs["Base Color"].default_value = DEMOLITION_COLOR
        principled.inputs["Alpha"].default_value = DEMOLITION_COLOR[3]
        principled.inputs["Roughness"].default_value = 0.72
    if hasattr(material, "surface_render_method"):
        material.surface_render_method = "DITHERED"
    return material


def create_demolition_references() -> tuple[int, list[str]]:
    remove_collection(DEMOLITION_COLLECTION)
    collection = bpy.data.collections.new(DEMOLITION_COLLECTION)
    bpy.context.scene.collection.children.link(collection)
    material = demolition_material()
    model = tool.Ifc.get()
    references: list[str] = []
    for entity in model.by_type("IfcWall"):
        if wall_status(entity) not in DEMOLITION_STATUSES:
            continue
        source = tool.Ifc.get_object(entity)
        if source is None or source.type != "MESH":
            raise RuntimeError(f"missing Blender mesh for demolition wall {entity.GlobalId}")
        source.hide_set(True)
        reference = bpy.data.objects.new(
            f"A102_DEMOLITION_REFERENCE_{entity.Tag or entity.GlobalId}",
            source.data.copy(),
        )
        reference.matrix_world = source.matrix_world.copy()
        reference.data.materials.clear()
        reference.data.materials.append(material)
        reference.color = DEMOLITION_COLOR
        reference.display_type = "SOLID"
        reference.show_in_front = False
        reference["source_ifc_global_id"] = entity.GlobalId
        reference["source_wall_status"] = wall_status(entity)
        collection.objects.link(reference)
        references.append(entity.GlobalId)
    collection.hide_viewport = True
    collection.hide_render = True
    return len(references), sorted(references)


def configure_viewport() -> None:
    for obj in bpy.context.scene.objects:
        obj.show_in_front = False
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
        space.overlay.show_relationship_lines = False


def main() -> dict:
    removed_options = remove_obsolete_route_options()
    fixed_position_labels_updated = label_fixed_equipment_positions()
    demolition_count, demolition_global_ids = create_demolition_references()
    configure_viewport()
    if demolition_count != 13:
        raise RuntimeError(f"expected 13 demolition walls, found {demolition_count}")
    return {
        "demolition_source_walls_hidden": demolition_count,
        "demolition_reference_collection": DEMOLITION_COLLECTION,
        "demolition_reference_collection_hidden_by_default": True,
        "demolition_reference_alpha": DEMOLITION_COLOR[3],
        "demolition_global_ids": demolition_global_ids,
        "obsolete_route_option_objects_removed": removed_options,
        "fixed_position_labels_updated": fixed_position_labels_updated,
        "show_in_front": False,
        "xray": False,
        "wireframes": False,
    }


RESULT = main()
print(RESULT)
