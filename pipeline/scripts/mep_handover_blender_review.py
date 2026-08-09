"""Create depth-aware Blender markers for developer handover MEP references."""

from __future__ import annotations

import json
from pathlib import Path

import bpy


ROOT_COLLECTION = "DEVELOPER_HANDOVER_MEP_REFERENCE"
GROUPS = {
    "electric": ("01_ELECTRIC_PURPLE", (0.46, 0.12, 1.00, 1.0)),
    "switch": ("02_SWITCH_ORANGE", (1.00, 0.28, 0.02, 1.0)),
    "plum": ("03_PLUM_BLUE", (0.02, 0.36, 1.00, 1.0)),
    "relocate": ("04_RELOCATION_REVIEW_RED", (1.00, 0.01, 0.01, 1.0)),
}
FURNITURE_YELLOW = (1.0, 0.65, 0.0, 1.0)
IFC_CONTEXT_GREY = (0.62, 0.62, 0.62, 1.0)


def project_root() -> Path:
    blend_path = Path(bpy.data.filepath).resolve()
    if blend_path.parent.name == "mep-positioning" and blend_path.parent.parent.name == "build":
        return blend_path.parents[2]
    raise RuntimeError(f"MEP review blend is outside the expected project build directory: {blend_path}")


def ensure_collection(name: str, parent: bpy.types.Collection) -> bpy.types.Collection:
    collection = bpy.data.collections.get(name)
    if collection is None:
        collection = bpy.data.collections.new(name)
    if collection not in parent.children.values():
        parent.children.link(collection)
    return collection


def clear_collection(collection: bpy.types.Collection) -> None:
    for obj in list(collection.objects):
        bpy.data.objects.remove(obj, do_unlink=True)


def marker(collection: bpy.types.Collection, record: dict, colour, default_z_mm: float) -> bpy.types.Object:
    x_mm, y_mm = record["ifc_position_mm"]
    height_mm = record.get("installation_height_mm")
    z_mm = height_mm if height_mm is not None else default_z_mm
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=2, radius=0.065, location=(x_mm / 1000.0, y_mm / 1000.0, z_mm / 1000.0))
    obj = bpy.context.object
    obj.name = record["candidate_id"]
    for existing in list(obj.users_collection):
        existing.objects.unlink(obj)
    collection.objects.link(obj)
    obj.color = colour
    obj.show_in_front = False
    obj["candidate_id"] = record["candidate_id"]
    obj["candidate_role"] = record["candidate_role"]
    obj["candidate_space"] = record["candidate_space"]["candidate_reference"]
    obj["source_status"] = "developer_handover_existing_reference"
    obj["installation_height_mm"] = height_mm if height_mm is not None else -1.0
    obj["automatic_ifc_write_allowed"] = False
    return obj


def configure_viewport() -> None:
    layout_workspace = bpy.data.workspaces.get("Layout")
    if layout_workspace is not None:
        for window in bpy.context.window_manager.windows:
            window.workspace = layout_workspace
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type == "CONSOLE":
                area.type = "VIEW_3D"
        if not any(area.type == "VIEW_3D" for area in screen.areas):
            largest = max(screen.areas, key=lambda area: area.width * area.height)
            largest.type = "VIEW_3D"
        for area in screen.areas:
            if area.type != "VIEW_3D":
                continue
            view_spaces = [space for space in area.spaces if space.type == "VIEW_3D"]
            if not view_spaces:
                continue
            space = view_spaces[0]
            space.shading.type = "SOLID"
            space.shading.color_type = "OBJECT"
            space.shading.show_xray = False
            space.overlay.show_wireframes = False
            space.overlay.show_relationship_lines = False
            space.region_3d.view_perspective = "ORTHO"
            space.region_3d.view_location = (0.0, 0.0, 0.8)
            space.region_3d.view_distance = 18.0


def build_review() -> dict:
    path = project_root() / "build/mep-positioning/mep-renovation-delta-candidate.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    root = ensure_collection(ROOT_COLLECTION, bpy.context.scene.collection)
    collections = {}
    for key, (name, _colour) in GROUPS.items():
        collection = ensure_collection(name, root)
        clear_collection(collection)
        collection.hide_viewport = False
        collection.hide_render = False
        collections[key] = collection

    counts = {key: 0 for key in GROUPS}
    for record in report["electrical"]:
        if record["disposition_candidate"] == "relocation_or_replacement_review":
            key = "relocate"
        else:
            key = "switch" if record["kind"] == "developer_switch_or_control" else "electric"
        marker(collections[key], record, GROUPS[key][1], 0.0)
        counts[key] += 1
    for record in report["plumbing"]:
        marker(collections["plum"], record, GROUPS["plum"][1], 0.0)
        counts["plum"] += 1

    previous_review = bpy.data.collections.get("MEP_POSITIONING_REVIEW")
    if previous_review is not None:
        previous_review.hide_viewport = True
        previous_review.hide_render = True
    # Keep the current IFC as neutral depth context. Only the imported handover
    # reference markers carry review colours; furniture remains yellow.
    for obj in bpy.context.scene.objects:
        if obj.name.startswith("Ifc"):
            obj.color = FURNITURE_YELLOW if obj.name.startswith("IfcFurniture") else IFC_CONTEXT_GREY
            obj.hide_set(False)
    for obj in bpy.context.scene.objects:
        obj.show_in_front = False
    configure_viewport()
    bpy.context.scene["developer_handover_mep_legend"] = "purple=electric; orange=switch/control; blue=plumbing; red=relocation review"
    bpy.context.scene["developer_handover_mep_is_current_design"] = False
    bpy.context.scene["developer_handover_mep_writes_ifc"] = False
    return {
        **counts,
        "show_in_front": False,
        "solid_view": True,
        "formal_ifc_write": False,
        "blend_save": False,
    }


RESULT = build_review()
print(RESULT)
