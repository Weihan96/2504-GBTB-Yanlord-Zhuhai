"""Create a depth-aware Blender review for first-round renovation electrical demands."""

from __future__ import annotations

import json
from pathlib import Path

import bpy


ROOT_COLLECTION = "RENOVATION_ELEC_ROUND1"
GROUPS = {
    "bedside": ("01_BEDSIDE_LIGHT_CYAN", (0.00, 0.84, 0.92, 1.0)),
    "socket": ("02_NEW_SOCKET_GREEN", (0.22, 0.70, 0.30, 1.0)),
    "cabinet": ("03_CABINET_POWER_AMBER", (1.00, 0.55, 0.00, 1.0)),
    "recheck": ("04_KITCHEN_SOCKET_RECHECK_RED", (1.00, 0.02, 0.02, 1.0)),
}
FURNITURE_YELLOW = (1.0, 0.65, 0.0, 1.0)
IFC_CONTEXT_GREY = (0.62, 0.62, 0.62, 1.0)


def project_root() -> Path:
    blend_path = Path(bpy.data.filepath).resolve()
    if blend_path.parent.name in {"elec", "mep-positioning"} and blend_path.parent.parent.name == "build":
        return blend_path.parents[2]
    raise RuntimeError(f"electrical review blend is outside the expected project build directory: {blend_path}")


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


def move_to_collection(obj: bpy.types.Object, collection: bpy.types.Collection) -> None:
    for existing in list(obj.users_collection):
        existing.objects.unlink(obj)
    collection.objects.link(obj)


def point_marker(collection: bpy.types.Collection, record: dict, colour, radius: float) -> bpy.types.Object:
    position = record["position_mm"]
    bpy.ops.mesh.primitive_ico_sphere_add(
        subdivisions=2,
        radius=radius,
        location=tuple(float(value) / 1000.0 for value in position),
    )
    obj = bpy.context.object
    obj.name = record["candidate_id"]
    move_to_collection(obj, collection)
    obj.color = colour
    obj.show_in_front = False
    obj["candidate_id"] = record["candidate_id"]
    obj["candidate_kind"] = record["kind"]
    obj["room_reference"] = record["room_reference"]
    obj["automatic_ifc_write_allowed"] = False
    return obj


def cabinet_marker(collection: bpy.types.Collection, record: dict, colour) -> bpy.types.Object:
    position = record["marker_position_mm"]
    bpy.ops.mesh.primitive_cube_add(size=0.18, location=tuple(float(value) / 1000.0 for value in position))
    obj = bpy.context.object
    obj.name = record["candidate_id"]
    move_to_collection(obj, collection)
    obj.color = colour
    obj.show_in_front = False
    obj["candidate_id"] = record["candidate_id"]
    obj["candidate_kind"] = record["kind"]
    obj["room_reference"] = record["room_reference"]
    obj["coordinate_status"] = "assembly_zone_only"
    obj["source_global_ids"] = ";".join(record["source_global_ids"])
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
            max(screen.areas, key=lambda area: area.width * area.height).type = "VIEW_3D"
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
            space.region_3d.view_location = (-1.5, -1.0, 0.8)
            space.region_3d.view_distance = 17.0


def build_review() -> dict:
    report_path = project_root() / "build/elec/elec-renovation-round1-candidate.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    root = ensure_collection(ROOT_COLLECTION, bpy.context.scene.collection)
    collections = {}
    for key, (name, _colour) in GROUPS.items():
        collection = ensure_collection(name, root)
        clear_collection(collection)
        collection.hide_viewport = False
        collection.hide_render = False
        collections[key] = collection

    for row in report["bedside_light_candidates"]:
        point_marker(collections["bedside"], row, GROUPS["bedside"][1], 0.075)
    for row in report["new_socket_candidates"]:
        point_marker(collections["socket"], row, GROUPS["socket"][1], 0.075)
    for row in report["cabinet_power_zones"]:
        cabinet_marker(collections["cabinet"], row, GROUPS["cabinet"][1])
    for row in report["kitchen_socket_rechecks"]:
        point_marker(collections["recheck"], row, GROUPS["recheck"][1], 0.060)

    for name in {"DEVELOPER_HANDOVER_MEP_REFERENCE", "MEP_POSITIONING_REVIEW"}:
        previous = bpy.data.collections.get(name)
        if previous is not None:
            previous.hide_viewport = True
            previous.hide_render = True
    for obj in bpy.context.scene.objects:
        if obj.name.startswith("Ifc"):
            obj.color = FURNITURE_YELLOW if obj.name.startswith("IfcFurniture") else IFC_CONTEXT_GREY
            obj.hide_set(False)
        obj.show_in_front = False
    configure_viewport()
    bpy.context.scene["renovation_elec_round1_legend"] = (
        "cyan=bedside lights; green=new sockets; amber=cabinet power zones; red=kitchen socket recheck"
    )
    bpy.context.scene["renovation_elec_round1_writes_ifc"] = False
    return {
        "bedside": len(report["bedside_light_candidates"]),
        "socket": len(report["new_socket_candidates"]),
        "cabinet": len(report["cabinet_power_zones"]),
        "recheck": len(report["kitchen_socket_rechecks"]),
        "show_in_front": False,
        "solid_view": True,
        "formal_ifc_write": False,
        "blend_save": False,
    }


RESULT = build_review()
print(RESULT)
