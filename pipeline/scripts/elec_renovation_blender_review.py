"""Create a depth-aware Blender review for first-round renovation electrical demands."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import bpy
import ifcopenshell


ROOT_COLLECTION = "RENOVATION_ELEC_ROUND1"
GROUPS = {
    "bedside": ("01_BEDSIDE_LIGHT_CYAN", (0.00, 0.84, 0.92, 1.0)),
    "socket": ("02_NEW_SOCKET_GREEN", (0.22, 0.70, 0.30, 1.0)),
    "cabinet": ("03_CABINET_POWER_AMBER", (1.00, 0.55, 0.00, 1.0)),
    "recheck": ("04_KITCHEN_SOCKET_RECHECK_RED", (1.00, 0.02, 0.02, 1.0)),
    "control": ("05_DOORWAY_CONTROL_MAGENTA", (0.84, 0.20, 0.52, 1.0)),
    "ap": ("06_BEDROOM_AP_BLUE", (0.13, 0.55, 0.90, 1.0)),
    "router": ("07_LIVING_STUDY_ROUTER_VIOLET", (0.44, 0.28, 0.91, 1.0)),
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


def object_ifc_definition_id(obj: bpy.types.Object) -> int | None:
    direct = obj.get("ifc_definition_id")
    if direct is not None:
        return int(direct)
    legacy_group = obj.get("BIMObjectProperties")
    if legacy_group is not None and "ifc_definition_id" in legacy_group:
        return int(legacy_group["ifc_definition_id"])
    return None


def hide_demolition_walls(root: Path) -> dict[str, int]:
    register_path = root / "pipeline/decisions/a102-demolition-review.csv"
    with register_path.open(encoding="utf-8-sig", newline="") as handle:
        global_ids = {
            row["candidate_global_id"]
            for row in csv.DictReader(handle)
            if row["ifc_status_candidate"] == "DEMOLISH"
        }
    model = ifcopenshell.open(str(root / "2504 GBTB Yanlord Zhuhai.ifc"))
    entity_ids = {model.by_guid(global_id).id() for global_id in global_ids}
    matches = []
    for obj in bpy.context.scene.objects:
        definition_id = object_ifc_definition_id(obj)
        if definition_id is None:
            continue
        if int(definition_id) in entity_ids:
            matches.append(obj)
            obj.hide_set(True)
            obj.hide_render = True
    demolition_reference = bpy.data.collections.get("A102_DEMOLITION_REFERENCE")
    if demolition_reference is not None:
        demolition_reference.hide_viewport = True
        demolition_reference.hide_render = True
    if len(matches) not in {0, 13}:
        raise RuntimeError(f"electrical review contains a partial DEMOLISH set: {len(matches)} of 13")
    visible_after = sum(not obj.hide_get() for obj in matches)
    if visible_after:
        raise RuntimeError(f"{visible_after} loaded DEMOLISH walls remain visible")
    return {"loaded": len(matches), "hidden": len(matches), "visible": visible_after}


def build_review() -> dict:
    root_path = project_root()
    report_path = root_path / "build/elec/elec-renovation-round1-candidate.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    control_network_path = root_path / "build/elec/elec-control-network-candidate.json"
    control_network = json.loads(control_network_path.read_text(encoding="utf-8"))
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
    for row in control_network["control_coordination_zones"]:
        point_marker(collections["control"], row, GROUPS["control"][1], 0.085)
    for row in control_network["network_coordination_zones"]:
        key = "ap" if row["network_role"] == "wireless_access_point" else "router"
        point_marker(collections[key], row, GROUPS[key][1], 0.085)

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
    demolition_visibility = hide_demolition_walls(root_path)
    configure_viewport()
    bpy.context.scene["renovation_elec_round1_legend"] = (
        "cyan=bedside lights; green=new sockets; amber=cabinet power zones; red=kitchen socket recheck; "
        "magenta=doorway control zones; blue=bedroom AP room zones; violet=living/study router room zones"
    )
    bpy.context.scene["renovation_elec_round1_writes_ifc"] = False
    return {
        "bedside": len(report["bedside_light_candidates"]),
        "socket": len(report["new_socket_candidates"]),
        "cabinet": len(report["cabinet_power_zones"]),
        "recheck": len(report["kitchen_socket_rechecks"]),
        "control": len(control_network["control_coordination_zones"]),
        "ap": sum(row["network_role"] == "wireless_access_point" for row in control_network["network_coordination_zones"]),
        "router": sum(row["network_role"] == "router_no_AP" for row in control_network["network_coordination_zones"]),
        "show_in_front": False,
        "solid_view": True,
        "formal_ifc_write": False,
        "demolition_walls": demolition_visibility,
        "blend_save": False,
    }


RESULT = build_review()
print(RESULT)
