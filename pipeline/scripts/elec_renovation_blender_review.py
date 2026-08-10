"""Create a depth-aware Blender review for first-round renovation electrical demands."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import bpy
import ifcopenshell
from mathutils import Vector

import bonsai.tool as bonsai_tool


ROOT_COLLECTION = "RENOVATION_ELEC_ROUND1"
REVIEW_LABEL_COLLECTION = "08_CURRENT_REVIEW_CALLOUTS"
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
REVIEW_LABEL_Z = 3.15
FOCUS_LABELS = {
    "NS-01": ("E303  NS-01  XY / Z=650", (0.18, 0.18)),
    "NS-02": ("E303  NS-02  XY / Z=300", (0.18, -0.20)),
    "A106-AP-R09": ("E304  AP-R09  XY / Z=2720", (0.18, 0.18)),
    "A106-AP-R14": ("E304  AP-R14  XY / Z=2720", (0.18, 0.18)),
    "NET-ROUTER-LIVING-STUDY": ("E304  ROUTER @ ENTRY CABINET / Z=TBD", (0.18, 0.18)),
}
CONTROL_JAMB_CLEARANCE_M = 0.15


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


def mark_review_overlay(obj: bpy.types.Object, candidate_id: str) -> None:
    obj.show_in_front = False
    obj["candidate_id"] = candidate_id
    obj["review_overlay_only"] = True
    obj["automatic_ifc_write_allowed"] = False


def router_cabinet_review_zone(collection: bpy.types.Collection, record: dict, colour) -> bpy.types.Object:
    minimum_x, minimum_y, maximum_x, maximum_y = [float(value) / 1000.0 for value in record["cabinet_bbox_ifc_mm"]]
    review_z = float(record["review_overlay_position_mm"][2]) / 1000.0
    bpy.ops.mesh.primitive_cube_add(
        size=1.0,
        location=((minimum_x + maximum_x) / 2.0, (minimum_y + maximum_y) / 2.0, review_z),
    )
    obj = bpy.context.object
    obj.name = f"PLAN-ZONE::{record['candidate_id']}"
    obj.dimensions = (maximum_x - minimum_x, maximum_y - minimum_y, 0.025)
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    move_to_collection(obj, collection)
    obj.color = colour
    mark_review_overlay(obj, record["candidate_id"])
    obj["coordinate_status"] = record["coordinate_status"]
    obj["installation_z_mm"] = "TBD"
    obj["weak_current_box_bottom_aff_mm"] = record["weak_current_box_bottom_aff_mm"]
    obj["evidence_note"] = "plan-only cabinet bay; weak-current-box H+350 is not the router installation height"
    return obj


def focus_callout(
    collection: bpy.types.Collection,
    candidate_id: str,
    label: str,
    position_mm,
    colour,
    offset,
) -> list[bpy.types.Object]:
    offset_x, offset_y = offset
    x, y, _z = (float(value) / 1000.0 for value in position_mm)
    label_x = x + offset_x
    label_y = y + offset_y

    bpy.ops.mesh.primitive_cylinder_add(
        vertices=32,
        radius=0.105,
        depth=0.025,
        location=(x, y, REVIEW_LABEL_Z),
    )
    anchor = bpy.context.object
    anchor.name = f"CALLOUT::{candidate_id}"
    move_to_collection(anchor, collection)
    anchor.color = colour
    mark_review_overlay(anchor, candidate_id)

    curve_data = bpy.data.curves.new(f"LEADER::{candidate_id}", "CURVE")
    curve_data.dimensions = "3D"
    curve_data.bevel_depth = 0.012
    curve_data.bevel_resolution = 2
    spline = curve_data.splines.new("POLY")
    spline.points.add(1)
    spline.points[0].co = (x, y, REVIEW_LABEL_Z + 0.02, 1.0)
    spline.points[1].co = (label_x, label_y, REVIEW_LABEL_Z + 0.02, 1.0)
    leader = bpy.data.objects.new(f"LEADER::{candidate_id}", curve_data)
    collection.objects.link(leader)
    leader.color = colour
    mark_review_overlay(leader, candidate_id)

    font_data = bpy.data.curves.new(f"LABEL::{candidate_id}", "FONT")
    font_data.body = label
    font_data.align_x = "RIGHT" if offset_x < 0 else "LEFT"
    font_data.align_y = "CENTER"
    font_data.size = 0.13
    font_data.extrude = 0.003
    text = bpy.data.objects.new(f"LABEL::{candidate_id}", font_data)
    text_padding = -0.04 if offset_x < 0 else 0.04
    text.location = (label_x + text_padding, label_y, REVIEW_LABEL_Z + 0.035)
    collection.objects.link(text)
    text.color = colour
    mark_review_overlay(text, candidate_id)
    return [anchor, leader, text]


def control_wall_side_options(
    collection: bpy.types.Collection,
    record: dict,
    colour,
) -> list[dict]:
    model = bonsai_tool.Ifc.get()
    door = model.by_guid(record["source_global_ids"][0])
    door_object = bonsai_tool.Ifc.get_object(door)
    if door_object is None:
        raise RuntimeError(f"door object is not loaded for {record['candidate_id']}")
    local_width_m = max(corner[0] for corner in door_object.bound_box) - min(
        corner[0] for corner in door_object.bound_box
    )
    if not 0.7 <= local_width_m <= 1.2:
        raise RuntimeError(f"unexpected door opening width for {record['candidate_id']}: {local_width_m:.3f} m")
    axis = door_object.matrix_world.to_quaternion() @ Vector((1.0, 0.0, 0.0))
    axis.z = 0.0
    axis.normalize()
    offset_m = local_width_m / 2.0 + CONTROL_JAMB_CLEARANCE_M
    centre = Vector(tuple(float(value) / 1000.0 for value in record["position_mm"]))
    options = []
    for suffix, direction in (("A", -1.0), ("B", 1.0)):
        position = centre + axis * offset_m * direction
        option_id = f"{record['candidate_id']}-{suffix}"
        option_record = {
            "candidate_id": option_id,
            "kind": "doorway_control_wall_side_option",
            "room_reference": record["room_reference"],
            "position_mm": [position.x * 1000.0, position.y * 1000.0, centre.z * 1000.0],
        }
        marker = point_marker(collection, option_record, colour, 0.055)
        marker["source_candidate_id"] = record["candidate_id"]
        marker["coordinate_status"] = "wall_side_option_not_final"
        marker["jamb_clearance_mm"] = CONTROL_JAMB_CLEARANCE_M * 1000.0
        options.append(option_record)
    return options


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
            space.region_3d.view_distance = 12.5


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
    callout_collection = ensure_collection(REVIEW_LABEL_COLLECTION, root)
    clear_collection(callout_collection)
    callout_collection.hide_viewport = False
    callout_collection.hide_render = False

    focus_records = []
    control_options = []

    for row in report["bedside_light_candidates"]:
        point_marker(collections["bedside"], row, GROUPS["bedside"][1], 0.075)
    for row in report["new_socket_candidates"]:
        point_marker(collections["socket"], row, GROUPS["socket"][1], 0.075)
        if row["candidate_id"] in FOCUS_LABELS:
            label, offset = FOCUS_LABELS[row["candidate_id"]]
            focus_records.append((row["candidate_id"], label, row["position_mm"], GROUPS["socket"][1], offset))
    for row in report["cabinet_power_zones"]:
        cabinet_marker(collections["cabinet"], row, GROUPS["cabinet"][1])
    for row in report["kitchen_socket_rechecks"]:
        point_marker(collections["recheck"], row, GROUPS["recheck"][1], 0.060)
    for row in control_network["control_coordination_zones"]:
        zone = point_marker(collections["control"], row, GROUPS["control"][1], 0.085)
        zone.name = f"ZONE::{row['candidate_id']}"
        zone["coordinate_status"] = "doorway_zone_only_not_wall_position"
        options = control_wall_side_options(collections["control"], row, GROUPS["control"][1])
        for index, option in enumerate(options):
            suffix = option["candidate_id"].rsplit("-", 1)[1]
            short_name = "ENTRY" if row["candidate_id"] == "CTRL-ENTRY" else "MASTER"
            label = f"E302  {short_name}-{suffix}  WALL SIDE / Z=1300"
            horizontal_offset = -0.18 if row["candidate_id"] == "CTRL-ENTRY" else 0.18
            label_offset = (horizontal_offset, -0.12 if index == 0 else 0.12)
            focus_records.append(
                (option["candidate_id"], label, option["position_mm"], GROUPS["control"][1], label_offset)
            )
            control_options.append(option)
    for row in control_network["network_coordination_zones"]:
        key = "ap" if row["network_role"] == "wireless_access_point" else "router"
        marker_position = row["position_mm"] if key == "ap" else row["review_overlay_position_mm"]
        if key == "ap":
            point_marker(collections[key], row, GROUPS[key][1], 0.085)
        else:
            router_cabinet_review_zone(collections[key], row, GROUPS[key][1])
        if row["candidate_id"] in FOCUS_LABELS:
            label, offset = FOCUS_LABELS[row["candidate_id"]]
            focus_records.append((row["candidate_id"], label, marker_position, GROUPS[key][1], offset))

    for candidate_id, label, position_mm, colour, offset in focus_records:
        focus_callout(callout_collection, candidate_id, label, position_mm, colour, offset)

    for key in {"bedside", "cabinet", "recheck"}:
        collections[key].hide_viewport = True
        collections[key].hide_render = True

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
        "magenta=doorway control zones; blue=bedroom AP mechanical positions; violet=entry-cabinet router plan evidence"
    )
    bpy.context.scene["renovation_elec_round1_writes_ifc"] = False
    bpy.context.scene["renovation_elec_current_review"] = (
        "E302 doorway controls show A/B wall-side options at Z=1300; E303 sockets display candidate Z; "
        "E304 APs display candidate Z=2720; the router is a plan-only entry-cabinet evidence zone with installation Z TBD. "
        "The weak-current-box H+350 datum is not used as router height. Callouts are review overlays at Z=3.15 m; "
        "AP source markers retain their true installation depth."
    )
    return {
        "bedside": len(report["bedside_light_candidates"]),
        "socket": len(report["new_socket_candidates"]),
        "cabinet": len(report["cabinet_power_zones"]),
        "recheck": len(report["kitchen_socket_rechecks"]),
        "control": len(control_network["control_coordination_zones"]),
        "ap": sum(row["network_role"] == "wireless_access_point" for row in control_network["network_coordination_zones"]),
        "router": sum(row["network_role"] == "router_no_AP" for row in control_network["network_coordination_zones"]),
        "focus_callouts": len(focus_records),
        "focus_ids": sorted(candidate_id for candidate_id, _label, _position, _colour, _offset in focus_records),
        "control_wall_side_options": len(control_options),
        "default_hidden_reference_groups": ["bedside", "cabinet", "recheck"],
        "review_label_z_m": REVIEW_LABEL_Z,
        "show_in_front": False,
        "solid_view": True,
        "formal_ifc_write": False,
        "demolition_walls": demolition_visibility,
        "blend_save": False,
    }


RESULT = build_review()
print(RESULT)
