#!/usr/bin/env python3
"""Render the approved Duofix type in its main-bathroom project context."""

from __future__ import annotations

import contextlib
import hashlib
import json
from collections import Counter
from pathlib import Path

import addon_utils
import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-duofix-sigma-224-212"
DERIVED_IFC = PRODUCT_DIR / "Geberit-Duofix-Sigma-224-212-derived-drawing.ifc"
DERIVED_REPORT = PRODUCT_DIR / "Geberit-Duofix-Sigma-224-212-derived-drawing-report.json"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
TARGET_GLOBAL_ID = "3pvAlH5C14v8uVEJ1LmK8M"
ROOM_GLOBAL_ID = "3a4COIs5X7lgDirMBDT4Vs"
ROOM_NAME = "主卫湿区"
CAMERA_PREFIX = "GEBERIT_DUOFIX_224_212_MAIN_BATH"
BLEND = PRODUCT_DIR / "Geberit-Duofix-Sigma-224-212-main-bathroom-context.blend"
MANIFEST = PRODUCT_DIR / "main-bathroom-bonsai-render-manifest.json"
RENDERS = {
    "PLAN": PRODUCT_DIR / "main-bathroom-bonsai-plan.png",
    "FRONT": PRODUCT_DIR / "main-bathroom-bonsai-front-elevation.png",
    "SIDE": PRODUCT_DIR / "main-bathroom-bonsai-side-elevation.png",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def object_bounds(obj):
    corners = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    minimum = Vector(tuple(min(point[axis] for point in corners) for axis in range(3)))
    maximum = Vector(tuple(max(point[axis] for point in corners) for axis in range(3)))
    return minimum, maximum


def overlaps(minimum, maximum, crop_minimum, crop_maximum):
    return all(maximum[axis] >= crop_minimum[axis] and minimum[axis] <= crop_maximum[axis] for axis in range(3))


def ifc_entity_for_object(tool, obj):
    with contextlib.suppress(Exception):
        return tool.Ifc.get_entity(obj)
    return None


def camera_axis(obj, local_axis):
    return (obj.matrix_world.to_3x3() @ Vector(local_axis)).normalized()


def projected_extent(bounds, center, axis):
    return max(abs((corner - center).dot(axis)) for corner in bounds) * 2.0


if sha256(FORMAL_IFC) != FORMAL_SHA256:
    raise RuntimeError("formal IFC hash mismatch before main-bathroom render")
report = json.loads(DERIVED_REPORT.read_text(encoding="utf-8"))
if not report.get("pass") or report.get("formal_ifc_sha256") != FORMAL_SHA256:
    raise RuntimeError("derived IFC report gate failed")
if report.get("derived_ifc_sha256") != sha256(DERIVED_IFC):
    raise RuntimeError("derived IFC hash does not match its verified report")

with contextlib.suppress(Exception):
    addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
addon_utils.enable("bonsai_bridge", default_set=False, persistent=False)
from bonsai import tool

tool.IfcGit.load_project(str(DERIVED_IFC))
ifc = tool.Ifc.get()
if not ifc:
    raise RuntimeError("Bonsai did not load the derived IFC")
target_entity = ifc.by_guid(TARGET_GLOBAL_ID)
room_entity = ifc.by_guid(ROOM_GLOBAL_ID)
if target_entity is None or room_entity is None or room_entity.LongName != ROOM_NAME:
    raise RuntimeError("main-bathroom target or room identity is missing")
target = tool.Ifc.get_object(target_entity)
room = tool.Ifc.get_object(room_entity)
if target is None or room is None:
    raise RuntimeError("main-bathroom target or room is not loaded in Blender")

representative = ifc.by_guid(report["representative_global_id"])
representation_identifiers = {
    representation.RepresentationIdentifier
    for representation in representative.Representation.Representations
}
expected_identifiers = set(report["representations"].values())
if not expected_identifiers.issubset(representation_identifiers):
    raise RuntimeError("approved derived drawing representations are not persisted in the loaded IFC")

room_minimum, room_maximum = object_bounds(room)
crop_minimum = room_minimum + Vector((-0.35, -0.35, -0.20))
crop_maximum = room_maximum + Vector((0.35, 0.35, 0.20))
allowed_classes = {
    "IfcBeam",
    "IfcBuildingElementProxy",
    "IfcColumn",
    "IfcCovering",
    "IfcCurtainWall",
    "IfcDoor",
    "IfcFlowFitting",
    "IfcFlowSegment",
    "IfcFlowTerminal",
    "IfcFurniture",
    "IfcMember",
    "IfcPlate",
    "IfcRailing",
    "IfcSanitaryTerminal",
    "IfcSlab",
    "IfcWall",
    "IfcWallStandardCase",
    "IfcWindow",
}
visible = []
visible_classes = Counter()
entity_classes = {}
visible_bounds = {}
for obj in bpy.context.scene.objects:
    if obj.type not in {"MESH", "CURVE"}:
        obj.hide_render = True
        continue
    entity = ifc_entity_for_object(tool, obj)
    if entity is None or entity.is_a() not in allowed_classes:
        obj.hide_render = True
        continue
    minimum, maximum = object_bounds(obj)
    obj.hide_render = not overlaps(minimum, maximum, crop_minimum, crop_maximum)
    if not obj.hide_render:
        visible.append(obj)
        visible_classes[entity.is_a()] += 1
        entity_classes[obj] = entity.is_a()
        visible_bounds[obj] = (minimum, maximum)
        obj.color = (0.78, 0.82, 0.86, 1.0)

target.hide_render = False
target.color = (0.01, 0.24, 0.92, 1.0)
if target not in visible:
    visible.append(target)
    visible_classes[target_entity.is_a()] += 1
    entity_classes[target] = target_entity.is_a()
    visible_bounds[target] = object_bounds(target)

scene = bpy.context.scene
for existing in list(bpy.data.objects):
    if existing.type == "CAMERA":
        bpy.data.objects.remove(existing, do_unlink=True)
scene.render.engine = "BLENDER_WORKBENCH"
scene.render.resolution_x = 1500
scene.render.resolution_y = 1100
scene.render.resolution_percentage = 100
scene.render.image_settings.file_format = "PNG"
scene.render.film_transparent = False
scene.display.shading.light = "STUDIO"
scene.display.shading.color_type = "OBJECT"
scene.display.shading.show_shadows = False
scene.display.shading.show_cavity = True
scene.display.shading.cavity_type = "WORLD"
if hasattr(scene.display.shading, "show_outline"):
    scene.display.shading.show_outline = True
scene.display.shading.show_xray = True
scene.display.shading.xray_alpha = 0.48
scene.display.shading.background_type = "VIEWPORT"
scene.display.shading.background_color = (0.985, 0.988, 0.992)

target_minimum, target_maximum = object_bounds(target)
target_center = (target_minimum + target_maximum) / 2.0
axis_x = camera_axis(target, (1.0, 0.0, 0.0))
axis_y = camera_axis(target, (0.0, 1.0, 0.0))
axis_z = camera_axis(target, (0.0, 0.0, 1.0))
room_center = (crop_minimum + crop_maximum) / 2.0
context_bounds = [
    Vector((x, y, z))
    for x in (crop_minimum.x, crop_maximum.x)
    for y in (crop_minimum.y, crop_maximum.y)
    for z in (crop_minimum.z, crop_maximum.z)
]
distance = 12.0
jobs = [
    ("PLAN", axis_z),
    ("FRONT", -axis_y),
    ("SIDE", -axis_x),
]
render_records = []
for label, direction in jobs:
    for obj in visible:
        obj.hide_render = False
    foreground_hidden = []
    if label in {"FRONT", "SIDE"}:
        structural_classes = {
            "IfcCovering", "IfcCurtainWall", "IfcDoor", "IfcSlab",
            "IfcWall", "IfcWallStandardCase", "IfcWindow",
        }
        target_depth = target_center.dot(direction)
        for obj in visible:
            if obj is target or entity_classes.get(obj) not in structural_classes:
                continue
            minimum, maximum = visible_bounds[obj]
            center = (minimum + maximum) / 2.0
            if center.dot(direction) > target_depth:
                obj.hide_render = True
                foreground_hidden.append(obj)
    camera_data = bpy.data.cameras.new(f"{CAMERA_PREFIX}_{label}")
    camera = bpy.data.objects.new(f"{CAMERA_PREFIX}_{label}", camera_data)
    scene.collection.objects.link(camera)
    camera.location = room_center + direction * distance
    camera.rotation_euler = (room_center - camera.location).to_track_quat("-Z", "Y").to_euler()
    camera.data.type = "ORTHO"
    bpy.context.view_layer.update()
    camera_right = camera.matrix_world.to_3x3() @ Vector((1.0, 0.0, 0.0))
    camera_up = camera.matrix_world.to_3x3() @ Vector((0.0, 1.0, 0.0))
    projected_width = projected_extent(context_bounds, room_center, camera_right)
    projected_height = projected_extent(context_bounds, room_center, camera_up)
    aspect = scene.render.resolution_x / scene.render.resolution_y
    camera.data.ortho_scale = max(projected_height, projected_width / aspect) * 1.04
    scene.camera = camera
    output = RENDERS[label]
    scene.render.filepath = str(output)
    bpy.ops.render.render(write_still=True)
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"missing main-bathroom render: {output}")
    render_records.append({
        "view": label.lower(),
        "path": relative(output),
        "sha256": sha256(output),
        "bytes": output.stat().st_size,
        "camera": camera.name,
        "camera_type": camera.data.type,
        "view_direction_world": [round(value, 6) for value in direction],
        "camera_axis_basis": "main-bathroom Duofix instance local axes",
        "resolution": [scene.render.resolution_x, scene.render.resolution_y],
        "render_operation": "bpy.ops.render.render(write_still=True)",
        "foreground_structural_objects_hidden": len(foreground_hidden),
    })

for obj in bpy.context.scene.objects:
    obj.select_set(False)
target.select_set(True)
bpy.context.view_layer.objects.active = target
bpy.ops.wm.save_as_mainfile(filepath=str(BLEND))
if not BLEND.is_file():
    raise RuntimeError("main-bathroom Bonsai blend was not saved")
if sha256(FORMAL_IFC) != FORMAL_SHA256:
    raise RuntimeError("formal IFC changed during main-bathroom render")

manifest = {
    "schema_version": 1,
    "generator": "pipeline/scripts/render_geberit_duofix_sigma_224_212_main_bathroom.py",
    "mode": "actual_bonsai_derived_ifc_main_bathroom_context_camera_render",
    "formal_ifc": relative(FORMAL_IFC),
    "formal_ifc_sha256": FORMAL_SHA256,
    "formal_ifc_bytes_unchanged": True,
    "derived_ifc": relative(DERIVED_IFC),
    "derived_ifc_sha256": sha256(DERIVED_IFC),
    "derived_report": relative(DERIVED_REPORT),
    "derived_report_sha256": sha256(DERIVED_REPORT),
    "persisted_drawing_representation_identifiers": sorted(expected_identifiers),
    "target_global_id": TARGET_GLOBAL_ID,
    "target_ifc_class": target_entity.is_a(),
    "target_type_name": next(relation.RelatingType.Name for relation in target_entity.IsTypedBy),
    "room_global_id": ROOM_GLOBAL_ID,
    "room_name": ROOM_NAME,
    "room_bbox_m": {
        "minimum": [round(value, 6) for value in room_minimum],
        "maximum": [round(value, 6) for value in room_maximum],
    },
    "project_context_retained": True,
    "visible_context_object_count": len(visible),
    "visible_context_ifc_classes": dict(sorted(visible_classes.items())),
    "target_highlight_color": "blue",
    "context_color": "light_grey",
    "saved_bonsai_session": {
        "path": relative(BLEND),
        "sha256": sha256(BLEND),
        "bytes": BLEND.stat().st_size,
        "camera_count": len(render_records),
    },
    "renders": render_records,
    "blender_version": bpy.app.version_string,
    "ifc_schema": ifc.schema,
    "pass": True,
}
MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps({
    "manifest": relative(MANIFEST),
    "renders": [record["path"] for record in render_records],
    "visible_context_object_count": len(visible),
    "pass": True,
}, indent=2, ensure_ascii=False))
