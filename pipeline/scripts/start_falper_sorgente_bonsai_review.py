#!/usr/bin/env python3
"""Start an isolated interactive Bonsai review session for Falper WFB."""

import contextlib
from pathlib import Path

import addon_utils
import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    ROOT
    / "output/review/highpoly-types/falper-sorgente/"
    "Falper-Sorgente-WFB-bonsai-isolated.ifc"
)
BLEND = (
    ROOT
    / "output/review/highpoly-types/falper-sorgente/"
    "Falper-Sorgente-WFB-bonsai-review.blend"
)
GLOBAL_ID = "350tdaubr8QP3Cu2YMQZIN"
CAMERA_RENDERS = {
    "PLAN": "bonsai-camera-plan.png",
    "FRONT": "bonsai-camera-front-elevation.png",
    "SIDE": "bonsai-camera-side-elevation.png",
    "ISO": "bonsai-camera-iso.png",
}


def render_body_cameras(obj: bpy.types.Object) -> None:
    for existing in list(bpy.data.objects):
        if existing.name.startswith("FALPER_CAM_"):
            bpy.data.objects.remove(existing, do_unlink=True)

    scene = bpy.context.scene
    scene.render.engine = "BLENDER_WORKBENCH"
    scene.render.resolution_x = 1200
    scene.render.resolution_y = 1200
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"
    scene.render.film_transparent = False
    scene.display.shading.light = "STUDIO"
    scene.display.shading.color_type = "SINGLE"
    scene.display.shading.single_color = (0.78, 0.80, 0.83)
    scene.display.shading.show_shadows = True
    scene.display.shading.show_cavity = True
    scene.display.shading.cavity_type = "WORLD"
    scene.display.shading.background_type = "VIEWPORT"
    scene.display.shading.background_color = (0.035, 0.040, 0.050)

    corners = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]
    minimum = Vector(
        (
            min(point.x for point in corners),
            min(point.y for point in corners),
            min(point.z for point in corners),
        )
    )
    maximum = Vector(
        (
            max(point.x for point in corners),
            max(point.y for point in corners),
            max(point.z for point in corners),
        )
    )
    center = (minimum + maximum) * 0.5
    extent = maximum - minimum
    jobs = [
        ("PLAN", center + Vector((0.0, 0.0, 3.0)), max(extent.x, extent.y) * 1.20),
        ("FRONT", center + Vector((0.0, -3.0, 0.0)), max(extent.x, extent.z) * 1.18),
        ("SIDE", center + Vector((3.0, 0.0, 0.0)), max(extent.y, extent.z) * 1.18),
        ("ISO", center + Vector((1.8, -1.8, 1.5)), max(extent) * 1.42),
    ]
    for label, location, scale in jobs:
        camera_data = bpy.data.cameras.new(f"FALPER_CAM_{label}")
        camera = bpy.data.objects.new(f"FALPER_CAM_{label}", camera_data)
        scene.collection.objects.link(camera)
        camera.location = location
        camera.rotation_euler = (center - location).to_track_quat("-Z", "Y").to_euler()
        camera.data.type = "ORTHO"
        camera.data.ortho_scale = scale
        scene.camera = camera
        scene.render.filepath = str(BLEND.parent / CAMERA_RENDERS[label])
        bpy.ops.render.render(write_still=True)


with contextlib.suppress(Exception):
    addon_utils.disable(
        "bl_ext.user_default.project_control", default_set=False, handle_error=None
    )
addon_utils.enable("bonsai_bridge", default_set=False, persistent=False)
import bonsai.tool as tool

tool.IfcGit.load_project(str(SOURCE))
ifc = tool.Ifc.get()
product = ifc.by_guid(GLOBAL_ID)
if product is None:
    raise RuntimeError("Falper representative is missing from Bonsai session")
obj = tool.Ifc.get_object(product)
if obj is None:
    raise RuntimeError("Bonsai did not create the Falper Blender object")
for other in list(bpy.context.scene.objects):
    if other != obj and other.type in {"MESH", "CURVE", "SURFACE", "META", "FONT"}:
        bpy.data.objects.remove(other, do_unlink=True)
bpy.ops.object.select_all(action="DESELECT")
obj.hide_set(False)
obj.select_set(True)
bpy.context.view_layer.objects.active = obj
representation_id = obj.data.BIMMeshProperties.ifc_definition_id
representation = ifc.by_id(representation_id)
if representation.ContextOfItems.ContextIdentifier != "Body":
    raise RuntimeError("camera renders must use the actual IFC Body representation")
render_body_cameras(obj)
bpy.ops.wm.save_as_mainfile(filepath=str(BLEND))
preferences = bpy.context.preferences.addons["bonsai_bridge"].preferences
preferences.allow_edits = True
bpy.ops.bonsai_mcp.start_bridge()
print(f"FALPER_BONSAI_READY {obj.name} {SOURCE}")
