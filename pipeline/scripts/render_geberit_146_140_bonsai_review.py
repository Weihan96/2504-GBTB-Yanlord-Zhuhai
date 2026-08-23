#!/usr/bin/env python3
"""Load the isolated Geberit IFC in Bonsai, add cameras, and execute Blender renders."""

import contextlib
import hashlib
import json
from pathlib import Path

import addon_utils
import bpy
from mathutils import Vector


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-146-140"
SOURCE = PRODUCT_DIR / "Geberit-146-140-bonsai-isolated.ifc"
BLEND = PRODUCT_DIR / "Geberit-146-140-bonsai-review.blend"
MANIFEST = PRODUCT_DIR / "bonsai-review-manifest.json"
GLOBAL_ID = "1rhZG98PPCSxaLeMFLTYb9"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
RENDERS = {
    "PLAN": "bonsai-camera-plan.png",
    "FRONT": "bonsai-camera-front-elevation.png",
    "SIDE": "bonsai-camera-side-elevation.png",
    "ISO": "bonsai-camera-iso.png",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def render_cameras(obj):
    scene = bpy.context.scene
    for existing in list(bpy.data.objects):
        if existing.type == "CAMERA":
            bpy.data.objects.remove(existing, do_unlink=True)
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
    minimum = Vector(tuple(min(point[axis] for point in corners) for axis in range(3)))
    maximum = Vector(tuple(max(point[axis] for point in corners) for axis in range(3)))
    center = (minimum + maximum) * 0.5
    extent = maximum - minimum
    jobs = [
        ("PLAN", center + Vector((0.0, 0.0, 3.0)), max(extent.x, extent.y) * 1.18),
        # The representative is placed with its product depth on world X. Its
        # rounded user-facing end is at minimum X, so these cameras preserve
        # the manufacturer's A/L front/side semantics in the saved session.
        ("FRONT", center + Vector((-3.0, 0.0, 0.0)), max(extent.y, extent.z) * 1.18),
        ("SIDE", center + Vector((0.0, -3.0, 0.0)), max(extent.x, extent.z) * 1.18),
        ("ISO", center + Vector((-1.8, -1.8, 1.5)), max(extent) * 1.42),
    ]
    results = []
    for label, location, scale in jobs:
        camera_data = bpy.data.cameras.new(f"GEBERIT_146_140_CAM_{label}")
        camera = bpy.data.objects.new(f"GEBERIT_146_140_CAM_{label}", camera_data)
        scene.collection.objects.link(camera)
        camera.location = location
        camera.rotation_euler = (center - location).to_track_quat("-Z", "Y").to_euler()
        camera.data.type = "ORTHO"
        camera.data.ortho_scale = scale
        scene.camera = camera
        target = PRODUCT_DIR / RENDERS[label]
        scene.render.filepath = str(target)
        bpy.ops.render.render(write_still=True)
        if not target.is_file():
            raise RuntimeError(f"Bonsai camera render missing: {target}")
        results.append({
            "view": label.lower(),
            "camera": camera.name,
            "camera_type": camera.data.type,
            "render_operation": "bpy.ops.render.render(write_still=True)",
            "path": relative(target),
            "sha256": sha256(target),
            "resolution": [scene.render.resolution_x, scene.render.resolution_y],
        })
    return results, minimum, maximum


if sha256(FORMAL_IFC) != FORMAL_SHA256:
    raise RuntimeError("formal IFC hash mismatch before Bonsai render")
with contextlib.suppress(Exception):
    addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
addon_utils.enable("bonsai_bridge", default_set=False, persistent=False)
import bonsai.tool as tool

tool.IfcGit.load_project(str(SOURCE))
ifc = tool.Ifc.get()
product = ifc.by_guid(GLOBAL_ID)
if product is None:
    raise RuntimeError("Geberit representative is missing from Bonsai session")
obj = tool.Ifc.get_object(product)
if obj is None:
    raise RuntimeError("Bonsai did not create the Geberit Blender object")
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
renders, minimum, maximum = render_cameras(obj)
bpy.ops.wm.save_as_mainfile(filepath=str(BLEND))
camera_names = sorted(item.name for item in bpy.data.objects if item.type == "CAMERA")
if len(camera_names) != 4:
    raise RuntimeError("saved Bonsai review must contain exactly four cameras")
if sha256(FORMAL_IFC) != FORMAL_SHA256:
    raise RuntimeError("formal IFC bytes changed during Bonsai render")
payload = {
    "schema_version": 1,
    "generator": "pipeline/scripts/render_geberit_146_140_bonsai_review.py",
    "mode": "actual_bonsai_ifc_body_camera_render",
    "representative_global_id": GLOBAL_ID,
    "geometry_product_count": 1,
    "whole_model_render": False,
    "formal_ifc": relative(FORMAL_IFC),
    "formal_ifc_sha256": FORMAL_SHA256,
    "formal_ifc_bytes_unchanged": True,
    "isolated_ifc": relative(SOURCE),
    "isolated_ifc_sha256": sha256(SOURCE),
    "bonsai_session": {
        "path": relative(BLEND),
        "sha256": sha256(BLEND),
        "saved_active_representation": "Body",
        "saved_camera_count": len(camera_names),
        "saved_camera_names": camera_names,
        "ifc_definition_id": representation_id,
        "ifc_context_identifier": representation.ContextOfItems.ContextIdentifier,
        "object_name": obj.name,
        "object_bounds_m": {
            "minimum": [round(value, 9) for value in minimum],
            "maximum": [round(value, 9) for value in maximum],
        },
    },
    "renders": renders,
    "pass": True,
}
MANIFEST.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"manifest": relative(MANIFEST), "renders": [item["path"] for item in renders], "pass": True}, indent=2))
