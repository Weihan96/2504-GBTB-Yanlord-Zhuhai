#!/usr/bin/env python3
"""Audit the saved Falper Bonsai scene from inside Blender."""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import bpy
import ifcopenshell


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/falper-sorgente"
BLEND = PRODUCT_DIR / "Falper-Sorgente-WFB-bonsai-review.blend"
ISOLATED_IFC = PRODUCT_DIR / "Falper-Sorgente-WFB-bonsai-isolated.ifc"
DEFAULT_OUTPUT = PRODUCT_DIR / "bonsai-blend-audit.json"
GLOBAL_ID = "350tdaubr8QP3Cu2YMQZIN"
EXPECTED_CAMERAS = {
    "FALPER_CAM_PLAN",
    "FALPER_CAM_FRONT",
    "FALPER_CAM_SIDE",
    "FALPER_CAM_ISO",
}
RENDERS = {
    "plan": PRODUCT_DIR / "bonsai-camera-plan.png",
    "front_elevation": PRODUCT_DIR / "bonsai-camera-front-elevation.png",
    "side_elevation": PRODUCT_DIR / "bonsai-camera-side-elevation.png",
    "isometric": PRODUCT_DIR / "bonsai-camera-iso.png",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parse_args() -> argparse.Namespace:
    arguments = sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else []
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args(arguments)


def main() -> None:
    args = parse_args()
    meshes = [obj for obj in bpy.context.scene.objects if obj.type == "MESH"]
    cameras = [obj for obj in bpy.context.scene.objects if obj.type == "CAMERA"]
    if len(meshes) != 1:
        raise RuntimeError(f"expected one IFC geometry object, found {len(meshes)}")
    camera_names = {camera.name for camera in cameras}
    if camera_names != EXPECTED_CAMERAS:
        raise RuntimeError(f"unexpected saved camera set: {sorted(camera_names)}")

    obj = meshes[0]
    product_id = int(obj.BIMObjectProperties.ifc_definition_id)
    representation_id = int(obj.data.BIMMeshProperties.ifc_definition_id)
    model = ifcopenshell.open(ISOLATED_IFC)
    product = model.by_id(product_id)
    representation = model.by_id(representation_id)
    context_identifier = representation.ContextOfItems.ContextIdentifier
    if product.GlobalId != GLOBAL_ID:
        raise RuntimeError(f"unexpected product GlobalId: {product.GlobalId}")
    if context_identifier != "Body":
        raise RuntimeError(f"unexpected active representation: {context_identifier}")

    payload = {
        "schema_version": 1,
        "generator": relative(Path(__file__)),
        "blend": {"path": relative(BLEND), "sha256": sha256(BLEND)},
        "isolated_ifc": {
            "path": relative(ISOLATED_IFC),
            "sha256": sha256(ISOLATED_IFC),
        },
        "geometry_product_count": len(meshes),
        "whole_model_render": False,
        "product": {
            "object_name": obj.name,
            "ifc_definition_id": product_id,
            "ifc_class": product.is_a(),
            "global_id": product.GlobalId,
        },
        "representation": {
            "ifc_definition_id": representation_id,
            "context_identifier": context_identifier,
        },
        "cameras": [
            {
                "name": camera.name,
                "type": camera.data.type,
                "orthographic_scale_m": round(float(camera.data.ortho_scale), 6),
            }
            for camera in sorted(cameras, key=lambda item: item.name)
        ],
        "saved_camera_count": len(cameras),
        "render_engine": bpy.context.scene.render.engine,
        "resolution_px": [
            bpy.context.scene.render.resolution_x,
            bpy.context.scene.render.resolution_y,
        ],
        "renders": {
            name: {"path": relative(path), "sha256": sha256(path)}
            for name, path in RENDERS.items()
        },
        "pass": True,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print("FALPER_BONSAI_BLEND_AUDIT_PASS " + json.dumps(payload, separators=(",", ":")))


if __name__ == "__main__":
    main()
