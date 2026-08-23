#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for the sxb010 kitchen proxy."""

import os
import runpy
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/sxb010",
    "HIGHPOLY_ISOLATED_IFC": "sxb010-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "sxb010-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "1O9JRXCI56VRUbpuLJy86Z",
    "HIGHPOLY_CAMERA_PREFIX": "SXB010",
    "HIGHPOLY_PRODUCT_LABEL": "sxb010 kitchen proxy",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_sxb010_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")

manifest_path = ROOT / "output/review/highpoly-types/sxb010/bonsai-review-manifest.json"
payload = json.loads(manifest_path.read_text(encoding="utf-8"))
product_dir = manifest_path.parent
semantic_aliases = {
    "plan": ("bonsai-camera-front-elevation.png", "bonsai-camera-semantic-plan.png"),
    "front": ("bonsai-camera-plan.png", "bonsai-camera-semantic-front-elevation.png"),
    "side": ("bonsai-camera-side-elevation.png", "bonsai-camera-semantic-side-elevation.png"),
}
for source_name, semantic_name in semantic_aliases.values():
    shutil.copy2(product_dir / source_name, product_dir / semantic_name)
payload["semantic_view_mapping"] = {
    "plan": {
        "render_view_record": "front",
        "source_camera_path": "output/review/highpoly-types/sxb010/bonsai-camera-front-elevation.png",
        "path": "output/review/highpoly-types/sxb010/bonsai-camera-semantic-plan.png",
        "project_axes": [0, 2],
    },
    "front": {
        "render_view_record": "plan",
        "source_camera_path": "output/review/highpoly-types/sxb010/bonsai-camera-plan.png",
        "path": "output/review/highpoly-types/sxb010/bonsai-camera-semantic-front-elevation.png",
        "project_axes": [0, 1],
    },
    "side": {
        "render_view_record": "side",
        "source_camera_path": "output/review/highpoly-types/sxb010/bonsai-camera-side-elevation.png",
        "path": "output/review/highpoly-types/sxb010/bonsai-camera-semantic-side-elevation.png",
        "project_axes": [2, 1],
    },
}
payload["semantic_view_mapping_reason"] = (
    "The legacy proxy local axes are X=width, Y=height and Z=depth; semantic drawing views "
    "must not inherit the shared renderer's generic local-axis filenames."
)
manifest_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
