#!/usr/bin/env python3
"""Render one unchanged actual-Body Side comparison without replacing review renders."""

import hashlib
import json
import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
MANIFEST = PRODUCT_DIR / "bonsai-review-manifest.json"
RENDERS = (
    PRODUCT_DIR / "bonsai-camera-plan.png",
    PRODUCT_DIR / "bonsai-camera-front-elevation.png",
    PRODUCT_DIR / "bonsai-camera-side-elevation.png",
    PRODUCT_DIR / "bonsai-camera-iso.png",
)
ACTUAL_SIDE = PRODUCT_DIR / "bonsai-camera-side-elevation-actual-body.png"
ACTUAL_MANIFEST = PRODUCT_DIR / "actual-body-comparison-manifest.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


saved = {path: path.read_bytes() for path in (*RENDERS, MANIFEST)}
try:
    os.environ.update({
        "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/gessi316-54294",
        "HIGHPOLY_ISOLATED_IFC": "Gessi316-45089-54294-bonsai-isolated.ifc",
        "HIGHPOLY_BLEND": "Gessi316-45089-54294-actual-body-comparison.blend",
        "HIGHPOLY_GLOBAL_ID": "2iKOL78$H0N9Yd9$ky3pW4",
        "HIGHPOLY_CAMERA_PREFIX": "GESSI316_45089_54294_ACTUAL",
        "HIGHPOLY_PRODUCT_LABEL": "Gessi316 45089_54294 unchanged actual Body",
        "HIGHPOLY_GENERATOR": "pipeline/scripts/render_gessi316_54294_actual_body_comparison.py",
    })
    runpy.run_path(
        str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"),
        run_name="__main__",
    )
    actual = json.loads(MANIFEST.read_text(encoding="utf-8"))
    actual_side = next(item for item in actual["renders"] if item["view"] == "side")
    ACTUAL_SIDE.write_bytes((PRODUCT_DIR / Path(actual_side["path"]).name).read_bytes())
    actual_side["path"] = ACTUAL_SIDE.relative_to(ROOT).as_posix()
    actual_side["sha256"] = sha256(ACTUAL_SIDE)
    ACTUAL_MANIFEST.write_text(json.dumps({
        "schema_version": 1,
        "mode": "unchanged_actual_project_ifc_body_side_comparison",
        "representative_global_id": actual["representative_global_id"],
        "formal_ifc_sha256": actual["formal_ifc_sha256"],
        "formal_ifc_bytes_unchanged": actual["formal_ifc_bytes_unchanged"],
        "isolated_ifc": actual["isolated_ifc"],
        "isolated_ifc_sha256": actual["isolated_ifc_sha256"],
        "bonsai_session": actual["bonsai_session"],
        "side_render": actual_side,
        "purpose": "Disclose the shorter adjustable project Body beside the review-only official maximum-reach proxy; this image must not be presented as the modified review proxy.",
        "pass": True,
    }, indent=2) + "\n", encoding="utf-8")
finally:
    for path, content in saved.items():
        path.write_bytes(content)

print(json.dumps({
    "side": ACTUAL_SIDE.relative_to(ROOT).as_posix(),
    "manifest": ACTUAL_MANIFEST.relative_to(ROOT).as_posix(),
    "pass": True,
}, indent=2))
