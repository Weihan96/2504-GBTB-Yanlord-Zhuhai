#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for RODA Orson 002 / CHA02."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/cha02",
    "HIGHPOLY_ISOLATED_IFC": "RODA-Orson-002-CHA02-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "RODA-Orson-002-CHA02-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "3YoxxZCgbF3Ap7gztAKkcs",
    "HIGHPOLY_CAMERA_PREFIX": "RODA_ORSON_002_CHA02",
    "HIGHPOLY_PRODUCT_LABEL": "RODA Orson 002 CHA02",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_cha02_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
