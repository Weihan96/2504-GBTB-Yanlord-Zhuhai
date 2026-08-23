#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Baxter Casablanca / BED01."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/bed01",
    "HIGHPOLY_ISOLATED_IFC": "Baxter-Casablanca-BED01-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Baxter-Casablanca-BED01-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "3IQBEqO5vDI8Z9k1Ltge_N",
    "HIGHPOLY_CAMERA_PREFIX": "BAXTER_CASABLANCA_BED01",
    "HIGHPOLY_PRODUCT_LABEL": "Baxter Casablanca BED01",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_bed01_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
