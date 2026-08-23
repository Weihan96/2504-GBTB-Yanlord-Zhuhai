#!/usr/bin/env python3
"""Run the shared actual-Bonsai renderer for Baxter Miami Soft I03."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/miamisoft-i03",
    "HIGHPOLY_ISOLATED_IFC": "Baxter-Miami-Soft-I03-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Baxter-Miami-Soft-I03-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "1hdNtfQRPCGfHLfqzdp6w0",
    "HIGHPOLY_CAMERA_PREFIX": "BAXTER_MIAMI_SOFT_I03",
    "HIGHPOLY_PRODUCT_LABEL": "Baxter Miami Soft I03 roll cushion",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_miamisoft_i03_bonsai_review.py",
    "HIGHPOLY_FRONT_DIRECTION_SIGN": "-1",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
