#!/usr/bin/env python3
"""Run the shared actual-Bonsai renderer for Baxter Marilyn 01."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/marilyn-01",
    "HIGHPOLY_ISOLATED_IFC": "Baxter-Marilyn-01-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Baxter-Marilyn-01-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "3l2Ji4k2H9oOTDGWYUq7uV",
    "HIGHPOLY_CAMERA_PREFIX": "BAXTER_MARILYN_01",
    "HIGHPOLY_PRODUCT_LABEL": "Baxter Marilyn bergere armchair 86 x 100 x 94 cm",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_marilyn_01_bonsai_review.py",
    "HIGHPOLY_FRONT_DIRECTION_SIGN": "-1",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
