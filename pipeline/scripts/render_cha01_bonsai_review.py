#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Baxter Colette / CHA01."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/cha01",
    "HIGHPOLY_ISOLATED_IFC": "Baxter-Colette-CHA01-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Baxter-Colette-CHA01-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "1luHljRzDAhPNXxTDNu7qB",
    "HIGHPOLY_CAMERA_PREFIX": "BAXTER_COLETTE_CHA01",
    "HIGHPOLY_PRODUCT_LABEL": "Baxter Colette CHA01",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_cha01_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
