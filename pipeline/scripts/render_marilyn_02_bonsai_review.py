#!/usr/bin/env python3
"""Run the shared actual-Bonsai renderer for Baxter Marilyn 02 pouf."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/marilyn-02",
    "HIGHPOLY_ISOLATED_IFC": "Baxter-Marilyn-02-pouf-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Baxter-Marilyn-02-pouf-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "1THxa7p7n97w$wtLn4THjz",
    "HIGHPOLY_CAMERA_PREFIX": "BAXTER_MARILYN_02_POUF",
    "HIGHPOLY_PRODUCT_LABEL": "Baxter Marilyn pouf with swivel base 80 x 62 x 45 cm",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_marilyn_02_bonsai_review.py",
    "HIGHPOLY_FRONT_DIRECTION_SIGN": "-1",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
