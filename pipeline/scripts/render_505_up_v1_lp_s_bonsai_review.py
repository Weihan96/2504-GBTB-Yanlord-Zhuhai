#!/usr/bin/env python3
"""Run the shared actual-Bonsai renderer for one Molteni 505 UP instance."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/505-up-v1-lp-s",
    "HIGHPOLY_ISOLATED_IFC": "Molteni-505-UP-V1-LP-S-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Molteni-505-UP-V1-LP-S-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "19MpdkWqXC7uhUNhLQgrce",
    "HIGHPOLY_CAMERA_PREFIX": "MOLTENI_505_UP_V1_LP_S",
    "HIGHPOLY_PRODUCT_LABEL": "Molteni 505 UP / project 505 UP V1.LP.S",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_505_up_v1_lp_s_bonsai_review.py",
    "HIGHPOLY_FRONT_DIRECTION_SIGN": "1",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
