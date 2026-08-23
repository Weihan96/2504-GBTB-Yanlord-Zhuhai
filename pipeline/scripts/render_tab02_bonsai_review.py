#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for RODA Bernardo 367 / TAB02."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/tab02",
    "HIGHPOLY_ISOLATED_IFC": "RODA-Bernardo-367-TAB02-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "RODA-Bernardo-367-TAB02-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "2eX84IyLr8_e34nuOoRPLQ",
    "HIGHPOLY_CAMERA_PREFIX": "RODA_BERNARDO_367_TAB02",
    "HIGHPOLY_PRODUCT_LABEL": "RODA Bernardo 367 side table / TAB02",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_tab02_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
