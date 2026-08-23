#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for project type FAU02."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/fau02",
    "HIGHPOLY_ISOLATED_IFC": "FAU02-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "FAU02-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "36ZX3QPyD7SvlXsDKMP8rY",
    "HIGHPOLY_CAMERA_PREFIX": "FAU02",
    "HIGHPOLY_PRODUCT_LABEL": "FAU02",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_fau02_bonsai_review.py",
})
runpy.run_path(
    str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"),
    run_name="__main__",
)
