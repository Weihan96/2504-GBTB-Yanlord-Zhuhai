#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Baxter Viktor / BED02."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/bed02",
    "HIGHPOLY_ISOLATED_IFC": "Baxter-Viktor-BED02-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Baxter-Viktor-BED02-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "1i_pqgLv9A7uuV7MjaArBW",
    "HIGHPOLY_CAMERA_PREFIX": "BAXTER_VIKTOR_BED02",
    "HIGHPOLY_PRODUCT_LABEL": "Baxter Viktor BED02",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_bed02_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
