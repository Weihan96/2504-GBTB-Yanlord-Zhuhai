#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Pivot + Senzafine / WD01."""

import os
import runpy
from pathlib import Path

import addon_utils


addon_utils.enable("bonsai_bridge", default_set=False, persistent=False)

ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/wd01",
    "HIGHPOLY_ISOLATED_IFC": "Poliform-Pivot-Senzafine-WD01-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Poliform-Pivot-Senzafine-WD01-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "2mPTt7$nvBSAn28I9NWt6T",
    "HIGHPOLY_CAMERA_PREFIX": "POLIFORM_PIVOT_SENZAFINE_WD01",
    "HIGHPOLY_PRODUCT_LABEL": "Poliform Pivot + Senzafine WD01",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_wd01_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
