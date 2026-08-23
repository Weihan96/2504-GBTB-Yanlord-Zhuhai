#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Geberit 154.154.00.1."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/geberit-154-154-00",
    "HIGHPOLY_ISOLATED_IFC": "Geberit-154-154-00-1-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Geberit-154-154-00-1-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "14EazrLgP8whYZWY_yCuKy",
    "HIGHPOLY_CAMERA_PREFIX": "GEBERIT_154_154_00_1",
    "HIGHPOLY_PRODUCT_LABEL": "Geberit CleanLine installation set 154.154.00.1",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_geberit_154_154_00_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
