#!/usr/bin/env python3
"""Render actual Geberit 154.154.00.1.F IFC Body with four Bonsai cameras."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/geberit-154-154-00-1-f",
    "HIGHPOLY_ISOLATED_IFC": "Geberit-154-154-00-1-F-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Geberit-154-154-00-1-F-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "2jjNIn9gHBYwNSWwlM5T_i",
    "HIGHPOLY_CAMERA_PREFIX": "GEBERIT_154_154_00_1_F",
    "HIGHPOLY_PRODUCT_LABEL": "Geberit CleanLine 154.154.00.1 project flange component .F",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_geberit_154_154_00_1_f_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
