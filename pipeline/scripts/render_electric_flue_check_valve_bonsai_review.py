#!/usr/bin/env python3
"""Run the shared actual-Bonsai renderer for the electric flue proxy."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/electric-flue-check-valve",
    "HIGHPOLY_ISOLATED_IFC": "Electric-flue-check-valve-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Electric-flue-check-valve-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "1faflkXXH6M9cnYPE9Liir",
    "HIGHPOLY_CAMERA_PREFIX": "ELECTRIC_FLUE_CHECK_VALVE",
    "HIGHPOLY_PRODUCT_LABEL": "Electric flue check valve",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_electric_flue_check_valve_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
