#!/usr/bin/env python3
"""Render the configured TRAP01 Body with cameras saved in an actual Bonsai session."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/trap01",
    "HIGHPOLY_ISOLATED_IFC": "Geberit-151.116.11.1-TRAP01-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Geberit-151.116.11.1-TRAP01-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "2Ak2ma0lvBEA49UpplzUqi",
    "HIGHPOLY_CAMERA_PREFIX": "GEBERIT_TRAP01_151_116_11_1",
    "HIGHPOLY_PRODUCT_LABEL": "Geberit space-saving dip tube trap 151.116.11.1 d32, configured project geometry",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_trap01_bonsai_review.py",
    "HIGHPOLY_FRONT_DIRECTION_SIGN": "-1",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
