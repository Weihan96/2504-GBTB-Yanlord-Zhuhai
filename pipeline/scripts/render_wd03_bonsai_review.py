#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Senzafine / WD03."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/wd03",
    "HIGHPOLY_ISOLATED_IFC": "Poliform-Senzafine-WD03-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Poliform-Senzafine-WD03-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "3cmikd9MTB$egM5KQaNgUf",
    "HIGHPOLY_CAMERA_PREFIX": "POLIFORM_SENZAFINE_WD03",
    "HIGHPOLY_PRODUCT_LABEL": "Poliform Senzafine WD03",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_wd03_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
