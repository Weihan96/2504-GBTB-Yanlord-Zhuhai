#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Gessi316 45089_54294."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/gessi316-54294",
    "HIGHPOLY_ISOLATED_IFC": "Gessi316-45089-54294-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Gessi316-45089-54294-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "2iKOL78$H0N9Yd9$ky3pW4",
    "HIGHPOLY_CAMERA_PREFIX": "GESSI316_45089_54294",
    "HIGHPOLY_PRODUCT_LABEL": "Gessi316 45089_54294",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_gessi316_54294_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
