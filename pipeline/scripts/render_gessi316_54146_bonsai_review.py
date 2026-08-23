#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Gessi316 54146."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/gessi316-54146",
    "HIGHPOLY_ISOLATED_IFC": "Gessi316-54146-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Gessi316-54146-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "04DLh1Jk9Dcu9ibcaE0id8",
    "HIGHPOLY_CAMERA_PREFIX": "GESSI316_54146",
    "HIGHPOLY_PRODUCT_LABEL": "Gessi316 54146",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_gessi316_54146_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
