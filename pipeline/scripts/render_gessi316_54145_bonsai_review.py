#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Gessi316 54145."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/gessi316-54145",
    "HIGHPOLY_ISOLATED_IFC": "Gessi316-54145-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Gessi316-54145-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "3jT4sCgpHC98VSIUdGUNYH",
    "HIGHPOLY_CAMERA_PREFIX": "GESSI316_54145",
    "HIGHPOLY_PRODUCT_LABEL": "Gessi316 54145",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_gessi316_54145_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
