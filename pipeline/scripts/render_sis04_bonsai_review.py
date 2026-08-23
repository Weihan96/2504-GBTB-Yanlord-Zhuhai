#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for SIS04."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/sis04",
    "HIGHPOLY_ISOLATED_IFC": "sis04-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "sis04-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "1MzM8Ms2vFo8KEm503j9w2",
    "HIGHPOLY_CAMERA_PREFIX": "SIS04",
    "HIGHPOLY_PRODUCT_LABEL": "Molteni Sistema 7 Wall Unit 4 Doors / SIS04",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_sis04_bonsai_review.py",
})
runpy.run_path(
    str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"),
    run_name="__main__",
)
