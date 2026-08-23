#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for domestic-custom STREET."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/street",
    "HIGHPOLY_ISOLATED_IFC": "antoniolupi-Street-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "antoniolupi-Street-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "1FgLPMw$5B4wBH2ySMkXE1",
    "HIGHPOLY_CAMERA_PREFIX": "ANTONIOLUPI_STREET",
    "HIGHPOLY_PRODUCT_LABEL": "antoniolupi Street domestic-custom top / project STREET",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_street_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
