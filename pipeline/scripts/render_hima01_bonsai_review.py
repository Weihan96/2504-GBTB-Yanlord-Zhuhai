#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Poliform Hima HIMA01."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/hima01",
    "HIGHPOLY_ISOLATED_IFC": "Poliform-Hima-HIMA01-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Poliform-Hima-HIMA01-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "2xmcLzu1rDTeMzRuNxPDyE",
    "HIGHPOLY_CAMERA_PREFIX": "HIMA01",
    "HIGHPOLY_PRODUCT_LABEL": "Poliform Hima HIMA01",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_hima01_bonsai_review.py",
})
runpy.run_path(
    str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"),
    run_name="__main__",
)
