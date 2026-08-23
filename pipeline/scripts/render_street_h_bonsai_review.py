#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for STREET-H."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/street-h",
    "HIGHPOLY_ISOLATED_IFC": "antoniolupi-Street-H-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "antoniolupi-Street-H-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "2ajpw0I9n1dBypfISg3ejX",
    "HIGHPOLY_CAMERA_PREFIX": "ANTONIOLUPI_STREET_H",
    "HIGHPOLY_PRODUCT_LABEL": "antoniolupi Street / project STREET-H",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_street_h_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
