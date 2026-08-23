#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Baxter Beside BST02."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/bst02",
    "HIGHPOLY_ISOLATED_IFC": "Baxter-Beside-BST02-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Baxter-Beside-BST02-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "3hA0vKpcn44u4Tsx4tqiUz",
    "HIGHPOLY_CAMERA_PREFIX": "BAXTER_BESIDE_BST02",
    "HIGHPOLY_PRODUCT_LABEL": "Baxter Beside bedside table / BST02",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_bst02_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
