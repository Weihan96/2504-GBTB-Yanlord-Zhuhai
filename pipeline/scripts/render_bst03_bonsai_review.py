#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Baxter Stone BST03."""

import os
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
os.environ.update({
    "HIGHPOLY_PRODUCT_DIR": "output/review/highpoly-types/bst03",
    "HIGHPOLY_ISOLATED_IFC": "Baxter-Stone-BST03-bonsai-isolated.ifc",
    "HIGHPOLY_BLEND": "Baxter-Stone-BST03-bonsai-review.blend",
    "HIGHPOLY_GLOBAL_ID": "1HZoxe$df4cBb4UXXH5J2S",
    "HIGHPOLY_CAMERA_PREFIX": "BAXTER_STONE_BST03",
    "HIGHPOLY_PRODUCT_LABEL": "Baxter Stone freestanding bedside table with L drawer / BST03",
    "HIGHPOLY_GENERATOR": "pipeline/scripts/render_bst03_bonsai_review.py",
})
runpy.run_path(str(ROOT / "pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"), run_name="__main__")
