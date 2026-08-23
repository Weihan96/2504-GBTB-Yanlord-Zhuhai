#!/usr/bin/env python3
"""Run the shared actual-Bonsai camera renderer for Baxter Ninfea BST01."""

import os,runpy
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
os.environ.update({"HIGHPOLY_PRODUCT_DIR":"output/review/highpoly-types/bst01","HIGHPOLY_ISOLATED_IFC":"Baxter-Ninfea-BST01-bonsai-isolated.ifc","HIGHPOLY_BLEND":"Baxter-Ninfea-BST01-bonsai-review.blend","HIGHPOLY_GLOBAL_ID":"3eic1dzkn5heTIn4PhF37v","HIGHPOLY_CAMERA_PREFIX":"BAXTER_NINFEA_BST01","HIGHPOLY_PRODUCT_LABEL":"Baxter Ninfea bedside table / BST01","HIGHPOLY_GENERATOR":"pipeline/scripts/render_bst01_bonsai_review.py"})
runpy.run_path(str(ROOT/"pipeline/scripts/render_geberit_duofix_sigma_224_212_bonsai_review.py"),run_name="__main__")
