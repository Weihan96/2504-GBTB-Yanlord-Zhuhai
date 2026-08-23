#!/usr/bin/env python3
"""Overlay exact official Gessi 54294 native-DWG paths on project drawings."""

from falper_sorgente_linework import ROOT
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "2iKOL78$H0N9Yd9$ky3pW4"
shared.OVERLAY_ID_PREFIX = "gessi316-54294"
shared.GENERATOR = "pipeline/scripts/render_gessi316_54294_project_context.py"
shared.SOURCE_KIND = "native_dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54294 原生 DWG 图纸表达"
shared.SOURCE_DWG_SHA256 = "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4"
shared.OFFICIAL_CAD_USED = True
shared.PATH_KEY = "official_native_dwg_paths_mm"
shared.CLOSE_PATHS = False
shared.ALIGNMENT_MODE = "centre"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.60
shared.REVIEW_HIDDEN_CSS = ".official-elevation-anchor{display:none}"
shared.PLAN_SOURCE = ROOT / "drawings/Sanitary Plan.svg"
shared.FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-06-20-R12-NY.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-06-22-R12-NX.svg"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-sanitary-plan.svg", (3.0, 3.0), 1),
    ("front", "side", shared.FRONT_SOURCE, "project-context-front-elevation.svg", (2.0, 2.0), 0),
    ("side", "front", shared.SIDE_SOURCE, "project-context-side-elevation.svg", (2.0, 2.0), 0),
)


if __name__ == "__main__":
    shared.main()
