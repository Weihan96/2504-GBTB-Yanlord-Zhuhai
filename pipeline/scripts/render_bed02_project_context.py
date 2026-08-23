#!/usr/bin/env python3
"""Overlay the BED02 geometry-derived plan proxy on the project furniture plan."""

from falper_sorgente_linework import ROOT
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/bed02"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "1i_pqgLv9A7uuV7MjaArBW"
shared.OVERLAY_ID_PREFIX = "bed02"
shared.GENERATOR = "pipeline/scripts/render_bed02_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.8
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (14.0, 14.0), 1),
)


if __name__ == "__main__":
    shared.main()
