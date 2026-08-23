#!/usr/bin/env python3
"""Overlay the BED01 geometry-derived plan proxy on the project furniture plan."""

from falper_sorgente_linework import ROOT
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/bed01"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "3IQBEqO5vDI8Z9k1Ltge_N"
shared.OVERLAY_ID_PREFIX = "bed01"
shared.GENERATOR = "pipeline/scripts/render_bed01_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 1.8
shared.ALIGNMENT_MODE = "centre"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (14.0, 14.0), 1),
)


if __name__ == "__main__":
    shared.main()
