#!/usr/bin/env python3
"""Overlay the electric flue proxy on its two real project elevations."""

import re

from falper_sorgente_linework import ROOT
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/electric-flue-check-valve"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "1faflkXXH6M9cnYPE9Liir"
shared.OVERLAY_ID_PREFIX = "electric-flue-check-valve"
shared.GENERATOR = "pipeline/scripts/render_electric_flue_check_valve_project_context.py"
shared.PX_SOURCE = ROOT / "drawings/elevations/native/EL-03-07-R04-PX.svg"
shared.NY_SOURCE = ROOT / "drawings/elevations/native/EL-03-08-R04-NY.svg"
shared.CONTEXT_VIEWS = (
    ("r04-px", "plan", shared.PX_SOURCE, "project-context-r04-px-elevation.svg", (8.0, 6.0), 1),
    ("r04-ny", "front", shared.NY_SOURCE, "project-context-r04-ny-elevation.svg", (8.0, 6.0), 1),
)
shared.ALIGNMENT_MODE = "centre"


original_add_overlay = shared.add_overlay


def add_overlay(source, target, view, path):
    original_add_overlay(source, target, view, path)
    content = target.read_text(encoding="utf-8")
    content = content.replace(
        "</defs>",
        '<style type="text/css">.official-elevation-anchor,#noninteger-highlights{display:none}</style></defs>',
        1,
    )
    content = re.sub(
        r'<g id="noninteger-highlights"[^>]*>.*?</g>',
        "",
        content,
        flags=re.DOTALL,
    )
    target.write_text(content, encoding="utf-8")


shared.add_overlay = add_overlay


if __name__ == "__main__":
    shared.main()
