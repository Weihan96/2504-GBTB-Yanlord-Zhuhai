#!/usr/bin/env python3
"""Overlay exact Gessi316 54093 G000 native DWG on its project plan."""

from falper_sorgente_linework import ROOT, load_json, write_json
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54093"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "1jM_suNMPAw8_nvZejQSnp"
shared.OVERLAY_ID_PREFIX = "gessi316-54093"
shared.GENERATOR = "pipeline/scripts/render_gessi316_54093_project_context.py"
shared.SOURCE_KIND = "native_dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54093 G000 原生 DWG 图纸表达"
shared.SOURCE_DWG_SHA256 = "a9308fc8498d34c8cf2f68fd28aaf90b11e62bfc59424e0a5a3ac7f0d28d5047"
shared.OFFICIAL_CAD_USED = True
shared.PATH_KEY = "official_native_dwg_paths_mm"
shared.CLOSE_PATHS = False
shared.BLUE = "#1677c8"
shared.ALIGNMENT_MODE = "centre"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.03
shared.PLAN_SOURCE = ROOT / "drawings/Sanitary Plan.svg"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-sanitary-plan.svg", (15.0, 15.0), 0),
)


if __name__ == "__main__":
    shared.main()
    path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(path)
    manifest["semantic_view_mapping"] = {
        "plan": {
            "candidate_axes": [0, 1],
            "source": "drawings/Sanitary Plan.svg",
            "rotate_quarter_turns": 0,
            "project_direction": "+Z"
        }
    }
    manifest["context_view_scope"] = {
        "included": ["plan"],
        "excluded": ["front", "side"],
        "reason": "The 54093 GlobalId occurs in the clean FFL and Sanitary plan SVGs, but in no native project elevation SVG. Front and side project contexts are therefore not invented; the marker-heavy P-202 candidate is excluded from visual evidence."
    }
    manifest["project_projection_note"] = {
        "ifc_body_local_xyz_mm": [48.974943, 190.281433, 273.203123],
        "official_api_width_depth_height_mm": [50.0, 190.0, 273.0],
        "official_native_dwg_plan_envelope_mm": [49.973479, 190.312335],
        "geometry_stretched": False,
        "native_project_elevation_found": False
    }
    write_json(path, manifest)
