#!/usr/bin/env python3
"""Overlay exact Baxter Marilyn pouf DWG views on complete project drawings."""

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from render_sis04_project_context import add_white_review_background, suppress_plan_elevation_markers, write_uncached_png_preview
import render_marilyn_01_project_context as shared


PRODUCT_DIR = ROOT / "output/review/highpoly-types/marilyn-02"
PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-09-35-R22-NY.svg"
SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-09-36-R22-NX.svg"

shared.PRODUCT_DIR = PRODUCT_DIR
shared.CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.PLAN_SOURCE = PLAN_SOURCE
shared.FRONT_SOURCE = FRONT_SOURCE
shared.SIDE_SOURCE = SIDE_SOURCE
shared.GLOBAL_ID = "1THxa7p7n97w$wtLn4THjz"
shared.SOURCE_KIND = "native_dwg"
shared.SOURCE_LABEL_ZH = "Baxter 官方原生 DWG"
shared.SOURCE_DWG_SHA256 = "cb9825eecab7c78ee28e7668e933c2fe413395a63bf75f3afab8cc4959864724"
shared.OVERLAY_ID_PREFIX = "marilyn-02"
shared.GENERATOR = "pipeline/scripts/render_marilyn_02_project_context.py"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", PLAN_SOURCE, "project-context-furniture-plan.svg", (12.0, 12.0), 0),
    ("front", "front", FRONT_SOURCE, "project-context-r22-front-elevation.svg", (9.0, 7.0), 0),
    ("side", "side", SIDE_SOURCE, "project-context-r22-side-elevation.svg", (10.0, 7.0), 0),
)


def postprocess_manifest():
    manifest_path = PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    for record in manifest["views"]:
        crop = ROOT / record["review_crop"]
        preview = ROOT / record["review_preview"]
        if record["view"] == "plan":
            suppress_plan_elevation_markers(crop)
        else:
            shared.suppress_noninteger_diagnostic_highlights(crop)
        add_white_review_background(crop)
        write_uncached_png_preview(crop, preview)
        record["review_crop_sha256"] = sha256(crop)
        record["review_preview_sha256"] = sha256(preview)
    manifest["generator"] = shared.GENERATOR
    manifest["semantic_view_mapping"] = {
        "plan": {"candidate_axes": [0, 1], "source": relative(PLAN_SOURCE), "rotate_quarter_turns": 0, "project_direction": "+Z"},
        "front": {"candidate_axes": [0, 2], "source": relative(FRONT_SOURCE), "rotate_quarter_turns": 0, "project_direction": "-Y", "complete_project_projection_group_only": True},
        "side": {"candidate_axes": [1, 2], "source": relative(SIDE_SOURCE), "rotate_quarter_turns": 0, "project_direction": "-X", "complete_project_projection_group_only": True},
    }
    manifest["context_view_scope"] = {
        "included": ["plan", "front", "side"],
        "excluded": [],
        "reason": "The Furniture Plan and selected R22 -Y/-X elevations contain mechanically complete Marilyn 02 projections.",
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "front": "noninteger residual-error diagnostic highlights only",
        "side": "noninteger residual-error diagnostic highlights only",
        "scope": "review crops only; complete project SVGs retain their annotations",
        "walls_furniture_and_ifc_geometry_removed": False,
    }
    manifest["review_preview_background"] = {
        "plan": "opaque white review-only paper background",
        "front": "opaque white review-only paper background",
        "side": "opaque white review-only paper background",
        "full_project_svgs_remain_transparent": True,
        "drawing_geometry_changed": False,
    }
    manifest["project_projection_note"] = {
        "official_nominal_xyz_mm": [800.0, 620.0, 450.0],
        "ifc_body_local_xyz_mm": [810.61496, 591.948944, 456.581987],
        "project_ifc_description": "Pouf with swivel base W80D62H45",
        "interpretation": "Exact official native-DWG pouf views are overlaid at native 1:50 scale and aligned by centre/floor baseline; no view is stretched to force-fit the project projection.",
        "geometry_stretched": False,
    }
    write_json(manifest_path, manifest)


if __name__ == "__main__":
    shared.main()
    postprocess_manifest()
