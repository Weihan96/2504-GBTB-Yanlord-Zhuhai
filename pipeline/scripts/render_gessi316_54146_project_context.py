#!/usr/bin/env python3
"""Overlay exact Gessi316 54146 G000 native-DWG views on project drawings."""

from falper_sorgente_linework import ROOT, load_json
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54146"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "04DLh1Jk9Dcu9ibcaE0id8"
shared.OVERLAY_ID_PREFIX = "gessi316-54146"
shared.GENERATOR = "pipeline/scripts/render_gessi316_54146_project_context.py"
shared.SOURCE_KIND = "native_dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54146 G000 原生 DWG 图纸表达"
shared.SOURCE_DWG_SHA256 = "c8ddb90f61565d5273a32574f777812bbb1d9ef33e2b98479f1e7c381e69950d"
shared.OFFICIAL_CAD_USED = True
shared.PATH_KEY = "official_native_dwg_paths_mm"
shared.CLOSE_PATHS = False
shared.BLUE = "#1677c8"
shared.REVIEW_HIDDEN_CSS = "#noninteger-highlights{display:none}"
shared.PLAN_SOURCE = ROOT / "drawings/FFL Plan.svg"
shared.FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-06-20-R12-NY.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-06-21-R12-PX.svg"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-ffl-plan.svg", (12.0, 12.0), 0),
    ("front", "front", shared.FRONT_SOURCE, "project-context-front-elevation.svg", (7.0, 6.0), 0),
    ("side", "side", shared.SIDE_SOURCE, "project-context-side-elevation.svg", (7.0, 6.0), 0),
)

original_group_points = shared.group_points


def group_points(source):
    if source != shared.PLAN_SOURCE:
        return original_group_points(source)
    candidate = load_json(shared.CANDIDATE)
    points = [
        tuple(point)
        for path in candidate["views"]["plan"][shared.PATH_KEY]
        for point in path
    ]
    minimum = (min(point[0] for point in points), min(point[1] for point in points))
    maximum = (max(point[0] for point in points), max(point[1] for point in points))
    size = (maximum[0] - minimum[0], maximum[1] - minimum[1])
    # FFL Plan uses x_svg = 200 + 0.02 * project_x and
    # y_svg = 160 - 0.02 * project_y. The product's local plan centre is zero.
    centre = (200.0 + 0.02 * -5300.66871643066, 160.0 - 0.02 * 950.0)
    half = (size[0] * shared.PROJECT_SCALE_SVG_UNITS_PER_MM / 2.0, size[1] * shared.PROJECT_SCALE_SVG_UNITS_PER_MM / 2.0)
    return [
        (centre[0] - half[0], centre[1] - half[1]),
        (centre[0] + half[0], centre[1] + half[1]),
    ]


shared.group_points = group_points


if __name__ == "__main__":
    shared.main()
