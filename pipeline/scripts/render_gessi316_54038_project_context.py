#!/usr/bin/env python3
"""Overlay exact Gessi316 54038 native-DWG views on complete project drawings."""

from falper_sorgente_linework import ROOT, load_json
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54038"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "245NU$zZL0d9tYTVwwBdk$"
shared.OVERLAY_ID_PREFIX = "gessi316-54038"
shared.GENERATOR = "pipeline/scripts/render_gessi316_54038_project_context.py"
shared.SOURCE_KIND = "native_dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54038 原生 DWG 图纸表达"
shared.SOURCE_DWG_SHA256 = "2c56b532bbbb78e5668050d1a7b52a7d545153a85a05efc377c7a651890897f1"
shared.OFFICIAL_CAD_USED = True
shared.PATH_KEY = "official_native_dwg_paths_mm"
shared.CLOSE_PATHS = False
shared.BLUE = "#1677c8"
shared.REVIEW_HIDDEN_CSS = ".official-elevation-anchor,#noninteger-highlights{display:none}"
shared.PLAN_SOURCE = ROOT / "drawings/FFL PLAN.svg"
shared.FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-08-31-R17-NY.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-08-32-R17-NX.svg"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-ffl-plan.svg", (8.0, 8.0), 1),
    ("front", "side", shared.FRONT_SOURCE, "project-context-front-elevation.svg", (6.0, 6.0), 0),
    ("side", "front", shared.SIDE_SOURCE, "project-context-side-elevation.svg", (6.0, 6.0), 0),
)

shared.ALIGNMENT_MODE = "centre"
original_group_points = shared.group_points


def candidate_size(view: str, rotate_quarter_turns: int = 0) -> tuple[float, float]:
    candidate = load_json(shared.CANDIDATE)
    points = [
        tuple(point)
        for path in candidate["views"][view][shared.PATH_KEY]
        for point in path
    ]
    minimum, maximum = shared.bbox(points)
    size = (maximum[0] - minimum[0], maximum[1] - minimum[1])
    return size[::-1] if rotate_quarter_turns % 2 else size


def centred_bbox(centre: tuple[float, float], size_mm: tuple[float, float]):
    half = tuple(value * shared.PROJECT_SCALE_SVG_UNITS_PER_MM / 2.0 for value in size_mm)
    return [
        (centre[0] - half[0], centre[1] - half[1]),
        (centre[0] + half[0], centre[1] + half[1]),
    ]


def group_points(source):
    if source == shared.PLAN_SOURCE:
        # Representative placement: local X -> world +Y and local Y -> world -X.
        # Place the quarter-turned Body proxy at that world position rather than
        # stretching it to the smaller authored PLAN_VIEW geometry.
        placement = (-1839.811680, 1840.0)
        local_plan_centre = (49.3475115, -0.000389)
        centre = (
            200.0 + 0.02 * (placement[0] + local_plan_centre[0]),
            160.0 - 0.02 * (placement[1] + local_plan_centre[1]),
        )
        return centred_bbox(centre, candidate_size("plan", 1))
    original = original_group_points(source)
    minimum, maximum = shared.bbox(original)
    centre_x = (minimum[0] + maximum[0]) / 2.0
    view = "side" if source == shared.FRONT_SOURCE else "front"
    width_mm, height_mm = candidate_size(view)
    half_width = width_mm * shared.PROJECT_SCALE_SVG_UNITS_PER_MM / 2.0
    height = height_mm * shared.PROJECT_SCALE_SVG_UNITS_PER_MM
    return [(centre_x - half_width, minimum[1]), (centre_x + half_width, minimum[1] + height)]


shared.group_points = group_points


if __name__ == "__main__":
    shared.main()
