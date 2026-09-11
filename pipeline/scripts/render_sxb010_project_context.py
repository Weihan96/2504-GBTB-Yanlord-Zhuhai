#!/usr/bin/env python3
"""Overlay sxb010 semantic views on mechanically located project drawings."""

import ast
import json
import xml.etree.ElementTree as ET
from itertools import product

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/sxb010"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "1O9JRXCI56VRUbpuLJy86Z"
shared.OVERLAY_ID_PREFIX = "sxb010"
shared.GENERATOR = "pipeline/scripts/render_sxb010_project_context.py"
shared.PLAN_SOURCE = ROOT / "drawings/FFL PLAN.svg"
FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-04-10-R07-PY.svg"
SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-04-11-R07-PX.svg"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.18
shared.ALIGNMENT_MODE = "centre"
shared.CLOSE_PATHS = False
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-ffl-plan.svg", (8.0, 8.0), 0),
    ("front", "front", FRONT_SOURCE, "project-context-r07-front-elevation.svg", (6.0, 6.0), 0),
    ("side", "side", SIDE_SOURCE, "project-context-r07-side-elevation.svg", (8.0, 6.0), 0),
)

# Mechanically measured from the formal IFC Body of 1O9JRXCI56VRUbpuLJy86Z
# with USE_WORLD_COORDS enabled. These coordinates are positioning evidence only;
# the formal IFC remains read-only and the candidate paths remain local-axis views.
PROJECT_WORLD_BBOX_M = (
    (-2.578700412869453, -4.9666923713684, 0.833773986816406),
    (-0.766700962185859, -4.90329906082153, 2.4298496093749997),
)


def body_group_points(source):
    """Choose the dense Body projection when Body and Curve3D share one GUID."""
    groups = []
    for element in ET.parse(source).getroot().iter():
        if not element.tag.endswith("g") or shared.guid(element) != shared.GLOBAL_ID:
            continue
        points = [
            point
            for child in element.iter()
            if child.tag.endswith("path")
            for point in shared.path_points(child)
        ]
        if points:
            groups.append(points)
    if len(groups) != 2:
        raise RuntimeError(f"expected Body and Curve3D projections for {shared.GLOBAL_ID}, found {len(groups)}")
    return max(groups, key=len)


def attribute_with_suffix(element, suffix):
    return next((value for key, value in element.attrib.items() if key.endswith(suffix)), None)


def projected_world_bbox_points(source):
    """Project the IFC world bbox through the native drawing plane and matrix."""
    drawing = next(
        (
            element
            for element in ET.parse(source).getroot().iter()
            if element.tag.endswith("g")
            and attribute_with_suffix(element, "plane")
            and attribute_with_suffix(element, "matrix3")
        ),
        None,
    )
    if drawing is None:
        raise RuntimeError(f"native drawing transform is missing from {source}")
    plane = ast.literal_eval(attribute_with_suffix(drawing, "plane"))
    matrix = ast.literal_eval(attribute_with_suffix(drawing, "matrix3"))
    origin = [plane[axis][3] for axis in range(3)]
    camera_x = [plane[axis][0] for axis in range(3)]
    camera_y = [plane[axis][1] for axis in range(3)]

    projected = []
    for world in product(*zip(*PROJECT_WORLD_BBOX_M)):
        relative = [world[axis] - origin[axis] for axis in range(3)]
        local_x = sum(relative[axis] * camera_x[axis] for axis in range(3))
        local_y = sum(relative[axis] * camera_y[axis] for axis in range(3))
        projected.append(
            (
                matrix[0][0] * local_x + matrix[0][2],
                matrix[1][2] - matrix[1][1] * local_y,
            )
        )
    return projected


def context_target_points(source):
    if source == shared.PLAN_SOURCE:
        return body_group_points(source)
    if source in {FRONT_SOURCE, SIDE_SOURCE}:
        return projected_world_bbox_points(source)
    raise RuntimeError(f"unregistered sxb010 context source: {source}")


shared.group_points = context_target_points


def suppress_plan_elevation_markers(path):
    """Remove only marker-anchor lines that obscure this blind in the review."""
    content = path.read_text(encoding="utf-8")
    lines = [line for line in content.splitlines() if 'class="official-elevation-anchor"' not in line]
    if len(lines) == len(content.splitlines()):
        raise RuntimeError(f"expected official elevation anchors in {path}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    plan_record = next(item for item in manifest["views"] if item["view"] == "plan")
    plan_output = ROOT / plan_record["output"]
    plan_crop = ROOT / plan_record["review_crop"]
    plan_preview = ROOT / plan_record["review_preview"]
    suppress_plan_elevation_markers(plan_output)
    suppress_plan_elevation_markers(plan_crop)
    shared.write_png_preview(plan_crop, plan_preview)
    plan_record["output_sha256"] = sha256(plan_output)
    plan_record["review_crop_sha256"] = sha256(plan_crop)
    plan_record["review_preview_sha256"] = sha256(plan_preview)
    manifest["semantic_view_mapping"] = {
        "plan": {"candidate_axes": [0, 2], "source": relative(shared.PLAN_SOURCE)},
        "front": {"candidate_axes": [0, 1], "source": relative(FRONT_SOURCE)},
        "side": {"candidate_axes": [2, 1], "source": relative(SIDE_SOURCE)},
    }
    manifest["project_world_bbox_m"] = [list(item) for item in PROJECT_WORLD_BBOX_M]
    manifest["context_target_derivation"] = {
        "plan": "densest complete Body projection for the matching IFC GlobalId in FFL PLAN.svg",
        "front": "formal IFC Body world bbox projected through the native R07 +Y drawing ifc:plane and ifc:matrix3",
        "side": "formal IFC Body world bbox projected through the native R07 +X drawing ifc:plane and ifc:matrix3",
        "no_nonuniform_fit_or_blank_space_guessing": True,
    }
    manifest["legacy_plan_source_group_count"] = 2
    manifest["selected_legacy_plan_source_group"] = "densest_Body_projection_not_Curve3D"
    manifest["rejected_context_evidence"] = {
        "source": "drawings/elevations/EL-P02-bonsai-public-long-section.svg",
        "reason": "the matching GUID contains only an approximately 2.113 x 0.662 SVG-unit clipped fragment, not a complete semantic elevation",
        "used_as_alignment_target": False,
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "reason": "the EL-04 marker cluster overlaps the blind location; walls, furniture, grids and IFC geometry remain retained",
        "geometry_removed": False,
    }
    write_json(manifest_path, manifest)
    print(json.dumps({"updated_manifest": relative(manifest_path), "pass": manifest["pass"]}, indent=2))
