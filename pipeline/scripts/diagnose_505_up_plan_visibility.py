#!/usr/bin/env python3
"""Build a diagnostic-only 505 UP Plan depth and component visibility audit.

This script does not change the candidate, approval, derived IFC, Bonsai
session, or any final Drawing. It derives visible horizontal regions from the
representative Body, from the Plan camera direction, and overlays the current
black candidate so every unresolved semantic boundary remains visible.
"""

import argparse
import json
import math
import sys
from hashlib import sha256 as hashlib_sha256
from pathlib import Path

import ifcopenshell
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

from int1_highpoly_type_review import mesh_for_one_product  # noqa: E402
from importlib import import_module  # noqa: E402

review = import_module("505_up_v1_lp_s_review")

FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
PRODUCT = ROOT / "output/review/highpoly-types/505-up-v1-lp-s"
GLOBAL_ID = "19MpdkWqXC7uhUNhLQgrce"
FRONT_SCENE_SHA256 = "85331ec1c29c42e053f06dce7fcf58a68f4cba8d594cb0044d28c0af1a6f10bf"
SIDE_SCENE_SHA256 = "e5c16902267f3ba56157230aca3ff7b234b446f7e47b5a7b7cecc248c3643c84"
FRONT_PNG_SHA256 = "8ff7aa9aab8ae33789f9b504ebbea63cbbf9910f75f08a8c581be08fd37361d9"
SIDE_PNG_SHA256 = "3879b3ae748a515d5d935192217afe30a23ab1c662c48e5eeb17e5857eff26f3"
PLAN_SCENE_SHA256 = "f67a0a3867e105e6ffc2a4507ba07f40d233c72142d46c5c06e2c9555ccf29d5"
DERIVED_IFC_SHA256 = "a49d1b9d2859fe596761ecc70c845cf888022f32c988a71946e8bfa4b0338417"
SESSION_IFC_SHA256 = "b8e9c5825225bd5f1df3a5ca5b5c9d0945e29acc9f1435c1e4dc98c3525b7687"


def file_sha256(path: Path) -> str:
    digest = hashlib_sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rounded_bbox(geometry):
    return [round(float(value), 3) for value in geometry.bounds]


def component_record(component):
    return {
        "root_vertex_index": component["root"],
        "minimum_local_mm": [round(float(value), 6) for value in component["minimum"]],
        "maximum_local_mm": [round(float(value), 6) for value in component["maximum"]],
        "size_mm": [round(float(value), 6) for value in component["size"]],
        "world_height_range_mm": [
            round(-float(component["maximum"][2]), 6),
            round(-float(component["minimum"][2]), 6),
        ],
    }


def horizontal_top_surfaces(vertices, faces, components, find):
    """Return the highest horizontal face union for every connected component."""
    by_root_height = {}
    for face in faces:
        points = [vertices[index] for index in face]
        local_z = [point[2] for point in points]
        if max(local_z) - min(local_z) > 0.05:
            continue
        polygon = Polygon([(point[0], point[1]) for point in points])
        if polygon.area <= 1e-6:
            continue
        world_height = round(-sum(local_z) / len(local_z), 3)
        root = find(face[0])
        by_root_height.setdefault((root, world_height), []).append(polygon)

    selected = []
    for component in components:
        heights = sorted(
            (height for root, height in by_root_height if root == component["root"]),
            reverse=True,
        )
        if not heights:
            continue
        top_height = heights[0]
        same_top = [
            polygon
            for (root, height), polygons in by_root_height.items()
            if root == component["root"] and top_height - height <= 0.1
            for polygon in polygons
        ]
        geometry = unary_union(same_top).buffer(0)
        if not geometry.is_empty:
            selected.append(
                {
                    "component": component,
                    "world_height_mm": top_height,
                    "geometry": geometry,
                }
            )

    visible = []
    covered = None
    heights = sorted({round(item["world_height_mm"], 1) for item in selected}, reverse=True)
    for height in heights:
        level = [item for item in selected if round(item["world_height_mm"], 1) == height]
        for item in level:
            geometry = item["geometry"]
            exposed = geometry if covered is None else geometry.difference(covered).buffer(0)
            if exposed.area > 1.0:
                visible.append({**item, "visible_geometry": exposed})
        level_union = unary_union([item["geometry"] for item in level]).buffer(0)
        covered = level_union if covered is None else unary_union([covered, level_union]).buffer(0)
    return visible


def semantic_key(item, main_root, back_root, slat_roots, display_root, lip_root, rail_roots):
    root = item["component"]["root"]
    if root in {main_root, back_root}:
        return "A_main_top"
    if root in slat_roots:
        return "B_right_slats"
    if root == display_root:
        return "C_front_left_display"
    if root == lip_root:
        return "D_mid_shelf_lip"
    if root in rail_roots:
        return "E_lower_front_rails"
    return "F_submillimetre_residual"


def polygons(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == "Polygon":
        return [geometry]
    if geometry.geom_type == "MultiPolygon":
        return list(geometry.geoms)
    return [part for part in getattr(geometry, "geoms", []) if part.geom_type == "Polygon"]


def svg_geometry_path(geometry, transform):
    commands = []
    for polygon in polygons(geometry):
        for ring in [polygon.exterior, *polygon.interiors]:
            points = [transform(float(x), float(y)) for x, y in ring.coords]
            if not points:
                continue
            commands.append("M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in points) + " Z")
    return " ".join(commands)


def svg_candidate_path(paths, transform):
    commands = []
    for path in paths:
        points = [transform(float(point[0]), float(point[1])) for point in path]
        if not points:
            continue
        command = "M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        if path[0] == path[-1]:
            command += " Z"
        commands.append(command)
    return " ".join(commands)


def candidate_segment_semantics(candidate_paths, region_geometries, probe_mm=2.0):
    """Name the visible region on both sides of every candidate segment."""
    semantic_keys = [
        "A_main_top",
        "B_right_slats",
        "C_front_left_display",
        "D_mid_shelf_lip",
        "E_lower_front_rails",
    ]
    probes = {
        key: region_geometries[key].buffer(0.35)
        for key in semantic_keys
    }

    def classify(x, y):
        point = Point(x, y)
        matches = [key for key in semantic_keys if probes[key].contains(point)]
        return matches[0] if matches else "background"

    records = []
    for path_index, path in enumerate(candidate_paths, 1):
        for segment_index, (start, end) in enumerate(zip(path, path[1:]), 1):
            dx, dy = end[0] - start[0], end[1] - start[1]
            length = math.hypot(dx, dy)
            if length <= 1e-9:
                continue
            midpoint = ((start[0] + end[0]) / 2.0, (start[1] + end[1]) / 2.0)
            normal = (-dy / length, dx / length)
            side_a = classify(midpoint[0] + normal[0] * probe_mm, midpoint[1] + normal[1] * probe_mm)
            side_b = classify(midpoint[0] - normal[0] * probe_mm, midpoint[1] - normal[1] * probe_mm)
            records.append(
                {
                    "segment_id": f"P{path_index:02d}-S{segment_index:02d}",
                    "path_id": f"P{path_index:02d}",
                    "start_mm": [round(float(value), 6) for value in start],
                    "end_mm": [round(float(value), 6) for value in end],
                    "length_mm": round(length, 6),
                    "side_a": side_a,
                    "side_b": side_b,
                    "semantic_separation_valid": side_a != side_b,
                }
            )
    return records


def render_svg(path, title, subtitle, region_geometries, candidate_paths, depth=False):
    width, height = 1600, 1050
    plot = (80.0, 175.0, 1120.0, 760.0)
    all_geometry = unary_union(list(region_geometries.values())).buffer(0)
    min_x, min_y, max_x, max_y = all_geometry.bounds
    # Project instance local Plan coordinates into the approved world-handed view:
    # screen horizontal = world +X = -local X; screen up = world +Y = local Y.
    screen_min_x, screen_max_x = -max_x, -min_x
    screen_min_y, screen_max_y = -max_y, -min_y
    span_x, span_y = screen_max_x - screen_min_x, screen_max_y - screen_min_y
    scale = min(plot[2] / span_x, plot[3] / span_y)
    offset_x = plot[0] + (plot[2] - span_x * scale) / 2.0
    offset_y = plot[1] + (plot[3] - span_y * scale) / 2.0

    def transform(local_x, local_y):
        return (
            offset_x + (-local_x - screen_min_x) * scale,
            offset_y + (-local_y - screen_min_y) * scale,
        )

    colours = {
        "A_main_top": "#4c78a8",
        "B_right_slats": "#17a398",
        "C_front_left_display": "#f28e2b",
        "D_mid_shelf_lip": "#edc948",
        "E_lower_front_rails": "#b279a2",
        "F_submillimetre_residual": "#b9b9b9",
    }
    heights = {
        "A_main_top": 2380.0,
        "B_right_slats": 1976.0,
        "C_front_left_display": 1610.705,
        "D_mid_shelf_lip": 470.0,
        "E_lower_front_rails": 440.0,
        "F_submillimetre_residual": 76.0,
    }

    fills = []
    for key, geometry in region_geometries.items():
        if geometry.is_empty:
            continue
        if depth:
            tone = round(38 + 200 * max(0.0, min(1.0, heights[key] / 2380.0)))
            colour = f"rgb({tone},{tone},{tone})"
        else:
            colour = colours[key]
        fills.append(
            f'<path d="{svg_geometry_path(geometry, transform)}" fill="{colour}" '
            'fill-rule="evenodd" stroke="#ffffff" stroke-width="1.2" vector-effect="non-scaling-stroke"/>'
        )

    semantic_boundaries = unary_union([geometry.boundary for geometry in region_geometries.values()])
    semantic_path = []
    for line in getattr(semantic_boundaries, "geoms", [semantic_boundaries]):
        if not hasattr(line, "coords"):
            continue
        points = [transform(float(x), float(y)) for x, y in line.coords]
        if len(points) >= 2:
            semantic_path.append("M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in points))

    labels = []
    for key, geometry in region_geometries.items():
        if geometry.is_empty or key == "F_submillimetre_residual":
            continue
        point = geometry.representative_point()
        x, y = transform(point.x, point.y)
        labels.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="15" fill="#ffffff" stroke="#111" stroke-width="1.5"/>'
            f'<text x="{x:.2f}" y="{y + 5:.2f}" text-anchor="middle" font-size="15" font-weight="700">{key[0]}</text>'
        )

    residual_boxes = []
    for index in range(1, 5):
        candidate = candidate_paths[index]
        xs = [point[0] for point in candidate]
        ys = [point[1] for point in candidate]
        x1, y1 = transform(max(xs), max(ys))
        x2, y2 = transform(min(xs), min(ys))
        residual_boxes.append(
            f'<rect x="{min(x1,x2)-6:.2f}" y="{min(y1,y2)-6:.2f}" width="{abs(x2-x1)+12:.2f}" height="{abs(y2-y1)+12:.2f}" '
            'fill="none" stroke="#d62728" stroke-width="2.5"/>'
        )

    legend_items = [
        ("A", "柜体顶面 + 背板顶面", "2380 mm", colours["A_main_top"]),
        ("B", "右侧 14 片格栅顶面", "1976 mm", colours["B_right_slats"]),
        ("C", "前方左侧 DISPLAY 顶面", "1610.705 mm", colours["C_front_left_display"]),
        ("D", "中部前沿构件可见窄带", "470 mm", colours["D_mid_shelf_lip"]),
        ("E", "更低前沿横向构件", "440 mm", colours["E_lower_front_rails"]),
    ]
    legend = []
    for idx, (letter, label, elevation, colour) in enumerate(legend_items):
        y = 235 + idx * 82
        fill = colour if not depth else "#f4f4f4"
        legend.append(
            f'<rect x="1265" y="{y-25}" width="34" height="34" rx="5" fill="{fill}" stroke="#333"/>'
            f'<text x="1282" y="{y-2}" text-anchor="middle" font-size="17" font-weight="700">{letter}</text>'
            f'<text x="1315" y="{y-8}" font-size="17" font-weight="600">{label}</text>'
            f'<text x="1315" y="{y+18}" font-size="15" fill="#555">世界顶面高程 {elevation}</text>'
        )

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="#f8f7f4"/>
<text x="80" y="65" font-family="Arial, sans-serif" font-size="31" font-weight="700">{title}</text>
<text x="80" y="105" font-family="Arial, sans-serif" font-size="18" fill="#444">{subtitle}</text>
<text x="80" y="135" font-family="Arial, sans-serif" font-size="15" fill="#666">画面方向：左→右 = 项目 world +X；上→下 = 项目 world +Y。黑线是当前候选，不是本轮定稿。</text>
<rect x="55" y="155" width="1170" height="810" rx="12" fill="#ffffff" stroke="#c8c8c8"/>
{''.join(fills)}
<path d="{' '.join(semantic_path)}" fill="none" stroke="#d62728" stroke-width="2.2" stroke-dasharray="9 6" vector-effect="non-scaling-stroke"/>
<path d="{svg_candidate_path(candidate_paths, transform)}" fill="none" stroke="#050505" stroke-width="2.4" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
{''.join(labels)}
{''.join(residual_boxes)}
<text x="1260" y="165" font-family="Arial, sans-serif" font-size="22" font-weight="700">可见区域语义</text>
{''.join(legend)}
<line x1="1265" y1="675" x2="1310" y2="675" stroke="#050505" stroke-width="3"/>
<text x="1320" y="681" font-family="Arial, sans-serif" font-size="16">当前候选黑线</text>
<line x1="1265" y1="720" x2="1310" y2="720" stroke="#d62728" stroke-width="3" stroke-dasharray="9 6"/>
<text x="1320" y="726" font-family="Arial, sans-serif" font-size="16">深度/组件区域边界</text>
<rect x="1265" y="755" width="45" height="28" fill="none" stroke="#d62728" stroke-width="2.5"/>
<text x="1320" y="776" font-family="Arial, sans-serif" font-size="16">P02–P05 无语义微小闭环</text>
<text x="1265" y="845" font-family="Arial, sans-serif" font-size="16" font-weight="700" fill="#9b1c1c">判定：当前 Plan 不可批准</text>
<text x="1265" y="875" font-family="Arial, sans-serif" font-size="15" fill="#444">红色虚线仍露出的地方，表示黑线没有</text>
<text x="1265" y="900" font-family="Arial, sans-serif" font-size="15" fill="#444">完成可见区域的语义分割。</text>
<text x="80" y="1010" font-family="Arial, sans-serif" font-size="14" fill="#666">来源：代表高模 Body / {GLOBAL_ID}；按连通组件、水平顶面与从上到下遮挡顺序机械提取。诊断专用，未写回 IFC。</text>
</svg>'''
    path.write_text(svg, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--output", type=Path, default=PRODUCT)
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    if file_sha256(source) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")

    candidate_path = output / "candidate-representations.json"
    candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
    candidate_paths = candidate["views"]["plan"]["proxy_paths_mm"]
    if len(candidate_paths) != 20:
        raise RuntimeError("expected the rejected 20-path Plan candidate")

    product, vertices, faces = mesh_for_one_product(ifcopenshell.open(source), GLOBAL_ID)
    components, find = review.mesh_components(vertices, faces)
    surfaces = horizontal_top_surfaces(vertices, faces, components, find)

    main = max(components, key=lambda item: item["size"][0] * item["size"][1] * item["size"][2])
    back = next(
        item for item in components
        if abs(item["size"][0] - 2912.0) <= 0.5 and abs(item["size"][1] - 16.0) <= 0.5
    )
    slats = {
        item["root"] for item in components
        if abs(item["size"][0] - 22.0) <= 0.5
        and abs(item["size"][1] - 25.0) <= 0.5
        and abs(item["size"][2] - 1528.0) <= 0.5
    }
    display = next(
        item for item in components
        if abs(item["size"][0] - 610.0) <= 0.75
        and abs(item["size"][1] - 420.462) <= 0.75
        and abs(item["size"][2] - 766.0) <= 0.75
    )
    lip = next(
        item for item in components
        if abs(item["size"][0] - 1248.0) <= 0.5
        and abs(item["size"][1] - 330.0) <= 0.5
        and abs(item["size"][2] - 42.0) <= 0.5
    )
    rails = {
        item["root"] for item in components
        if abs(item["size"][1] - 32.0) <= 0.5
        and abs(item["size"][2] - 376.0) <= 0.5
    }
    if len(slats) != 14 or len(rails) != 2:
        raise RuntimeError(f"component classification drifted: slats={len(slats)} rails={len(rails)}")

    region_parts = {key: [] for key in (
        "A_main_top", "B_right_slats", "C_front_left_display",
        "D_mid_shelf_lip", "E_lower_front_rails", "F_submillimetre_residual",
    )}
    surface_records = []
    for item in surfaces:
        key = semantic_key(
            item, main["root"], back["root"], slats, display["root"], lip["root"], rails
        )
        visible = item["visible_geometry"]
        if key == "F_submillimetre_residual" and visible.area < 500.0:
            region_parts[key].append(visible)
        elif key != "F_submillimetre_residual" and visible.area >= 10.0:
            region_parts[key].append(visible)
        else:
            continue
        surface_records.append(
            {
                "semantic_region": key,
                "component": component_record(item["component"]),
                "top_world_height_mm": round(item["world_height_mm"], 6),
                "visible_area_mm2": round(visible.area, 3),
                "visible_bbox_local_xy_mm": rounded_bbox(visible),
                "visible_part_count": len(polygons(visible)),
            }
        )

    region_geometries = {
        key: unary_union(parts).buffer(0) if parts else Polygon()
        for key, parts in region_parts.items()
    }
    component_svg = output / "505-up-plan-component-visibility-diagnostic.svg"
    depth_svg = output / "505-up-plan-depth-visibility-diagnostic.svg"
    render_svg(
        component_svg,
        "505 UP Plan：连通组件分色 + 当前黑线叠加",
        "按上方 Plan 相机方向提取实际可见顶面；红色虚线是组件/深度区域边界。",
        region_geometries,
        candidate_paths,
        depth=False,
    )
    render_svg(
        depth_svg,
        "505 UP Plan：真实顶视深度 + 当前黑线叠加",
        "亮色更高、暗色更低；相同投影位置仅保留从上向下首先命中的可见面。",
        region_geometries,
        candidate_paths,
        depth=True,
    )

    path_semantics = [
        {
            "path_id": "P01",
            "classification": "exterior_union_silhouette_only",
            "two_sides": ["visible geometry from several components", "background"],
            "verdict": "outer boundary is useful, but it merges A/C/D/E and cannot define their internal regions",
        },
        *[
            {
                "path_id": f"P{index:02d}",
                "classification": "unexplained_numeric_residual_loop",
                "two_sides": ["sub-millimetre tessellation residue near local y=0", "same top assembly"],
                "verdict": "invalid; no product or visible-face semantic boundary",
            }
            for index in range(2, 6)
        ],
        *[
            {
                "path_id": f"P{index:02d}",
                "classification": "right_slat_top_outline",
                "two_sides": ["B slat top at world 1976 mm", "background or lower front rail"],
                "verdict": "valid visible component boundary",
            }
            for index in range(6, 20)
        ],
        {
            "path_id": "P20",
            "classification": "main_top_to_display_depth_boundary",
            "two_sides": ["A main top at world 2380 mm", "C display top at world 1610.705 mm"],
            "verdict": "valid but incomplete; it does not segment A/D, C/E, or D/E",
        },
    ]
    segment_semantics = candidate_segment_semantics(candidate_paths, region_geometries)
    audit = {
        "schema_version": 1,
        "profile_key": "505-up-v1-lp-s",
        "representative_global_id": GLOBAL_ID,
        "status": "diagnostic_only_plan_revision_required",
        "user_rejection_evidence": "这里不行啊。我真是服了你就。哪怕按照3d模型从上到下它的深度可见面的来看，你都知道你这个图画里有问题。或者说，你反过来思考一下你的每个线分割的区域，它的语义分割对吗？",
        "plan_camera_visibility_rule": {
            "camera_direction_world": [0.0, 0.0, -1.0],
            "local_to_world_z_sign": -1,
            "finding": "local z near 0 is world floor level, not the product top",
            "occlusion": "at every Plan XY location only the highest world-Z horizontal surface remains visible",
        },
        "component_count": len(components),
        "semantic_regions": {
            "A_main_top": "continuous cabinet top and rear-panel top at world 2380 mm",
            "B_right_slats": "14 individually visible slat tops at world 1976 mm",
            "C_front_left_display": "DISPLAY protruding top visible beyond local y≈320.96 mm at world 1610.705 mm",
            "D_mid_shelf_lip": "narrow visible strip at world 470 mm",
            "E_lower_front_rails": "two lower front members, partly visible at world 440 mm",
            "F_submillimetre_residual": "tessellation/bevel residue without drawing semantics",
        },
        "visible_surface_records": surface_records,
        "candidate_path_semantics": path_semantics,
        "candidate_segment_semantics": segment_semantics,
        "candidate_segment_summary": {
            "segment_count": len(segment_semantics),
            "valid_two_side_semantic_count": sum(
                item["semantic_separation_valid"] for item in segment_semantics
            ),
            "invalid_same_side_or_unexplained_count": sum(
                not item["semantic_separation_valid"] for item in segment_semantics
            ),
            "invalid_segment_ids": [
                item["segment_id"]
                for item in segment_semantics
                if not item["semantic_separation_valid"]
            ],
        },
        "missing_semantic_separation": [
            "A main top versus D mid-shelf lip where the lower component overlaps the top footprint in projection",
            "C display versus E lower front rails where the display occludes the rails",
            "D mid-shelf lip versus E lower front rails at their depth transition",
            "P02-P05 must be removed because neither side names a distinct visible component or face region",
        ],
        "decision": {
            "current_plan_candidate_pass": False,
            "final_plan_updated": False,
            "derived_ifc_updated": False,
            "bonsai_session_updated": False,
            "create_drawing_called": False,
            "reason": "diagnostic regions are now explicit, but the simplified Plan boundary policy still requires owner review before a new candidate is written",
        },
        "frozen_views": {
            "plan_rejected_scene_svg_sha256": file_sha256(output / "project-drawings/MOLTENI-505-UP-ENTRANCE-PLAN.svg"),
            "front_scene_svg_sha256": file_sha256(output / "project-drawings/MOLTENI-505-UP-ENTRANCE-FRONT.svg"),
            "side_scene_svg_sha256": file_sha256(output / "project-drawings/MOLTENI-505-UP-ENTRANCE-SIDE.svg"),
            "front_rendered_png_sha256": file_sha256(output / "project-drawings/rendered/pdf-page-2.png"),
            "side_rendered_png_sha256": file_sha256(output / "project-drawings/rendered/pdf-page-3.png"),
            "derived_ifc_sha256": file_sha256(output / "505-up-v1-lp-s-derived-drawing.ifc"),
            "bonsai_session_ifc_sha256": file_sha256(output / "505-up-v1-lp-s-bonsai-drawing-session.ifc"),
        },
        "formal_ifc_sha256": file_sha256(source),
        "outputs": {
            "component_visibility_svg": str(component_svg.relative_to(ROOT)),
            "depth_visibility_svg": str(depth_svg.relative_to(ROOT)),
        },
        "course_evidence": {
            "lesson": "085000 Introduction to Drawings",
            "course_fact": "Plan camera direction, camera Depth, and Underlay/Linework/Annotation layers determine the generated SVG; the SVG must be visually inspected after Create Drawing.",
            "current_version_inference": "Depth-visible Body surfaces are diagnosed mechanically before changing the persisted LINEWORK annotation in Blender 4.5.3 / Bonsai 0.8.4.",
        },
        "pass": True,
    }
    expected_frozen = {
        "plan_rejected_scene_svg_sha256": PLAN_SCENE_SHA256,
        "front_scene_svg_sha256": FRONT_SCENE_SHA256,
        "side_scene_svg_sha256": SIDE_SCENE_SHA256,
        "front_rendered_png_sha256": FRONT_PNG_SHA256,
        "side_rendered_png_sha256": SIDE_PNG_SHA256,
        "derived_ifc_sha256": DERIVED_IFC_SHA256,
        "bonsai_session_ifc_sha256": SESSION_IFC_SHA256,
    }
    if audit["frozen_views"] != expected_frozen:
        raise RuntimeError(f"Front/Side freeze drifted: {audit['frozen_views']}")
    audit_path = output / "505-up-plan-visibility-semantics.json"
    audit["outputs"]["component_visibility_svg_sha256"] = file_sha256(component_svg)
    audit["outputs"]["depth_visibility_svg_sha256"] = file_sha256(depth_svg)
    for key, name in (
        ("component_visibility_png", "505-up-plan-component-visibility-diagnostic.png"),
        ("depth_visibility_png", "505-up-plan-depth-visibility-diagnostic.png"),
    ):
        raster = output / name
        if raster.exists():
            audit["outputs"][key] = str(raster.relative_to(ROOT))
            audit["outputs"][f"{key}_sha256"] = file_sha256(raster)
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    markdown = output / "505-up-plan-visibility-audit.md"
    markdown.write_text(
        "# 505 UP Plan 深度可见性诊断\n\n"
        "本轮只做诊断，不改最终 Plan、不写派生 IFC、不刷新 Bonsai Drawing。Front 与 Side 已冻结。\n\n"
        "## 结论\n\n"
        "上一轮把局部 `z≈0` 的 76 mm 板误认作顶盖；实例放置矩阵反转了局部 Z，所以它们实际位于世界高程 0–76 mm。真正顶面为 A（2380 mm）。当前黑线 P01 把 A/C/D/E 合成外轮廓，P20 只补了 A/C；D/E 及 C/E 的深度区域没有完整语义分割，P02–P05 更是没有产品语义的微小闭环，因此当前 Plan 不可批准。\n\n"
        "## 诊断图\n\n"
        "![组件分色](./505-up-plan-component-visibility-diagnostic.svg)\n\n"
        "![深度图](./505-up-plan-depth-visibility-diagnostic.svg)\n\n"
        "红色虚线为实际可见区域边界，黑线为当前候选。红线未被黑线覆盖的位置就是缺失的语义分割。\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "product": product.Name,
        "component_count": len(components),
        "visible_surface_record_count": len(surface_records),
        "component_svg": str(component_svg),
        "depth_svg": str(depth_svg),
        "audit": str(audit_path),
        "final_plan_updated": False,
        "formal_ifc_sha256": file_sha256(source),
        "pass": True,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
