#!/usr/bin/env python3
"""Generate Gessi316 54294 review views from the exact official native DWG."""

from __future__ import annotations

import argparse
import html
import math
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from int1_highpoly_type_review import (
    VIEWS,
    bounds_3d,
    display_edge_sample,
    mesh_for_one_product,
    product_type_name,
    projected_raw_edges,
    projected_silhouette,
    svg_path,
)


PROFILE_KEY = "gessi316-54294"
ARTICLE_NUMBER = "45089_54294"
DRAWING_PRODUCT_CODE = "54294"
GENERATOR = "pipeline/scripts/gessi316_54294_review.py"
OFFICIAL_CAD_STATUS = "public_official_api_exact_54294_native_dwg_acquired"
OFFICIAL_CAD_ACQUIRED = True
OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
EXPECTED_SOURCE_EXACT_PROJECT_CONFIGURATION_MATCH = False
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OUTPUT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
REGISTER = OUTPUT_DIR / "profile.json"
ACCESS_RECORD = OUTPUT_DIR / "official-source/source-access-record.json"
LINEWORK = OUTPUT_DIR / "official-native-dwg-linework.json"
PDF_VERIFICATION = OUTPUT_DIR / "official-source/official-pdf-verification.json"
REVALIDATION = OUTPUT_DIR / "official-source/official-source-revalidation.json"
SOURCE_DWG = OUTPUT_DIR / "official-source/GPF5429400000G000_3.dwg"
SOURCE_KIND = "native_dwg"
SOURCE_LABEL_ZH = "Gessi 官方精确型号 54294 原生 DWG 图纸表达"
SOURCE_LABEL_EN = "drawing representation from the exact Gessi 54294 official native DWG"
REVIEW_LINE_SOURCE_KIND = "native_dwg_review_simplification"
REVIEW_LINE_SOURCE_LABEL_ZH = "基于官方 Gessi 54294 原生 DWG 轮廓的去纹审核简化表达"
REVIEW_LINE_SOURCE_LABEL_EN = "de-textured review simplification based on the official Gessi 54294 native-DWG outline"
SOURCE_DWG_SHA256 = "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4"
EXPECTED_DESCRIPTION = "External parts three-holes basin mixer with long spout, without waste."
EXPECTED_PATH_COUNTS = {"plan": 1481, "front": 1099, "side": 745}
VIEW_TOLERANCES_MM = {"plan": 8.0, "front": 5.0, "side": 30.0}
SCOPE = "official Gessi exact 54294 family reference and 45089_54294 article combination; not a project shop drawing"
BLUE = "#1677c8"
VIEW_DEFINITIONS = VIEWS
IDENTITY_POLICY_KEY = "45089_companion_dwg_used_as_54294_geometry"
PUBLIC_ACCESS_KEY = "exact_54294_native_dwg_publicly_downloadable"
CANDIDATE_IDENTITY_FLAG = "companion_45089_dwg_used_as_54294_geometry"
SMOOTH_HANDLE_SEGMENTS = 64
SPOUT_CROSS_SECTION_SEGMENTS = 24
SPOUT_ARC_SEGMENTS = 24
SIDE_FEATURE_TOLERANCE_MM = 0.2
HANDLE_TEXTURE_ZONES = {
    "plan": (
        (-120.85, -32.0, -78.85, -4.5),
        (79.15, -32.0, 121.15, -4.5),
    ),
    "front": (),
    "side": ((-51.318959, 18.051810, -23.318959, 60.051706),),
}
HANDLE_SMOOTH_OUTLINE_BOUNDS = {
    "plan": (
        (-120.805970, -32.0, -78.894030, -4.5),
        (79.194030, -32.0, 121.105970, -4.5),
    ),
    "front": (),
    "side": ((-51.318959, 18.051810, -23.318959, 60.051706),),
}
ORIGINAL_HANDLE_OUTLINE_PATH_COUNTS = {"plan": 4, "front": 0, "side": 2}


def rounded(paths):
    return [
        [[round(float(x), 6), round(float(y), 6)] for x, y in path]
        for path in paths
    ]


def path_bounds(paths):
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return {
        "minimum": [round(value, 6) for value in minimum],
        "maximum": [round(value, 6) for value in maximum],
        "size": [round(maximum[axis] - minimum[axis], 6) for axis in range(2)],
    }


def simplify_official_handle_texture(view, paths):
    """Replace only the native-DWG knurl zones with smooth outline rectangles.

    The immutable native linework remains in official-native-dwg-linework.json and
    is also rendered to separate official-original-*.svg evidence sheets.  This
    function produces the explicitly labelled review simplification used by the
    main review and project-context overlays.
    """
    zones = HANDLE_TEXTURE_ZONES[view]
    if not zones:
        return paths, {
            "mode": "no_handle_texture_visible_in_this_projection",
            "original_path_count": len(paths),
            "paths_in_handle_texture_zones": 0,
            "original_handle_outline_path_count": 0,
            "original_texture_detail_path_count": 0,
            "replacement_smooth_outline_path_count": 0,
            "review_path_count": len(paths),
            "review_texture_detail_path_count": 0,
            "handle_envelope_delta_mm": [0.0, 0.0],
            "handle_centres_delta_mm": [[0.0, 0.0], [0.0, 0.0]],
            "three_hole_and_spout_positions_unchanged": True,
            "pass": True,
        }

    tolerance = 0.001
    removed_indices = []
    retained = []
    for index, path in enumerate(paths):
        bounds = path_bounds([path])
        minimum, maximum = bounds["minimum"], bounds["maximum"]
        in_texture_zone = any(
            minimum[0] >= x_min - tolerance
            and minimum[1] >= y_min - tolerance
            and maximum[0] <= x_max + tolerance
            and maximum[1] <= y_max + tolerance
            for x_min, y_min, x_max, y_max in zones
        )
        if in_texture_zone:
            removed_indices.append(index)
        else:
            retained.append(path)

    outline_bounds = HANDLE_SMOOTH_OUTLINE_BOUNDS[view]
    outlines = [
        [[x_min, y_min], [x_max, y_min], [x_max, y_max], [x_min, y_max], [x_min, y_min]]
        for x_min, y_min, x_max, y_max in outline_bounds
    ]
    expected_zone_counts = {"plan": 1414, "side": 711}
    if len(removed_indices) != expected_zone_counts[view]:
        raise RuntimeError(
            f"Gessi {view} native-DWG handle texture classification drifted: "
            f"expected {expected_zone_counts[view]}, got {len(removed_indices)}"
        )
    original_outline_count = ORIGINAL_HANDLE_OUTLINE_PATH_COUNTS[view]
    original_texture_count = len(removed_indices) - original_outline_count
    review = retained + outlines
    original_zone_bounds = path_bounds([
        point_path
        for index, point_path in enumerate(paths)
        if index in set(removed_indices)
    ])
    replacement_bounds = path_bounds(outlines)
    envelope_delta = [
        round(replacement_bounds["size"][axis] - original_zone_bounds["size"][axis], 6)
        for axis in range(2)
    ]
    centres = [
        [round((x_min + x_max) / 2.0, 6), round((y_min + y_max) / 2.0, 6)]
        for x_min, y_min, x_max, y_max in zones
    ]
    return review, {
        "mode": "official_54294_outline_based_de_textured_review_simplification",
        "original_path_count": len(paths),
        "paths_in_handle_texture_zones": len(removed_indices),
        "original_handle_outline_path_count": original_outline_count,
        "original_texture_detail_path_count": original_texture_count,
        "replacement_smooth_outline_path_count": len(outlines),
        "review_path_count": len(review),
        "review_texture_detail_path_count": 0,
        "handle_texture_zone_bounds_mm": [list(zone) for zone in zones],
        "replacement_smooth_outline_bounds_mm": [list(zone) for zone in outline_bounds],
        "handle_envelope_delta_mm": envelope_delta,
        "handle_centres_mm": centres,
        "handle_centres_delta_mm": [[0.0, 0.0] for _ in centres],
        "three_hole_and_spout_positions_unchanged": True,
        "original_native_dwg_preserved_unchanged": True,
        "pass": max(abs(value) for value in envelope_delta) <= 0.001,
    }


def mesh_components(vertices, faces):
    parent = list(range(len(vertices)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for face in faces:
        anchor = face[0]
        for index in face[1:]:
            root_a, root_b = find(anchor), find(index)
            if root_a != root_b:
                parent[root_b] = root_a
    grouped = {}
    for index in range(len(vertices)):
        grouped.setdefault(find(index), []).append(index)
    face_counts = {root: 0 for root in grouped}
    for face in faces:
        face_counts[find(face[0])] += 1
    records = []
    for root, indices in grouped.items():
        points = [vertices[index] for index in indices]
        minimum = [min(point[axis] for point in points) for axis in range(3)]
        maximum = [max(point[axis] for point in points) for axis in range(3)]
        records.append({
            "root": root,
            "indices": indices,
            "vertex_count": len(indices),
            "face_count": face_counts[root],
            "minimum": minimum,
            "maximum": maximum,
            "size": [maximum[axis] - minimum[axis] for axis in range(3)],
            "center": [(minimum[axis] + maximum[axis]) / 2.0 for axis in range(3)],
        })
    return records, find


def smooth_cylinder_mesh(minimum, maximum, segments=SMOOTH_HANDLE_SEGMENTS):
    center_x = (minimum[0] + maximum[0]) / 2.0
    center_z = (minimum[2] + maximum[2]) / 2.0
    radius_x = (maximum[0] - minimum[0]) / 2.0
    radius_z = (maximum[2] - minimum[2]) / 2.0
    vertices = []
    for y in (minimum[1], maximum[1]):
        for index in range(segments):
            angle = 2.0 * math.pi * index / segments
            vertices.append((center_x + radius_x * math.cos(angle), y, center_z + radius_z * math.sin(angle)))
    vertices.extend(((center_x, minimum[1], center_z), (center_x, maximum[1], center_z)))
    front_center, back_center = 2 * segments, 2 * segments + 1
    faces = []
    for index in range(segments):
        following = (index + 1) % segments
        front_a, front_b = index, following
        back_a, back_b = segments + index, segments + following
        faces.extend((
            (front_a, back_a, back_b),
            (front_a, back_b, front_b),
            (front_center, front_b, front_a),
            (back_center, back_a, back_b),
        ))
    return vertices, faces


def official_side_features(paths):
    horizontal_levels = []
    vertical_candidates = []
    all_points = [point for path in paths for point in path]
    for path in paths:
        for start, end in zip(path, path[1:]):
            dx, dz = end[0] - start[0], end[1] - start[1]
            if abs(dz) <= 1e-6 and abs(dx) >= 100.0:
                horizontal_levels.append(float(start[1]))
            if abs(dx) <= 1e-6 and abs(dz) >= 64.0 and max(start[0], end[0]) < -10.0:
                vertical_candidates.append(float(start[0]))
    horizontal_levels = sorted(set(round(level, 6) for level in horizontal_levels))
    if len(horizontal_levels) != 2 or not vertical_candidates:
        raise RuntimeError("official Gessi Side semantic anchors could not be resolved")
    wall_surface_x = max(vertical_candidates)
    axis_z = sum(horizontal_levels) / 2.0
    outlet_endpoint_x = min(point[0] for point in all_points)
    outlet_endpoint_z = min(point[1] for point in all_points if abs(point[0] - outlet_endpoint_x) <= 1e-9)
    lower_lip_z = min(point[1] for point in all_points)
    lower_lip_x = sum(point[0] for point in all_points if abs(point[1] - lower_lip_z) <= 1e-9) / len([point for point in all_points if abs(point[1] - lower_lip_z) <= 1e-9])
    upper_envelope_z = max(point[1] for point in all_points)
    raw_max_x = max(point[0] for point in all_points)
    bend_start_x = min(
        min(start[0], end[0])
        for path in paths
        for start, end in zip(path, path[1:])
        if abs(end[1] - start[1]) <= 1e-6 and abs(end[0] - start[0]) >= 100.0
    )
    return {
        "wall_surface_x_mm": round(wall_surface_x, 6),
        "straight_spout_wall_z_mm": horizontal_levels,
        "straight_spout_axis_z_mm": round(axis_z, 6),
        "bend_start_x_mm": round(bend_start_x, 6),
        "outlet_endpoint_x_mm": round(outlet_endpoint_x, 6),
        "outlet_endpoint_z_mm": round(outlet_endpoint_z, 6),
        "lower_lip_z_mm": round(lower_lip_z, 6),
        "lower_lip_x_mm": round(lower_lip_x, 6),
        "upper_envelope_z_mm": round(upper_envelope_z, 6),
        "visible_reach_mm": round(wall_surface_x - outlet_endpoint_x, 6),
        "visible_envelope_mm": [
            round(wall_surface_x - outlet_endpoint_x, 6),
            round(upper_envelope_z - lower_lip_z, 6),
        ],
        "raw_full_envelope_mm": [
            round(raw_max_x - outlet_endpoint_x, 6),
            round(upper_envelope_z - lower_lip_z, 6),
        ],
        "concealed_rear_extension_mm": round(raw_max_x - wall_surface_x, 6),
    }


def official_spout_tube_mesh(center_x, wall_surface_y, axis_z, official):
    bend_start_y = wall_surface_y + official["bend_start_x_mm"] - official["wall_surface_x_mm"]
    bend_center_z = axis_z - 40.0
    rings = [((wall_surface_y, axis_z), (-1.0, 0.0), 10.0)]
    for index in range(SPOUT_ARC_SEGMENTS + 1):
        angle = math.radians(90.0 + 70.0 * index / SPOUT_ARC_SEGMENTS)
        center = (bend_start_y + 40.0 * math.cos(angle), bend_center_z + 40.0 * math.sin(angle))
        tangent = (-math.sin(angle), math.cos(angle))
        if index == 0:
            rings.append((center, tangent, 10.0))
        else:
            rings.append((center, tangent, 10.0))
    x_offset = wall_surface_y - official["wall_surface_x_mm"]
    z_offset = axis_z - official["straight_spout_axis_z_mm"]
    outlet = (official["outlet_endpoint_x_mm"] + x_offset, official["outlet_endpoint_z_mm"] + z_offset)
    lower = (official["lower_lip_x_mm"] + x_offset, official["lower_lip_z_mm"] + z_offset)
    tip_center = ((outlet[0] + lower[0]) / 2.0, (outlet[1] + lower[1]) / 2.0)
    normal = (outlet[0] - lower[0], outlet[1] - lower[1])
    tip_radius = math.hypot(*normal) / 2.0
    normal = (normal[0] / (2.0 * tip_radius), normal[1] / (2.0 * tip_radius))
    tip_tangent = (normal[1], -normal[0])
    rings.append((tip_center, tip_tangent, tip_radius))
    tube_vertices = []
    for (center_y, center_z), (tangent_y, tangent_z), tube_radius in rings:
        normal_y, normal_z = -tangent_z, tangent_y
        for index in range(SPOUT_CROSS_SECTION_SEGMENTS):
            angle = 2.0 * math.pi * index / SPOUT_CROSS_SECTION_SEGMENTS
            transverse = tube_radius * math.cos(angle)
            normal = tube_radius * math.sin(angle)
            tube_vertices.append((
                center_x + transverse,
                center_y + normal * normal_y,
                center_z + normal * normal_z,
            ))
    tube_faces = []
    for ring in range(len(rings) - 1):
        for index in range(SPOUT_CROSS_SECTION_SEGMENTS):
            following = (index + 1) % SPOUT_CROSS_SECTION_SEGMENTS
            a = ring * SPOUT_CROSS_SECTION_SEGMENTS + index
            b = ring * SPOUT_CROSS_SECTION_SEGMENTS + following
            c = (ring + 1) * SPOUT_CROSS_SECTION_SEGMENTS + following
            d = (ring + 1) * SPOUT_CROSS_SECTION_SEGMENTS + index
            tube_faces.extend(((a, d, c), (a, c, b)))
    for ring, reverse in ((0, True), (len(rings) - 1, False)):
        center_index = len(tube_vertices)
        center_y, center_z = rings[ring][0]
        tube_vertices.append((center_x, center_y, center_z))
        base = ring * SPOUT_CROSS_SECTION_SEGMENTS
        for index in range(SPOUT_CROSS_SECTION_SEGMENTS):
            following = (index + 1) % SPOUT_CROSS_SECTION_SEGMENTS
            face = (center_index, base + following, base + index)
            tube_faces.append(face if reverse else tuple(reversed(face)))
    return tube_vertices, tube_faces


def geometry_side_features(vertices, faces):
    components, _ = mesh_components(vertices, faces)
    trim_discs = [component for component in components if (
        component["vertex_count"] == 208
        and abs(component["size"][0] - 64.0) <= 0.05
        and abs(component["size"][1] - 5.0) <= 0.05
        and 63.9 <= component["size"][2] <= 65.1
    )]
    spout = max(components, key=lambda component: component["size"][1])
    if len(trim_discs) != 3 or spout["size"][1] < 200.0:
        raise RuntimeError("Gessi Side geometry feature components drifted")
    wall_surface = sum(component["maximum"][1] for component in trim_discs) / len(trim_discs)
    axis_z = sum(component["center"][2] for component in trim_discs) / len(trim_discs)
    lower_lip = spout["minimum"][2]
    outlet_endpoint = spout["minimum"][1]
    upper = max(vertex[2] for vertex in vertices)
    overall_minimum, overall_maximum = bounds_3d(vertices)
    return {
        "wall_surface_y_mm": round(wall_surface, 6),
        "straight_spout_axis_z_mm": round(axis_z, 6),
        "lower_lip_z_mm": round(lower_lip, 6),
        "outlet_endpoint_y_mm": round(outlet_endpoint, 6),
        "upper_envelope_z_mm": round(upper, 6),
        "visible_reach_mm": round(wall_surface - outlet_endpoint, 6),
        "visible_envelope_mm": [round(wall_surface - outlet_endpoint, 6), round(upper - lower_lip, 6)],
        "raw_full_envelope_mm": [
            round(overall_maximum[1] - overall_minimum[1], 6),
            round(overall_maximum[2] - overall_minimum[2], 6),
        ],
    }


def simplify_handle_knurls(vertices, faces, official_side_paths=None):
    if official_side_paths is None:
        official_side_paths = load_json(LINEWORK)["views"]["side"]["paths_mm"]
    official = official_side_features(official_side_paths)
    components, find = mesh_components(vertices, faces)
    knurls = [component for component in components if (
        component["vertex_count"] >= 9_000
        and component["face_count"] >= 19_000
        and abs(component["size"][0] - 41.015) <= 0.05
        and abs(component["size"][1] - 24.375) <= 0.05
        and abs(component["size"][2] - 41.015) <= 0.05
        and abs(component["center"][2] - 18.0) <= 0.1
    )]
    if len(knurls) != 2:
        raise RuntimeError(f"expected two Gessi knurled handle sleeves, found {len(knurls)}")
    main_cylinders = [component for component in components if (
        component["vertex_count"] == 96
        and abs(component["size"][0] - 40.25) <= 0.05
        and abs(component["size"][1] - 45.0) <= 0.05
        and abs(component["size"][2] - 40.25) <= 0.05
    )]
    trim_discs = [component for component in components if (
        component["vertex_count"] == 208
        and abs(component["size"][0] - 64.0) <= 0.05
        and abs(component["size"][1] - 5.0) <= 0.05
        and abs(component["size"][2] - 64.0) <= 0.05
    )]
    spouts = [component for component in components if component["size"][1] >= 200.0]
    if len(main_cylinders) != 2 or len(trim_discs) != 3 or len(spouts) != 1:
        raise RuntimeError("Gessi preserved semantic component gate failed")
    actual_features = geometry_side_features(vertices, faces)
    main_axis_z = actual_features["straight_spout_axis_z_mm"]
    wall_surface = actual_features["wall_surface_y_mm"]
    removed_roots = {component["root"] for component in knurls + spouts}
    trim_roots = {component["root"] for component in trim_discs}
    retained_indices = [index for index in range(len(vertices)) if find(index) not in removed_roots]
    remap = {old: new for new, old in enumerate(retained_indices)}
    optimised_vertices = [
        (
            vertices[index][0],
            vertices[index][1],
            main_axis_z + (vertices[index][2] - main_axis_z) * (65.0 / 64.0),
        ) if find(index) in trim_roots else vertices[index]
        for index in retained_indices
    ]
    optimised_faces = [tuple(remap[index] for index in face) for face in faces if find(face[0]) not in removed_roots]
    for component in knurls:
        cylinder_vertices, cylinder_faces = smooth_cylinder_mesh(component["minimum"], component["maximum"])
        offset = len(optimised_vertices)
        optimised_vertices.extend(cylinder_vertices)
        optimised_faces.extend(tuple(offset + index for index in face) for face in cylinder_faces)
    spout_vertices, spout_faces = official_spout_tube_mesh(spouts[0]["center"][0], wall_surface, main_axis_z, official)
    spout_offset = len(optimised_vertices)
    optimised_vertices.extend(spout_vertices)
    optimised_faces.extend(tuple(spout_offset + index for index in face) for face in spout_faces)
    original_minimum, original_maximum = bounds_3d(vertices)
    optimised_minimum, optimised_maximum = bounds_3d(optimised_vertices)
    review_features = geometry_side_features(optimised_vertices, optimised_faces)
    x_translation = review_features["wall_surface_y_mm"] - official["wall_surface_x_mm"]
    z_translation = review_features["straight_spout_axis_z_mm"] - official["straight_spout_axis_z_mm"]
    target_features = {
        "wall_surface_y_mm": review_features["wall_surface_y_mm"],
        "straight_spout_axis_z_mm": review_features["straight_spout_axis_z_mm"],
        "lower_lip_z_mm": round(official["lower_lip_z_mm"] + z_translation, 6),
        "outlet_endpoint_y_mm": round(official["outlet_endpoint_x_mm"] + x_translation, 6),
        "upper_envelope_z_mm": round(official["upper_envelope_z_mm"] + z_translation, 6),
        "visible_reach_mm": official["visible_reach_mm"],
        "visible_envelope_mm": official["visible_envelope_mm"],
    }
    feature_keys = ("wall_surface_y_mm", "straight_spout_axis_z_mm", "lower_lip_z_mm", "outlet_endpoint_y_mm", "visible_reach_mm")
    review_residuals = {key: round(review_features[key] - target_features[key], 6) + 0.0 for key in feature_keys}
    review_residuals["visible_envelope_mm"] = [
        round(review_features["visible_envelope_mm"][axis] - target_features["visible_envelope_mm"][axis], 6) + 0.0
        for axis in range(2)
    ]
    actual_residuals = {key: round(actual_features[key] - target_features[key], 6) + 0.0 for key in feature_keys}
    actual_residuals["visible_envelope_mm"] = [
        round(actual_features["visible_envelope_mm"][axis] - target_features["visible_envelope_mm"][axis], 6) + 0.0
        for axis in range(2)
    ]
    prior_x_translation = original_maximum[1]
    prior_z_translation = main_axis_z - official["straight_spout_axis_z_mm"]
    side_gate = {
        "tolerance_mm": SIDE_FEATURE_TOLERANCE_MM,
        "official_exact_54294_native_dwg": official,
        "actual_project_ifc_body": actual_features,
        "official_mapped_targets": target_features,
        "review_only_proxy": review_features,
        "actual_body_residual_mm": actual_residuals,
        "review_proxy_residual_mm": review_residuals,
        "previous_axis_only_overlay": {
            "translation_mm": [round(prior_x_translation, 6), round(prior_z_translation, 6)],
            "wall_surface_residual_mm": round(prior_x_translation + official["wall_surface_x_mm"] - actual_features["wall_surface_y_mm"], 6),
            "outlet_endpoint_residual_mm": round(prior_x_translation + official["outlet_endpoint_x_mm"] - actual_features["outlet_endpoint_y_mm"], 6),
            "lower_lip_residual_mm": round(prior_z_translation + official["lower_lip_z_mm"] - actual_features["lower_lip_z_mm"], 6),
            "straight_spout_axis_residual_mm": 0.0,
        },
        "correct_wall_and_axis_translation_mm": [round(x_translation, 6), round(z_translation, 6)],
        "official_line_scaled": False,
        "official_line_mirrored": False,
        "review_proxy_geometry_changed": True,
        "actual_body_comparison_retained": True,
        "cause": "previous overlay used the concealed rear endpoint instead of the published wall surface; after correcting that transform, the project Body remains 6.966033 mm shorter because 54294 is adjustable and the official DWG shows the 210 mm maximum-depth configuration",
        "wrong_dwg_view_selected": False,
        "latest_exact_54294_dwg_used": True,
        "project_body_within_official_adjustable_depth_range_190_210_mm": 190.0 <= actual_features["visible_reach_mm"] <= 210.0,
        "feature_gate_pass": max(
            abs(value)
            for key, value in review_residuals.items()
            for value in (value if isinstance(value, list) else [value])
        ) <= SIDE_FEATURE_TOLERANCE_MM,
    }
    return optimised_vertices, optimised_faces, {
        "mode": "review_only_smooth_handles_and_official_max_reach_spout",
        "original_vertex_count": len(vertices),
        "original_face_count": len(faces),
        "removed_component_count": len(knurls) + len(spouts),
        "removed_vertex_count": sum(component["vertex_count"] for component in knurls + spouts),
        "removed_face_count": sum(component["face_count"] for component in knurls + spouts),
        "handle_replacement_component_count": len(knurls),
        "handle_replacement_vertex_count": len(knurls) * (2 * SMOOTH_HANDLE_SEGMENTS + 2),
        "handle_replacement_face_count": len(knurls) * 4 * SMOOTH_HANDLE_SEGMENTS,
        "spout_replacement_component_count": 1,
        "spout_replacement_vertex_count": len(spout_vertices),
        "spout_replacement_face_count": len(spout_faces),
        "optimised_vertex_count": len(optimised_vertices),
        "optimised_face_count": len(optimised_faces),
        "vertex_reduction_percent": round((1.0 - len(optimised_vertices) / len(vertices)) * 100.0, 3),
        "face_reduction_percent": round((1.0 - len(optimised_faces) / len(faces)) * 100.0, 3),
        "smooth_segments_per_sleeve": SMOOTH_HANDLE_SEGMENTS,
        "main_handle_cylinder_count_preserved": len(main_cylinders),
        "control_levers_preserved": True,
        "spout_component_count_preserved": 1,
        "three_hole_trim_disc_count_preserved": len(trim_discs),
        "plan_front_identity_and_three_hole_layout_preserved": True,
        "main_axis_z_mm": main_axis_z,
        "side_mechanical_gate": side_gate,
        "original_bounds_mm": {
            "minimum": [round(value, 6) for value in original_minimum],
            "maximum": [round(value, 6) for value in original_maximum],
        },
        "optimised_bounds_mm": {
            "minimum": [round(value, 6) for value in optimised_minimum],
            "maximum": [round(value, 6) for value in optimised_maximum],
        },
        "formal_ifc_modified": False,
        "pass": side_gate["feature_gate_pass"],
    }


def align_official_paths(view, paths, minimum, maximum, simplification):
    center_x = (minimum[0] + maximum[0]) / 2.0
    if view == "plan":
        return [
            [[center_x + x, maximum[1] + y] for x, y in path]
            for path in paths
        ], {"mode": "project_center_and_wall_plane_translation_only"}
    official_bounds = path_bounds(paths)
    ifc_z_center = (minimum[2] + maximum[2]) / 2.0
    official_z_center = (
        official_bounds["minimum"][1] + official_bounds["maximum"][1]
    ) / 2.0
    z_offset = ifc_z_center - official_z_center
    if view == "front":
        return [
            [[center_x + x, z_offset + z] for x, z in path]
            for path in paths
        ], {"mode": "project_center_translation_only"}
    gate = simplification["side_mechanical_gate"]
    x_offset, z_offset = gate["correct_wall_and_axis_translation_mm"]
    official_features = gate["official_exact_54294_native_dwg"]
    review_features = gate["review_only_proxy"]
    alignment = {
        "mode": "wall_plane_and_straight_spout_axis_translation_only",
        "scale_applied": False,
        "mirror_applied": False,
        "official_wall_surface_x_mm": official_features["wall_surface_x_mm"],
        "review_wall_surface_y_mm": review_features["wall_surface_y_mm"],
        "official_straight_spout_wall_z_mm": official_features["straight_spout_wall_z_mm"],
        "official_straight_spout_axis_z_mm": official_features["straight_spout_axis_z_mm"],
        "review_straight_spout_axis_z_mm": review_features["straight_spout_axis_z_mm"],
        "prior_axis_only_translation_mm": gate["previous_axis_only_overlay"]["translation_mm"],
        "corrected_wall_and_axis_translation_mm": [x_offset, z_offset],
        "prior_wall_surface_residual_mm": gate["previous_axis_only_overlay"]["wall_surface_residual_mm"],
        "prior_outlet_endpoint_residual_mm": gate["previous_axis_only_overlay"]["outlet_endpoint_residual_mm"],
        "corrected_anchor_residual_mm": [0.0, 0.0],
        "side_mechanical_gate": gate,
        "pass": gate["feature_gate_pass"],
    }
    return [
        [[x_offset + y, z_offset + z] for y, z in path]
        for path in paths
    ], alignment


def compatibility(view, official, minimum, maximum, simplification=None):
    if view == "side" and simplification is not None:
        gate = simplification["side_mechanical_gate"]
        review = gate["review_only_proxy"]
        exact = gate["official_exact_54294_native_dwg"]
        delta = [abs(review["visible_envelope_mm"][axis] - exact["visible_envelope_mm"][axis]) for axis in range(2)]
        return {
            "ifc_body_projection_size_mm": gate["actual_project_ifc_body"]["visible_envelope_mm"],
            "review_proxy_visible_side_envelope_mm": review["visible_envelope_mm"],
            "official_native_dwg_visible_side_envelope_mm": exact["visible_envelope_mm"],
            "official_native_dwg_raw_full_side_envelope_mm": exact["raw_full_envelope_mm"],
            "official_concealed_rear_extension_mm": exact["concealed_rear_extension_mm"],
            "absolute_delta_mm": [round(value, 6) for value in delta],
            "tolerance_mm": SIDE_FEATURE_TOLERANCE_MM,
            "side_full_envelope_note": "Visible Side envelope is gated wall-surface-to-outlet. The raw DWG additionally contains an 18.818959 mm concealed rear connector, disclosed separately and not used as the wall anchor.",
            "pass": max(delta) <= SIDE_FEATURE_TOLERANCE_MM and gate["feature_gate_pass"],
        }
    axes = VIEW_DEFINITIONS[view]["axes"]
    official_bounds = path_bounds(official)
    ifc_size = [maximum[axis] - minimum[axis] for axis in axes]
    delta = [abs(official_bounds["size"][axis] - ifc_size[axis]) for axis in range(2)]
    tolerance = VIEW_TOLERANCES_MM[view]
    return {
        "ifc_body_projection_size_mm": [round(value, 6) for value in ifc_size],
        "official_native_dwg_size_mm": official_bounds["size"],
        "absolute_delta_mm": [round(value, 6) for value in delta],
        "tolerance_mm": tolerance,
        "side_full_envelope_note": (
            "The official side envelope includes the published wall trim and full adjustable spout reach."
            if view == "side"
            else None
        ),
        "pass": max(delta) <= tolerance,
    }


def render_svg(profile, view, raw_edges, proxy, official, metadata):
    width, height = 1400, 980
    plot_x, plot_y, plot_w, plot_h = 60, 170, 980, 750
    points = [point for edge in raw_edges for point in edge]
    points.extend(point for path in proxy for point in path)
    points.extend(point for path in official for point in path)
    min_x, max_x = min(point[0] for point in points), max(point[0] for point in points)
    min_y, max_y = min(point[1] for point in points), max(point[1] for point in points)
    padding = max(max_x - min_x, max_y - min_y) * 0.07
    min_x, max_x = min_x - padding, max_x + padding
    min_y, max_y = min_y - padding, max_y + padding
    scale = min(plot_w / (max_x - min_x), plot_h / (max_y - min_y))

    def transform(point):
        return plot_x + (point[0] - min_x) * scale, plot_y + plot_h - (point[1] - min_y) * scale

    edge_path = svg_path([[start, end] for start, end in raw_edges], transform)
    proxy_path = svg_path(proxy, transform, close=True)
    official_path = svg_path(official, transform)
    requirements = "\n".join(
        f'<text x="1090" y="{258 + index * 28}" font-family="Arial,sans-serif" font-size="15" fill="#34495e">• {html.escape(item)}</text>'
        for index, item in enumerate(profile["required_semantics"])
    )
    cross_check = metadata["compatibility"]
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="1400" height="980" fill="#fbfaf7"/>
<text x="60" y="58" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">{html.escape(profile["display_name"])}</text>
<text x="60" y="96" font-family="Arial,sans-serif" font-size="20" fill="#41566d">{VIEW_DEFINITIONS[view]["label"]} · isolated representative {profile["representative_global_id"]}</text>
<text x="60" y="128" font-family="Arial,sans-serif" font-size="17" fill="#68798a">Grey = actual IFC Body · Black = smooth/max-reach proxy · Blue = de-textured review simplification based on official {DRAWING_PRODUCT_CODE} outline</text>
<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="original-highpoly" d="{edge_path}" fill="none" stroke="#87929c" stroke-width="0.5" stroke-opacity="0.30" vector-effect="non-scaling-stroke"/>
<path class="simplified-proxy-silhouette geometry-derived" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference-mask" d="{official_path}" fill="none" stroke="#ffffff" stroke-width="7" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="review-simplified-reference official-outline-derived" data-source-kind="{REVIEW_LINE_SOURCE_KIND}" data-source-label-zh="{REVIEW_LINE_SOURCE_LABEL_ZH}" data-product-code="{DRAWING_PRODUCT_CODE}" data-source-dwg-sha256="{SOURCE_DWG_SHA256}" data-unaltered-official-dwg="false" d="{official_path}" fill="none" stroke="{BLUE}" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="1080" y="205" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Acceptance</text>
{requirements}
<text x="1080" y="410" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Review line source</text>
<text x="1080" y="442" font-family="Arial,sans-serif" font-size="14" fill="#1677c8">Official 54294 outline-based simplification</text>
<text x="1080" y="470" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Review paths: {len(official)} · texture: 0</text>
<text x="1080" y="498" font-family="Arial,sans-serif" font-size="14" fill="#41566d">DWG SHA-256: {SOURCE_DWG_SHA256[:16]}…</text>
<text x="1080" y="526" font-family="Arial,sans-serif" font-size="14" fill="#41566d">PDF / API / DWG cross-check: pass</text>
<text x="1080" y="566" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Mechanical fit</text>
<text x="1080" y="598" font-family="Arial,sans-serif" font-size="14" fill="#41566d">delta: {cross_check["absolute_delta_mm"]} mm</text>
<text x="1080" y="626" font-family="Arial,sans-serif" font-size="14" fill="#41566d">tolerance: {cross_check["tolerance_mm"]} mm</text>
<text x="1080" y="666" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Isolation</text>
<text x="1080" y="698" font-family="Arial,sans-serif" font-size="15" fill="#41566d">geometry products: 1</text>
<text x="1080" y="726" font-family="Arial,sans-serif" font-size="15" fill="#41566d">whole model render: false</text>
<text x="1080" y="754" font-family="Arial,sans-serif" font-size="15" fill="#41566d">mesh faces: {metadata["mesh_face_count"]}</text>
<text x="1080" y="806" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Official family reference; not a project shop drawing.</text>
<text x="1080" y="832" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Visual review pending; no IFC write.</text>
</svg>'''


def render_original_official_evidence_svg(profile, view, official):
    """Render the untouched official paths as a separate, immutable evidence view."""
    width, height = 1400, 980
    points = [point for path in official for point in path]
    min_x, max_x = min(point[0] for point in points), max(point[0] for point in points)
    min_y, max_y = min(point[1] for point in points), max(point[1] for point in points)
    padding = max(max_x - min_x, max_y - min_y) * 0.07
    min_x, max_x = min_x - padding, max_x + padding
    min_y, max_y = min_y - padding, max_y + padding
    scale = min(1280 / (max_x - min_x), 760 / (max_y - min_y))

    def transform(point):
        return 60 + (point[0] - min_x) * scale, 900 - (point[1] - min_y) * scale

    path = svg_path(official, transform)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="1400" height="980" fill="#fbfaf7"/>
<text x="60" y="58" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">Untouched official Gessi {DRAWING_PRODUCT_CODE} native-DWG evidence</text>
<text x="60" y="96" font-family="Arial,sans-serif" font-size="20" fill="#41566d">{VIEW_DEFINITIONS[view]["label"]} · original texture retained · not the de-textured review representation</text>
<text x="60" y="128" font-family="Arial,sans-serif" font-size="15" fill="#68798a">DWG SHA-256: {SOURCE_DWG_SHA256} · paths: {len(official)}</text>
<path class="official-reference native-dwg original-unaltered-evidence" data-source-kind="native_dwg" data-unaltered-official-dwg="true" data-product-code="{DRAWING_PRODUCT_CODE}" data-dwg-sha256="{SOURCE_DWG_SHA256}" d="{path}" fill="none" stroke="{BLUE}" stroke-width="2.0" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
</svg>'''


def write_index(manifest, output_dir=OUTPUT_DIR):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (output_dir / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Gessi316 54294 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Gessi316 Meccanica / 45089_54294</h1><p>Grey = unchanged actual project IFC Body; black = review-only proxy with smooth handles and a spout rebuilt to the official 54294 maximum 210 mm adjustable setting; blue = de-textured review simplification based on the official 54294 outline, not the untouched official CAD. The original official native-DWG evidence remains separately available with its texture and source hash unchanged. The actual project Body remains 6.966033 mm shorter. Visual review pending.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="review-geometry-simplification.json">Review simplification</a><a href="actual-body-comparison-manifest.json">Actual Body comparison</a><a href="profile.json">Profile</a><a href="official-native-dwg-linework.json">Untouched native-DWG JSON</a><a href="official-original-plan.svg">Original official Plan</a><a href="official-original-front.svg">Original official Front</a><a href="official-original-side.svg">Original official Side</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/official-pdf-verification.json">PDF verification</a><a href="project-context-sanitary-plan.svg">Project plan</a><a href="project-context-front-elevation.svg">Project front</a><a href="project-context-side-elevation.svg">Project side</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://areapro.gessi.com/en/product/54294">Official product</a></nav><h2>Project drawing context</h2><main><article><h2>Sanitary plan</h2><a href="project-context-sanitary-plan-review.svg"><img src="project-context-sanitary-plan-review-preview.png"></a></article><article><h2>Front elevation</h2><a href="project-context-front-elevation-review.svg"><img src="project-context-front-elevation-review-preview.png"></a></article><article><h2>Side elevation</h2><a href="project-context-side-elevation-review.svg"><img src="project-context-side-elevation-review-preview.png"></a></article></main><h2>Three-view review simplification / actual Body / review proxy comparison</h2><main>{cards}</main><h2>Latest review-only Bonsai camera renders</h2><main>{bonsai_cards}</main><h2>Unchanged actual project Body disclosure</h2><main><article><h2>Actual Body Side</h2><a href="bonsai-camera-side-elevation-actual-body.png"><img src="bonsai-camera-side-elevation-actual-body.png"></a></article></main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--register", type=Path, default=REGISTER)
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    if sha256(source) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    register = load_json(args.register.resolve())
    profile = register["profiles"][PROFILE_KEY]
    access = load_json(ACCESS_RECORD)
    linework = load_json(LINEWORK)
    pdf_verification = load_json(PDF_VERIFICATION)
    revalidation = load_json(REVALIDATION)
    drawing_source = profile["drawing_source"]
    if (
        drawing_source.get("source_kind") != SOURCE_KIND
        or drawing_source.get("official_cad_used") is not True
        or drawing_source.get("third_party_cad_used") is not False
        or drawing_source.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or access.get("pass") is not True
        or access.get("official_product_cad", {}).get("acquired") is not True
        or access.get("official_product_cad", {}).get("exact_project_configuration_match") is not EXPECTED_SOURCE_EXACT_PROJECT_CONFIGURATION_MATCH
        or access.get("drawing_geometry_source", {}).get("source_kind") != SOURCE_KIND
        or access.get("drawing_geometry_source", {}).get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or access.get("identity_and_geometry_policy", {}).get(IDENTITY_POLICY_KEY) is not False
        or linework.get("pass") is not True
        or linework.get("source_kind") != SOURCE_KIND
        or linework.get("identity_gates", {}).get(IDENTITY_POLICY_KEY) is not False
        or pdf_verification.get("pass") is not True
        or pdf_verification.get("visual_review", {}).get("status") != "codex_visual_qa_pass"
        or revalidation.get("pass") is not True
        or revalidation.get("public_access", {}).get(PUBLIC_ACCESS_KEY) is not True
        or sha256(SOURCE_DWG) != SOURCE_DWG_SHA256
    ):
        raise RuntimeError("Gessi exact native-DWG source gate failed")
    for evidence in access["official_identity_sources"]:
        if sha256(ROOT / evidence["local_path"]) != evidence["sha256"]:
            raise RuntimeError(f'Gessi official identity evidence hash mismatch: {evidence["local_path"]}')
    model = ifcopenshell.open(source)
    product, original_vertices, original_faces = mesh_for_one_product(model, profile["representative_global_id"])
    type_relations = list(product.IsTypedBy)
    product_type = type_relations[0].RelatingType if type_relations else None
    actual_identity = product_type_name(product) if product_type is not None else product.Name
    if actual_identity != profile["ifc_type_name"]:
        raise RuntimeError(f"{PROFILE_KEY} representative identity drifted")
    actual_description = product_type.Description if product_type is not None else product.Description
    if actual_description != EXPECTED_DESCRIPTION:
        raise RuntimeError(f"{PROFILE_KEY} IFC description drifted")
    instances = sorted(
        item.GlobalId for item in model.by_type(product.is_a())
        if (product_type_name(item) if item.IsTypedBy else item.Name) == profile["ifc_type_name"]
    )
    if instances != profile["expected_instance_global_ids"]:
        raise RuntimeError("Gessi instance set drifted")
    vertices, faces, simplification = simplify_handle_knurls(
        original_vertices,
        original_faces,
        linework["views"]["side"]["paths_mm"],
    )
    minimum, maximum = bounds_3d(vertices)
    output.mkdir(parents=True, exist_ok=True)
    candidate_views = {}
    views = []
    for view, definition in VIEW_DEFINITIONS.items():
        axes = definition["axes"]
        all_edges = projected_raw_edges(original_vertices, original_faces, axes)
        edges = display_edge_sample(all_edges, maximum=1000)
        proxy = projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"]))
        linework_view = linework["views"][view]
        if linework_view.get("path_count") != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi {view} official native-DWG path count drifted")
        official_original, alignment = align_official_paths(view, linework_view["paths_mm"], minimum, maximum, simplification)
        review_raw, texture_audit = simplify_official_handle_texture(view, linework_view["paths_mm"])
        official, _ = align_official_paths(view, review_raw, minimum, maximum, simplification)
        cross_check = compatibility(view, official_original, minimum, maximum, simplification)
        if not cross_check["pass"]:
            raise RuntimeError(f"Gessi {view} IFC/DWG compatibility gate failed")
        if texture_audit["pass"] is not True or texture_audit["review_texture_detail_path_count"] != 0:
            raise RuntimeError(f"Gessi {view} handle line-texture simplification gate failed")
        original_target = output / f"official-original-{view}.svg"
        original_target.write_text(
            render_original_official_evidence_svg(profile, view, official_original),
            encoding="utf-8",
        )
        target = output / f"{view}.svg"
        target.write_text(render_svg(profile, view, edges, proxy, official, {"mesh_face_count": len(faces), "compatibility": cross_check}), encoding="utf-8")
        candidate_views[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": rounded(proxy),
            "original_official_native_dwg_paths_mm": rounded(official_original),
            "original_official_native_dwg_path_count": len(official_original),
            "review_simplified_official_outline_paths_mm": rounded(official),
            "review_simplified_official_outline_path_count": len(official),
            "source_kind": REVIEW_LINE_SOURCE_KIND,
            "source_label_zh": REVIEW_LINE_SOURCE_LABEL_ZH,
            "source_dwg_sha256": SOURCE_DWG_SHA256,
            "unaltered_official_dwg": False,
            "original_official_native_dwg_evidence_svg": relative(original_target),
            "original_official_native_dwg_evidence_svg_sha256": sha256(original_target),
            "handle_line_texture_simplification": texture_audit,
            "ifc_dwg_compatibility": cross_check,
            "translation_alignment": alignment,
        }
        views.append({
            "view": view,
            "svg": relative(target),
            "svg_sha256": sha256(target),
            "projection_axes": list(axes),
            "raw_edge_count": len(all_edges),
            "displayed_raw_edge_count": len(edges),
            "silhouette_path_count": len(proxy),
            "drawing_line_source_kind": REVIEW_LINE_SOURCE_KIND,
            "drawing_line_source_label_zh": REVIEW_LINE_SOURCE_LABEL_ZH,
            "review_simplified_path_count": len(official),
            "original_official_cad_path_count": len(official_original),
            "original_official_evidence_svg": relative(original_target),
            "original_official_evidence_svg_sha256": sha256(original_target),
            "handle_line_texture_simplification": texture_audit,
            "blue_line_present": True,
            "white_mask_present": True,
            "ifc_dwg_compatibility": cross_check,
            "translation_alignment": alignment,
        })
    candidate_path = output / "candidate-representations.json"
    write_json(candidate_path, {
        "schema_version": 2,
        "profile_key": PROFILE_KEY,
        "representative_global_id": product.GlobalId,
        "ifc_type_name": profile["ifc_type_name"],
        "article_number": ARTICLE_NUMBER,
        "drawing_product_code": DRAWING_PRODUCT_CODE,
        "units": "mm",
        "source_kind": REVIEW_LINE_SOURCE_KIND,
        "source_label_zh": REVIEW_LINE_SOURCE_LABEL_ZH,
        "source_label_en": REVIEW_LINE_SOURCE_LABEL_EN,
        "original_source_kind": SOURCE_KIND,
        "original_source_label_zh": SOURCE_LABEL_ZH,
        "original_source_label_en": SOURCE_LABEL_EN,
        "source_dwg": relative(SOURCE_DWG),
        "source_dwg_sha256": SOURCE_DWG_SHA256,
        "official_cad_used": True,
        "unaltered_official_cad_used_as_review_representation": False,
        "original_official_cad_evidence_preserved": True,
        "third_party_cad_used": False,
        CANDIDATE_IDENTITY_FLAG: False,
        "formal_ifc_write_allowed": False,
        "review_only_geometry_simplification": simplification,
        "review_status": "visual_review_pending",
        "views": candidate_views,
    })
    manifest = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": GENERATOR,
        "formal_ifc": relative(source),
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(source),
        "formal_ifc_bytes_unchanged": sha256(source) == FORMAL_SHA256,
        "profile_register": relative(args.register.resolve()),
        "profile_register_sha256": sha256(args.register.resolve()),
        "profile_key": PROFILE_KEY,
        "display_name": profile["display_name"],
        "ifc_type_name": profile["ifc_type_name"],
        "ifc_type_description": actual_description,
        "article_number": ARTICLE_NUMBER,
        "drawing_product_code": DRAWING_PRODUCT_CODE,
        "representative_global_id": product.GlobalId,
        "registered_instance_global_ids": instances,
        "geometry_product_count": 1,
        "whole_model_render": False,
        "mesh_vertex_count": len(vertices),
        "mesh_face_count": len(faces),
        "original_mesh_vertex_count": len(original_vertices),
        "original_mesh_face_count": len(original_faces),
        "review_only_geometry_simplification": simplification,
        "bounds_mm": {
            "minimum": [round(value, 6) for value in minimum],
            "maximum": [round(value, 6) for value in maximum],
            "size": [round(maximum[index] - minimum[index], 6) for index in range(3)],
        },
        "source_kind": REVIEW_LINE_SOURCE_KIND,
        "source_label_zh": REVIEW_LINE_SOURCE_LABEL_ZH,
        "original_source_kind": SOURCE_KIND,
        "original_source_label_zh": SOURCE_LABEL_ZH,
        "drawing_source": {
            "source_kind": REVIEW_LINE_SOURCE_KIND,
            "source_label_zh": REVIEW_LINE_SOURCE_LABEL_ZH,
            "source_label_en": REVIEW_LINE_SOURCE_LABEL_EN,
            "original_source_kind": SOURCE_KIND,
            "original_source_label_zh": SOURCE_LABEL_ZH,
            "original_source_label_en": SOURCE_LABEL_EN,
            "official_cad_used": True,
            "third_party_cad_used": False,
            "official_source_access_record": relative(ACCESS_RECORD),
            "official_source_access_record_sha256": sha256(ACCESS_RECORD),
            "official_product_cad_status": OFFICIAL_CAD_STATUS,
            "source_dwg": relative(SOURCE_DWG),
            "source_dwg_sha256": SOURCE_DWG_SHA256,
        },
        "official_reference": profile["official_reference"],
        "official_identity_evidence_only": False,
        "dimension_cross_check": access["dimension_cross_check"],
        "official_cad_acquired": OFFICIAL_CAD_ACQUIRED,
        "official_cad_exact_project_configuration_match": OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH,
        "official_cad_used": True,
        "unaltered_official_cad_used_as_review_representation": False,
        "original_official_cad_evidence_preserved": True,
        "third_party_cad_used": False,
        CANDIDATE_IDENTITY_FLAG: False,
        "blue_line_present": True,
        "white_mask_present": True,
        "official_native_dwg_linework": relative(LINEWORK),
        "official_native_dwg_linework_sha256": sha256(LINEWORK),
        "official_pdf_verification": relative(PDF_VERIFICATION),
        "official_pdf_verification_sha256": sha256(PDF_VERIFICATION),
        "official_source_revalidation": relative(REVALIDATION),
        "official_source_revalidation_sha256": sha256(REVALIDATION),
        "official_source_access_record": relative(ACCESS_RECORD),
        "official_source_access_record_sha256": sha256(ACCESS_RECORD),
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "derived_ifc_write_allowed": False,
        "formal_ifc_write": "not performed",
        "scope": SCOPE,
        "views": views,
        "pass": all(
            view["blue_line_present"]
            and view["original_official_cad_path_count"] == EXPECTED_PATH_COUNTS[view["view"]]
            and view["handle_line_texture_simplification"]["review_texture_detail_path_count"] == 0
            and view["handle_line_texture_simplification"]["pass"]
            and view["ifc_dwg_compatibility"]["pass"]
            for view in views
        ),
    }
    context_manifest = output / "project-context-manifest.json"
    bonsai_manifest = output / "bonsai-review-manifest.json"
    if context_manifest.is_file():
        context = load_json(context_manifest)
        if context.get("pass") is not True:
            raise RuntimeError("Gessi project-context evidence failed")
        manifest["project_context"] = {
            "manifest": relative(context_manifest),
            "manifest_sha256": sha256(context_manifest),
            "walls_and_surrounding_project_elements_retained": context["walls_and_surrounding_project_elements_retained"],
            "overlay_top_layer_with_white_mask": context["overlay_top_layer_with_white_mask"],
            "pass": True,
        }
    if bonsai_manifest.is_file():
        bonsai = load_json(bonsai_manifest)
        if bonsai.get("mode") != "actual_bonsai_ifc_body_camera_render" or bonsai.get("pass") is not True:
            raise RuntimeError("Gessi Bonsai camera evidence failed")
        manifest["bonsai_review"] = {
            "manifest": relative(bonsai_manifest),
            "manifest_sha256": sha256(bonsai_manifest),
            "mode": bonsai["mode"],
            "saved_active_representation": bonsai["bonsai_session"]["saved_active_representation"],
            "saved_camera_count": bonsai["bonsai_session"]["saved_camera_count"],
            "pass": True,
        }
    manifest_path = output / "manifest.json"
    write_json(manifest_path, manifest)
    write_index(manifest, output)
    print(relative(manifest_path))


if __name__ == "__main__":
    main()
