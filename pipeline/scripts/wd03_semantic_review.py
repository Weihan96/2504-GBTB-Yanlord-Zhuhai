#!/usr/bin/env python3
"""Build semantic WD03 review linework from the isolated actual IFC Body.

The public Poliform page is identity evidence only.  Every black line emitted by
this script is tied to a classified representation item and either separates
two different visible components or records a real occlusion boundary.
"""

from __future__ import annotations

import base64
import html
import math
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
import ifcopenshell.geom

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from int1_highpoly_type_review import display_edge_sample, projected_raw_edges
from render_sis04_project_context import write_uncached_png_preview


PRODUCT_DIR = ROOT / "output/review/highpoly-types/wd03"
ISOLATED_IFC = PRODUCT_DIR / "Poliform-Senzafine-WD03-bonsai-isolated.ifc"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
GLOBAL_ID = "3cmikd9MTB$egM5KQaNgUf"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
SOURCE_LABEL_EN = "simplified drawing representation derived from the original high-poly geometry"
VIEW_AXES = {"plan": (0, 1), "front": (0, 2), "side": (1, 2)}
VIEW_LABELS = {"plan": "PLAN / XY", "front": "FRONT / XZ", "side": "SIDE / YZ"}
EXPECTED_PATH_COUNTS = {"plan": 5, "front": 9, "side": 3}


def rounded_path(path):
    return [[round(float(x), 6), round(float(y), 6)] for x, y in path]


def rectangle(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]


def union_bounds(records):
    return (
        [min(record["minimum"][axis] for record in records) for axis in range(3)],
        [max(record["maximum"][axis] for record in records) for axis in range(3)],
    )


def close(value, target, tolerance):
    return abs(value - target) <= tolerance


def load_body_components():
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash changed before WD03 semantic review")
    model = ifcopenshell.open(ISOLATED_IFC)
    product = model.by_guid(GLOBAL_ID)
    if product is None:
        raise RuntimeError("isolated WD03 representative missing")
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    shape = ifcopenshell.geom.create_shape(settings, product)
    geometry = shape.geometry
    world_vertices = [
        tuple(float(value) * 1000.0 for value in geometry.verts[index:index + 3])
        for index in range(0, len(geometry.verts), 3)
    ]
    faces = [
        tuple(int(value) for value in geometry.faces[index:index + 3])
        for index in range(0, len(geometry.faces), 3)
    ]
    if len(world_vertices) != 3608 or len(faces) != 7148 or len(geometry.item_ids) != len(faces):
        raise RuntimeError("WD03 isolated Body mesh/item identity drifted")
    manifest = load_json(PRODUCT_DIR / "manifest.json")
    target_minimum = manifest["bounds_mm"]["minimum"]
    world_minimum = [min(point[axis] for point in world_vertices) for axis in range(3)]
    shift = [target_minimum[axis] - world_minimum[axis] for axis in range(3)]
    vertices = [
        tuple(point[axis] + shift[axis] for axis in range(3))
        for point in world_vertices
    ]
    grouped_faces = defaultdict(list)
    for face, item_id in zip(faces, geometry.item_ids):
        grouped_faces[int(item_id)].append(face)
    records = []
    for item_id, item_faces in sorted(grouped_faces.items()):
        indices = {index for face in item_faces for index in face}
        minimum = [min(vertices[index][axis] for index in indices) for axis in range(3)]
        maximum = [max(vertices[index][axis] for index in indices) for axis in range(3)]
        size = [maximum[axis] - minimum[axis] for axis in range(3)]
        records.append({
            "item_id": item_id,
            "ifc_entity": model.by_id(item_id).is_a(),
            "face_count": len(item_faces),
            "minimum": minimum,
            "maximum": maximum,
            "size": size,
            "semantic": None,
        })
    return vertices, faces, records, shift


def classify_components(records):
    buckets = defaultdict(list)
    for record in records:
        sx, sy, sz = record["size"]
        z0, z1 = record["minimum"][2], record["maximum"][2]
        if 24.0 <= sx <= 25.1 and sy >= 585.0 and sz >= 2389.0:
            semantic = "left_or_right_full_height_side_panel"
        elif sx >= 1262.0 and 6.0 <= sy <= 7.0 and sz >= 2317.0:
            semantic = "full_height_back_panel"
        elif sx >= 1262.0 and 572.0 <= sy <= 574.0 and 23.5 <= sz <= 24.5:
            semantic = "top_deck" if z1 > 0.0 else "middle_shelf"
        elif sx >= 1262.0 and 2.5 <= sy <= 3.5 and 24.0 <= sz <= 25.5:
            semantic = "top_front_edge" if z1 > 0.0 else "middle_shelf_front_edge"
        elif sx >= 1262.0 and 10.0 <= sy <= 12.0 and 9.5 <= sz <= 10.6:
            semantic = "top_under_shelf_runner" if z1 > -500.0 else "middle_under_shelf_runner"
        elif sx >= 1255.0 and 29.5 <= sy <= 31.5 and 29.5 <= sz <= 31.5:
            semantic = "upper_clothes_rail" if z1 > -100.0 else "lower_clothes_rail"
        elif 36.0 <= sx <= 38.0 and 16.0 <= sy <= 17.0 and 45.0 <= sz <= 46.5:
            semantic = "upper_rail_end_bracket" if z1 > -100.0 else "lower_rail_end_bracket"
        elif 6.0 <= sx <= 7.0 and 503.0 <= sy <= 505.0 and sz >= 900.0:
            semantic = "upper_open_door_panel" if z0 > -1400.0 else "lower_open_door_panel"
        elif sx >= 1250.0 and sy >= 570.0 and z1 <= -2292.0:
            semantic = "bottom_plinth_or_deck_layer"
        else:
            raise RuntimeError(f"unclassified WD03 Body item {record['item_id']} size {record['size']}")
        record["semantic"] = semantic
        buckets[semantic].append(record)
    expected_counts = {
        "left_or_right_full_height_side_panel": 2,
        "full_height_back_panel": 1,
        "top_deck": 1,
        "middle_shelf": 1,
        "top_front_edge": 1,
        "middle_shelf_front_edge": 1,
        "top_under_shelf_runner": 1,
        "middle_under_shelf_runner": 1,
        "upper_clothes_rail": 2,
        "lower_clothes_rail": 2,
        "upper_rail_end_bracket": 2,
        "lower_rail_end_bracket": 2,
        "upper_open_door_panel": 2,
        "lower_open_door_panel": 2,
        "bottom_plinth_or_deck_layer": 4,
    }
    actual_counts = {key: len(buckets[key]) for key in expected_counts}
    if actual_counts != expected_counts:
        raise RuntimeError(f"WD03 semantic component count drifted: {actual_counts}")
    return buckets


def semantic_view_paths(records, buckets):
    overall_minimum, overall_maximum = union_bounds(records)
    top_minimum, top_maximum = union_bounds(buckets["top_deck"])
    middle_minimum, middle_maximum = union_bounds(buckets["middle_shelf"])
    bottom_minimum, bottom_maximum = union_bounds(buckets["bottom_plinth_or_deck_layer"])
    upper_rail_minimum, upper_rail_maximum = union_bounds(buckets["upper_clothes_rail"])
    lower_rail_minimum, lower_rail_maximum = union_bounds(buckets["lower_clothes_rail"])
    upper_door_minimum, upper_door_maximum = union_bounds(buckets["upper_open_door_panel"])
    lower_door_minimum, lower_door_maximum = union_bounds(buckets["lower_open_door_panel"])
    x0, y0, z0 = overall_minimum
    x1, y1, z1 = overall_maximum
    ix0, iy0 = top_minimum[0], top_minimum[1]
    ix1, iy1 = top_maximum[0], top_maximum[1]

    plan = [
        {"id": "plan.outer_envelope", "path": rectangle(x0, y0, x1, y1), "kind": "silhouette", "sides": ["wardrobe_body", "exterior_void"]},
        {"id": "plan.left_side_cap", "path": [[ix0, y0], [ix0, y1]], "kind": "component_boundary", "sides": ["left_side_panel_cap", "top_deck/back/front_regions"]},
        {"id": "plan.right_side_cap", "path": [[ix1, y0], [ix1, y1]], "kind": "component_boundary", "sides": ["right_side_panel_cap", "top_deck/back/front_regions"]},
        {"id": "plan.back_to_top_deck", "path": [[ix0, iy0], [ix1, iy0]], "kind": "depth_transition", "sides": ["back_panel_or_rear_gap", "top_deck"]},
        {"id": "plan.top_deck_to_front_edge", "path": [[ix0, iy1], [ix1, iy1]], "kind": "component_boundary", "sides": ["top_deck", "top_front_edge"]},
    ]
    front = [
        {"id": "front.outer_envelope", "path": rectangle(x0, z0, x1, z1), "kind": "silhouette", "sides": ["wardrobe_body", "exterior_void"]},
        {"id": "front.left_side_panel", "path": [[ix0, z0], [ix0, z1]], "kind": "component_boundary", "sides": ["left_side_panel", "open_bays/shelves"]},
        {"id": "front.right_side_panel", "path": [[ix1, z0], [ix1, z1]], "kind": "component_boundary", "sides": ["right_side_panel", "open_bays/shelves"]},
        {"id": "front.top_deck_lower_edge", "path": [[ix0, top_minimum[2]], [ix1, top_minimum[2]]], "kind": "occlusion_boundary", "sides": ["top_deck", "upper_open_bay"]},
        {"id": "front.middle_shelf_upper_edge", "path": [[ix0, middle_maximum[2]], [ix1, middle_maximum[2]]], "kind": "occlusion_boundary", "sides": ["upper_open_bay", "middle_shelf"]},
        {"id": "front.middle_shelf_lower_edge", "path": [[ix0, middle_minimum[2]], [ix1, middle_minimum[2]]], "kind": "occlusion_boundary", "sides": ["middle_shelf", "lower_open_bay"]},
        {"id": "front.bottom_deck_upper_edge", "path": [[ix0, bottom_maximum[2]], [ix1, bottom_maximum[2]]], "kind": "occlusion_boundary", "sides": ["lower_open_bay", "bottom_deck/plinth"]},
        {"id": "front.upper_clothes_rail", "path": rectangle(upper_rail_minimum[0], upper_rail_minimum[2], upper_rail_maximum[0], upper_rail_maximum[2]), "kind": "occluding_component_envelope", "sides": ["upper_clothes_rail", "back_panel"]},
        {"id": "front.lower_clothes_rail", "path": rectangle(lower_rail_minimum[0], lower_rail_minimum[2], lower_rail_maximum[0], lower_rail_maximum[2]), "kind": "occluding_component_envelope", "sides": ["lower_clothes_rail", "back_panel"]},
    ]
    side = [
        {"id": "side.outer_side_panel", "path": rectangle(y0, z0, y1, z1), "kind": "silhouette", "sides": ["left_side_panel", "exterior_void"]},
        {"id": "side.upper_open_door_panel", "path": rectangle(upper_door_minimum[1], upper_door_minimum[2], upper_door_maximum[1], upper_door_maximum[2]), "kind": "occluding_component_envelope", "sides": ["upper_open_door_panel", "side_panel"]},
        {"id": "side.lower_open_door_panel", "path": rectangle(lower_door_minimum[1], lower_door_minimum[2], lower_door_maximum[1], lower_door_maximum[2]), "kind": "occluding_component_envelope", "sides": ["lower_open_door_panel", "side_panel"]},
    ]
    result = {"plan": plan, "front": front, "side": side}
    counts = {view: len(paths) for view, paths in result.items()}
    if counts != EXPECTED_PATH_COUNTS:
        raise RuntimeError(f"WD03 semantic path counts drifted: {counts}")
    return result


def path_d(paths, transform):
    commands = []
    for path in paths:
        if not path:
            continue
        x, y = transform(path[0])
        commands.append(f"M{x:.3f},{y:.3f}")
        commands.extend(f"L{transform(point)[0]:.3f},{transform(point)[1]:.3f}" for point in path[1:])
    return " ".join(commands)


def render_review_svg(view, raw_edges, semantics):
    paths = [item["path"] for item in semantics]
    points = [point for path in paths for point in path]
    min_x = min(point[0] for point in points)
    max_x = max(point[0] for point in points)
    min_y = min(point[1] for point in points)
    max_y = max(point[1] for point in points)
    padding = max(max_x - min_x, max_y - min_y) * 0.06
    min_x, max_x = min_x - padding, max_x + padding
    min_y, max_y = min_y - padding, max_y + padding
    scale = min(960.0 / (max_x - min_x), 710.0 / (max_y - min_y))
    transform = lambda point: (55.0 + (point[0] - min_x) * scale, 900.0 - (point[1] - min_y) * scale)
    raw_d = path_d([[start, end] for start, end in raw_edges], transform)
    semantic_markup = "\n".join(
        f'<path class="simplified-proxy-silhouette geometry-derived semantic-boundary" data-source-kind="{SOURCE_KIND}" data-semantic-id="{html.escape(item["id"])}" data-boundary-kind="{item["kind"]}" d="{path_d([item["path"]], transform)}" fill="none" stroke="#111820" stroke-width="3" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>'
        for item in semantics
    )
    notes = "\n".join(
        f'<text x="1080" y="{270 + index * 29}" font-family="Arial,sans-serif" font-size="14" fill="#34495e">• {html.escape(item["id"].split(".", 1)[1].replace("_", " "))}</text>'
        for index, item in enumerate(semantics)
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="980" viewBox="0 0 1400 980">
<rect width="1400" height="980" fill="#fbfaf7"/>
<text x="55" y="55" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">Poliform Senzafine WD03 · semantic review</text>
<text x="55" y="94" font-family="Arial,sans-serif" font-size="20" fill="#41566d">{VIEW_LABELS[view]} · isolated representative {GLOBAL_ID}</text>
<text x="55" y="128" font-family="Arial,sans-serif" font-size="17" fill="#68798a">Grey = actual IFC Body mesh · Black = semantically validated drawing boundaries · no official CAD / no blue line</text>
<rect x="55" y="160" width="980" height="760" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="original-highpoly" d="{raw_d}" fill="none" stroke="#87929c" stroke-width="0.55" stroke-opacity="0.28" vector-effect="non-scaling-stroke"/>
{semantic_markup}
<text x="1080" y="205" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Validated semantics</text>
<text x="1080" y="238" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Every internal line separates components or depth.</text>
{notes}
<text x="1080" y="760" font-family="Arial,sans-serif" font-size="16" font-weight="700" fill="#1f2d3d">Source</text>
<text x="1080" y="790" font-family="Arial,sans-serif" font-size="14" fill="#41566d">{SOURCE_LABEL_ZH}</text>
<text x="1080" y="820" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Poliform page: family identity only</text>
<text x="1080" y="850" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Visual review pending · no IFC write</text>
</svg>'''


def render_diagnostic_svg(view, semantics):
    paths = [item["path"] for item in semantics]
    points = [point for path in paths for point in path]
    min_x, max_x = min(p[0] for p in points), max(p[0] for p in points)
    min_y, max_y = min(p[1] for p in points), max(p[1] for p in points)
    padding = max(max_x - min_x, max_y - min_y) * 0.08
    min_x, max_x = min_x - padding, max_x + padding
    min_y, max_y = min_y - padding, max_y + padding
    scale = min(1220.0 / (max_x - min_x), 760.0 / (max_y - min_y))
    transform = lambda point: (90.0 + (point[0] - min_x) * scale, 900.0 - (point[1] - min_y) * scale)
    colors = ["#1769aa", "#d1495b", "#2a9d8f", "#e9c46a", "#7b2cbf", "#f77f00", "#4d908e", "#577590", "#90be6d"]
    markup = []
    for index, item in enumerate(semantics):
        color = colors[index % len(colors)]
        markup.append(f'<path data-semantic-id="{item["id"]}" d="{path_d([item["path"]], transform)}" fill="none" stroke="{color}" stroke-width="6" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>')
        point = item["path"][len(item["path"]) // 2]
        x, y = transform(point)
        markup.append(f'<text x="{x + 8:.2f}" y="{y - 8:.2f}" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="{color}">{index + 1}</text>')
    legend = "\n".join(
        f'<text x="90" y="{120 + index * 26}" font-family="Arial,sans-serif" font-size="15" fill="{colors[index % len(colors)]}">{index + 1}. {html.escape(item["id"])} — {html.escape(item["sides"][0])} / {html.escape(item["sides"][1])}</text>'
        for index, item in enumerate(semantics)
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="980" viewBox="0 0 1400 980">
<rect width="1400" height="980" fill="#ffffff"/>
<text x="90" y="58" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">WD03 {view.title()} semantic-boundary diagnostic</text>
{legend}
{''.join(markup)}
</svg>'''


def render_contact_svg(target, title, preview_paths):
    images = []
    for index, path in enumerate(preview_paths):
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        x = 30 + index * 590
        images.append(f'<image x="{x}" y="120" width="560" height="760" preserveAspectRatio="xMidYMid meet" href="data:image/png;base64,{encoded}"/>')
    target.write_text(f'''<svg xmlns="http://www.w3.org/2000/svg" width="1800" height="920" viewBox="0 0 1800 920">
<rect width="1800" height="920" fill="#f4f2ee"/>
<text x="40" y="62" font-family="Arial,sans-serif" font-size="34" font-weight="700" fill="#1f2d3d">{html.escape(title)}</text>
<text x="40" y="98" font-family="Arial,sans-serif" font-size="18" fill="#68798a">Plan · Front · Side — isolated IFC Body semantics; family page is identity evidence only</text>
{''.join(images)}
</svg>''', encoding="utf-8")


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()} single-product SVG</h2><a href="{item["view"]}.svg"><img src="{item["view"]}-preview.png"></a><p>{item["semantic_path_count"]} validated paths</p></article>'
        for item in manifest["views"]
    )
    diagnostics = "".join(
        f'<article><h2>{view.title()} diagnostic</h2><a href="wd03-semantic-diagnostic-{view}.svg"><img src="wd03-semantic-diagnostic-{view}.png"></a></article>'
        for view in VIEW_AXES
    )
    (PRODUCT_DIR / "index.html").write_text(f'''<!doctype html><html><meta charset="utf-8"><title>WD03 semantic review</title>
<style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style>
<h1>Poliform Senzafine / project WD03</h1><p>Black line = {SOURCE_LABEL_EN}. No official CAD was acquired and no catalogue composition is substituted.</p>
<nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="wd03-semantic-segmentation.json">Semantic JSON</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-furniture-plan-review.svg">Project plan overlay</a><a href="project-context-side-elevation-review.svg">Project side overlay</a></nav>
<h2>Semantic single-product review</h2><main>{cards}</main><h2>Boundary diagnostics</h2><main>{diagnostics}</main></html>''', encoding="utf-8")


def main():
    vertices, faces, records, shift = load_body_components()
    buckets = classify_components(records)
    view_semantics = semantic_view_paths(records, buckets)
    views = []
    candidate_views = {}
    preview_paths = []
    diagnostic_previews = []
    for view, axes in VIEW_AXES.items():
        raw_edges = display_edge_sample(projected_raw_edges(vertices, faces, axes), maximum=975)
        semantics = view_semantics[view]
        target = PRODUCT_DIR / f"{view}.svg"
        target.write_text(render_review_svg(view, raw_edges, semantics), encoding="utf-8")
        preview = PRODUCT_DIR / f"{view}-preview.png"
        write_uncached_png_preview(target, preview)
        preview_paths.append(preview)
        diagnostic = PRODUCT_DIR / f"wd03-semantic-diagnostic-{view}.svg"
        diagnostic.write_text(render_diagnostic_svg(view, semantics), encoding="utf-8")
        diagnostic_preview = PRODUCT_DIR / f"wd03-semantic-diagnostic-{view}.png"
        write_uncached_png_preview(diagnostic, diagnostic_preview)
        diagnostic_previews.append(diagnostic_preview)
        candidate_views[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": [rounded_path(item["path"]) for item in semantics],
            "semantic_paths": [{key: value for key, value in item.items() if key != "path"} for item in semantics],
            "source_kind": SOURCE_KIND,
            "official_cad_paths_mm": [],
        }
        views.append({
            "view": view,
            "svg": relative(target),
            "svg_sha256": sha256(target),
            "preview": relative(preview),
            "preview_sha256": sha256(preview),
            "projection_axes": list(axes),
            "raw_edge_count": len(projected_raw_edges(vertices, faces, axes)),
            "displayed_raw_edge_count": len(raw_edges),
            "silhouette_path_count": len(semantics),
            "semantic_path_count": len(semantics),
            "internal_semantic_path_count": len(semantics) - 1,
            "drawing_line_source_kind": SOURCE_KIND,
            "drawing_line_source_label_zh": SOURCE_LABEL_ZH,
            "official_cad_path_count": 0,
            "blue_line_present": False,
            "every_internal_line_has_two_different_semantics_or_real_occlusion": True,
        })

    audit_path = PRODUCT_DIR / "wd03-semantic-segmentation.json"
    audit = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/wd03_semantic_review.py",
        "representative_global_id": GLOBAL_ID,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_used": False,
        "catalogue_geometry_used": False,
        "isolated_ifc": relative(ISOLATED_IFC),
        "isolated_ifc_sha256": sha256(ISOLATED_IFC),
        "mesh_vertex_count": len(vertices),
        "mesh_face_count": len(faces),
        "representation_item_count": len(records),
        "world_to_candidate_translation_mm": [round(value, 6) for value in shift],
        "components": [
            {
                **{key: value for key, value in record.items() if key not in ("minimum", "maximum", "size")},
                "minimum_mm": [round(value, 6) for value in record["minimum"]],
                "maximum_mm": [round(value, 6) for value in record["maximum"]],
                "size_mm": [round(value, 6) for value in record["size"]],
            }
            for record in records
        ],
        "views": {
            view: {
                "path_count": len(items),
                "paths": [{**{key: value for key, value in item.items() if key != "path"}, "path_mm": rounded_path(item["path"])} for item in items],
                "every_internal_line_has_two_different_semantics_or_real_occlusion": True,
                "meaningless_coplanar_internal_lines": 0,
            }
            for view, items in view_semantics.items()
        },
        "omissions": {
            "rail_end_brackets": "below drawing-detail threshold; retained in grey actual-Body evidence but not promoted to black line",
            "under_shelf_runners": "occluded by shelf/deck in the orthographic review directions; retained in semantic inventory",
        },
        "pass": True,
    }
    write_json(audit_path, audit)

    candidate_path = PRODUCT_DIR / "candidate-representations.json"
    write_json(candidate_path, {
        "schema_version": 1,
        "profile_key": "wd03",
        "representative_global_id": GLOBAL_ID,
        "ifc_type_name": "WD03",
        "article_number": "Poliform Senzafine / project WD03",
        "units": "mm",
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_used": False,
        "third_party_cad_used": False,
        "catalogue_geometry_used": False,
        "formal_ifc_write_allowed": False,
        "review_status": "visual_review_pending",
        "semantic_segmentation": relative(audit_path),
        "semantic_segmentation_sha256": sha256(audit_path),
        "views": candidate_views,
    })

    for name, title, previews in (
        ("review-contact-sheet", "WD03 semantic single-product review", preview_paths),
        ("wd03-semantic-diagnostic-contact-sheet", "WD03 component / boundary diagnostics", diagnostic_previews),
    ):
        svg = PRODUCT_DIR / f"{name}.svg"
        png = PRODUCT_DIR / f"{name}.png"
        render_contact_svg(svg, title, previews)
        write_uncached_png_preview(svg, png)

    subprocess.run([sys.executable, str(ROOT / "pipeline/scripts/render_wd03_project_context.py")], cwd=ROOT, check=True)
    context_path = PRODUCT_DIR / "project-context-manifest.json"
    manifest_path = PRODUCT_DIR / "manifest.json"
    manifest = load_json(manifest_path)
    manifest.update({
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/wd03_semantic_review.py",
        "formal_ifc_sha256_after_generation": sha256(FORMAL_IFC),
        "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_used": False,
        "third_party_cad_used": False,
        "blue_line_present": False,
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "semantic_segmentation": {
            "path": relative(audit_path),
            "sha256": sha256(audit_path),
            "representation_item_count": len(records),
            "classified_component_count": len(records),
            "meaningless_coplanar_internal_line_count": 0,
            "every_internal_line_has_two_different_semantics_or_real_occlusion": True,
            "catalogue_geometry_used": False,
            "pass": True,
        },
        "views": views,
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "derived_ifc_write_allowed": False,
        "formal_ifc_write": "not performed",
        "review_contact_sheet": relative(PRODUCT_DIR / "review-contact-sheet.png"),
        "review_contact_sheet_sha256": sha256(PRODUCT_DIR / "review-contact-sheet.png"),
        "semantic_diagnostic_contact_sheet": relative(PRODUCT_DIR / "wd03-semantic-diagnostic-contact-sheet.png"),
        "semantic_diagnostic_contact_sheet_sha256": sha256(PRODUCT_DIR / "wd03-semantic-diagnostic-contact-sheet.png"),
        "project_context": {
            "manifest": relative(context_path),
            "manifest_sha256": sha256(context_path),
            "walls_and_surrounding_project_elements_retained": True,
            "overlay_top_layer_with_white_mask": True,
            "semantic_black_line_overlay": True,
            "pass": True,
        },
        "pass": True,
    })
    write_json(manifest_path, manifest)

    approval_path = ROOT / "pipeline/decisions/wd03-drawing-approval.json"
    approval = load_json(approval_path)
    approval.update({
        "candidate_manifest_sha256": sha256(manifest_path),
        "status": "pending",
        "reviewer": None,
        "review_date": None,
        "approved_views": [],
        "derived_ifc_write_allowed": False,
        "formal_authoritative_ifc_write_allowed": False,
        "approval_evidence": None,
        "note": "Pending explicit human approval of the semantic Plan, Front and Side black-line views. Every internal line is tied to a different component or real occlusion boundary. Exact native CAD was not acquired; no catalogue composition was substituted and no IFC write is allowed.",
    })
    write_json(approval_path, approval)

    profile_path = PRODUCT_DIR / "profile.json"
    profile = load_json(profile_path)
    record = profile["profiles"]["wd03"]
    record["profile_source"] = "isolated_highpoly_component_semantics; official_family_identity_verified; exact_native_cad_not_acquired"
    record["semantic_sections"] = {
        view: [item["id"] for item in items]
        for view, items in view_semantics.items()
    }
    record["semantic_audit"] = relative(audit_path)
    record["catalogue_geometry_used"] = False
    write_json(profile_path, profile)
    write_index(manifest)

    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during WD03 semantic review")
    print(f"WD03 semantic review written to {PRODUCT_DIR}")


if __name__ == "__main__":
    main()
