#!/usr/bin/env python3
"""Generate a three-view candidate from exactly one high-poly IFC instance."""

from __future__ import annotations

import argparse
import html
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
import ifcopenshell.geom
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon
from shapely.ops import linemerge, unary_union

from falper_sorgente_linework import EXPECTED, ROOT, SCOPE, load_json, relative, sha256, write_json


DEFAULT_REGISTER = ROOT / "pipeline/decisions/highpoly-drawing-profile-register.json"
DEFAULT_OUTPUT = ROOT / "output/review/highpoly-types"
AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
VIEWS = {
    "plan": {"axes": (0, 1), "label": "PLAN / XY"},
    "front": {"axes": (0, 2), "label": "FRONT / XZ"},
    "side": {"axes": (1, 2), "label": "SIDE / YZ"},
}


def product_type_name(product) -> str | None:
    return next((relation.RelatingType.Name for relation in product.IsTypedBy), None)


def mesh_for_one_product(model, global_id: str):
    product = model.by_guid(global_id)
    if product is None:
        raise RuntimeError(f"representative IFC product not found: {global_id}")
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, False)
    # Isolation invariant: this is the only create_shape call in this script.
    shape = ifcopenshell.geom.create_shape(settings, product)
    raw_vertices = shape.geometry.verts
    raw_faces = shape.geometry.faces
    vertices = [
        (
            float(raw_vertices[index]) * 1000.0,
            float(raw_vertices[index + 1]) * 1000.0,
            float(raw_vertices[index + 2]) * 1000.0,
        )
        for index in range(0, len(raw_vertices), 3)
    ]
    faces = [
        (int(raw_faces[index]), int(raw_faces[index + 1]), int(raw_faces[index + 2]))
        for index in range(0, len(raw_faces), 3)
    ]
    return product, vertices, faces


def bounds_3d(vertices):
    return (
        [min(point[axis] for point in vertices) for axis in range(3)],
        [max(point[axis] for point in vertices) for axis in range(3)],
    )


def projected(point, axes):
    return point[axes[0]], point[axes[1]]


def unique_mesh_edges(faces):
    edges = set()
    for face in faces:
        for start, end in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edges.add(tuple(sorted((start, end))))
    return sorted(edges)


def projected_raw_edges(vertices, faces, axes):
    return [
        (projected(vertices[start], axes), projected(vertices[end], axes))
        for start, end in unique_mesh_edges(faces)
    ]


def display_edge_sample(edges, maximum: int = 600):
    """Keep a deterministic, render-safe sample of original mesh edges."""
    if len(edges) <= maximum:
        return edges
    step = math.ceil(len(edges) / maximum)
    return edges[::step]


def display_path_sample(paths, maximum_points: int = 97):
    """Reduce display vertices without changing the registered source paths."""
    sampled_paths = []
    for path in paths:
        if len(path) <= maximum_points:
            sampled_paths.append(path)
            continue
        step = math.ceil((len(path) - 1) / (maximum_points - 1))
        sampled = path[::step]
        if sampled[-1] != path[-1]:
            sampled.append(path[-1])
        sampled_paths.append(sampled)
    return sampled_paths


def projected_silhouette(vertices, faces, axes, tolerance_mm: float):
    triangles = []
    for face in faces:
        polygon = Polygon([projected(vertices[index], axes) for index in face])
        if polygon.is_valid and polygon.area > 1e-6:
            triangles.append(polygon)
    if not triangles:
        raise RuntimeError("representative mesh produced no projected triangles")
    merged = unary_union(triangles).buffer(0).simplify(tolerance_mm, preserve_topology=True)
    polygons = list(merged.geoms) if isinstance(merged, MultiPolygon) else [merged]
    paths = []
    for polygon in polygons:
        if not isinstance(polygon, Polygon) or polygon.is_empty:
            continue
        paths.append(list(polygon.exterior.coords))
        paths.extend(list(interior.coords) for interior in polygon.interiors)
    return paths


def triangle_plane_segment(triangle, axis: int, value: float, tolerance: float = 1e-7):
    points = []
    for start, end in ((triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0])):
        a = start[axis] - value
        b = end[axis] - value
        if abs(a) <= tolerance and abs(b) <= tolerance:
            points.extend((start, end))
        elif abs(a) <= tolerance:
            points.append(start)
        elif abs(b) <= tolerance:
            points.append(end)
        elif a * b < 0:
            ratio = a / (a - b)
            points.append(
                tuple(start[index] + ratio * (end[index] - start[index]) for index in range(3))
            )
    unique = []
    for point in points:
        if not any(sum((point[i] - existing[i]) ** 2 for i in range(3)) < 1e-10 for existing in unique):
            unique.append(point)
    if len(unique) < 2:
        return None
    if len(unique) == 2:
        return unique[0], unique[1]
    return max(
        ((a, b) for index, a in enumerate(unique) for b in unique[index + 1 :]),
        key=lambda pair: sum((pair[0][i] - pair[1][i]) ** 2 for i in range(3)),
    )


def semantic_section_paths(vertices, faces, axes, section, minimum, maximum):
    axis = AXIS_INDEX[section["axis"]]
    plane = minimum[axis] + (maximum[axis] - minimum[axis]) * float(section["fraction"])
    lines = []
    for face in faces:
        segment = triangle_plane_segment([vertices[index] for index in face], axis, plane)
        if not segment:
            continue
        start, end = projected(segment[0], axes), projected(segment[1], axes)
        if math.dist(start, end) >= 0.5:
            lines.append(LineString((start, end)))
    if not lines:
        return [], plane
    merged = unary_union(lines)
    try:
        merged = linemerge(merged)
    except ValueError:
        pass
    geometries = list(merged.geoms) if isinstance(merged, MultiLineString) else [merged]
    return (
        [list(geometry.coords) for geometry in geometries if isinstance(geometry, LineString) and geometry.length >= 2.0],
        plane,
    )


def native_wfb_paths(profile: dict, view: str, minimum, maximum):
    reference = profile["official_reference"]
    if reference["source_kind"] != "native_dwg" or reference["model_code"] != "WFB":
        raise RuntimeError("review blue line must be native WFB DWG")
    linework_path = ROOT / reference["linework_path"]
    verification_path = ROOT / reference["verification_path"]
    linework = load_json(linework_path)
    verification = load_json(verification_path)
    if linework["source_kind"] != "native_dwg" or not verification["pass"]:
        raise RuntimeError("native DWG linework is not mechanically verified")
    variant = linework["variants"]["WFB"]
    if variant["source_dwg_sha256"] != EXPECTED["wfb_2d"]:
        raise RuntimeError("WFB DWG hash mismatch")
    center_x = (minimum[0] + maximum[0]) / 2.0
    center_y = (minimum[1] + maximum[1]) / 2.0
    base_z = minimum[2]
    if view == "plan":
        paths = [
            [(center_x + point[0], center_y + point[1]) for point in path]
            for path in variant["views"]["plan"]["paths_mm"]
        ]
    else:
        center = center_x if view == "front" else center_y
        paths = [
            [(center + point[0], base_z + point[1]) for point in path]
            for path in variant["views"]["elevation"]["paths_mm"]
        ]
    metadata = {
        **reference,
        "source_kind": "native_dwg",
        "linework_register": relative(linework_path),
        "linework_register_sha256": sha256(linework_path),
        "verification_register": relative(verification_path),
        "verification_register_sha256": sha256(verification_path),
        "verification_pass": True,
        "pdf_cross_check_sha256": EXPECTED["pdf"],
    }
    return paths, metadata


def svg_path(paths, transform, close: bool = False):
    commands = []
    for path in paths:
        if len(path) < 2:
            continue
        first = transform(path[0])
        commands.append(f"M {first[0]:.3f} {first[1]:.3f}")
        for point in path[1:]:
            value = transform(point)
            commands.append(f"L {value[0]:.3f} {value[1]:.3f}")
        if close:
            commands.append("Z")
    return " ".join(commands)


def render_svg(profile, view, raw_edges, silhouette, semantic, official, metadata):
    width, height = 1400, 980
    plot_x, plot_y, plot_w, plot_h = 60, 170, 980, 750
    points = []
    for start, end in raw_edges:
        points.extend((start, end))
    for group in (silhouette, semantic, official):
        for path in group:
            points.extend(path)
    min_x, max_x = min(point[0] for point in points), max(point[0] for point in points)
    min_y, max_y = min(point[1] for point in points), max(point[1] for point in points)
    padding = max(max_x - min_x, max_y - min_y) * 0.06
    min_x, max_x, min_y, max_y = min_x - padding, max_x + padding, min_y - padding, max_y + padding
    scale = min(plot_w / (max_x - min_x), plot_h / (max_y - min_y))

    def transform(point):
        return plot_x + (point[0] - min_x) * scale, plot_y + plot_h - (point[1] - min_y) * scale

    raw_path = svg_path([[start, end] for start, end in raw_edges], transform)
    silhouette_path = svg_path(silhouette, transform, close=True)
    semantic_path = svg_path(semantic, transform)
    official_path = svg_path(official, transform)
    required = "\n".join(
        f'  <text x="1092" y="{255 + index * 28}" font-family="Arial, sans-serif" font-size="15" fill="#34495e">• {html.escape(item)}</text>'
        for index, item in enumerate(profile["required_semantics"])
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="1400" height="980" fill="#fbfaf7"/>
  <text x="60" y="58" font-family="Arial, sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">{html.escape(profile["display_name"])}</text>
  <text x="60" y="96" font-family="Arial, sans-serif" font-size="20" fill="#41566d">{VIEWS[view]["label"]} · isolated representative {profile["representative_global_id"]}</text>
  <text x="60" y="128" font-family="Arial, sans-serif" font-size="17" fill="#68798a">Grey = original high-poly projection sample · Black = simplified proxy · Blue = official WFB native DWG</text>
  <rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#ffffff" stroke="#cad2d9" stroke-width="2"/>
  <path class="original-highpoly" d="{raw_path}" fill="none" stroke="#87929c" stroke-width="0.55" stroke-opacity="0.28" vector-effect="non-scaling-stroke"/>
  <path class="simplified-proxy-silhouette" d="{silhouette_path}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
  <path class="simplified-proxy-semantics" d="{semantic_path}" fill="none" stroke="#111820" stroke-width="1.6" stroke-linecap="round" vector-effect="non-scaling-stroke"/>
  <path class="official-reference-mask" d="{official_path}" fill="none" stroke="#ffffff" stroke-width="7" stroke-linecap="round" vector-effect="non-scaling-stroke"/>
  <path class="official-reference native-dwg" data-source-kind="native_dwg" data-model-code="WFB" d="{official_path}" fill="none" stroke="#1677c8" stroke-width="2.8" stroke-linecap="round" vector-effect="non-scaling-stroke"/>
  <text x="1080" y="205" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Acceptance</text>
{required}
  <text x="1080" y="510" font-family="Arial, sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Isolation proof</text>
  <text x="1080" y="542" font-family="Arial, sans-serif" font-size="15" fill="#41566d">geometry products: 1</text>
  <text x="1080" y="568" font-family="Arial, sans-serif" font-size="15" fill="#41566d">whole model render: false</text>
  <text x="1080" y="594" font-family="Arial, sans-serif" font-size="15" fill="#41566d">mesh faces: {metadata["mesh_face_count"]}</text>
  <text x="1080" y="620" font-family="Arial, sans-serif" font-size="15" fill="#41566d">raw edges: {metadata["raw_edge_count"]}</text>
  <text x="1080" y="646" font-family="Arial, sans-serif" font-size="15" fill="#41566d">displayed grey edges: {metadata["displayed_raw_edge_count"]}</text>
  <text x="1080" y="672" font-family="Arial, sans-serif" font-size="15" fill="#41566d">proxy paths: {metadata["proxy_path_count"]}</text>
  <text x="1080" y="700" font-family="Arial, sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Official reference</text>
  <text x="1080" y="732" font-family="Arial, sans-serif" font-size="15" fill="#1677c8">WFB native DWG · PDF cross-check passed</text>
  <text x="1080" y="780" font-family="Arial, sans-serif" font-size="14" fill="#68798a">Family reference only.</text>
  <text x="1080" y="804" font-family="Arial, sans-serif" font-size="14" fill="#68798a">Visual review pending; no IFC write.</text>
</svg>'''


def render_index(output_dir: Path, manifest: dict) -> None:
    cards = "\n".join(
        f'''<article><h2>{view["view"].title()}</h2><a href="{view["view"]}.svg"><img src="{view["view"]}.svg" alt="{view["view"]} review"></a></article>'''
        for view in manifest["views"]
    )
    content = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Falper Sorgente review</title>
<style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}nav a{{margin-right:20px}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:white;padding:14px;border-radius:10px}}img{{width:100%;height:auto}}code{{font-size:13px}}</style>
<h1>Falper Sorgente WFB native-DWG review</h1><p>Grey = original high-poly; black = proxy; blue = native WFB.dwg. Status: visual review pending.</p>
<nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate paths</a><a href="../../../../pipeline/decisions/falper-sorgente-dwg-pdf-verification.json">DWG/PDF verification</a></nav>
<main>{cards}</main><p><code>{html.escape(manifest["formal_ifc_sha256"])}</code></p></html>'''
    (output_dir / "index.html").write_text(content, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--register", type=Path, default=DEFAULT_REGISTER)
    parser.add_argument("--profile-key", default="falper-sorgente")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source = args.input.resolve()
    register_path = args.register.resolve()
    output_root = args.output_root.resolve()
    formal_before = sha256(source)
    register = load_json(register_path)
    if formal_before != register["formal_ifc_sha256"] or formal_before != EXPECTED["formal_ifc"]:
        raise RuntimeError("formal IFC hash mismatch")
    profile = register["profiles"][args.profile_key]
    model = ifcopenshell.open(source)
    representative, vertices, faces = mesh_for_one_product(model, profile["representative_global_id"])
    actual_type = product_type_name(representative)
    if actual_type != profile["ifc_type_name"]:
        raise RuntimeError(f"profile type mismatch: expected {profile['ifc_type_name']}, got {actual_type}")
    instances = sorted(
        product.GlobalId
        for product in model.by_type(representative.is_a())
        if product_type_name(product) == actual_type
    )
    if instances != sorted(profile["expected_instance_global_ids"]):
        raise RuntimeError("registered IFC type instance set drifted")
    minimum, maximum = bounds_3d(vertices)
    output_dir = output_root / args.profile_key
    output_dir.mkdir(parents=True, exist_ok=True)
    view_records = []
    candidate_views = {}
    resolved_reference = None
    for view, definition in VIEWS.items():
        axes = definition["axes"]
        all_raw_edges = projected_raw_edges(vertices, faces, axes)
        raw_edges = display_edge_sample(all_raw_edges)
        silhouette = projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"]))
        semantic = []
        sections = []
        for section in profile["semantic_sections"].get(view, []):
            paths, plane = semantic_section_paths(vertices, faces, axes, section, minimum, maximum)
            semantic.extend(paths)
            sections.append({**section, "plane_mm": round(plane, 6), "path_count": len(paths)})
        official, reference = native_wfb_paths(profile, view, minimum, maximum)
        displayed_official = display_path_sample(official)
        resolved_reference = reference
        metadata = {
            "mesh_face_count": len(faces),
            "raw_edge_count": len(all_raw_edges),
            "displayed_raw_edge_count": len(raw_edges),
            "proxy_path_count": len(silhouette) + len(semantic),
        }
        target = output_dir / f"{view}.svg"
        target.write_text(
            render_svg(profile, view, raw_edges, silhouette, semantic, displayed_official, metadata),
            encoding="utf-8",
        )
        candidate_views[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": [
                [[round(float(x), 6), round(float(y), 6)] for x, y in path]
                for path in silhouette + semantic
            ],
            "official_native_dwg_paths_mm": [
                [[round(float(x), 6), round(float(y), 6)] for x, y in path]
                for path in official
            ],
        }
        view_records.append(
            {
                "view": view,
                "projection_axes": list(axes),
                "svg": relative(target),
                "svg_sha256": sha256(target),
                "raw_edge_count": len(all_raw_edges),
                "displayed_raw_edge_count": len(raw_edges),
                "silhouette_path_count": len(silhouette),
                "semantic_path_count": len(semantic),
                "official_reference_path_count": len(official),
                "official_reference_source_point_count": sum(len(path) for path in official),
                "official_reference_display_point_count": sum(len(path) for path in displayed_official),
                "blue_line_source_kind": "native_dwg",
                "blue_line_model_code": "WFB",
                "semantic_sections": sections,
            }
        )
    candidate_path = output_dir / "candidate-representations.json"
    candidate = {
        "schema_version": 2,
        "profile_key": args.profile_key,
        "representative_global_id": representative.GlobalId,
        "ifc_type_name": actual_type,
        "units": "mm",
        "scope": SCOPE,
        "source_kind": "geometry_derived_simplified_proxy",
        "formal_ifc_write_allowed": False,
        "review_status": "visual_review_pending",
        "views": candidate_views,
    }
    write_json(candidate_path, candidate)
    manifest = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/int1_highpoly_type_review.py",
        "formal_ifc": relative(source),
        "formal_ifc_sha256": formal_before,
        "formal_ifc_sha256_after_generation": sha256(source),
        "formal_ifc_bytes_unchanged": sha256(source) == formal_before,
        "profile_register": relative(register_path),
        "profile_register_sha256": sha256(register_path),
        "profile_key": args.profile_key,
        "display_name": profile["display_name"],
        "ifc_type_name": actual_type,
        "representative_global_id": representative.GlobalId,
        "registered_instance_global_ids": instances,
        "geometry_product_count": 1,
        "whole_model_render": False,
        "formal_ifc_write": False,
        "mesh_vertex_count": len(vertices),
        "mesh_face_count": len(faces),
        "local_bounds_mm": {
            "minimum": [round(value, 6) for value in minimum],
            "maximum": [round(value, 6) for value in maximum],
        },
        "official_reference": resolved_reference,
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "views": view_records,
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "pass": len(view_records) == 3 and all(record["blue_line_source_kind"] == "native_dwg" for record in view_records),
    }
    manifest_path = output_dir / "manifest.json"
    write_json(manifest_path, manifest)
    render_index(output_dir, manifest)
    if sha256(source) != formal_before:
        raise RuntimeError("formal IFC bytes changed during review generation")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
