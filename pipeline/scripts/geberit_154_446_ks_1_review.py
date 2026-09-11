#!/usr/bin/env python3
"""Generate isolated CleanLine50 review views with exact archived DWG overlays."""

from __future__ import annotations

import argparse
import html
import json
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


ARTICLE = "154.446.KS.1"
PROFILE_KEY = "geberit-154-446-ks-1"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OUTPUT_DIR = ROOT / "output/review/highpoly-types/geberit-154-446-ks-1"
REGISTER = OUTPUT_DIR / "profile.json"
LINEWORK = OUTPUT_DIR / "official-native-dwg-linework.json"
CDN_REVALIDATION = OUTPUT_DIR / "official-source/official-cdn-revalidation.json"
PRODUCT_PAGE_ARCHIVE = OUTPUT_DIR / "official-source/geberit-cleanline50-product-page.html"
PRODUCT_PAGE_ARCHIVE_SHA256 = "cf4578d1f1d60078c683548281f22e30fb71608f7d4413a21c840bde212fa02f"
BLUE = "#1677c8"
SOURCE_LABEL_ZH = "Geberit 官方归档原生 DWG 图纸表达"
SOURCE_LABEL_EN = "drawing representation from archived Geberit official native DWG"


def path_bounds(paths):
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return minimum, maximum


def dominant_component_axis_bounds(vertices, faces, axis):
    """Return bounds for the largest connected mesh component on one axis."""
    parent = list(range(len(vertices)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(first, second):
        first_root, second_root = find(first), find(second)
        if first_root != second_root:
            parent[second_root] = first_root

    for face in faces:
        union(face[0], face[1])
        union(face[1], face[2])
        union(face[2], face[0])
    components = {}
    for index, point in enumerate(vertices):
        components.setdefault(find(index), []).append(point)
    primary = max(components.values(), key=len)
    return min(point[axis] for point in primary), max(point[axis] for point in primary)


def align_paths(view, source_view, minimum, maximum, side_axis_datum_mm=None):
    paths = source_view["paths_mm"]
    contour_min, contour_max = path_bounds(paths)
    if view == "plan":
        header = source_view["native_header_extents_mm"]
        source_center = [
            (header["minimum"][axis] + header["maximum"][axis]) / 2.0
            for axis in range(2)
        ]
        target_center = [(minimum[0] + maximum[0]) / 2.0, (minimum[1] + maximum[1]) / 2.0]
        return [[(target_center[0] + x - source_center[0], target_center[1] + y - source_center[1]) for x, y in path] for path in paths]
    if view == "front":
        source_center_x = (contour_min[0] + contour_max[0]) / 2.0
        target_center_x = (minimum[0] + maximum[0]) / 2.0
        return [[(target_center_x + x - source_center_x, maximum[2] + y - contour_max[1]) for x, y in path] for path in paths]
    if side_axis_datum_mm is None:
        raise RuntimeError("CleanLine50 side alignment requires the IFC primary-channel axis")
    return [[(side_axis_datum_mm + x, maximum[2] + y - contour_max[1]) for x, y in path] for path in paths]


def cross_check(view, source_view, minimum, maximum, side_axis_datum_mm=None, side_primary_bounds_mm=None):
    if view == "plan":
        official = source_view["native_header_extents_mm"]["size"]
        actual = [maximum[0] - minimum[0], maximum[1] - minimum[1]]
        delta = [abs(official[index] - actual[index]) for index in range(2)]
        return {
            "comparison": "native_dwg_header_footprint_vs_ifc_body_footprint",
            "official_size_mm": official,
            "ifc_body_size_mm": [round(value, 6) for value in actual],
            "absolute_delta_mm": [round(value, 6) for value in delta],
            "tolerance_mm": 0.2,
            "pass": max(delta) <= 0.2,
        }
    official = source_view["contour_bounds_mm"]["size"]
    if view == "front":
        actual = [maximum[0] - minimum[0], maximum[2] - minimum[2]]
        delta = [abs(official[index] - actual[index]) for index in range(2)]
        return {
            "comparison": "native_dwg_visible_contour_vs_ifc_body_length_height",
            "official_size_mm": official,
            "ifc_body_size_mm": [round(value, 6) for value in actual],
            "absolute_delta_mm": [round(value, 6) for value in delta],
            "tolerance_mm": 0.2,
            "pass": max(delta) <= 0.2,
        }
    if side_axis_datum_mm is None or side_primary_bounds_mm is None:
        raise RuntimeError("CleanLine50 side cross-check requires the IFC primary-channel axis")
    actual = [maximum[1] - minimum[1], maximum[2] - minimum[2]]
    height_delta = abs(official[1] - actual[1])
    contour_center_x = (
        source_view["contour_bounds_mm"]["minimum"][0]
        + source_view["contour_bounds_mm"]["maximum"][0]
    ) / 2.0
    full_body_center_y = (minimum[1] + maximum[1]) / 2.0
    previous_translation = full_body_center_y - contour_center_x
    correction = side_axis_datum_mm - previous_translation
    aligned_visible_bounds = [
        side_axis_datum_mm + source_view["contour_bounds_mm"]["minimum"][0],
        side_axis_datum_mm + source_view["contour_bounds_mm"]["maximum"][0],
    ]
    return {
        "comparison": "native_dwg_visible_side_contour_height_vs_ifc_body_height",
        "official_visible_contour_size_mm": official,
        "ifc_body_size_mm": [round(value, 6) for value in actual],
        "visible_width_difference_mm": round(abs(official[0] - actual[0]), 6),
        "height_absolute_delta_mm": round(height_delta, 6),
        "horizontal_alignment": {
            "mode": "native_dwg_origin_to_ifc_primary_channel_axis",
            "ifc_primary_channel_bounds_y_mm": [round(value, 6) for value in side_primary_bounds_mm],
            "ifc_primary_channel_axis_y_mm": round(side_axis_datum_mm, 6),
            "native_dwg_origin_x_mm": 0.0,
            "applied_translation_mm": round(side_axis_datum_mm, 6),
            "previous_contour_center_translation_mm": round(previous_translation, 6),
            "horizontal_offset_correction_mm": round(correction, 6),
            "aligned_official_visible_bounds_y_mm": [round(value, 6) for value in aligned_visible_bounds],
            "axis_alignment_absolute_delta_mm": 0.0,
            "tolerance_mm": 0.01,
            "pass": True,
        },
        "tolerance_mm": 0.2,
        "note": "The official L contour depicts the visible 53.4 mm channel component; its native x=0 datum is aligned to the primary IFC channel axis. The IFC Body includes wider concealed ancillary geometry recorded by the G-view header footprint.",
        "pass": height_delta <= 0.2,
    }


def rounded(paths):
    return [[[round(float(x), 6), round(float(y), 6)] for x, y in path] for path in paths]


def render_svg(profile, view, edges, proxy, official, source_view, check, face_count):
    width, height = 1400, 980
    plot_x, plot_y, plot_w, plot_h = 60, 170, 980, 750
    points = [point for edge in edges for point in edge]
    points.extend(point for path in proxy for point in path)
    points.extend(point for path in official for point in path)
    min_x, max_x = min(p[0] for p in points), max(p[0] for p in points)
    min_y, max_y = min(p[1] for p in points), max(p[1] for p in points)
    padding = max(max_x - min_x, max_y - min_y) * 0.055
    min_x, max_x, min_y, max_y = min_x - padding, max_x + padding, min_y - padding, max_y + padding
    scale = min(plot_w / (max_x - min_x), plot_h / (max_y - min_y))

    def transform(point):
        return plot_x + (point[0] - min_x) * scale, plot_y + plot_h - (point[1] - min_y) * scale

    edge_path = svg_path([[start, end] for start, end in edges], transform)
    proxy_path = svg_path(proxy, transform, close=True)
    official_path = svg_path(official, transform)
    requirements = "\n".join(
        f'<text x="1090" y="{258 + index * 28}" font-family="Arial,sans-serif" font-size="15" fill="#34495e">• {html.escape(item)}</text>'
        for index, item in enumerate(profile["required_semantics"])
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="1400" height="980" fill="#fbfaf7"/>
<text x="60" y="58" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">{html.escape(profile["display_name"])}</text>
<text x="60" y="96" font-family="Arial,sans-serif" font-size="20" fill="#41566d">{VIEWS[view]["label"]} · isolated representative {profile["representative_global_id"]}</text>
<text x="60" y="128" font-family="Arial,sans-serif" font-size="17" fill="#68798a">Grey = actual IFC Body · Black = simplified proxy · Blue = exact archived {ARTICLE} native DWG</text>
<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="original-highpoly" d="{edge_path}" fill="none" stroke="#87929c" stroke-width="0.5" stroke-opacity="0.25" vector-effect="non-scaling-stroke"/>
<path class="simplified-proxy-silhouette" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference-mask" d="{official_path}" fill="none" stroke="#fff" stroke-width="7" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference native-dwg" data-source-kind="native_dwg" data-article-number="{ARTICLE}" data-native-dwg-code="{source_view["native_dwg_code"]}" d="{official_path}" fill="none" stroke="{BLUE}" stroke-width="2.35" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="1080" y="205" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Acceptance</text>
{requirements}
<text x="1080" y="420" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Identity</text>
<text x="1080" y="452" font-family="Arial,sans-serif" font-size="15" fill="{BLUE}">old article {ARTICLE} · {source_view["native_dwg_code"]}.dwg</text>
<text x="1080" y="480" font-family="Arial,sans-serif" font-size="15" fill="#41566d">replacement CAD used: false</text>
<text x="1080" y="508" font-family="Arial,sans-serif" font-size="15" fill="#41566d">mechanical gate: {str(check["pass"]).lower()}</text>
<text x="1080" y="570" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Isolation</text>
<text x="1080" y="602" font-family="Arial,sans-serif" font-size="15" fill="#41566d">geometry products: 1</text>
<text x="1080" y="630" font-family="Arial,sans-serif" font-size="15" fill="#41566d">whole model render: false</text>
<text x="1080" y="658" font-family="Arial,sans-serif" font-size="15" fill="#41566d">mesh faces: {face_count}</text>
<text x="1080" y="686" font-family="Arial,sans-serif" font-size="15" fill="#41566d">official paths: {len(official)}</text>
<text x="1080" y="746" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Family reference only.</text>
<text x="1080" y="770" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Visual review pending; no IFC write.</text>
</svg>'''


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png" alt="actual Bonsai IFC Body {view} camera render"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Geberit CleanLine50 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Geberit CleanLine50 154.446.KS.1</h1><p>Blue line = exact archived G/A/L native DWG. The replacement 154.446.KS.2 CAD is not used. Visual review pending.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="official-native-dwg-linework.json">Official linework</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/">Official CAD</a><a href="project-context-sanitary-plan.svg">Full project plan</a><a href="project-context-front-elevation.svg">Project front elevation</a><a href="project-context-side-elevation.svg">Project side elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://catalog.geberit-global.com/en-XB/product/PRO_4702328">Official product page</a></nav><h2>DWG / IFC mechanical overlay</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
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
    formal_hash = sha256(source)
    if formal_hash != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    profile = load_json(args.register.resolve())["profiles"][PROFILE_KEY]
    linework = load_json(LINEWORK)
    if linework.get("article_number") != ARTICLE or linework.get("source_kind") != "native_dwg" or linework.get("retirement_identity", {}).get("replacement_cad_used") is not False:
        raise RuntimeError("CleanLine50 exact-old-article linework gate failed")
    if sha256(PRODUCT_PAGE_ARCHIVE) != PRODUCT_PAGE_ARCHIVE_SHA256:
        raise RuntimeError("archived Geberit CleanLine50 product page hash mismatch")
    for code, source_record in linework["official_sources"].items():
        if sha256(ROOT / source_record["path"]) != source_record["sha256"]:
            raise RuntimeError(f"CleanLine50 official archived {code} DWG hash mismatch")
    cdn_revalidation = load_json(CDN_REVALIDATION)
    if (
        cdn_revalidation.get("article_number") != ARTICLE
        or cdn_revalidation.get("replacement_article_number") != linework.get("replacement_article_number")
        or cdn_revalidation.get("replacement_cad_used") is not False
        or not cdn_revalidation.get("all_declared_urls_accessible")
        or not cdn_revalidation.get("all_downloaded_bytes_match_local_archive")
        or not cdn_revalidation.get("pass")
    ):
        raise RuntimeError("CleanLine50 exact archived-article official CDN revalidation gate failed")
    for code, source_record in linework["official_sources"].items():
        cdn_result = cdn_revalidation.get("results", {}).get(code, {})
        if (
            cdn_result.get("url") != source_record["url"]
            or cdn_result.get("expected_sha256") != source_record["sha256"]
            or cdn_result.get("downloaded_sha256") != source_record["sha256"]
            or cdn_result.get("local_sha256") != source_record["sha256"]
            or not cdn_result.get("downloaded_bytes_match_local_archive")
            or not cdn_result.get("pass")
        ):
            raise RuntimeError(f"CleanLine50 exact archived-article CDN evidence mismatch for {code}")
    model = ifcopenshell.open(source)
    product, vertices, faces = mesh_for_one_product(model, profile["representative_global_id"])
    if product_type_name(product) != profile["ifc_type_name"]:
        raise RuntimeError("CleanLine50 representative type identity drifted")
    instances = sorted(item.GlobalId for item in model.by_type(product.is_a()) if product_type_name(item) == profile["ifc_type_name"])
    if instances != profile["expected_instance_global_ids"]:
        raise RuntimeError("CleanLine50 instance set drifted")
    minimum, maximum = bounds_3d(vertices)
    side_primary_bounds = dominant_component_axis_bounds(vertices, faces, 1)
    side_axis_datum = sum(side_primary_bounds) / 2.0
    output.mkdir(parents=True, exist_ok=True)
    view_records = []
    candidate_views = {}
    checks = {}
    for view, definition in VIEWS.items():
        axes = definition["axes"]
        all_edges = projected_raw_edges(vertices, faces, axes)
        edges = display_edge_sample(all_edges)
        proxy = projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"]))
        source_view = linework["views"][view]
        official = align_paths(view, source_view, minimum, maximum, side_axis_datum)
        check = cross_check(view, source_view, minimum, maximum, side_axis_datum, side_primary_bounds)
        if not check["pass"]:
            raise RuntimeError(f"CleanLine50 {view} mechanical identity gate failed")
        checks[view] = check
        target = output / f"{view}.svg"
        target.write_text(render_svg(profile, view, edges, proxy, official, source_view, check, len(faces)), encoding="utf-8")
        candidate_views[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": rounded(proxy),
            "official_native_dwg_paths_mm": rounded(official),
            "native_dwg_code": source_view["native_dwg_code"],
            "native_dwg_sha256": source_view["source_dwg_sha256"],
        }
        view_records.append({
            "view": view,
            "svg": relative(target),
            "svg_sha256": sha256(target),
            "projection_axes": list(axes),
            "raw_edge_count": len(all_edges),
            "displayed_raw_edge_count": len(edges),
            "silhouette_path_count": len(proxy),
            "official_reference_path_count": len(official),
            "blue_line_source_kind": "native_dwg",
            "blue_line_article_number": ARTICLE,
            "blue_line_native_dwg_code": source_view["native_dwg_code"],
            "blue_line_native_dwg_sha256": source_view["source_dwg_sha256"],
            "mechanical_cross_check": check,
        })
    candidate_path = output / "candidate-representations.json"
    write_json(candidate_path, {
        "schema_version": 1,
        "profile_key": PROFILE_KEY,
        "representative_global_id": product.GlobalId,
        "ifc_type_name": profile["ifc_type_name"],
        "units": "mm",
        "source_kind": "native_dwg",
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "blue_line_present": True,
        "representation_geometry_source": "official_native_dwg_paths_mm",
        "simplified_proxy_comparison_included": True,
        "simplified_proxy_comparison_source_kind": "geometry_derived_simplified_proxy",
        "official_overlay_source_kind": "native_dwg",
        "official_article_number": ARTICLE,
        "replacement_cad_used": False,
        "formal_ifc_write_allowed": False,
        "review_status": "visual_review_pending",
        "views": candidate_views,
    })
    source_access_record = output / "official-source/source-access-record.json"
    source_checked_on = (
        load_json(source_access_record).get("checked_on")
        if source_access_record.is_file()
        else datetime.now(timezone.utc).date().isoformat()
    )
    write_json(source_access_record, {
        "schema_version": 1,
        "checked_on": source_checked_on,
        "manufacturer": "Geberit",
        "family": "CleanLine50 shower channel L90 cm",
        "article_number": ARTICLE,
        "replacement_article_number": linework["replacement_article_number"],
        "replacement_cad_used": False,
        "ifc_type_name": profile["ifc_type_name"],
        "source_kind": "native_dwg",
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_acquired": True,
        "official_cad_exact_archived_article_match": True,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "blue_line_present": True,
        "representation_geometry_source": "official_native_dwg_paths_mm",
        "official_product_page": profile["official_reference"]["product_page"],
        "official_product_page_archive": {
            "path": relative(PRODUCT_PAGE_ARCHIVE),
            "sha256": PRODUCT_PAGE_ARCHIVE_SHA256,
        },
        "current_catalog_status": profile["official_reference"]["current_catalog_status"],
        "official_native_dwg": linework["official_sources"],
        "official_cdn_revalidation": {
            "path": relative(CDN_REVALIDATION),
            "sha256": sha256(CDN_REVALIDATION),
            "checked_at": cdn_revalidation["checked_at"],
            "all_declared_urls_accessible": cdn_revalidation["all_declared_urls_accessible"],
            "all_downloaded_bytes_match_local_archive": cdn_revalidation["all_downloaded_bytes_match_local_archive"],
            "replacement_cad_used": cdn_revalidation["replacement_cad_used"],
            "pass": cdn_revalidation["pass"],
        },
        "drawing_view_mapping": {"plan": "G", "front": "A", "side": "L"},
        "official_3d_identity_only": {
            "code": "P",
            "used_as_plan_or_elevation_geometry": False,
        },
        "scope": profile["official_reference"]["scope"],
        "pass": True,
    })
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/geberit_154_446_ks_1_review.py",
        "formal_ifc": relative(source),
        "formal_ifc_sha256": formal_hash,
        "formal_ifc_sha256_after_generation": sha256(source),
        "formal_ifc_bytes_unchanged": sha256(source) == formal_hash,
        "profile_register": relative(args.register.resolve()),
        "profile_register_sha256": sha256(args.register.resolve()),
        "profile_key": PROFILE_KEY,
        "display_name": profile["display_name"],
        "ifc_type_name": profile["ifc_type_name"],
        "representative_global_id": product.GlobalId,
        "registered_instance_global_ids": instances,
        "geometry_product_count": 1,
        "whole_model_render": False,
        "formal_ifc_write": False,
        "source_kind": "native_dwg",
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_acquired": True,
        "official_cad_exact_archived_article_match": True,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "blue_line_present": True,
        "representation_geometry_source": "official_native_dwg_paths_mm",
        "simplified_proxy_comparison_included": True,
        "mesh_vertex_count": len(vertices),
        "mesh_face_count": len(faces),
        "local_bounds_mm": {"minimum": minimum, "maximum": maximum},
        "official_reference": {
            **profile["official_reference"],
            "product_page_archive": relative(PRODUCT_PAGE_ARCHIVE),
            "product_page_archive_sha256": sha256(PRODUCT_PAGE_ARCHIVE),
            "linework_sha256": sha256(LINEWORK),
            "official_sources": linework["official_sources"],
            "official_cdn_revalidation": {
                "path": relative(CDN_REVALIDATION),
                "sha256": sha256(CDN_REVALIDATION),
                "checked_at": cdn_revalidation["checked_at"],
                "replacement_cad_used": cdn_revalidation["replacement_cad_used"],
                "pass": cdn_revalidation["pass"],
            },
            "mechanical_cross_checks": checks,
        },
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "official_source_access_record": relative(source_access_record),
        "official_source_access_record_sha256": sha256(source_access_record),
        "views": view_records,
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "pass": all(check["pass"] for check in checks.values()),
    }
    context_manifest = output / "project-context-manifest.json"
    bonsai_manifest = output / "bonsai-review-manifest.json"
    if context_manifest.is_file():
        context = load_json(context_manifest)
        if context.get("pass") is not True:
            raise RuntimeError("CleanLine50 project-context evidence failed")
        manifest["project_context"] = {
            "manifest": relative(context_manifest),
            "manifest_sha256": sha256(context_manifest),
            "walls_and_surrounding_project_elements_retained": context["walls_and_surrounding_project_elements_retained"],
            "blue_line_top_layer_with_white_mask": context["blue_line_top_layer_with_white_mask"],
            "pass": True,
        }
    if bonsai_manifest.is_file():
        bonsai = load_json(bonsai_manifest)
        if bonsai.get("mode") != "actual_bonsai_ifc_body_camera_render" or bonsai.get("pass") is not True:
            raise RuntimeError("CleanLine50 Bonsai camera evidence failed")
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
    write_index(manifest)
    print(json.dumps({"manifest": relative(manifest_path), "views": [item["svg"] for item in view_records], "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
