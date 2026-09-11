#!/usr/bin/env python3
"""Generate an isolated Duofix Sigma review with exact 224.212.00.2 DWGs."""

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


ARTICLE = "224.212.00.2"
PROFILE_KEY = "geberit-duofix-sigma-224-212"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OUTPUT_DIR = ROOT / "output/review/highpoly-types/geberit-duofix-sigma-224-212"
REGISTER = OUTPUT_DIR / "profile.json"
LINEWORK = OUTPUT_DIR / "official-native-dwg-linework.json"
SOURCE_REVALIDATION = OUTPUT_DIR / "official-source/official-source-revalidation.json"
CATALOGUE_EXTRACT = OUTPUT_DIR / "official-source/geberit-224212-catalogue-page11.pdf"
CATALOGUE_EXTRACT_SHA256 = "8c3a17569e6fa53d9a769ab584fed79a865118c4a27018facd619b4a5cc5cdd3"
BOX_LABEL = OUTPUT_DIR / "official-source/project-received-224212-box-label.jpg"
BOX_LABEL_SHA256 = "522952db10076af64be712e2843033dd69747428c15bf1c5a9085a82db840fa0"
BLUE = "#1677c8"
SOURCE_LABEL_ZH = "基于 Geberit 精确型号原生 DWG 的官方图纸表达"
SOURCE_LABEL_EN = "official drawing representation from exact Geberit article native DWG"


def path_bounds(paths):
    points = [point for path in paths for point in path]
    return (
        [min(point[axis] for point in points) for axis in range(2)],
        [max(point[axis] for point in points) for axis in range(2)],
    )


def align_paths(view, source_view, minimum, maximum):
    paths = source_view["paths_mm"]
    source_min, source_max = path_bounds(paths)
    source_center = [(source_min[axis] + source_max[axis]) / 2.0 for axis in range(2)]
    if view == "plan":
        target_center = [(minimum[0] + maximum[0]) / 2.0, (minimum[1] + maximum[1]) / 2.0]
        return [[(target_center[0] + x - source_center[0], target_center[1] + y - source_center[1]) for x, y in path] for path in paths]
    if view == "front":
        target_center_x = (minimum[0] + maximum[0]) / 2.0
        return [[(target_center_x + x - source_center[0], maximum[2] + y - source_max[1]) for x, y in path] for path in paths]
    target_center_y = (minimum[1] + maximum[1]) / 2.0
    return [[(target_center_y + x - source_center[0], maximum[2] + y - source_max[1]) for x, y in path] for path in paths]


def orient_ifc_projection(view, paths, minimum, maximum):
    """Match IFC projections to the manufacturer drawing-view direction."""
    if view != "side":
        return paths
    # Geberit's L drawing is a left-side elevation. The generic Y/Z IFC
    # projection is viewed from local +X, so reflect its horizontal axis to
    # express the same local -X viewing direction as the official L view.
    horizontal_axis_sum = minimum[1] + maximum[1]
    return [
        [(horizontal_axis_sum - x, y) for x, y in path]
        for path in paths
    ]


def side_handedness_check(proxy, official, minimum, maximum):
    alignment_center = (minimum[1] + maximum[1]) / 2.0

    def horizontal_centroid(paths):
        points = [point for path in paths for point in path]
        return sum(point[0] for point in points) / len(points)

    proxy_offset = horizontal_centroid(proxy) - alignment_center
    official_offset = horizontal_centroid(official) - alignment_center
    return {
        "comparison": "asymmetric_path_centroids_share_left_view_handedness",
        "view_direction": "left_side_local_negative_x",
        "proxy_horizontal_centroid_offset_mm": round(proxy_offset, 6),
        "official_horizontal_centroid_offset_mm": round(official_offset, 6),
        "same_horizontal_side": proxy_offset * official_offset > 0.0,
        "pass": proxy_offset * official_offset > 0.0,
    }


def cross_check(view, source_view, minimum, maximum):
    official = source_view["contour_bounds_mm"]["size"]
    if view == "plan":
        actual = [maximum[0] - minimum[0], maximum[1] - minimum[1]]
        delta = [abs(official[index] - actual[index]) for index in range(2)]
        return {
            "comparison": "native_dwg_full_plan_contour_vs_ifc_body_footprint",
            "official_size_mm": official,
            "ifc_body_size_mm": [round(value, 6) for value in actual],
            "absolute_delta_mm": [round(value, 6) for value in delta],
            "width_tolerance_mm": 1.0,
            "depth_tolerance_mm": 20.0,
            "note": "The exact official G contour includes approximately 14.8 mm more installation depth than the current IFC Body envelope; width identity remains within 1 mm.",
            "pass": delta[0] <= 1.0 and delta[1] <= 20.0,
        }
    if view == "front":
        actual = [maximum[0] - minimum[0], maximum[2] - minimum[2]]
        delta = [abs(official[index] - actual[index]) for index in range(2)]
        return {
            "comparison": "native_dwg_front_contour_vs_ifc_body_width_height",
            "official_size_mm": official,
            "ifc_body_size_mm": [round(value, 6) for value in actual],
            "absolute_delta_mm": [round(value, 6) for value in delta],
            "tolerance_mm": 1.0,
            "pass": max(delta) <= 1.0,
        }
    actual = [maximum[1] - minimum[1], maximum[2] - minimum[2]]
    delta = [abs(official[index] - actual[index]) for index in range(2)]
    return {
        "comparison": "native_dwg_side_contour_vs_ifc_body_depth_height",
        "official_size_mm": official,
        "ifc_body_size_mm": [round(value, 6) for value in actual],
        "absolute_delta_mm": [round(value, 6) for value in delta],
        "depth_tolerance_mm": 20.0,
        "height_tolerance_mm": 1.0,
        "note": "The exact official L contour includes the same installation-depth extent recorded in the G view; height identity remains within 1 mm.",
        "pass": delta[0] <= 20.0 and delta[1] <= 1.0,
    }


def rounded(paths):
    return [[[round(float(x), 6), round(float(y), 6)] for x, y in path] for path in paths]


def render_svg(profile, view, edges, proxy, official, source_view, check, face_count):
    width, height = 1400, 980
    plot_x, plot_y, plot_w, plot_h = 60, 170, 980, 750
    points = [point for edge in edges for point in edge]
    points.extend(point for path in proxy for point in path)
    points.extend(point for path in official for point in path)
    min_x, max_x = min(point[0] for point in points), max(point[0] for point in points)
    min_y, max_y = min(point[1] for point in points), max(point[1] for point in points)
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
<text x="60" y="128" font-family="Arial,sans-serif" font-size="17" fill="#68798a">Grey = actual IFC Body · Black = simplified proxy · Blue = exact official {ARTICLE} native DWG</text>
<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="original-highpoly" d="{edge_path}" fill="none" stroke="#87929c" stroke-width="0.5" stroke-opacity="0.25" vector-effect="non-scaling-stroke"/>
<path class="simplified-proxy-silhouette" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference-mask" d="{official_path}" fill="none" stroke="#fff" stroke-width="7" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference native-dwg" data-source-kind="native_dwg" data-article-number="{ARTICLE}" data-native-dwg-code="{source_view["native_dwg_code"]}" data-source-sha256="{source_view["source_dwg_sha256"]}" d="{official_path}" fill="none" stroke="{BLUE}" stroke-width="2.1" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="1080" y="205" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Acceptance</text>
{requirements}
<text x="1080" y="420" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Identity</text>
<text x="1080" y="452" font-family="Arial,sans-serif" font-size="15" fill="{BLUE}">{ARTICLE} · {source_view["native_dwg_code"]}.dwg</text>
<text x="1080" y="480" font-family="Arial,sans-serif" font-size="15" fill="#41566d">all official contour paths shown: {len(official)}</text>
<text x="1080" y="508" font-family="Arial,sans-serif" font-size="15" fill="#41566d">mechanical gate: {str(check["pass"]).lower()}</text>
<text x="1080" y="570" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Isolation</text>
<text x="1080" y="602" font-family="Arial,sans-serif" font-size="15" fill="#41566d">geometry products: 1</text>
<text x="1080" y="630" font-family="Arial,sans-serif" font-size="15" fill="#41566d">whole model render: false</text>
<text x="1080" y="658" font-family="Arial,sans-serif" font-size="15" fill="#41566d">mesh faces: {face_count}</text>
<text x="1080" y="686" font-family="Arial,sans-serif" font-size="15" fill="#41566d">official paths: {len(official)}</text>
<text x="1080" y="746" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Exact article family reference; not a shop drawing.</text>
<text x="1080" y="770" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Visual review pending; no IFC write.</text>
</svg>'''


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Geberit Duofix Sigma review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Geberit Duofix Sigma / {ARTICLE}</h1><p>Blue line = exact official G/A/L native DWG. Visual review pending.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="official-native-dwg-linework.json">Official linework</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/">Official CAD</a><a href="project-context-sanitary-plan.svg">Full project plan</a><a href="project-context-front-elevation.svg">Project front elevation</a><a href="project-context-side-elevation.svg">Project side elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a></nav><h2>DWG / IFC mechanical overlay</h2><main>{cards}</main><h2>Actual Bonsai camera renders</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
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
    if linework.get("article_number") != ARTICLE or linework.get("source_kind") != "native_dwg" or not linework.get("pass"):
        raise RuntimeError("Duofix exact-article native-DWG gate failed")
    if sha256(CATALOGUE_EXTRACT) != CATALOGUE_EXTRACT_SHA256:
        raise RuntimeError("Duofix official catalogue extract hash mismatch")
    if sha256(BOX_LABEL) != BOX_LABEL_SHA256:
        raise RuntimeError("Duofix project received-box label hash mismatch")
    for code, source_record in linework["official_sources"].items():
        source_path = ROOT / source_record["path"]
        if sha256(source_path) != source_record["sha256"]:
            raise RuntimeError(f"Duofix official {code} DWG hash mismatch")
    source_revalidation = load_json(SOURCE_REVALIDATION)
    if (
        source_revalidation.get("article_number") != ARTICLE
        or not source_revalidation.get("all_declared_dwg_urls_accessible")
        or not source_revalidation.get("all_downloaded_dwg_bytes_match_local_archive")
        or not source_revalidation.get("catalogue", {}).get("downloaded_article_page_render_matches_local_extract")
        or not source_revalidation.get("catalogue", {}).get("pass")
        or not source_revalidation.get("pass")
    ):
        raise RuntimeError("Duofix exact-article official-source revalidation gate failed")
    for code, source_record in linework["official_sources"].items():
        result = source_revalidation.get("dwg_results", {}).get(code, {})
        if (
            result.get("url") != source_record["url"]
            or result.get("expected_sha256") != source_record["sha256"]
            or result.get("downloaded_sha256") != source_record["sha256"]
            or result.get("local_sha256") != source_record["sha256"]
            or not result.get("downloaded_bytes_match_local_archive")
            or not result.get("pass")
        ):
            raise RuntimeError(f"Duofix exact-article official CDN evidence mismatch for {code}")
    model = ifcopenshell.open(source)
    product, vertices, faces = mesh_for_one_product(model, profile["representative_global_id"])
    if product_type_name(product) != profile["ifc_type_name"]:
        raise RuntimeError("Duofix representative type identity drifted")
    instances = sorted(item.GlobalId for item in model.by_type(product.is_a()) if product_type_name(item) == profile["ifc_type_name"])
    if instances != sorted(profile["expected_instance_global_ids"]):
        raise RuntimeError("Duofix instance set drifted")
    minimum, maximum = bounds_3d(vertices)
    output.mkdir(parents=True, exist_ok=True)
    records, candidates, checks = [], {}, {}
    for view, definition in VIEWS.items():
        axes = definition["axes"]
        all_edges = projected_raw_edges(vertices, faces, axes)
        edges = orient_ifc_projection(view, display_edge_sample(all_edges), minimum, maximum)
        proxy = orient_ifc_projection(
            view,
            projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"])),
            minimum,
            maximum,
        )
        source_view = linework["views"][view]
        official = align_paths(view, source_view, minimum, maximum)
        check = cross_check(view, source_view, minimum, maximum)
        if view == "side":
            check["handedness"] = side_handedness_check(proxy, official, minimum, maximum)
            check["pass"] = check["pass"] and check["handedness"]["pass"]
        if not check["pass"]:
            raise RuntimeError(f"Duofix {view} mechanical identity gate failed")
        checks[view] = check
        target = output / f"{view}.svg"
        target.write_text(render_svg(profile, view, edges, proxy, official, source_view, check, len(faces)), encoding="utf-8")
        candidates[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": rounded(proxy),
            "official_native_dwg_paths_mm": rounded(official),
            "native_dwg_code": source_view["native_dwg_code"],
            "native_dwg_sha256": source_view["source_dwg_sha256"],
            "view_direction": "left_side_local_negative_x" if view == "side" else definition["label"],
        }
        records.append({
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
            "view_direction": "left_side_local_negative_x" if view == "side" else definition["label"],
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
        "formal_ifc_write_allowed": False,
        "review_status": "visual_review_pending",
        "views": candidates,
    })
    source_access_record = output / "official-source/source-access-record.json"
    write_json(source_access_record, {
        "schema_version": 1,
        "checked_on": datetime.now(timezone.utc).date().isoformat(),
        "manufacturer": "Geberit",
        "family": "Duofix element for wall-hung WC with Sigma concealed cistern 12 cm",
        "article_number": ARTICLE,
        "ifc_type_name": profile["ifc_type_name"],
        "source_kind": "native_dwg",
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_acquired": True,
        "official_cad_exact_project_configuration_match": True,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "blue_line_present": True,
        "representation_geometry_source": "official_native_dwg_paths_mm",
        "official_product_catalogue": profile["official_reference"]["official_catalogue"],
        "official_catalogue_extract": {
            "path": relative(CATALOGUE_EXTRACT),
            "sha256": CATALOGUE_EXTRACT_SHA256,
        },
        "project_received_identity_evidence": {
            "path": relative(BOX_LABEL),
            "sha256": BOX_LABEL_SHA256,
        },
        "official_native_dwg": linework["official_sources"],
        "official_source_revalidation": {
            "path": relative(SOURCE_REVALIDATION),
            "sha256": sha256(SOURCE_REVALIDATION),
            "checked_at": source_revalidation["checked_at"],
            "all_declared_dwg_urls_accessible": source_revalidation["all_declared_dwg_urls_accessible"],
            "all_downloaded_dwg_bytes_match_local_archive": source_revalidation["all_downloaded_dwg_bytes_match_local_archive"],
            "catalogue_page_render_matches_local_extract": source_revalidation["catalogue"]["downloaded_article_page_render_matches_local_extract"],
            "pass": source_revalidation["pass"],
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
        "generator": "pipeline/scripts/geberit_duofix_sigma_224_212_review.py",
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
        "official_cad_exact_project_configuration_match": True,
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
            "linework_sha256": sha256(LINEWORK),
            "official_sources": linework["official_sources"],
            "official_source_revalidation": {
                "path": relative(SOURCE_REVALIDATION),
                "sha256": sha256(SOURCE_REVALIDATION),
                "checked_at": source_revalidation["checked_at"],
                "catalogue_downloaded_pdf_sha256": source_revalidation["catalogue"]["downloaded_pdf_sha256"],
                "catalogue_page_render_matches_local_extract": source_revalidation["catalogue"]["downloaded_article_page_render_matches_local_extract"],
                "pass": source_revalidation["pass"],
            },
            "mechanical_cross_checks": checks,
            "received_box_label": linework["received_box_label"],
            "nominal_dimensions_mm": linework["nominal_dimensions_mm"],
        },
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "official_source_access_record": relative(source_access_record),
        "official_source_access_record_sha256": sha256(source_access_record),
        "views": records,
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "pass": all(check["pass"] for check in checks.values()),
    }
    context_manifest = output / "project-context-manifest.json"
    bonsai_manifest = output / "bonsai-review-manifest.json"
    if context_manifest.is_file():
        context = load_json(context_manifest)
        if context.get("pass") is not True:
            raise RuntimeError("Duofix project-context evidence failed")
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
            raise RuntimeError("Duofix Bonsai camera evidence failed")
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
    print(json.dumps({"manifest": relative(manifest_path), "views": [item["svg"] for item in records], "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
