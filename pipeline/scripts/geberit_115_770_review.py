#!/usr/bin/env python3
"""Generate isolated Geberit Sigma01 three-view review with exact official DWG overlays."""

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
    semantic_section_paths,
    svg_path,
)


OUTPUT_DIR = ROOT / "output/review/highpoly-types/geberit-115-770"
REGISTER = OUTPUT_DIR / "profile.json"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
LINEWORK = OUTPUT_DIR / "official-native-dwg-linework.json"
PRODUCT_PAGE_ARCHIVE = OUTPUT_DIR / "official-source/geberit-sigma01-product-page.html"
PRODUCT_PAGE_ARCHIVE_SHA256 = "46ed8825fcc0dae55300ecd458bd4969163f026d68e0ed81ba37c39c3180a0d3"
CDN_REVALIDATION = OUTPUT_DIR / "official-source/official-cdn-revalidation.json"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
PROFILE_KEY = "geberit-115-770"
ARTICLE = "115.770.11.5"
BLUE = "#1677c8"
SOURCE_LABEL_ZH = "基于 Geberit 精确型号原生 DWG 的官方图纸表达"
SOURCE_LABEL_EN = "official drawing representation from exact Geberit article native DWG"


def round_paths(paths):
    return [
        [[round(float(x), 6), round(float(y), 6)] for x, y in path]
        for path in paths
    ]


def align_official_paths(view: str, paths, minimum, maximum):
    """Place official G/A/L native-DWG coordinates in the IFC Body local frame."""
    if view == "plan":
        offset_x = (minimum[0] + maximum[0]) / 2.0
        wall_y = maximum[1]
        return [[(offset_x + x, wall_y + y) for x, y in path] for path in paths]
    if view == "front":
        offset_x = (minimum[0] + maximum[0]) / 2.0
        offset_z = (minimum[2] + maximum[2]) / 2.0
        return [[(offset_x + x, offset_z + y) for x, y in path] for path in paths]
    wall_y = maximum[1]
    offset_z = (minimum[2] + maximum[2]) / 2.0
    return [[(wall_y - x, offset_z + y) for x, y in path] for path in paths]


def bounds(paths):
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return minimum, maximum


def dimension_cross_check(view: str, official, minimum, maximum):
    axes = VIEWS[view]["axes"]
    official_min, official_max = bounds(official)
    official_size = [official_max[i] - official_min[i] for i in range(2)]
    ifc_size = [maximum[axis] - minimum[axis] for axis in axes]
    delta = [abs(official_size[i] - ifc_size[i]) for i in range(2)]
    # The project Body is a legacy authoring approximation (254 x 170.18 x
    # 12.45 mm); the exact current manufacturer article is 245.12 x 164.23 x
    # 12.01 mm in native DWG.  Keep both geometries unscaled and gate the
    # documented modelling delta instead of stretching official linework.
    tolerance = 10.0
    return {
        "official_native_dwg_size_mm": [round(value, 6) for value in official_size],
        "ifc_body_size_mm": [round(value, 6) for value in ifc_size],
        "absolute_delta_mm": [round(value, 6) for value in delta],
        "tolerance_mm": tolerance,
        "pass": max(delta) <= tolerance,
    }


def render_svg(profile, view, raw_edges, silhouette, semantic, official, metadata):
    width, height = 1400, 980
    plot_x, plot_y, plot_w, plot_h = 60, 170, 980, 750
    points = [point for edge in raw_edges for point in edge]
    for group in (silhouette, semantic, official):
        points.extend(point for path in group for point in path)
    min_x, max_x = min(p[0] for p in points), max(p[0] for p in points)
    min_y, max_y = min(p[1] for p in points), max(p[1] for p in points)
    padding = max(max_x - min_x, max_y - min_y) * 0.06
    min_x, max_x = min_x - padding, max_x + padding
    min_y, max_y = min_y - padding, max_y + padding
    scale = min(plot_w / (max_x - min_x), plot_h / (max_y - min_y))

    def transform(point):
        return plot_x + (point[0] - min_x) * scale, plot_y + plot_h - (point[1] - min_y) * scale

    raw_path = svg_path([[start, end] for start, end in raw_edges], transform)
    proxy_path = svg_path(silhouette, transform, close=True)
    semantic_path = svg_path(semantic, transform)
    official_path = svg_path(official, transform)
    requirements = "\n".join(
        f'  <text x="1092" y="{255 + index * 28}" font-family="Arial, sans-serif" font-size="15" fill="#34495e">• {html.escape(item)}</text>'
        for index, item in enumerate(profile["required_semantics"])
    )
    check = metadata["dimension_cross_check"]
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <rect width="1400" height="980" fill="#fbfaf7"/>
  <text x="60" y="58" font-family="Arial, sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">{html.escape(profile["display_name"])}</text>
  <text x="60" y="96" font-family="Arial, sans-serif" font-size="20" fill="#41566d">{VIEWS[view]["label"]} · isolated representative {profile["representative_global_id"]}</text>
  <text x="60" y="128" font-family="Arial, sans-serif" font-size="17" fill="#68798a">Grey = original IFC Body · Black = simplified proxy · Blue = official {ARTICLE} native DWG</text>
  <rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#ffffff" stroke="#cad2d9" stroke-width="2"/>
  <path class="original-highpoly" d="{raw_path}" fill="none" stroke="#87929c" stroke-width="0.55" stroke-opacity="0.28" vector-effect="non-scaling-stroke"/>
  <path class="simplified-proxy-silhouette" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
  <path class="simplified-proxy-semantics" d="{semantic_path}" fill="none" stroke="#111820" stroke-width="1.6" stroke-linecap="round" vector-effect="non-scaling-stroke"/>
  <path class="official-reference-mask" d="{official_path}" fill="none" stroke="#ffffff" stroke-width="7" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
  <path class="official-reference native-dwg" data-source-kind="native_dwg" data-article-number="{ARTICLE}" data-native-dwg-code="{metadata["native_dwg_code"]}" d="{official_path}" fill="none" stroke="{BLUE}" stroke-width="2.35" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
  <text x="1080" y="205" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Acceptance</text>
{requirements}
  <text x="1080" y="410" font-family="Arial, sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Identity check</text>
  <text x="1080" y="442" font-family="Arial, sans-serif" font-size="15" fill="{BLUE}">exact article {ARTICLE} · {metadata["native_dwg_code"]}.dwg</text>
  <text x="1080" y="470" font-family="Arial, sans-serif" font-size="15" fill="#41566d">DWG/IFC size delta: {max(check["absolute_delta_mm"]):.3f} mm</text>
  <text x="1080" y="498" font-family="Arial, sans-serif" font-size="15" fill="#41566d">mechanical gate: {str(check["pass"]).lower()}</text>
  <text x="1080" y="550" font-family="Arial, sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Isolation proof</text>
  <text x="1080" y="582" font-family="Arial, sans-serif" font-size="15" fill="#41566d">geometry products: 1</text>
  <text x="1080" y="610" font-family="Arial, sans-serif" font-size="15" fill="#41566d">whole model render: false</text>
  <text x="1080" y="638" font-family="Arial, sans-serif" font-size="15" fill="#41566d">mesh faces: {metadata["mesh_face_count"]}</text>
  <text x="1080" y="666" font-family="Arial, sans-serif" font-size="15" fill="#41566d">official paths: {metadata["official_path_count"]}</text>
  <text x="1080" y="720" font-family="Arial, sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Review gate</text>
  <text x="1080" y="752" font-family="Arial, sans-serif" font-size="14" fill="#68798a">Official family reference only.</text>
  <text x="1080" y="776" font-family="Arial, sans-serif" font-size="14" fill="#68798a">Visual review pending; no IFC write.</text>
</svg>'''


def render_index(manifest):
    cards = "\n".join(
        f'<article><h2>{view["view"].title()}</h2><a href="{view["view"]}.svg"><img src="{view["view"]}.svg" alt="{view["view"]} review"></a></article>'
        for view in manifest["views"]
    )
    bonsai_cards = "\n".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png" alt="actual Bonsai IFC Body {view} camera render"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    context_cards = "\n".join(
        f'<article><h2>{label}</h2><a href="{svg}"><img src="{preview}" alt="{label} project context review"></a></article>'
        for label, svg, preview in (
            ("Plan · instance 2gFgc", "project-context-sanitary-plan-2gFgcOYEXEaQWAzcKulTFt-review.svg", "project-context-sanitary-plan-2gFgcOYEXEaQWAzcKulTFt-review-preview.png"),
            ("Plan · instance 2lDPs", "project-context-sanitary-plan-2lDPsdQevFSfeOThtjSlPG-review.svg", "project-context-sanitary-plan-2lDPsdQevFSfeOThtjSlPG-review-preview.png"),
            ("Front elevation", "project-context-front-elevation-review.svg", "project-context-front-elevation-review-preview.png"),
            ("Side elevation", "project-context-side-elevation-review.svg", "project-context-side-elevation-review-preview.png"),
        )
    )
    content = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Geberit Sigma01 115.770.11.5 review</title>
<style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}nav a{{margin-right:20px}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:white;padding:14px;border-radius:10px}}img{{width:100%;height:auto}}code{{font-size:13px}}</style>
<h1>Geberit Sigma01 115.770.11.5 native-DWG review</h1><p>Grey = actual IFC Body; black = proxy; blue = exact official 115.770.11.5 G/A/L native DWG. The formal IFC's WCSEAT classification is documented as incorrect; this product is a flush actuator plate. Status: visual review pending; derived IFC write is gated.</p>
<nav><a href="review-contact-sheet.png">Review contact sheet</a><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate paths</a><a href="profile.json">Profile</a><a href="official-native-dwg-linework.json">Official linework</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/">Official source folder</a><a href="project-context-sanitary-plan.svg">Full project plan</a><a href="project-context-front-elevation.svg">Project front elevation</a><a href="project-context-side-elevation.svg">Project side elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://catalog.kz.geberit.com/ru-KZ/product/PRO_100555">Official product page</a></nav>
<h2>Project drawing context</h2><main>{context_cards}</main><h2>DWG / IFC mechanical overlay</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><p><code>{manifest["formal_ifc_sha256"]}</code></p></html>'''
    (OUTPUT_DIR / "index.html").write_text(content, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--register", type=Path, default=REGISTER)
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    formal_before = sha256(source)
    if formal_before != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    register = load_json(args.register.resolve())
    profile = register["profiles"][PROFILE_KEY]
    linework = load_json(LINEWORK)
    if linework.get("source_kind") != "native_dwg" or linework.get("article_number") != ARTICLE or not linework.get("pass"):
        raise RuntimeError("official Geberit native-DWG identity gate failed")
    if sha256(PRODUCT_PAGE_ARCHIVE) != PRODUCT_PAGE_ARCHIVE_SHA256:
        raise RuntimeError("archived Geberit official product page hash mismatch")
    for code, source_record in linework["official_sources"].items():
        if sha256(ROOT / source_record["path"]) != source_record["sha256"]:
            raise RuntimeError(f"Geberit Sigma01 official {code} DWG hash mismatch")
    cdn_revalidation = load_json(CDN_REVALIDATION)
    if (
        cdn_revalidation.get("article_number") != ARTICLE
        or not cdn_revalidation.get("all_declared_urls_accessible")
        or not cdn_revalidation.get("all_downloaded_bytes_match_local_archive")
        or not cdn_revalidation.get("pass")
    ):
        raise RuntimeError("Geberit 115.770 official CDN revalidation gate failed")
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
            raise RuntimeError(f"Geberit 115.770 official CDN evidence mismatch for {code}")
    model = ifcopenshell.open(source)
    product, vertices, faces = mesh_for_one_product(model, profile["representative_global_id"])
    if product_type_name(product) != profile["ifc_type_name"]:
        raise RuntimeError("Geberit representative type identity drifted")
    instances = sorted(
        item.GlobalId for item in model.by_type(product.is_a())
        if product_type_name(item) == profile["ifc_type_name"]
    )
    if instances != sorted(profile["expected_instance_global_ids"]):
        raise RuntimeError("Geberit instance set drifted")
    minimum, maximum = bounds_3d(vertices)
    output.mkdir(parents=True, exist_ok=True)
    views = []
    candidate_views = {}
    cross_checks = {}
    for view, definition in VIEWS.items():
        axes = definition["axes"]
        all_edges = projected_raw_edges(vertices, faces, axes)
        display_edges = display_edge_sample(all_edges)
        silhouette = projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"]))
        semantic = []
        sections = []
        for section in profile["semantic_sections"].get(view, []):
            paths, plane = semantic_section_paths(vertices, faces, axes, section, minimum, maximum)
            semantic.extend(paths)
            sections.append({**section, "plane_mm": round(plane, 6), "path_count": len(paths)})
        source_view = linework["views"][view]
        official = align_official_paths(view, source_view["paths_mm"], minimum, maximum)
        check = dimension_cross_check(view, official, minimum, maximum)
        if not check["pass"]:
            raise RuntimeError(f"official native DWG vs IFC dimension gate failed for {view}")
        cross_checks[view] = check
        metadata = {
            "mesh_face_count": len(faces),
            "official_path_count": len(official),
            "native_dwg_code": source_view["native_dwg_code"],
            "dimension_cross_check": check,
        }
        target = output / f"{view}.svg"
        target.write_text(
            render_svg(profile, view, display_edges, silhouette, semantic, official, metadata),
            encoding="utf-8",
        )
        candidate_views[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": round_paths(silhouette + semantic),
            "official_native_dwg_paths_mm": round_paths(official),
            "native_dwg_code": source_view["native_dwg_code"],
            "native_dwg_sha256": source_view["source_dwg_sha256"],
        }
        views.append({
            "view": view,
            "svg": relative(target),
            "svg_sha256": sha256(target),
            "projection_axes": list(axes),
            "raw_edge_count": len(all_edges),
            "displayed_raw_edge_count": len(display_edges),
            "silhouette_path_count": len(silhouette),
            "semantic_path_count": len(semantic),
            "semantic_sections": sections,
            "official_reference_path_count": len(official),
            "official_reference_source_point_count": sum(len(path) for path in official),
            "blue_line_source_kind": "native_dwg",
            "blue_line_article_number": ARTICLE,
            "blue_line_native_dwg_code": source_view["native_dwg_code"],
            "blue_line_native_dwg_sha256": source_view["source_dwg_sha256"],
            "dimension_cross_check": check,
        })
    candidate_path = output / "candidate-representations.json"
    candidate = {
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
        "views": candidate_views,
    }
    write_json(candidate_path, candidate)
    source_access_record = output / "official-source/source-access-record.json"
    write_json(source_access_record, {
        "schema_version": 1,
        "checked_on": datetime.now(timezone.utc).date().isoformat(),
        "manufacturer": "Geberit",
        "family": "Sigma01 dual-flush actuator plate",
        "article_number": ARTICLE,
        "ifc_type_name": profile["ifc_type_name"],
        "source_kind": "native_dwg",
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_acquired": True,
        "official_cad_exact_article_match": True,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "blue_line_present": True,
        "representation_geometry_source": "official_native_dwg_paths_mm",
        "official_product_page": profile["official_reference"]["product_page"],
        "official_product_page_archive": {
            "path": relative(PRODUCT_PAGE_ARCHIVE),
            "sha256": PRODUCT_PAGE_ARCHIVE_SHA256,
        },
        "official_pdf": {
            "path": profile["official_reference"]["official_pdf"],
            "sha256": profile["official_reference"]["official_pdf_sha256"],
        },
        "official_native_dwg": linework["official_sources"],
        "official_cdn_revalidation": {
            "path": relative(CDN_REVALIDATION),
            "sha256": sha256(CDN_REVALIDATION),
            "checked_at": cdn_revalidation["checked_at"],
            "all_declared_urls_accessible": cdn_revalidation["all_declared_urls_accessible"],
            "all_downloaded_bytes_match_local_archive": cdn_revalidation["all_downloaded_bytes_match_local_archive"],
            "pass": cdn_revalidation["pass"],
        },
        "drawing_view_mapping": {"plan": "G", "front": "A", "side": "L"},
        "official_3d_identity_only": {
            "code": "P",
            "used_as_plan_or_elevation_geometry": False,
        },
        "scope": profile["official_reference"]["scope"],
        "formal_ifc_semantic_issue": profile["formal_ifc_semantic_issue"],
        "pass": True,
    })
    manifest = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/geberit_115_770_review.py",
        "formal_ifc": relative(source),
        "formal_ifc_sha256": formal_before,
        "formal_ifc_sha256_after_generation": sha256(source),
        "formal_ifc_bytes_unchanged": sha256(source) == formal_before,
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
        "official_cad_exact_article_match": True,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "blue_line_present": True,
        "representation_geometry_source": "official_native_dwg_paths_mm",
        "simplified_proxy_comparison_included": True,
        "mesh_vertex_count": len(vertices),
        "mesh_face_count": len(faces),
        "local_bounds_mm": {"minimum": minimum, "maximum": maximum},
        "required_semantics": profile["required_semantics"],
        "formal_ifc_semantic_issue": profile["formal_ifc_semantic_issue"],
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
                "pass": cdn_revalidation["pass"],
            },
            "dimension_cross_checks": cross_checks,
        },
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "official_source_access_record": relative(source_access_record),
        "official_source_access_record_sha256": sha256(source_access_record),
        "project_context": relative(output / "project-context-manifest.json") if (output / "project-context-manifest.json").is_file() else None,
        "bonsai_review": relative(output / "bonsai-review-manifest.json") if (output / "bonsai-review-manifest.json").is_file() else None,
        "views": views,
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "pass": all(check["pass"] for check in cross_checks.values()),
    }
    manifest_path = output / "manifest.json"
    write_json(manifest_path, manifest)
    render_index(manifest)
    print(json.dumps({"manifest": relative(manifest_path), "views": [item["svg"] for item in views], "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
