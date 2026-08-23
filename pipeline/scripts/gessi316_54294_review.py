#!/usr/bin/env python3
"""Generate Gessi316 54294 review views from the exact official native DWG."""

from __future__ import annotations

import argparse
import html
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
OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = True
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


def align_official_paths(view, paths, minimum, maximum):
    center_x = (minimum[0] + maximum[0]) / 2.0
    if view == "plan":
        return [
            [[center_x + x, maximum[1] + y] for x, y in path]
            for path in paths
        ]
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
        ]
    return [
        [[maximum[1] + y, z_offset + z] for y, z in path]
        for path in paths
    ]


def compatibility(view, official, minimum, maximum):
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
<text x="60" y="128" font-family="Arial,sans-serif" font-size="17" fill="#68798a">Grey = actual IFC Body · Black = simplified proxy · Blue = exact official Gessi {DRAWING_PRODUCT_CODE} native DWG</text>
<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="original-highpoly" d="{edge_path}" fill="none" stroke="#87929c" stroke-width="0.5" stroke-opacity="0.30" vector-effect="non-scaling-stroke"/>
<path class="simplified-proxy-silhouette geometry-derived" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference-mask" d="{official_path}" fill="none" stroke="#ffffff" stroke-width="7" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference native-dwg" data-source-kind="native_dwg" data-product-code="{DRAWING_PRODUCT_CODE}" data-dwg-sha256="{SOURCE_DWG_SHA256}" d="{official_path}" fill="none" stroke="{BLUE}" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="1080" y="205" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Acceptance</text>
{requirements}
<text x="1080" y="410" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Official source</text>
<text x="1080" y="442" font-family="Arial,sans-serif" font-size="14" fill="#1677c8">Gessi {DRAWING_PRODUCT_CODE} native DWG</text>
<text x="1080" y="470" font-family="Arial,sans-serif" font-size="14" fill="#41566d">DWG paths: {len(official)}</text>
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
        f'''<!doctype html><html><meta charset="utf-8"><title>Gessi316 54294 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Gessi316 Meccanica / 45089_54294</h1><p>Grey = actual IFC Body; black = simplified proxy; blue = exact official Gessi 54294 native DWG. Visual review pending.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-native-dwg-linework.json">Native DWG linework</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/official-pdf-verification.json">PDF verification</a><a href="project-context-sanitary-plan.svg">Project plan</a><a href="project-context-front-elevation.svg">Project front</a><a href="project-context-side-elevation.svg">Project side</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://areapro.gessi.com/en/product/54294">Official product</a></nav><h2>Project drawing context</h2><main><article><h2>Sanitary plan</h2><a href="project-context-sanitary-plan-review.svg"><img src="project-context-sanitary-plan-review-preview.png"></a></article><article><h2>Front elevation</h2><a href="project-context-front-elevation-review.svg"><img src="project-context-front-elevation-review-preview.png"></a></article><article><h2>Side elevation</h2><a href="project-context-side-elevation-review.svg"><img src="project-context-side-elevation-review-preview.png"></a></article></main><h2>Three-view native DWG / IFC comparison</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
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
        or access.get("official_product_cad", {}).get("exact_project_configuration_match") is not True
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
    product, vertices, faces = mesh_for_one_product(model, profile["representative_global_id"])
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
    minimum, maximum = bounds_3d(vertices)
    output.mkdir(parents=True, exist_ok=True)
    candidate_views = {}
    views = []
    for view, definition in VIEW_DEFINITIONS.items():
        axes = definition["axes"]
        all_edges = projected_raw_edges(vertices, faces, axes)
        edges = display_edge_sample(all_edges, maximum=1000)
        proxy = projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"]))
        linework_view = linework["views"][view]
        if linework_view.get("path_count") != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi {view} official native-DWG path count drifted")
        official = align_official_paths(view, linework_view["paths_mm"], minimum, maximum)
        cross_check = compatibility(view, official, minimum, maximum)
        if not cross_check["pass"]:
            raise RuntimeError(f"Gessi {view} IFC/DWG compatibility gate failed")
        target = output / f"{view}.svg"
        target.write_text(render_svg(profile, view, edges, proxy, official, {"mesh_face_count": len(faces), "compatibility": cross_check}), encoding="utf-8")
        candidate_views[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": rounded(proxy),
            "official_native_dwg_paths_mm": rounded(official),
            "source_kind": SOURCE_KIND,
            "source_dwg_sha256": SOURCE_DWG_SHA256,
            "official_native_dwg_path_count": len(official),
            "ifc_dwg_compatibility": cross_check,
        }
        views.append({
            "view": view,
            "svg": relative(target),
            "svg_sha256": sha256(target),
            "projection_axes": list(axes),
            "raw_edge_count": len(all_edges),
            "displayed_raw_edge_count": len(edges),
            "silhouette_path_count": len(proxy),
            "drawing_line_source_kind": SOURCE_KIND,
            "drawing_line_source_label_zh": SOURCE_LABEL_ZH,
            "official_cad_path_count": len(official),
            "blue_line_present": True,
            "white_mask_present": True,
            "ifc_dwg_compatibility": cross_check,
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
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "source_dwg": relative(SOURCE_DWG),
        "source_dwg_sha256": SOURCE_DWG_SHA256,
        "official_cad_used": True,
        "third_party_cad_used": False,
        CANDIDATE_IDENTITY_FLAG: False,
        "formal_ifc_write_allowed": False,
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
        "bounds_mm": {
            "minimum": [round(value, 6) for value in minimum],
            "maximum": [round(value, 6) for value in maximum],
            "size": [round(maximum[index] - minimum[index], 6) for index in range(3)],
        },
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "drawing_source": {
            "source_kind": SOURCE_KIND,
            "source_label_zh": SOURCE_LABEL_ZH,
            "source_label_en": SOURCE_LABEL_EN,
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
        "pass": all(view["blue_line_present"] and view["official_cad_path_count"] == EXPECTED_PATH_COUNTS[view["view"]] and view["ifc_dwg_compatibility"]["pass"] for view in views),
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
    write_index(manifest)
    print(relative(manifest_path))


if __name__ == "__main__":
    main()
