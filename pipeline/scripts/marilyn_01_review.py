#!/usr/bin/env python3
"""Generate the Baxter Marilyn 01 native-DWG/high-poly three-view review."""

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
from render_sis04_project_context import write_uncached_png_preview


PROFILE_KEY = "marilyn-01"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OUTPUT_DIR = ROOT / "output/review/highpoly-types/marilyn-01"
REGISTER = OUTPUT_DIR / "profile.json"
ACCESS_RECORD = OUTPUT_DIR / "official-source/source-access-record.json"
SOURCE_REVALIDATION = OUTPUT_DIR / "official-source/official-source-revalidation.json"
SOURCE_REVALIDATION_TYPE = "Marilyn 01"
ADJACENT_VARIANT_GEOMETRY_FIELD = "marilyn_02_geometry_used"
EXACT_3DS_EVIDENCE_FIELD = "bergere_3ds_member_sha256"
LINEWORK = OUTPUT_DIR / "official-native-dwg-linework.json"
SOURCE_LABEL_ZH = "基于 Baxter 精确型号原生 DWG 的官方图纸表达"
SOURCE_LABEL_EN = "official drawing representation from exact Baxter model native DWG"
BLUE = "#1677c8"
EXPECTED_DESCRIPTION = "Bergère armchair with swivel base W86D100H78"
EXPECTED_PATH_COUNTS = {"plan": 24, "front": 68, "side": 54}
SCOPE = "exact Baxter Marilyn bergere 86 x 100 x 94 cm family CAD reference; not a project shop drawing"
ARTICLE_LABEL = "Baxter Marilyn bergere armchair 86 x 100 x 94 cm"
COMPARISON_TOLERANCE_MM = 60.0
GENERATOR = "pipeline/scripts/marilyn_01_review.py"
INDEX_TITLE = "Baxter Marilyn bergere / project Marilyn 01"
INDEX_DESCRIPTION = "Blue line = exact Baxter native DWG bergere cluster. Grey = actual IFC Body. Black = simplified proxy. The project IFC Description height 78 cm is stale; official current sources and Body geometry identify the 86 x 100 x 94 cm bergere."
DISCLOSURE_LINES = (
    "IFC Description height 78 cm is stale.",
    "Official/current and Body identify H94 cm.",
)


def rounded(paths):
    return [[[round(float(x), 6), round(float(y), 6)] for x, y in path] for path in paths]


def path_bounds(paths):
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return {
        "minimum": [round(value, 6) for value in minimum],
        "maximum": [round(value, 6) for value in maximum],
        "size": [round(maximum[axis] - minimum[axis], 6) for axis in range(2)],
    }


def render_svg(profile, view, raw_edges, proxy, official, metadata):
    width, height = 1400, 980
    plot_x, plot_y, plot_w, plot_h = 60, 170, 980, 750
    points = [point for edge in raw_edges for point in edge]
    points.extend(point for path in proxy for point in path)
    points.extend(point for path in official for point in path)
    min_x, max_x = min(point[0] for point in points), max(point[0] for point in points)
    min_y, max_y = min(point[1] for point in points), max(point[1] for point in points)
    padding = max(max_x - min_x, max_y - min_y) * 0.055
    min_x, max_x = min_x - padding, max_x + padding
    min_y, max_y = min_y - padding, max_y + padding
    scale = min(plot_w / (max_x - min_x), plot_h / (max_y - min_y))

    def transform(point):
        return plot_x + (point[0] - min_x) * scale, plot_y + plot_h - (point[1] - min_y) * scale

    edge_path = svg_path([[start, end] for start, end in raw_edges], transform)
    proxy_path = svg_path(proxy, transform, close=True)
    official_path = svg_path(official, transform)
    requirements = "\n".join(
        f'<text x="1080" y="{254 + index * 27}" font-family="Arial,sans-serif" font-size="14" fill="#34495e">• {html.escape(item)}</text>'
        for index, item in enumerate(profile["required_semantics"])
    )
    comparison = metadata["comparison"]
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="1400" height="980" fill="#fbfaf7"/>
<text x="60" y="58" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">{html.escape(profile["display_name"])}</text>
<text x="60" y="96" font-family="Arial,sans-serif" font-size="20" fill="#41566d">{VIEWS[view]["label"]} · isolated representative {profile["representative_global_id"]}</text>
<text x="60" y="128" font-family="Arial,sans-serif" font-size="17" fill="#68798a">Grey = actual IFC Body · Black = simplified proxy · Blue = Baxter native DWG</text>
<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="original-highpoly" d="{edge_path}" fill="none" stroke="#87929c" stroke-width="0.5" stroke-opacity="0.25" vector-effect="non-scaling-stroke"/>
<path class="simplified-proxy-silhouette geometry-derived" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference-mask" d="{official_path}" fill="none" stroke="#ffffff" stroke-width="7.2" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference native-dwg" data-source-kind="native_dwg" data-source-sha256="{metadata["dwg_sha256"]}" d="{official_path}" fill="none" stroke="{BLUE}" stroke-width="2.8" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="1080" y="205" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Acceptance</text>
{requirements}
<text x="1080" y="410" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Mechanical comparison</text>
<text x="1080" y="442" font-family="Arial,sans-serif" font-size="14" fill="#41566d">DWG visible: {comparison["official_visible_size_mm"]}</text>
<text x="1080" y="470" font-family="Arial,sans-serif" font-size="14" fill="#41566d">IFC Body: {comparison["ifc_projection_size_mm"]}</text>
<text x="1080" y="498" font-family="Arial,sans-serif" font-size="14" fill="#41566d">max delta: {comparison["maximum_absolute_delta_mm"]:.3f} mm</text>
<text x="1080" y="526" font-family="Arial,sans-serif" font-size="14" fill="#41566d">tolerance: {comparison["tolerance_mm"]:.1f} mm · pass</text>
<text x="1080" y="586" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Source</text>
<text x="1080" y="618" font-family="Arial,sans-serif" font-size="14" fill="{BLUE}">{SOURCE_LABEL_ZH}</text>
<text x="1080" y="646" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Marilyn_Abaco.dwg · AC1032 · mm</text>
<text x="1080" y="674" font-family="Arial,sans-serif" font-size="14" fill="#41566d">third-party CAD used: false</text>
<text x="1080" y="730" font-family="Arial,sans-serif" font-size="14" fill="#68798a">{html.escape(DISCLOSURE_LINES[0])}</text>
<text x="1080" y="756" font-family="Arial,sans-serif" font-size="14" fill="#68798a">{html.escape(DISCLOSURE_LINES[1])}</text>
<text x="1080" y="812" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Visual review pending; no IFC write.</text>
</svg>'''


def write_index(manifest):
    cards = "".join(f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>' for item in manifest["views"])
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>{html.escape(INDEX_TITLE)} review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>{html.escape(INDEX_TITLE)}</h1><p>{html.escape(INDEX_DESCRIPTION)}</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-native-dwg-linework.json">Native DWG linework</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-furniture-plan.svg">Project furniture plan</a><a href="project-context-r22-front-elevation.svg">R22 front</a><a href="project-context-r22-side-elevation.svg">R22 side</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.baxter.it/en/products/marilyn-sofas-and-armchairs">Official page</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture plan</h2><img src="project-context-furniture-plan-review-preview.png"></article><article><h2>R22 front</h2><img src="project-context-r22-front-elevation-review-preview.png"></article><article><h2>R22 side</h2><img src="project-context-r22-side-elevation-review-preview.png"></article></main><h2>Native DWG / IFC comparison</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--register", type=Path, default=REGISTER)
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    source, output = args.input.resolve(), args.output.resolve()
    if sha256(source) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    profile = load_json(args.register.resolve())["profiles"][PROFILE_KEY]
    access = load_json(ACCESS_RECORD)
    source_revalidation = load_json(SOURCE_REVALIDATION)
    linework = load_json(LINEWORK)
    if (
        profile["drawing_source"].get("source_kind") != "native_dwg"
        or profile["drawing_source"].get("official_cad_used") is not True
        or profile["drawing_source"].get("third_party_cad_used") is not False
        or access.get("drawing_geometry_source", {}).get("source_kind") != "native_dwg"
        or access.get("native_cad_selection", {}).get("pass") is not True
        or access.get("dimension_cross_check", {}).get("pass") is not True
        or linework.get("source_kind") != "native_dwg"
        or linework.get("source_dwg_sha256") != access["native_cad_selection"]["source_dwg_sha256"]
        or linework.get("pass") is not True
        or source_revalidation.get("project_ifc_type_name") != SOURCE_REVALIDATION_TYPE
        or source_revalidation.get(ADJACENT_VARIANT_GEOMETRY_FIELD) is not False
        or source_revalidation.get("product_page", {}).get("all_required_identity_tokens_found") is not True
        or source_revalidation.get("native_zip", {}).get("dwg_member_sha256") != linework.get("source_dwg_sha256")
        or source_revalidation.get("native_zip", {}).get(EXACT_3DS_EVIDENCE_FIELD) != linework.get("exact_model_3ds_sha256")
        or source_revalidation.get("technical_pdf", {}).get("current_identity_page_render_matches_archive") is not True
        or source_revalidation.get("pass") is not True
    ):
        raise RuntimeError("Marilyn native-DWG source gate failed")
    model = ifcopenshell.open(source)
    product, vertices, faces = mesh_for_one_product(model, profile["representative_global_id"])
    if product_type_name(product) != profile["ifc_type_name"]:
        raise RuntimeError("Marilyn representative type identity drifted")
    product_type = next(relation.RelatingType for relation in product.IsTypedBy)
    if product_type.Description != EXPECTED_DESCRIPTION:
        raise RuntimeError("Marilyn IFC Description drifted; update the disclosed source conflict")
    minimum, maximum = bounds_3d(vertices)
    size = [round(maximum[index] - minimum[index], 6) for index in range(3)]
    expected_size = access["dimension_cross_check"]["project_ifc_body_local_xyz_mm"]
    if any(abs(size[index] - expected_size[index]) > 0.001 for index in range(3)):
        raise RuntimeError("Marilyn Body dimensions drifted")
    output.mkdir(parents=True, exist_ok=True)
    candidate_views, views, comparisons = {}, [], {}
    for view, definition in VIEWS.items():
        axes = definition["axes"]
        all_edges = projected_raw_edges(vertices, faces, axes)
        edges = display_edge_sample(all_edges)
        proxy = projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"]))
        native_official = linework["views"][view]["paths_mm"]
        # The catalogue side view looks toward the opposite local-Y direction
        # from the project IFC YZ convention. Mirror only X in the 2D review /
        # representation coordinate system; no scale or shape fitting occurs.
        official = (
            [[[-float(x), float(y)] for x, y in path] for path in native_official]
            if view == "side"
            else native_official
        )
        if len(official) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Marilyn {view} native path count drifted")
        official_size = path_bounds(official)["size"]
        ifc_size = [size[axes[0]], size[axes[1]]]
        deltas = [abs(ifc_size[index] - official_size[index]) for index in range(2)]
        comparison = {
            "official_visible_size_mm": official_size,
            "ifc_projection_size_mm": ifc_size,
            "absolute_delta_mm": [round(value, 6) for value in deltas],
            "maximum_absolute_delta_mm": round(max(deltas), 6),
            "tolerance_mm": COMPARISON_TOLERANCE_MM,
            "geometry_scaled_to_match": False,
            "pass": max(deltas) <= COMPARISON_TOLERANCE_MM,
        }
        comparisons[view] = comparison
        target = output / f"{view}.svg"
        target.write_text(render_svg(profile, view, edges, proxy, official, {"comparison": comparison, "dwg_sha256": linework["source_dwg_sha256"]}), encoding="utf-8")
        preview = output / f"{view}-preview.png"
        write_uncached_png_preview(target, preview)
        candidate_views[view] = {
            "projection_axes": list(axes),
            "official_native_dwg_paths_mm": rounded(official),
            "proxy_paths_mm": rounded(proxy),
            "source_kind": "native_dwg",
            "source_dwg_sha256": linework["source_dwg_sha256"],
            "view_orientation_transform": "mirror_x_for_ifc_yz_direction" if view == "side" else "none",
            "geometry_scaled_to_match_ifc": False,
            "third_party_cad_used": False,
        }
        views.append({
            "view": view,
            "svg": relative(target),
            "svg_sha256": sha256(target),
            "preview": relative(preview),
            "preview_sha256": sha256(preview),
            "projection_axes": list(axes),
            "raw_edge_count": len(all_edges),
            "displayed_raw_edge_count": len(edges),
            "silhouette_path_count": len(proxy),
            "official_cad_path_count": len(official),
            "blue_line_source_kind": "native_dwg",
            "blue_line_native_dwg_sha256": linework["source_dwg_sha256"],
            "white_mask_below_blue": True,
            "view_orientation_transform": "mirror_x_for_ifc_yz_direction" if view == "side" else "none",
            "geometry_scaled_to_match_ifc": False,
            "blue_line_present": True,
            "comparison": comparison,
        })
    if not all(item["pass"] for item in comparisons.values()):
        raise RuntimeError("Marilyn native DWG / IFC visible-bounds comparison failed")
    candidate_path = output / "candidate-representations.json"
    write_json(candidate_path, {
        "schema_version": 1,
        "profile_key": PROFILE_KEY,
        "representative_global_id": product.GlobalId,
        "ifc_type_name": profile["ifc_type_name"],
        "article_number": ARTICLE_LABEL,
        "units": "mm",
        "source_kind": "native_dwg",
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "drawing_source": {
            "source_kind": "native_dwg",
            "source_label_zh": SOURCE_LABEL_ZH,
            "source_label_en": SOURCE_LABEL_EN,
            "official_cad_used": True,
            "third_party_cad_used": False,
        },
        "source_dwg": linework["source_dwg"],
        "source_dwg_sha256": linework["source_dwg_sha256"],
        "formal_ifc_write_allowed": False,
        "review_status": "visual_review_pending",
        "views": candidate_views,
    })
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": GENERATOR,
        "formal_ifc": relative(source),
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(source),
        "formal_ifc_bytes_unchanged": sha256(source) == FORMAL_SHA256,
        "profile_key": PROFILE_KEY,
        "display_name": profile["display_name"],
        "ifc_type_name": profile["ifc_type_name"],
        "ifc_type_description": product_type.Description,
        "representative_global_id": product.GlobalId,
        "geometry_product_count": 1,
        "whole_model_render": False,
        "mesh_vertex_count": len(vertices),
        "mesh_face_count": len(faces),
        "bounds_mm": {"minimum": [round(value, 6) for value in minimum], "maximum": [round(value, 6) for value in maximum], "size": size},
        "source_kind": "native_dwg",
        "source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_acquired": True,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "drawing_source": {
            "source_kind": "native_dwg",
            "source_label_zh": SOURCE_LABEL_ZH,
            "source_label_en": SOURCE_LABEL_EN,
            "official_cad_used": True,
            "third_party_cad_used": False,
        },
        "source_dwg": linework["source_dwg"],
        "source_dwg_sha256": linework["source_dwg_sha256"],
        "official_native_dwg_linework": relative(LINEWORK),
        "official_native_dwg_linework_sha256": sha256(LINEWORK),
        "official_source_access_record": relative(ACCESS_RECORD),
        "official_source_access_record_sha256": sha256(ACCESS_RECORD),
        "official_source_revalidation": {
            "path": relative(SOURCE_REVALIDATION),
            "sha256": sha256(SOURCE_REVALIDATION),
            "checked_at": source_revalidation["checked_at"],
            "product_page_identity_tokens_pass": source_revalidation["product_page"]["all_required_identity_tokens_found"],
            "native_zip_bytes_match_archive": source_revalidation["native_zip"]["downloaded_bytes_match_local_archive"],
            "technical_pdf_identity_page_render_matches_archive": source_revalidation["technical_pdf"]["current_identity_page_render_matches_archive"],
            "measurement_svgs_pass": all(item["pass"] for item in source_revalidation["measurement_svgs"].values()),
            ADJACENT_VARIANT_GEOMETRY_FIELD: source_revalidation[ADJACENT_VARIANT_GEOMETRY_FIELD],
            "pass": source_revalidation["pass"],
        },
        "dimension_cross_check": access["dimension_cross_check"],
        "native_dwg_ifc_view_comparisons": comparisons,
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "derived_ifc_write_allowed": False,
        "formal_ifc_write": "not performed",
        "scope": SCOPE,
        "views": views,
        "pass": sha256(source) == FORMAL_SHA256 and all(item["comparison"]["pass"] for item in views),
    }
    context_path = output / "project-context-manifest.json"
    if context_path.is_file():
        context = load_json(context_path)
        manifest["project_context"] = {"manifest": relative(context_path), "manifest_sha256": sha256(context_path), "walls_and_surrounding_project_elements_retained": context.get("walls_and_surrounding_project_elements_retained"), "overlay_top_layer_with_white_mask": context.get("overlay_top_layer_with_white_mask"), "pass": context.get("pass")}
        manifest["pass"] = manifest["pass"] and context.get("pass") is True
    bonsai_path = output / "bonsai-review-manifest.json"
    if bonsai_path.is_file():
        bonsai = load_json(bonsai_path)
        manifest["bonsai_review"] = {"manifest": relative(bonsai_path), "manifest_sha256": sha256(bonsai_path), "mode": bonsai.get("mode"), "saved_active_representation": bonsai.get("bonsai_session", {}).get("saved_active_representation"), "saved_camera_count": bonsai.get("bonsai_session", {}).get("saved_camera_count"), "pass": bonsai.get("pass")}
        manifest["pass"] = manifest["pass"] and bonsai.get("pass") is True
    manifest_path = output / "manifest.json"
    write_json(manifest_path, manifest)
    write_index(manifest)
    print(manifest_path)


if __name__ == "__main__":
    main()
