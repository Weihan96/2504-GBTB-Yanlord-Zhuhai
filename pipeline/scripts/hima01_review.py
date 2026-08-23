#!/usr/bin/env python3
"""Generate a single-instance Hima screen review without claiming official CAD."""

from __future__ import annotations

import argparse
import html
import json
from collections import defaultdict
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


PROFILE_KEY = "hima01"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OUTPUT_DIR = ROOT / "output/review/highpoly-types/hima01"
REGISTER = OUTPUT_DIR / "profile.json"
ACCESS_RECORD = OUTPUT_DIR / "official-source/source-access-record.json"
PAGE_EVIDENCE = OUTPUT_DIR / "official-source/official-product-page-evidence.json"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
SOURCE_LABEL_EN = "simplified drawing representation derived from the original high-poly geometry"


def rounded(paths):
    return [
        [[round(float(x), 6), round(float(y), 6)] for x, y in path]
        for path in paths
    ]


def screen_frame_components(vertices, faces):
    parent = list(range(len(vertices)))

    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(first, second):
        first, second = find(first), find(second)
        if first != second:
            parent[second] = first

    for first, second, third in faces:
        union(first, second)
        union(second, third)
    groups = defaultdict(list)
    for index, point in enumerate(vertices):
        groups[find(index)].append(point)
    records = []
    for points in groups.values():
        minimum = [min(point[axis] for point in points) for axis in range(3)]
        maximum = [max(point[axis] for point in points) for axis in range(3)]
        size = [maximum[axis] - minimum[axis] for axis in range(3)]
        horizontal = sorted(size[:2])
        if size[2] > 950.0 and horizontal[1] > 700.0 and horizontal[0] < 40.0:
            records.append({
                "vertex_count": len(points),
                "minimum_mm": [round(value, 6) for value in minimum],
                "maximum_mm": [round(value, 6) for value in maximum],
                "size_mm": [round(value, 6) for value in size],
                "centre_mm": [round((minimum[axis] + maximum[axis]) / 2.0, 6) for axis in range(3)],
            })
    return sorted(records, key=lambda item: item["centre_mm"][:2])


def render_svg(profile, view, raw_edges, proxy, metadata):
    width, height = 1400, 980
    plot_x, plot_y, plot_w, plot_h = 60, 170, 980, 750
    points = [point for edge in raw_edges for point in edge]
    points.extend(point for path in proxy for point in path)
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
    requirements = "\n".join(
        f'<text x="1090" y="{258 + index * 28}" font-family="Arial,sans-serif" font-size="15" fill="#34495e">• {html.escape(item)}</text>'
        for index, item in enumerate(profile["required_semantics"])
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="1400" height="980" fill="#fbfaf7"/>
<text x="60" y="58" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">{html.escape(profile["display_name"])}</text>
<text x="60" y="96" font-family="Arial,sans-serif" font-size="20" fill="#41566d">{VIEWS[view]["label"]} · isolated representative {profile["representative_global_id"]}</text>
<text x="60" y="128" font-family="Arial,sans-serif" font-size="17" fill="#68798a">Grey = actual IFC Body · Black = geometry-derived simplified proxy · Blue = none (official DWG not acquired)</text>
<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="original-highpoly" d="{edge_path}" fill="none" stroke="#87929c" stroke-width="0.5" stroke-opacity="0.27" vector-effect="non-scaling-stroke"/>
<path class="simplified-proxy-silhouette geometry-derived" data-source-kind="geometry_derived_simplified_proxy" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="3" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="1080" y="205" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Acceptance</text>
{requirements}
<text x="1080" y="410" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Source status</text>
<text x="1080" y="442" font-family="Arial,sans-serif" font-size="14" fill="#111820">{SOURCE_LABEL_ZH}</text>
<text x="1080" y="470" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Official 2D DWG published, not acquired</text>
<text x="1080" y="498" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Official CAD geometry used: false</text>
<text x="1080" y="526" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Third-party CAD used: false</text>
<text x="1080" y="586" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Isolation</text>
<text x="1080" y="618" font-family="Arial,sans-serif" font-size="15" fill="#41566d">geometry products: 1</text>
<text x="1080" y="646" font-family="Arial,sans-serif" font-size="15" fill="#41566d">whole model render: false</text>
<text x="1080" y="674" font-family="Arial,sans-serif" font-size="15" fill="#41566d">mesh faces: {metadata["mesh_face_count"]}</text>
<text x="1080" y="702" font-family="Arial,sans-serif" font-size="15" fill="#41566d">derived proxy paths: {len(proxy)}</text>
<text x="1080" y="758" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Official page is identity evidence only.</text>
<text x="1080" y="784" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Visual review pending; no IFC write.</text>
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
        f'''<!doctype html><html><meta charset="utf-8"><title>Poliform Hima review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Poliform Hima / HIMA01</h1><p>Black line = {SOURCE_LABEL_EN}. No blue line is shown because the official DWG has not been acquired. Visual review pending.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Official-source access record</a><a href="official-source/official-product-page-evidence.json">Official-page evidence</a><a href="project-context-furniture-plan.svg">Full project plan</a><a href="project-context-front-elevation.svg">Project front elevation</a><a href="project-context-side-elevation.svg">Project side elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.poliform.it/en/products/hima/">Official product page</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article><article><h2>Front elevation</h2><a href="project-context-front-elevation-review.svg"><img src="project-context-front-elevation-review-preview.png"></a></article><article><h2>Side elevation</h2><a href="project-context-side-elevation-review.svg"><img src="project-context-side-elevation-review-preview.png"></a></article></main><h2>High-poly / simplified proxy review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
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
    page_evidence = load_json(PAGE_EVIDENCE)
    drawing_source = profile["drawing_source"]
    if (
        drawing_source.get("source_kind") != "geometry_derived_simplified_proxy"
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or access.get("official_2d_dwg", {}).get("acquired") is not False
        or access.get("drawing_geometry_source", {}).get("source_kind") != "geometry_derived_simplified_proxy"
        or page_evidence.get("drawing_geometry_conclusion", {}).get("official_native_dwg_acquired") is not False
        or page_evidence.get("drawing_geometry_conclusion", {}).get("third_party_cad_used") is not False
        or sha256(PAGE_EVIDENCE) != access.get("official_product_page_evidence", {}).get("sha256")
        or access.get("official_dimension_cross_check", {}).get("pass") is not True
    ):
        raise RuntimeError("Hima no-official-CAD source gate failed")
    model = ifcopenshell.open(source)
    product, vertices, faces = mesh_for_one_product(model, profile["representative_global_id"])
    if product_type_name(product) != profile["ifc_type_name"]:
        raise RuntimeError("Hima representative type identity drifted")
    instances = sorted(
        item.GlobalId
        for item in model.by_type(product.is_a())
        if product_type_name(item) == profile["ifc_type_name"]
    )
    if instances != profile["expected_instance_global_ids"]:
        raise RuntimeError("Hima instance set drifted")
    minimum, maximum = bounds_3d(vertices)
    frames = screen_frame_components(vertices, faces)
    height_mm = maximum[2] - minimum[2]
    identity_checks = {
        "ifc_type_description": {
            "expected": "Poliform HIMA 3 elements 1m",
            "actual": next(relation.RelatingType.Description for relation in product.IsTypedBy),
            "pass": next(relation.RelatingType.Description for relation in product.IsTypedBy) == "Poliform HIMA 3 elements 1m",
        },
        "screen_element_count": {
            "method": "connected high-poly frame components with height >950 mm, long horizontal span >700 mm and thickness <40 mm",
            "expected": 3,
            "actual": len(frames),
            "components": frames,
            "pass": len(frames) == 3,
        },
        "nominal_height": {
            "official_nominal_mm": 1000.0,
            "ifc_body_height_mm": round(height_mm, 6),
            "absolute_delta_mm": round(abs(height_mm - 1000.0), 6),
            "tolerance_mm": 2.0,
            "pass": abs(height_mm - 1000.0) <= 2.0,
        },
    }
    if not all(check["pass"] for check in identity_checks.values()):
        raise RuntimeError("Hima identity checks failed")
    output.mkdir(parents=True, exist_ok=True)
    candidate_views = {}
    views = []
    for view, definition in VIEWS.items():
        axes = definition["axes"]
        all_edges = projected_raw_edges(vertices, faces, axes)
        edges = display_edge_sample(all_edges)
        proxy = projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"]))
        target = output / f"{view}.svg"
        target.write_text(render_svg(profile, view, edges, proxy, {"mesh_face_count": len(faces)}), encoding="utf-8")
        candidate_views[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": rounded(proxy),
            "source_kind": "geometry_derived_simplified_proxy",
            "official_cad_paths_mm": [],
        }
        views.append({
            "view": view,
            "svg": relative(target),
            "svg_sha256": sha256(target),
            "projection_axes": list(axes),
            "raw_edge_count": len(all_edges),
            "displayed_raw_edge_count": len(edges),
            "silhouette_path_count": len(proxy),
            "drawing_line_source_kind": "geometry_derived_simplified_proxy",
            "drawing_line_source_label_zh": SOURCE_LABEL_ZH,
            "official_cad_path_count": 0,
            "blue_line_present": False,
        })
    candidate_path = output / "candidate-representations.json"
    write_json(candidate_path, {
        "schema_version": 1,
        "profile_key": PROFILE_KEY,
        "representative_global_id": product.GlobalId,
        "ifc_type_name": profile["ifc_type_name"],
        "units": "mm",
        "source_kind": "geometry_derived_simplified_proxy",
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_used": False,
        "third_party_cad_used": False,
        "formal_ifc_write_allowed": False,
        "review_status": "visual_review_pending",
        "views": candidate_views,
    })
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/hima01_review.py",
        "formal_ifc": relative(source),
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(source),
        "formal_ifc_bytes_unchanged": sha256(source) == FORMAL_SHA256,
        "profile_register": relative(args.register.resolve()),
        "profile_register_sha256": sha256(args.register.resolve()),
        "profile_key": PROFILE_KEY,
        "display_name": profile["display_name"],
        "ifc_type_name": profile["ifc_type_name"],
        "ifc_type_description": "Poliform HIMA 3 elements 1m",
        "representative_global_id": product.GlobalId,
        "registered_instance_global_ids": instances,
        "geometry_product_count": 1,
        "whole_model_render": False,
        "formal_ifc_write": False,
        "source_kind": "geometry_derived_simplified_proxy",
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_used": False,
        "third_party_cad_used": False,
        "blue_line_present": False,
        "mesh_vertex_count": len(vertices),
        "mesh_face_count": len(faces),
        "local_bounds_mm": {"minimum": minimum, "maximum": maximum},
        "identity_checks": identity_checks,
        "drawing_source": {
            **drawing_source,
            "official_source_access_record": relative(ACCESS_RECORD),
            "official_source_access_record_sha256": sha256(ACCESS_RECORD),
            "official_product_page_evidence": relative(PAGE_EVIDENCE),
            "official_product_page_evidence_sha256": sha256(PAGE_EVIDENCE),
            "official_2d_dwg_status": access["official_2d_dwg"]["access_status"],
        },
        "official_reference": profile["official_reference"],
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "views": views,
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "pass": True,
    }
    context_manifest = output / "project-context-manifest.json"
    bonsai_manifest = output / "bonsai-review-manifest.json"
    if context_manifest.is_file():
        context = load_json(context_manifest)
        if context.get("pass") is not True:
            raise RuntimeError("Hima project-context evidence failed")
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
            raise RuntimeError("Hima Bonsai camera evidence failed")
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
    print(json.dumps({"manifest": relative(manifest_path), "views": [item["svg"] for item in views], "pass": True}, indent=2))


if __name__ == "__main__":
    main()
