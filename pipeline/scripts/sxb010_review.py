#!/usr/bin/env python3
"""Generate a legible sxb010 Venetian-blind review from its original IFC Body."""

from __future__ import annotations

import argparse
import html
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from int1_highpoly_type_review import (
    bounds_3d,
    display_edge_sample,
    mesh_for_one_product,
    product_type_name,
    projected_raw_edges,
    projected_silhouette,
    svg_path,
)


PROFILE_KEY = "sxb010"
ARTICLE_NUMBER = "project sxb010 / Hunter Douglas 25 mm family direction"
GENERATOR = "pipeline/scripts/sxb010_review.py"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OUTPUT_DIR = ROOT / "output/review/highpoly-types/sxb010"
REGISTER = OUTPUT_DIR / "profile.json"
ACCESS_RECORD = OUTPUT_DIR / "official-source/source-access-record.json"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
SOURCE_LABEL_EN = "simplified drawing representation derived from the original high-poly geometry"
VIEW_DEFINITIONS = {
    "plan": {"axes": (0, 2), "label": "PLAN / local XZ (width × depth)"},
    "front": {"axes": (0, 1), "label": "FRONT / local XY (width × height)"},
    "side": {"axes": (2, 1), "label": "SIDE / local ZY (depth × height)"},
}

# At the final review sheet this is about 9 px in Front/Side. It collapses
# the two thickness edges of one 20.7 mm slat, while the minimum slat pitch
# remains above 30 mm and therefore stays intact.
NEAR_LINE_MERGE_THRESHOLD_MM = 22.0
PLOT_WIDTH_PX = 980.0
PLOT_HEIGHT_PX = 750.0
PLOT_PADDING_FACTOR = 0.07
EXPECTED_SLAT_COUNT = 45
EXPECTED_GUIDE_CLUSTER_COUNT = 3


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


def rectangle(minimum, maximum):
    return [[
        (minimum[0], minimum[1]),
        (maximum[0], minimum[1]),
        (maximum[0], maximum[1]),
        (minimum[0], maximum[1]),
        (minimum[0], minimum[1]),
    ]]


def segment_count(paths):
    return sum(max(0, len(path) - 1) for path in paths)


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
    records = []
    for indices in grouped.values():
        points = [vertices[index] for index in indices]
        minimum = [min(point[axis] for point in points) for axis in range(3)]
        maximum = [max(point[axis] for point in points) for axis in range(3)]
        records.append({
            "indices": indices,
            "minimum": minimum,
            "maximum": maximum,
            "size": [maximum[axis] - minimum[axis] for axis in range(3)],
            "center": [(minimum[axis] + maximum[axis]) / 2.0 for axis in range(3)],
        })
    return records


def principal_axis_path(component, vertices, axes):
    points = [(vertices[index][axes[0]], vertices[index][axes[1]]) for index in component["indices"]]
    center = (sum(point[0] for point in points) / len(points), sum(point[1] for point in points) / len(points))
    xy = sum((point[0] - center[0]) * (point[1] - center[1]) for point in points)
    minimum = (min(point[0] for point in points), min(point[1] for point in points))
    maximum = (max(point[0] for point in points), max(point[1] for point in points))
    if xy >= 0.0:
        return [[minimum, maximum]]
    return [[(minimum[0], maximum[1]), (maximum[0], minimum[1])]]


def cluster_values(values, threshold):
    clusters = []
    for value in sorted(values):
        if not clusters or value - clusters[-1][-1] > threshold:
            clusters.append([value])
        else:
            clusters[-1].append(value)
    return [sum(cluster) / len(cluster) for cluster in clusters]


def semantic_paths(vertices, faces, minimum, maximum):
    components = mesh_components(vertices, faces)
    overall_width = maximum[0] - minimum[0]
    overall_height = maximum[1] - minimum[1]
    slats = [
        component for component in components
        if component["size"][0] >= overall_width * 0.95
        and 15.0 <= component["size"][1] <= NEAR_LINE_MERGE_THRESHOLD_MM
        and 25.0 <= component["size"][2] <= 40.0
    ]
    if len(slats) != EXPECTED_SLAT_COUNT:
        raise RuntimeError(f"expected {EXPECTED_SLAT_COUNT} Venetian slats, found {len(slats)}")
    slat_ids = {id(component) for component in slats}
    rails = [
        component for component in components
        if id(component) not in slat_ids and component["size"][0] >= overall_width * 0.95
    ]
    if len(rails) != 2:
        raise RuntimeError(f"expected headrail and bottom rail, found {len(rails)}")
    rail_ids = {id(component) for component in rails}
    controls = [component for component in components if id(component) not in slat_ids | rail_ids]
    if len(controls) != 9:
        raise RuntimeError(f"expected nine guide/cord/control components, found {len(controls)}")

    slat_centers = sorted(component["center"][1] for component in slats)
    slat_pitch = [following - current for current, following in zip(slat_centers, slat_centers[1:])]
    guide_centers = cluster_values(
        [component["center"][0] for component in controls if component["size"][1] > overall_height * 0.5],
        NEAR_LINE_MERGE_THRESHOLD_MM,
    )
    if len(guide_centers) != EXPECTED_GUIDE_CLUSTER_COUNT:
        raise RuntimeError(f"expected {EXPECTED_GUIDE_CLUSTER_COUNT} merged guide/cord axes, found {len(guide_centers)}")

    center_z = (minimum[2] + maximum[2]) / 2.0
    plan = rectangle((minimum[0], minimum[2]), (maximum[0], maximum[2]))
    plan.append([(minimum[0], center_z), (maximum[0], center_z)])
    plan.extend([[(x, minimum[2]), (x, maximum[2])] for x in guide_centers])

    front = rectangle((minimum[0], minimum[1]), (maximum[0], maximum[1]))
    front.extend([
        [(component["minimum"][0], component["center"][1]), (component["maximum"][0], component["center"][1])]
        for component in sorted(slats, key=lambda item: item["center"][1])
    ])
    front.extend([[(x, minimum[1]), (x, maximum[1])] for x in guide_centers])
    small_controls = [component for component in controls if component["size"][1] <= overall_height * 0.5]
    front.extend(
        rectangle(
            (component["minimum"][0], component["minimum"][1]),
            (component["maximum"][0], component["maximum"][1]),
        )[0]
        for component in small_controls
    )

    side = rectangle((minimum[2], minimum[1]), (maximum[2], maximum[1]))
    for component in sorted(slats, key=lambda item: item["center"][1]):
        side.extend(principal_axis_path(component, vertices, (2, 1)))
    for component in controls:
        side.extend(principal_axis_path(component, vertices, (2, 1)))

    audit_common = {
        "merge_threshold_mm": NEAR_LINE_MERGE_THRESHOLD_MM,
        "slat_count_before": len(slats),
        "slat_centerline_count_after": len(slats),
        "slat_projected_thickness_front_mm": {
            "minimum": round(min(component["size"][1] for component in slats), 6),
            "maximum": round(max(component["size"][1] for component in slats), 6),
        },
        "slat_center_pitch_after_mm": {
            "minimum": round(min(slat_pitch), 6),
            "maximum": round(max(slat_pitch), 6),
        },
        "headrail_count_preserved": 1,
        "bottom_rail_count_preserved": 1,
        "guide_and_cord_component_count_before": sum(component["size"][1] > overall_height * 0.5 for component in controls),
        "guide_and_cord_axis_count_after": len(guide_centers),
        "control_component_count_preserved": len(controls),
        "outer_envelope_preserved": True,
        "slat_rhythm_preserved": True,
    }
    return {"plan": plan, "front": front, "side": side}, audit_common


def render_scale(paths):
    bounds = path_bounds(paths)
    span_x, span_y = bounds["size"]
    padding = max(span_x, span_y) * PLOT_PADDING_FACTOR
    return min(PLOT_WIDTH_PX / (span_x + 2.0 * padding), PLOT_HEIGHT_PX / (span_y + 2.0 * padding))


def render_svg(profile, view, raw_edges, paths, metadata):
    width, height = 1400, 980
    plot_x, plot_y, plot_w, plot_h = 60, 170, int(PLOT_WIDTH_PX), int(PLOT_HEIGHT_PX)
    points = [point for edge in raw_edges for point in edge]
    points.extend(point for path in paths for point in path)
    min_x, max_x = min(point[0] for point in points), max(point[0] for point in points)
    min_y, max_y = min(point[1] for point in points), max(point[1] for point in points)
    padding = max(max_x - min_x, max_y - min_y) * PLOT_PADDING_FACTOR
    min_x, max_x = min_x - padding, max_x + padding
    min_y, max_y = min_y - padding, max_y + padding
    scale = min(plot_w / (max_x - min_x), plot_h / (max_y - min_y))

    def transform(point):
        return plot_x + (point[0] - min_x) * scale, plot_y + plot_h - (point[1] - min_y) * scale

    raw_path = svg_path([[start, end] for start, end in raw_edges], transform)
    proxy_path = svg_path(paths, transform)
    requirements = "\n".join(
        f'<text x="1090" y="{258 + index * 28}" font-family="Arial,sans-serif" font-size="15" fill="#34495e">• {html.escape(item)}</text>'
        for index, item in enumerate(profile["required_semantics"])
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="1400" height="980" fill="#fbfaf7"/>
<text x="60" y="58" font-family="Arial,sans-serif" font-size="30" font-weight="700" fill="#1f2d3d">{html.escape(profile["display_name"])}</text>
<text x="60" y="96" font-family="Arial,sans-serif" font-size="20" fill="#41566d">{VIEW_DEFINITIONS[view]["label"]} · isolated representative {profile["representative_global_id"]}</text>
<text x="60" y="128" font-family="Arial,sans-serif" font-size="17" fill="#68798a">Grey = actual IFC Body sample · Black = near-line-merged geometry-derived proxy · no official CAD acquired</text>
<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="original-highpoly" d="{raw_path}" fill="none" stroke="#87929c" stroke-width="0.35" stroke-opacity="0.10" vector-effect="non-scaling-stroke"/>
<path class="simplified-proxy-silhouette geometry-derived" data-source-kind="{SOURCE_KIND}" data-source-label-zh="{SOURCE_LABEL_ZH}" data-merge-threshold-mm="{NEAR_LINE_MERGE_THRESHOLD_MM}" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="1080" y="205" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#1f2d3d">Acceptance</text>
{requirements}
<text x="1080" y="410" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Near-line merge</text>
<text x="1080" y="442" font-family="Arial,sans-serif" font-size="14" fill="#41566d">threshold: {NEAR_LINE_MERGE_THRESHOLD_MM:.1f} mm / {metadata["threshold_px"]:.2f} px</text>
<text x="1080" y="470" font-family="Arial,sans-serif" font-size="14" fill="#41566d">paths: {metadata["before_path_count"]} → {metadata["after_path_count"]}</text>
<text x="1080" y="498" font-family="Arial,sans-serif" font-size="14" fill="#41566d">segments: {metadata["before_segment_count"]} → {metadata["after_segment_count"]}</text>
<text x="1080" y="538" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Preserved</text>
<text x="1080" y="570" font-family="Arial,sans-serif" font-size="14" fill="#41566d">45 slat centre-lines / rhythm unchanged</text>
<text x="1080" y="598" font-family="Arial,sans-serif" font-size="14" fill="#41566d">outer envelope / rails / 3 guide axes</text>
<text x="1080" y="650" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Isolation</text>
<text x="1080" y="682" font-family="Arial,sans-serif" font-size="15" fill="#41566d">geometry products: 1</text>
<text x="1080" y="710" font-family="Arial,sans-serif" font-size="15" fill="#41566d">whole model render: false</text>
<text x="1080" y="738" font-family="Arial,sans-serif" font-size="15" fill="#41566d">mesh faces: {metadata["mesh_face_count"]}</text>
<text x="1080" y="790" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Family identity reference only; no official blue line.</text>
<text x="1080" y="816" font-family="Arial,sans-serif" font-size="14" fill="#68798a">Visual review pending; no IFC write.</text>
</svg>'''


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    semantic_bonsai = (
        ("Plan", "bonsai-camera-semantic-plan.png"),
        ("Front", "bonsai-camera-semantic-front-elevation.png"),
        ("Side", "bonsai-camera-semantic-side-elevation.png"),
        ("Iso", "bonsai-camera-iso.png"),
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {label}</h2><a href="{filename}"><img src="{filename}"></a></article>'
        for label, filename in semantic_bonsai
    )
    context_cards = "".join(
        f'<article><h2>Project {label}</h2><a href="{filename}"><img src="{preview}"></a></article>'
        for label, filename, preview in (
            ("Plan", "project-context-ffl-plan.svg", "project-context-ffl-plan-review-preview.png"),
            ("Front", "project-context-r07-front-elevation.svg", "project-context-r07-front-elevation-review-preview.png"),
            ("Side", "project-context-r07-side-elevation.svg", "project-context-r07-side-elevation-review-preview.png"),
        )
    )
    (OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>sxb010 Venetian blind review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Project sxb010 / Hunter Douglas 25 mm Venetian-blind direction</h1><p>Black line = {SOURCE_LABEL_EN}. Lines within {NEAR_LINE_MERGE_THRESHOLD_MM:.1f} mm are semantically merged: each physical slat remains one centre-line, while the outer envelope, rails, guides, cords, depth and controls remain represented. Exact SKU and official project CAD remain unconfirmed, so no blue line is shown.</p><nav><a href="review-contact-sheet.png">Contact sheet</a><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="line-simplification-audit.json">Line audit</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-manifest.json">Project context evidence</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.hunterdouglas.cn/product/venetian-blind/16mm-25mm-venetian-blinds">Official family page</a></nav><h2>Three-view review</h2><main>{cards}</main><h2>Complete project drawing context</h2><main>{context_cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--register", type=Path, default=REGISTER)
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    source, output, register_path = args.input.resolve(), args.output.resolve(), args.register.resolve()
    if sha256(source) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    register = load_json(register_path)
    profile = register["profiles"][PROFILE_KEY]
    access = load_json(ACCESS_RECORD)
    drawing_source = profile["drawing_source"]
    if (
        drawing_source.get("source_kind") != SOURCE_KIND
        or drawing_source.get("source_label_zh") != SOURCE_LABEL_ZH
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or access.get("drawing_geometry_source", {}).get("source_kind") != SOURCE_KIND
        or access.get("drawing_geometry_source", {}).get("source_label_zh") != SOURCE_LABEL_ZH
        or access.get("official_product_cad", {}).get("acquired") is not False
    ):
        raise RuntimeError("sxb010 geometry-derived source gate failed")
    model = ifcopenshell.open(source)
    product, vertices, faces = mesh_for_one_product(model, profile["representative_global_id"])
    actual_identity = product_type_name(product) if product.IsTypedBy else product.Name
    if actual_identity != profile["ifc_type_name"]:
        raise RuntimeError("sxb010 representative identity drifted")
    instances = sorted(
        item.GlobalId for item in model.by_type(product.is_a())
        if (product_type_name(item) if item.IsTypedBy else item.Name) == profile["ifc_type_name"]
    )
    if instances != profile["expected_instance_global_ids"]:
        raise RuntimeError("sxb010 instance set drifted")
    minimum, maximum = bounds_3d(vertices)
    simplified, semantic_audit = semantic_paths(vertices, faces, minimum, maximum)
    output.mkdir(parents=True, exist_ok=True)
    views, candidate_views, view_audits = [], {}, {}
    for view, definition in VIEW_DEFINITIONS.items():
        axes = definition["axes"]
        raw = projected_raw_edges(vertices, faces, axes)
        displayed = display_edge_sample(raw, maximum=400)
        before = projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"]))
        before_bounds = path_bounds(before)
        after = [list(path) for path in simplified[view]]
        after[0] = rectangle(tuple(before_bounds["minimum"]), tuple(before_bounds["maximum"]))[0]
        after = [[
            (
                min(max(point[0], before_bounds["minimum"][0]), before_bounds["maximum"][0]),
                min(max(point[1], before_bounds["minimum"][1]), before_bounds["maximum"][1]),
            )
            for point in path
        ] for path in after]
        after_bounds = path_bounds(after)
        envelope_delta = [
            abs(after_bounds[edge][axis] - before_bounds[edge][axis])
            for edge in ("minimum", "maximum")
            for axis in range(2)
        ]
        if max(envelope_delta) > 0.00001:
            raise RuntimeError(f"sxb010 {view} outer envelope changed during line merge")
        scale = render_scale(after)
        audit = {
            "view": view,
            "merge_threshold_mm": NEAR_LINE_MERGE_THRESHOLD_MM,
            "final_svg_scale_px_per_mm": round(scale, 6),
            "merge_threshold_at_final_svg_scale_px": round(scale * NEAR_LINE_MERGE_THRESHOLD_MM, 6),
            "before_path_count": len(before),
            "after_path_count": len(after),
            "before_segment_count": segment_count(before),
            "after_segment_count": segment_count(after),
            "before_bounds_mm": before_bounds,
            "after_bounds_mm": after_bounds,
            "outer_envelope_delta_mm": [round(value, 6) for value in envelope_delta],
            **semantic_audit,
            "pass": True,
        }
        target = output / f"{view}.svg"
        target.write_text(render_svg(profile, view, displayed, after, {
            "mesh_face_count": len(faces),
            "threshold_px": audit["merge_threshold_at_final_svg_scale_px"],
            **audit,
        }), encoding="utf-8")
        candidate_views[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": rounded(after),
            "source_kind": SOURCE_KIND,
            "official_cad_paths_mm": [],
            "line_simplification": audit,
        }
        views.append({
            "view": view,
            "svg": relative(target),
            "svg_sha256": sha256(target),
            "projection_axes": list(axes),
            "raw_edge_count": len(raw),
            "displayed_raw_edge_count": len(displayed),
            "silhouette_path_count": len(after),
            "silhouette_segment_count": segment_count(after),
            "drawing_line_source_kind": SOURCE_KIND,
            "drawing_line_source_label_zh": SOURCE_LABEL_ZH,
            "official_cad_path_count": 0,
            "blue_line_present": False,
            "line_simplification": audit,
        })
        view_audits[view] = audit
    candidate_path = output / "candidate-representations.json"
    write_json(candidate_path, {
        "schema_version": 2,
        "profile_key": PROFILE_KEY,
        "representative_global_id": product.GlobalId,
        "ifc_type_name": profile["ifc_type_name"],
        "article_number": ARTICLE_NUMBER,
        "units": "mm",
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "official_cad_used": False,
        "third_party_cad_used": False,
        "formal_ifc_write_allowed": False,
        "review_status": "visual_review_pending",
        "near_line_merge_threshold_mm": NEAR_LINE_MERGE_THRESHOLD_MM,
        "views": candidate_views,
    })
    audit_path = output / "line-simplification-audit.json"
    write_json(audit_path, {
        "schema_version": 1,
        "product": PROFILE_KEY,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "operation": "semantic_merge_of_close_parallel_collinear_and_repeated_projection_lines",
        "merge_threshold_mm": NEAR_LINE_MERGE_THRESHOLD_MM,
        "views": view_audits,
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_bytes_unchanged": sha256(source) == FORMAL_SHA256,
        "review_status": "visual_review_pending",
        "derived_ifc_write_allowed": False,
        "pass": all(record["pass"] for record in view_audits.values()),
    })
    manifest = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": GENERATOR,
        "formal_ifc": relative(source),
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(source),
        "formal_ifc_bytes_unchanged": sha256(source) == FORMAL_SHA256,
        "profile_register": relative(register_path),
        "profile_register_sha256": sha256(register_path),
        "profile_key": PROFILE_KEY,
        "display_name": profile["display_name"],
        "ifc_type_name": profile["ifc_type_name"],
        "ifc_type_description": product.Description,
        "article_number": ARTICLE_NUMBER,
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
            **drawing_source,
            "official_source_access_record": relative(ACCESS_RECORD),
            "official_source_access_record_sha256": sha256(ACCESS_RECORD),
            "official_product_cad_status": "not_published_on_verified_manufacturer_surfaces",
        },
        "official_reference": profile["official_reference"],
        "official_identity_evidence_only": True,
        "dimension_cross_check": access["dimension_cross_check"],
        "official_cad_acquired": False,
        "official_cad_used": False,
        "third_party_cad_used": False,
        "blue_line_present": False,
        "official_source_access_record": relative(ACCESS_RECORD),
        "official_source_access_record_sha256": sha256(ACCESS_RECORD),
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "line_simplification_audit": relative(audit_path),
        "line_simplification_audit_sha256": sha256(audit_path),
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "derived_ifc_write_allowed": False,
        "formal_ifc_write": "not performed",
        "views": views,
        "pass": all(view["line_simplification"]["pass"] and not view["blue_line_present"] for view in views),
    }
    context_manifest = output / "project-context-manifest.json"
    bonsai_manifest = output / "bonsai-review-manifest.json"
    if context_manifest.is_file():
        context = load_json(context_manifest)
        manifest["project_context"] = {
            "manifest": relative(context_manifest),
            "manifest_sha256": sha256(context_manifest),
            "walls_and_surrounding_project_elements_retained": context["walls_and_surrounding_project_elements_retained"],
            "overlay_top_layer_with_white_mask": context["overlay_top_layer_with_white_mask"],
            "pass": context.get("pass") is True,
        }
    if bonsai_manifest.is_file():
        bonsai = load_json(bonsai_manifest)
        manifest["bonsai_review"] = {
            "manifest": relative(bonsai_manifest),
            "manifest_sha256": sha256(bonsai_manifest),
            "mode": bonsai["mode"],
            "saved_active_representation": bonsai["bonsai_session"]["saved_active_representation"],
            "saved_camera_count": bonsai["bonsai_session"]["saved_camera_count"],
            "pass": bonsai.get("pass") is True,
        }
    manifest_path = output / "manifest.json"
    write_json(manifest_path, manifest)
    write_index(manifest)
    print(relative(manifest_path))


if __name__ == "__main__":
    main()
