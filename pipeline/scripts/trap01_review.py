#!/usr/bin/env python3
"""Build TRAP01 configured-project and exact official-family three-view review."""

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
from render_sis04_project_context import write_uncached_png_preview


PROFILE_KEY = "trap01"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OUTPUT_DIR = ROOT / "output/review/highpoly-types/trap01"
REGISTER = OUTPUT_DIR / "profile.json"
ACCESS_RECORD = OUTPUT_DIR / "official-source/source-access-record.json"
LINEWORK = OUTPUT_DIR / "official-native-dwg-linework.json"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
SOURCE_LABEL_EN = "simplified drawing representation derived from the original high-poly geometry"
CONFIGURATION_NOTE = "The proxy preserves the shortened project installation configuration; it does not force-fit the official default family extension."
ARTICLE = "151.116.11.1"
TYPE_DESCRIPTION = "Space Saving Dip Tube Trap"
SCOPE = "exact Geberit 151.116.11.1 adjustable family reference; project instance is a shortened installation configuration; not a project shop drawing"
BLUE = "#1677c8"
EXPECTED_OFFICIAL_PATH_COUNTS = {"plan": 92, "front": 77, "side": 117}
ADJUSTMENT_RANGES_MM = {"horizontal_extension": [0.0, 252.0], "h": [85.0, 334.0]}


def rounded(paths):
    return [[[round(float(x), 6), round(float(y), 6)] for x, y in path] for path in paths]


def bounds_2d(paths):
    points = [point for path in paths for point in path]
    return (
        (min(point[0] for point in points), min(point[1] for point in points)),
        (max(point[0] for point in points), max(point[1] for point in points)),
    )


def panel_transform(paths, x, y, width, height):
    minimum, maximum = bounds_2d(paths)
    span_x, span_y = maximum[0] - minimum[0], maximum[1] - minimum[1]
    padding = max(span_x, span_y) * 0.055
    minimum = (minimum[0] - padding, minimum[1] - padding)
    maximum = (maximum[0] + padding, maximum[1] + padding)
    scale = min(width / (maximum[0] - minimum[0]), height / (maximum[1] - minimum[1]))
    return lambda point: (
        x + (width - (maximum[0] - minimum[0]) * scale) / 2.0 + (point[0] - minimum[0]) * scale,
        y + height - (height - (maximum[1] - minimum[1]) * scale) / 2.0 - (point[1] - minimum[1]) * scale,
    )


def path_length(paths):
    return round(sum(
        ((second[0] - first[0]) ** 2 + (second[1] - first[1]) ** 2) ** 0.5
        for path in paths
        for first, second in zip(path, path[1:])
    ), 6)


def interpolate(first, second, ratio):
    return [
        first[axis] + (second[axis] - first[axis]) * ratio
        for axis in range(2)
    ]


def split_path_at_configured_limits(path, limits):
    """Split one official path into retained and clipped fragments without scaling."""
    fragments = {"solid": [], "clipped": []}
    if len(path) < 2:
        return fragments

    current_kind = None
    current = []
    for first, second in zip(path, path[1:]):
        ratios = [0.0, 1.0]
        for axis, limit in limits:
            delta = second[axis] - first[axis]
            if abs(delta) > 1e-9:
                ratio = (limit - first[axis]) / delta
                if 1e-9 < ratio < 1.0 - 1e-9:
                    ratios.append(ratio)
        ratios = sorted(set(round(value, 12) for value in ratios))
        for start_ratio, end_ratio in zip(ratios, ratios[1:]):
            start = interpolate(first, second, start_ratio)
            end = interpolate(first, second, end_ratio)
            middle = interpolate(first, second, (start_ratio + end_ratio) / 2.0)
            kind = "clipped" if any(middle[axis] > limit + 1e-6 for axis, limit in limits) else "solid"
            if kind != current_kind:
                if current_kind is not None and len(current) >= 2:
                    fragments[current_kind].append(current)
                current_kind = kind
                current = [start, end]
            else:
                if current[-1] != start:
                    current.append(start)
                current.append(end)
    if current_kind is not None and len(current) >= 2:
        fragments[current_kind].append(current)
    return fragments


def align_and_segment_official(view, proxy, official):
    """Align exact native DWG by fixed-body datums, then mark only excess extension dashed."""
    proxy_minimum, proxy_maximum = bounds_2d(proxy)
    official_minimum, official_maximum = bounds_2d(official)
    if view == "plan":
        translation = [
            proxy_minimum[0] - official_minimum[0],
            (proxy_minimum[1] + proxy_maximum[1] - official_minimum[1] - official_maximum[1]) / 2.0,
        ]
        anchor = "left outlet datum plus fixed-body centre axis"
    elif view == "front":
        translation = [
            proxy_minimum[0] - official_minimum[0],
            proxy_minimum[1] - official_minimum[1],
        ]
        anchor = "lower outlet datum in X/Z"
    else:
        translation = [
            (proxy_minimum[0] + proxy_maximum[0] - official_minimum[0] - official_maximum[0]) / 2.0,
            proxy_minimum[1] - official_minimum[1],
        ]
        anchor = "fixed-body centre axis plus lower outlet datum"
    aligned = [
        [[point[0] + translation[0], point[1] + translation[1]] for point in path]
        for path in official
    ]
    limits = []
    if view in {"plan", "front"}:
        limits.append((0, proxy_maximum[0]))
    if view in {"front", "side"}:
        limits.append((1, proxy_maximum[1]))
    retained, clipped = [], []
    for path in aligned:
        fragments = split_path_at_configured_limits(path, limits)
        retained.extend(fragments["solid"])
        clipped.extend(fragments["clipped"])
    if not retained or not clipped:
        raise RuntimeError(f"TRAP01 {view} aligned overlay did not produce both retained and clipped linework")
    aligned_minimum, aligned_maximum = bounds_2d(aligned)
    return aligned, retained, clipped, {
        "alignment_method": "translation_only_from_fixed_body_and_connection_axis_datums",
        "anchor": anchor,
        "translation_mm": [round(value, 6) for value in translation],
        "scale": 1.0,
        "reflection": False,
        "configured_bounds_mm": {
            "minimum": [round(value, 6) for value in proxy_minimum],
            "maximum": [round(value, 6) for value in proxy_maximum],
        },
        "aligned_official_bounds_mm": {
            "minimum": [round(value, 6) for value in aligned_minimum],
            "maximum": [round(value, 6) for value in aligned_maximum],
        },
        "independent_configured_limits": {
            "horizontal_x_max_mm": round(proxy_maximum[0], 6) if view in {"plan", "front"} else None,
            "vertical_z_max_mm": round(proxy_maximum[1], 6) if view in {"front", "side"} else None,
        },
        "clipped_extension_mm": {
            "horizontal": round(max(0.0, aligned_maximum[0] - proxy_maximum[0]), 6) if view in {"plan", "front"} else 0.0,
            "vertical": round(max(0.0, aligned_maximum[1] - proxy_maximum[1]), 6) if view in {"front", "side"} else 0.0,
        },
        "source_path_count": len(official),
        "retained_solid_fragment_count": len(retained),
        "clipped_dashed_fragment_count": len(clipped),
        "retained_solid_length_mm": path_length(retained),
        "clipped_dashed_length_mm": path_length(clipped),
        "horizontal_and_vertical_classified_independently": True,
        "geometry_scaled_or_stretched": False,
        "pass": True,
    }


def render_svg(profile, view, raw_edges, proxy, aligned, retained, clipped, comparison, overlay_audit):
    width, height = 1600, 980
    panel = (50, 185, 1500, 700)
    transform = panel_transform(proxy + aligned, *panel)
    edge_path = svg_path([[start, end] for start, end in raw_edges], transform)
    proxy_path = svg_path(proxy, transform, close=True)
    retained_path = svg_path(retained, transform)
    clipped_path = svg_path(clipped, transform)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="1600" height="980" fill="#fbfaf7"/>
<text x="50" y="54" font-family="Arial,sans-serif" font-size="29" font-weight="700" fill="#1f2d3d">{html.escape(profile["display_name"])}</text>
<text x="50" y="92" font-family="Arial,sans-serif" font-size="20" fill="#41566d">{VIEWS[view]["label"]} · isolated representative {profile["representative_global_id"]}</text>
<text x="50" y="125" font-family="Arial,sans-serif" font-size="16" fill="#68798a">Official 151.116.11.1 native DWG translated onto the current shortened project outline; scale remains exactly 1:1</text>
<text x="50" y="158" font-family="Arial,sans-serif" font-size="16" fill="#111820">Grey = actual IFC Body · Black = current configured proxy · Blue solid = retained official range · Blue dashed = official adjustable excess clipped by this project configuration</text>
<rect x="{panel[0]}" y="{panel[1]}" width="{panel[2]}" height="{panel[3]}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<text x="70" y="220" font-family="Arial,sans-serif" font-size="17" font-weight="700" fill="#111820">ALIGNED CONFIGURATION OVERLAY · OFFICIAL DWG {comparison["native_dwg_code"]}</text>
<path class="original-highpoly" d="{edge_path}" fill="none" stroke="#87929c" stroke-width="0.5" stroke-opacity="0.28" vector-effect="non-scaling-stroke"/>
<path class="official-family-native-dwg clipped-adjustable-extension" data-source-kind="native_dwg" data-source-sha256="{comparison["source_dwg_sha256"]}" data-used-as-project-representation="false" data-installation-status="not-project-actual" d="{clipped_path}" fill="none" stroke="{BLUE}" stroke-width="3" stroke-dasharray="12 8" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-family-native-dwg retained-configured-range" data-source-kind="native_dwg" data-source-sha256="{comparison["source_dwg_sha256"]}" data-used-as-project-representation="false" d="{retained_path}" fill="none" stroke="{BLUE}" stroke-width="4.5" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="configured-simplified-proxy geometry-derived" data-source-kind="{SOURCE_KIND}" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="2.2" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="50" y="930" font-family="Arial,sans-serif" font-size="15" fill="#41566d">Project: {comparison["project_configured_size_mm"]} mm · Official default: {comparison["official_default_size_mm"]} mm</text>
<text x="50" y="955" font-family="Arial,sans-serif" font-size="15" fill="#68798a">Anchor: {html.escape(overlay_audit["anchor"])} · clipped excess H {overlay_audit["clipped_extension_mm"]["horizontal"]} mm / V {overlay_audit["clipped_extension_mm"]["vertical"]} mm · dashed blue is not an installed project line.</text>
</svg>'''


def write_index(manifest):
    cards = "".join(f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>' for item in manifest["views"])
    bonsai = "".join(f'<article><h2>Bonsai {label}</h2><img src="bonsai-camera-{name}.png"></article>' for label, name in (("Plan", "plan"), ("Front", "front-elevation"), ("Side", "side-elevation"), ("ISO", "iso")))
    (OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>TRAP01 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:16px}}</style><h1>Geberit 151.116.11.1 / TRAP01</h1><p>Black = current shortened project proxy. Official G/A/L/P DWGs are translated 1:1 onto the fixed-body datums: retained official range is solid blue; the published adjustable excess beyond the current horizontal or vertical endpoint is dashed blue and is not an installed project line.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="adjustable-overlay-audit.json">Adjustable overlay audit</a><a href="official-native-dwg-linework.json">Official DWG linework</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-sanitary-plan.svg">Project sanitary plan</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://catalog.geberit.us/en-US/product/PRO_185224">Official page</a></nav><h2>Configured candidate / aligned official family CAD</h2><main>{cards}</main><h2>Project context</h2><main><article><h2>Sanitary Plan</h2><img src="project-context-sanitary-plan-review-preview.png"></article></main><h2>Actual Bonsai camera renders</h2><main>{bonsai}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
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
    access, linework = load_json(ACCESS_RECORD), load_json(LINEWORK)
    drawing_source = profile["drawing_source"]
    if (
        drawing_source.get("source_kind") != SOURCE_KIND
        or drawing_source.get("official_cad_acquired") is not True
        or drawing_source.get("official_cad_used_as_representation") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or access.get("resolved_article") != ARTICLE
        or access.get("native_cad_selection", {}).get("pass") is not True
        or access.get("configuration_cross_check", {}).get("official_default_family_paths_used_as_project_representation") is not False
        or access.get("drawing_geometry_source", {}).get("source_kind") != SOURCE_KIND
        or linework.get("configuration_cross_check", {}).get("pass") is not True
        or linework.get("pass") is not True
    ):
        raise RuntimeError("TRAP01 configured-source gate failed")
    model = ifcopenshell.open(source)
    product, vertices, faces = mesh_for_one_product(model, profile["representative_global_id"])
    if product_type_name(product) != profile["ifc_type_name"]:
        raise RuntimeError("TRAP01 representative type identity drifted")
    product_type = next(relation.RelatingType for relation in product.IsTypedBy)
    if product_type.Description != TYPE_DESCRIPTION:
        raise RuntimeError("TRAP01 IFC type description drifted")
    instances = sorted(item.GlobalId for item in model.by_type(product.is_a()) if product_type_name(item) == profile["ifc_type_name"])
    if instances != profile["expected_instance_global_ids"]:
        raise RuntimeError("TRAP01 instance set drifted")
    minimum, maximum = bounds_3d(vertices)
    size = [round(maximum[index] - minimum[index], 6) for index in range(3)]
    if any(abs(size[index] - access["article_resolution"]["project_ifc_body_local_xyz_mm"][index]) > 0.001 for index in range(3)):
        raise RuntimeError("TRAP01 configured Body dimensions drifted")
    output.mkdir(parents=True, exist_ok=True)
    candidate_views, views, overlay_views = {}, [], {}
    configured_sizes = {"plan": [size[0], size[1]], "front": [size[0], size[2]], "side": [size[1], size[2]]}
    for view, definition in VIEWS.items():
        axes = definition["axes"]
        all_edges = projected_raw_edges(vertices, faces, axes)
        edges = display_edge_sample(all_edges)
        proxy = projected_silhouette(vertices, faces, axes, float(profile["silhouette_simplify_mm"]))
        official = linework["views"][view]["paths_mm"]
        if len(official) != EXPECTED_OFFICIAL_PATH_COUNTS[view]:
            raise RuntimeError(f"TRAP01 official {view} path count drifted")
        aligned, retained, clipped, overlay_audit = align_and_segment_official(view, proxy, official)
        overlay_views[view] = overlay_audit
        comparison = {
            "native_dwg_code": linework["views"][view]["native_dwg_code"],
            "source_dwg_sha256": linework["views"][view]["source_dwg_sha256"],
            "project_configured_size_mm": configured_sizes[view],
            "official_default_size_mm": linework["views"][view]["bounds_mm"]["size"],
            "geometry_scaled_or_stretched": False,
            "official_default_family_paths_used_as_project_representation": False,
        }
        target = output / f"{view}.svg"
        target.write_text(render_svg(profile, view, edges, proxy, aligned, retained, clipped, comparison, overlay_audit), encoding="utf-8")
        preview = output / f"{view}-preview.png"
        write_uncached_png_preview(target, preview)
        candidate_views[view] = {
            "projection_axes": list(axes),
            "proxy_paths_mm": rounded(proxy),
            "source_kind": SOURCE_KIND,
            "official_cad_paths_mm": [],
            "official_family_native_dwg_paths_mm": official,
            "adjustable_overlay": overlay_audit,
            "official_family_native_dwg_code": comparison["native_dwg_code"],
            "official_family_native_dwg_sha256": comparison["source_dwg_sha256"],
            "official_family_paths_used_as_project_representation": False,
            "geometry_scaled_or_stretched": False,
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
            "candidate_drawing_line_source_kind": SOURCE_KIND,
            "official_family_native_dwg_path_count": len(official),
            "official_family_retained_solid_fragment_count": len(retained),
            "official_family_clipped_dashed_fragment_count": len(clipped),
            "blue_line_present_as_separate_family_reference": True,
            "blue_line_present": True,
            "blue_line_used_as_project_representation": False,
            "comparison": comparison,
        })
    overlay_audit_path = output / "adjustable-overlay-audit.json"
    write_json(overlay_audit_path, {
        "schema_version": 1,
        "profile_key": PROFILE_KEY,
        "article_number": ARTICLE,
        "method": "translate official fixed-body and connection-axis datums onto the current black proxy; split horizontal and vertical excess independently",
        "line_semantics": {
            "black": "current shortened project configuration",
            "blue_solid": "official native DWG linework retained within the current configured horizontal and vertical limits",
            "blue_dashed": "official adjustable-family default/extended linework beyond the current configured endpoint; not installed project geometry and not an exact project shop drawing",
        },
        "official_adjustment_ranges_mm": ADJUSTMENT_RANGES_MM,
        "views": overlay_views,
        "geometry_scaled_or_stretched": False,
        "formal_ifc_write": "not performed",
        "pass": all(item["pass"] for item in overlay_views.values()),
    })
    candidate_path = output / "candidate-representations.json"
    write_json(candidate_path, {
        "schema_version": 1,
        "profile_key": PROFILE_KEY,
        "representative_global_id": product.GlobalId,
        "ifc_type_name": profile["ifc_type_name"],
        "article_number": ARTICLE,
        "units": "mm",
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "configuration_note": CONFIGURATION_NOTE,
        "official_cad_acquired": True,
        "official_cad_used": False,
        "official_cad_used_as_representation": False,
        "official_cad_family_reference_archived": True,
        "third_party_cad_used": False,
        "formal_ifc_write_allowed": False,
        "blue_line_present": True,
        "blue_line_role": "1:1 official adjustable-family review overlay only; dashed excess is not project actual geometry",
        "review_status": "visual_review_pending",
        "scope": SCOPE,
        "views": candidate_views,
    })
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/trap01_review.py",
        "formal_ifc": relative(source),
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(source),
        "formal_ifc_bytes_unchanged": sha256(source) == FORMAL_SHA256,
        "profile_key": PROFILE_KEY,
        "display_name": profile["display_name"],
        "ifc_type_name": profile["ifc_type_name"],
        "ifc_type_description": product_type.Description,
        "representative_global_id": product.GlobalId,
        "registered_instance_global_ids": instances,
        "geometry_product_count": 1,
        "whole_model_render": False,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "source_label_en": SOURCE_LABEL_EN,
        "configuration_note": CONFIGURATION_NOTE,
        "official_cad_acquired": True,
        "official_cad_used": False,
        "official_cad_used_as_representation": False,
        "third_party_cad_used": False,
        "blue_line_present": True,
        "mesh_vertex_count": len(vertices),
        "mesh_face_count": len(faces),
        "bounds_mm": {"minimum": [round(value, 6) for value in minimum], "maximum": [round(value, 6) for value in maximum], "size": size},
        "drawing_source": {**drawing_source, "official_source_access_record": relative(ACCESS_RECORD), "official_source_access_record_sha256": sha256(ACCESS_RECORD), "official_native_dwg_linework": relative(LINEWORK), "official_native_dwg_linework_sha256": sha256(LINEWORK)},
        "article_resolution": access["article_resolution"],
        "configuration_cross_check": access["configuration_cross_check"],
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "adjustable_overlay_audit": relative(overlay_audit_path),
        "adjustable_overlay_audit_sha256": sha256(overlay_audit_path),
        "views": views,
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "derived_ifc_write_allowed": False,
        "formal_ifc_write": "not performed",
        "scope": SCOPE,
        "pass": sha256(source) == FORMAL_SHA256,
    }
    for key, filename in (("project_context", "project-context-manifest.json"), ("bonsai_review", "bonsai-review-manifest.json")):
        path = output / filename
        if path.is_file():
            payload = load_json(path)
            manifest[key] = {"manifest": relative(path), "manifest_sha256": sha256(path), "mode": payload.get("mode"), "pass": payload.get("pass")}
            manifest["pass"] = manifest["pass"] and payload.get("pass") is True
    manifest_path = output / "manifest.json"
    write_json(manifest_path, manifest)
    write_index(manifest)
    print(manifest_path)


if __name__ == "__main__":
    main()
