#!/usr/bin/env python3
"""Place the configured TRAP01 Body proxy in the complete sanitary plan."""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import ifcopenshell

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from int1_highpoly_type_review import mesh_for_one_product, projected_raw_edges
from render_hima01_project_context import bbox, path_points, svg_path, write_review_crop
from render_sis04_project_context import add_white_review_background, suppress_plan_elevation_markers, write_uncached_png_preview


PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
PLAN_SOURCE = ROOT / "drawings/Sanitary Plan.svg"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
GLOBAL_ID = "2Ak2ma0lvBEA49UpplzUqi"
SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
PROJECT_SCALE = 0.02
TRANSLATION = (71.64, 210.999997)
GENERATOR = "pipeline/scripts/render_trap01_project_context.py"
NUMBER = re.compile(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?")


def guid(element):
    return next((value for key, value in element.attrib.items() if key.endswith("guid")), None)


def project_group_segments():
    groups = []
    for element in ET.parse(PLAN_SOURCE).getroot().iter():
        if not element.tag.endswith("g") or guid(element) != GLOBAL_ID:
            continue
        segments = []
        for child in element.iter():
            if not child.tag.endswith("path"):
                continue
            points = path_points(child)
            segments.extend(zip(points, points[1:]))
        if segments:
            groups.append(segments)
    if len(groups) != 1:
        raise RuntimeError(f"expected one TRAP01 projection group, found {len(groups)}")
    return groups[0]


def transform(point):
    return TRANSLATION[0] + PROJECT_SCALE * point[0], TRANSLATION[1] - PROJECT_SCALE * point[1]


def mechanical_translation_match():
    model = ifcopenshell.open(FORMAL_IFC)
    _, vertices, faces = mesh_for_one_product(model, GLOBAL_ID)
    raw_edges = projected_raw_edges(vertices, faces, (0, 1))
    lookup = defaultdict(list)
    for start, end in raw_edges:
        first, second = transform(start), transform(end)
        delta = (round(second[0] - first[0], 6), round(second[1] - first[1], 6))
        lookup[delta].append((first, second))
        lookup[(-delta[0], -delta[1])].append((second, first))
    translations = Counter()
    segments = project_group_segments()
    for start, end in segments:
        delta = (round(end[0] - start[0], 6), round(end[1] - start[1], 6))
        for first, _ in lookup.get(delta, []):
            translations[(round(start[0] - (first[0] - TRANSLATION[0]), 6), round(start[1] - (first[1] - TRANSLATION[1]), 6))] += 1
    best_translation, count = translations.most_common(1)[0]
    if max(abs(best_translation[index] - TRANSLATION[index]) for index in range(2)) > 0.000002 or count < 600:
        raise RuntimeError(f"TRAP01 project translation match drifted: {best_translation} / {count}")
    points = [point for segment in segments for point in segment]
    return {
        "method": "exact_quantized_project_hidden-line_edge_vectors_matched_to_isolated_IFC_Body_XY_edges",
        "project_visible_segment_count": len(segments),
        "matching_edge_vote_count": count,
        "translation_svg_units": list(TRANSLATION),
        "scale_svg_units_per_mm": PROJECT_SCALE,
        "flip_projected_y": True,
        "original_visible_projection_bbox": [[round(value, 6) for value in item] for item in bbox(points)],
        "original_projection_is_occluded": True,
        "pass": True,
    }


def add_overlay(source: Path, target: Path, path: str):
    content = source.read_text(encoding="utf-8")
    pattern = re.compile(rf'<g(?=[^>]*\bifc:guid="{re.escape(GLOBAL_ID)}")')
    content, hidden_count = pattern.subn('<g style="display:none" data-review-replaced-by="configured-body-proxy"', content)
    if hidden_count != 1:
        raise RuntimeError(f"expected one TRAP01 project group, hid {hidden_count}")
    group = f'''<g id="trap01-plan-{GLOBAL_ID}" class="configured-simplified-proxy project-context-overlay" data-source-kind="{SOURCE_KIND}" data-source-label-zh="{SOURCE_LABEL_ZH}" data-official-cad-acquired="true" data-official-cad-used-as-representation="false" data-third-party-cad-used="false" data-ifc-guid="{GLOBAL_ID}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <path class="configured-proxy-mask" d="{path}" fill="#ffffff" fill-rule="evenodd" stroke="#ffffff" stroke-width="0.16"/>
  <path class="configured-proxy" d="{path}" stroke="#111820" stroke-width="0.058"/>
</g>'''
    target.write_text(content.rsplit("</svg>", 1)[0] + group + "\n</svg>\n", encoding="utf-8")


def main():
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    candidate = load_json(CANDIDATE)
    if (
        candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("official_cad_acquired") is not True
        or candidate.get("official_cad_used_as_representation") is not False
        or candidate.get("third_party_cad_used") is not False
    ):
        raise RuntimeError("TRAP01 configured candidate source gate failed")
    translation = mechanical_translation_match()
    paths = candidate["views"]["plan"]["proxy_paths_mm"]
    transformed_points = [transform(point) for path in paths for point in path]
    full_bbox = bbox(transformed_points)
    overlay_path = svg_path(paths, transform)
    target = PRODUCT_DIR / "project-context-sanitary-plan.svg"
    add_overlay(PLAN_SOURCE, target, overlay_path)
    crop = PRODUCT_DIR / "project-context-sanitary-plan-review.svg"
    write_review_crop(target, crop, full_bbox, 12.0, 12.0)
    suppress_plan_elevation_markers(crop)
    add_white_review_background(crop)
    preview = PRODUCT_DIR / "project-context-sanitary-plan-review-preview.png"
    write_uncached_png_preview(crop, preview)
    record = {
        "view": "plan",
        "candidate_view": "plan",
        "source": relative(PLAN_SOURCE),
        "output": relative(target),
        "output_sha256": sha256(target),
        "review_crop": relative(crop),
        "review_crop_sha256": sha256(crop),
        "review_preview": relative(preview),
        "review_preview_sha256": sha256(preview),
        "overlay": {
            "ifc_guid": GLOBAL_ID,
            "source_kind": SOURCE_KIND,
            "source_label_zh": SOURCE_LABEL_ZH,
            "official_cad_acquired": True,
            "official_cad_used_as_representation": False,
            "third_party_cad_used": False,
            "path_count": len(paths),
            "full_configured_proxy_bbox": [[round(value, 6) for value in item] for item in full_bbox],
            "fixed_project_scale_preserved": True,
            "geometry_scaled_or_stretched": False,
            "mechanical_translation": translation,
        },
    }
    manifest = {
        "schema_version": 1,
        "generator": GENERATOR,
        "source_kind": SOURCE_KIND,
        "source_label_zh": SOURCE_LABEL_ZH,
        "official_cad_used": False,
        "official_cad_acquired": True,
        "official_cad_used_as_representation": False,
        "blue_line_present": False,
        "third_party_cad_used": False,
        "project_context_retained": True,
        "walls_and_surrounding_project_elements_retained": True,
        "overlay_top_layer_with_white_mask": True,
        "blue_product_cad_line_present": False,
        "project_scale_svg_units_per_mm": PROJECT_SCALE,
        "semantic_view_mapping": {"plan": {"candidate_axes": [0, 1], "source": relative(PLAN_SOURCE), "project_direction": "+Z"}},
        "context_view_scope": {
            "included": ["plan"],
            "excluded": ["front", "side"],
            "reason": "TRAP01 appears in the sanitary plan, but no native project elevation SVG contains its GlobalId. Project elevations are not invented.",
        },
        "review_annotation_suppression": {"plan": "official-elevation-anchor groups only", "scope": "review crop only", "walls_furniture_and_ifc_geometry_removed": False},
        "views": [record],
        "pass": translation["pass"],
    }
    path = PRODUCT_DIR / "project-context-manifest.json"
    write_json(path, manifest)
    print(json.dumps({"manifest": relative(path), "translation": translation, "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
