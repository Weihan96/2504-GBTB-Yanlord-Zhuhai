#!/usr/bin/env python3
"""Overlay the de-textured official-outline review simplification on drawings."""

from falper_sorgente_linework import ROOT
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.GLOBAL_ID = "2iKOL78$H0N9Yd9$ky3pW4"
shared.OVERLAY_ID_PREFIX = "gessi316-54294"
shared.GENERATOR = "pipeline/scripts/render_gessi316_54294_project_context.py"
shared.SOURCE_KIND = "native_dwg_review_simplification"
shared.SOURCE_LABEL_ZH = "基于官方 Gessi 54294 原生 DWG 轮廓的去纹审核简化表达"
shared.SOURCE_DWG_SHA256 = "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4"
shared.OFFICIAL_CAD_USED = True
shared.PATH_KEY = "review_simplified_official_outline_paths_mm"
shared.CLOSE_PATHS = False
shared.ALIGNMENT_MODE = "centre"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.60
shared.REVIEW_HIDDEN_CSS = ".official-elevation-anchor{display:none}"
shared.PLAN_SOURCE = ROOT / "drawings/Sanitary Plan.svg"
shared.FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-06-20-R12-NY.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-06-22-R12-NX.svg"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-sanitary-plan.svg", (3.0, 3.0), 1),
    ("front", "side", shared.FRONT_SOURCE, "project-context-front-elevation.svg", (2.0, 2.0), 0),
    ("side", "front", shared.SIDE_SOURCE, "project-context-side-elevation.svg", (2.0, 2.0), 0),
)


def add_review_simplified_overlay(source, target, view, path, overlay_bbox):
    content = source.read_text(encoding="utf-8")
    minimum, maximum = overlay_bbox
    padding = 0.18
    mask_x = minimum[0] - padding
    mask_y = minimum[1] - padding
    mask_width = maximum[0] - minimum[0] + 2 * padding
    mask_height = maximum[1] - minimum[1] + 2 * padding
    group = f'''<g id="{shared.OVERLAY_ID_PREFIX}-{view}-{shared.GLOBAL_ID}" class="review-simplified-reference official-outline-derived project-context-overlay" data-source-kind="{shared.SOURCE_KIND}" data-source-label-zh="{shared.SOURCE_LABEL_ZH}" data-source-dwg-sha256="{shared.SOURCE_DWG_SHA256}" data-official-cad-used="true" data-unaltered-official-dwg="false" data-ifc-guid="{shared.GLOBAL_ID}" fill="none" stroke-linecap="round" stroke-linejoin="round">
  <rect class="official-reference-envelope-mask" x="{mask_x:.6f}" y="{mask_y:.6f}" width="{mask_width:.6f}" height="{mask_height:.6f}" fill="#ffffff" stroke="none"/>
  <path class="official-reference-mask" d="{path}" stroke="#ffffff" stroke-width="0.15"/>
  <path class="review-simplified-reference official-outline-derived" d="{path}" stroke="{shared.BLUE}" stroke-width="0.058"/>
</g>'''
    if "</svg>" not in content:
        raise RuntimeError(f"invalid project SVG: {source}")
    content = content.rsplit("</svg>", 1)[0] + group + "\n</svg>\n"
    target.write_text(shared.rebase_external_resources(content, source, target), encoding="utf-8")


shared.add_overlay = add_review_simplified_overlay


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = shared.load_json(manifest_path)
    manifest["unaltered_official_cad_used_as_review_representation"] = False
    manifest["original_official_cad_evidence_preserved"] = True
    manifest["review_texture_detail_path_count"] = 0
    shared.write_json(manifest_path, manifest)
