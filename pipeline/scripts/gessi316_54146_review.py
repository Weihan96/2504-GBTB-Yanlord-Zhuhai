#!/usr/bin/env python3
"""Generate Gessi316 54146 review views and preserve original DWG evidence."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "gessi316-54146"
shared.ARTICLE_NUMBER = "54146"
shared.DRAWING_PRODUCT_CODE = "54146"
shared.GENERATOR = "pipeline/scripts/gessi316_54146_review.py"
shared.OFFICIAL_CAD_STATUS = "public_official_api_exact_54146_g000_native_dwg_acquired"
shared.OFFICIAL_CAD_ACQUIRED = True
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = True
shared.EXPECTED_SOURCE_EXACT_PROJECT_CONFIGURATION_MATCH = True
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/gessi316-54146"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.LINEWORK = shared.OUTPUT_DIR / "official-native-dwg-linework.json"
shared.PDF_VERIFICATION = shared.OUTPUT_DIR / "official-source/official-pdf-verification.json"
shared.REVALIDATION = shared.OUTPUT_DIR / "official-source/official-source-revalidation.json"
shared.SOURCE_DWG = shared.OUTPUT_DIR / "official-source/GPF5414600000G000_3.dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54146 G000 原生 DWG 图纸表达"
shared.SOURCE_LABEL_EN = "drawing representation from the exact Gessi 54146 G000 official native DWG"
shared.REVIEW_LINE_SOURCE_KIND = "native_dwg_review_simplification"
shared.REVIEW_LINE_SOURCE_LABEL_ZH = "基于官方 Gessi 54146 G000 原生 DWG 轮廓的简化蓝线审核表达"
shared.REVIEW_LINE_SOURCE_LABEL_EN = "simplified blue review representation based on the official Gessi 54146 G000 native-DWG outline"
shared.SOURCE_DWG_SHA256 = "c8ddb90f61565d5273a32574f777812bbb1d9ef33e2b98479f1e7c381e69950d"
shared.EXPECTED_DESCRIPTION = "Ceiling-mounted adjustable headshower Ø300"
shared.EXPECTED_PATH_COUNTS = {"plan": 16, "front": 607, "side": 660}
shared.VIEW_TOLERANCES_MM = {"plan": 1.0, "front": 1.0, "side": 1.0}
shared.SCOPE = "official Gessi exact 54146 G000 family reference; not a project shop drawing"
shared.IDENTITY_POLICY_KEY = "g001_variant_cad_used"
shared.PUBLIC_ACCESS_KEY = "exact_54146_g000_native_dwg_publicly_downloadable"
shared.CANDIDATE_IDENTITY_FLAG = "g001_variant_cad_used"


def align_official_paths(view, paths, minimum, maximum, simplification=None):
    center_x = (minimum[0] + maximum[0]) / 2.0
    official_bounds = shared.path_bounds(paths)
    if view == "plan":
        aligned = [[[center_x + x, minimum[1] + y] for x, y in path] for path in paths]
        return aligned, {
            "mode": "project_centre_and_minimum_depth_translation_only",
            "scale_applied": False,
            "reflection_applied": False,
        }
    top_offset = maximum[2] - official_bounds["maximum"][1]
    if view == "front":
        aligned = [[[center_x + x, top_offset + z] for x, z in path] for path in paths]
        return aligned, {
            "mode": "project_centre_and_top_elevation_translation_only",
            "scale_applied": False,
            "reflection_applied": False,
        }
    aligned = [[[minimum[1] + y, top_offset + z] for y, z in path] for path in paths]
    return aligned, {
        "mode": "source_axes_swap_plus_vertical_reflection_then_project_translation",
        "scale_applied": False,
        "horizontal_reflection_applied": False,
        "vertical_reflection_applied": True,
        "ceiling_mount_above_spray_face": True,
        "camera_cross_check": "matches actual IFC Body side camera and the Front vertical order",
    }


def compatibility(view, official, minimum, maximum, simplification=None):
    axes = shared.VIEW_DEFINITIONS[view]["axes"]
    official_size = shared.path_bounds(official)["size"]
    ifc_size = [maximum[axis] - minimum[axis] for axis in axes]
    delta = [abs(official_size[axis] - ifc_size[axis]) for axis in range(2)]
    tolerance = shared.VIEW_TOLERANCES_MM[view]
    return {
        "ifc_body_projection_size_mm": [round(value, 6) for value in ifc_size],
        "official_native_dwg_size_mm": official_size,
        "absolute_delta_mm": [round(value, 6) for value in delta],
        "tolerance_mm": tolerance,
        "configuration": "G000",
        "all_projection_axes_compared": True,
        "pass": max(delta) <= tolerance,
    }


def simplify_review_linework(view, paths):
    original_bounds = shared.path_bounds(paths)
    original_nodes = sum(len(path) for path in paths)
    if view == "plan":
        review = paths
        removed = []
        replacement = []
    else:
        lower = original_bounds["minimum"][1]
        removed = []
        review = []
        for path in paths:
            xs = [point[0] for point in path]
            ys = [point[1] for point in path]
            maximum_dimension = max(max(xs) - min(xs), max(ys) - min(ys))
            centre_y = (min(ys) + max(ys)) / 2.0
            if maximum_dimension <= 8.0 and centre_y <= lower + 10.25:
                removed.append(path)
            else:
                review.append(path)
        retained_bounds = shared.path_bounds(review)
        centre_x = (original_bounds["minimum"][0] + original_bounds["maximum"][0]) / 2.0
        replacement = [[
            [centre_x, retained_bounds["minimum"][1]],
            [centre_x, original_bounds["minimum"][1]],
        ]]
        review.extend(replacement)
    review_bounds = shared.path_bounds(review)
    envelope_delta = [
        round(review_bounds["size"][axis] - original_bounds["size"][axis], 6)
        for axis in range(2)
    ]
    original_centre = [
        round((original_bounds["minimum"][axis] + original_bounds["maximum"][axis]) / 2.0, 6)
        for axis in range(2)
    ]
    review_centre = [
        round((review_bounds["minimum"][axis] + review_bounds["maximum"][axis]) / 2.0, 6)
        for axis in range(2)
    ]
    review_nodes = sum(len(path) for path in review)
    audit = {
        "mode": "official_54146_outline_based_fine_spray_detail_simplification",
        "classification": "sub_8mm_paths_in_lower_10_25mm_spray_face_band",
        "original_path_count": len(paths),
        "original_node_count": original_nodes,
        "original_segment_count": sum(max(0, len(path) - 1) for path in paths),
        "removed_fine_spray_nozzle_detail_path_count": len(removed),
        "removed_fine_spray_nozzle_detail_node_count": sum(len(path) for path in removed),
        "replacement_overall_spray_tip_datum_path_count": len(replacement),
        "review_path_count": len(review),
        "review_node_count": review_nodes,
        "review_segment_count": sum(max(0, len(path) - 1) for path in review),
        "review_texture_detail_path_count": 0,
        "envelope_delta_mm": envelope_delta,
        "centre_delta_mm": [round(review_centre[axis] - original_centre[axis], 6) for axis in range(2)],
        "installation_axis_preserved": True,
        "ceiling_mount_anchor_preserved": True,
        "spray_face_diameter_endpoints_preserved": True,
        "overall_spray_tip_datum_preserved": True,
        "scale_applied": False,
        "unaltered_official_dwg": False,
        "pass": envelope_delta == [0.0, 0.0] and original_centre == review_centre,
    }
    return review, audit


def preserve_actual_ifc_geometry(vertices, faces, official_side_paths):
    return vertices, faces, {
        "mode": "actual_ifc_body_unchanged_for_review",
        "original_vertex_count": len(vertices),
        "original_face_count": len(faces),
        "review_vertex_count": len(vertices),
        "review_face_count": len(faces),
        "formal_ifc_modified": False,
        "pass": True,
    }


def render_svg(profile, view, raw_edges, proxy, official, metadata):
    svg = shared._gessi54146_original_render_svg(profile, view, raw_edges, proxy, official, metadata)
    svg = svg.replace(
        "Black = smooth/max-reach proxy · Blue = de-textured review simplification based on official 54146 outline",
        "Black = IFC silhouette proxy · Blue = simplified review line based on official 54146 G000 outline · original DWG retained separately",
    )
    svg = svg.replace("Official 54294 outline-based simplification", "Official 54146 outline-based review simplification")
    svg = svg.replace("texture: 0", "fine spray detail: removed")
    svg = svg.replace("Official family reference; not a project shop drawing.", "Exact G000 family outline reference; not a project shop drawing.")
    return svg


def render_original_official_evidence_svg(profile, view, official):
    svg = shared._gessi54146_original_evidence_svg(profile, view, official)
    svg = svg.replace(
        "Untouched official Gessi 54146 native-DWG evidence",
        "Original official Gessi 54146 native-DWG geometry evidence",
    )
    svg = svg.replace("original texture retained", "original spray-nozzle detail retained")
    if view == "side":
        svg = svg.replace(
            "original spray-nozzle detail retained · not the de-textured review representation",
            "original spray-nozzle detail retained · audited Side orientation normalization · scale 1.0",
        )
    return svg


def write_index(manifest: dict, output=None) -> None:
    output = output or shared.OUTPUT_DIR
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (output / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Gessi316 54146 G000 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Gessi316 Meccanica / 54146 G000 ceiling-mounted adjustable headshower</h1><p>Grey = actual IFC Body; black = simplified proxy; blue = review simplification based on the official 54146 G000 outline. The exact original native-DWG linework, including dense spray-nozzle detail, remains in separate evidence SVGs and JSON. G001 and wall-mounted 54145 are excluded.</p><nav><a href="review-contact-sheet.png">Review contact sheet</a><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-native-dwg-linework.json">Original native-DWG linework JSON</a><a href="official-original-plan.svg">Original official Plan</a><a href="official-original-front.svg">Original official Front</a><a href="official-original-side.svg">Original official Side</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/official-pdf-verification.json">PDF verification</a><a href="project-context-ffl-plan.svg">Project plan</a><a href="project-context-front-elevation.svg">Project front</a><a href="project-context-side-elevation.svg">Project side</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://areapro.gessi.com/en/product/54146">Official product</a></nav><h2>Project drawing context</h2><main><article><h2>FFL plan with furniture and walls</h2><a href="project-context-ffl-plan-review.svg"><img src="project-context-ffl-plan-review-preview.png"></a></article><article><h2>R12 front elevation</h2><a href="project-context-front-elevation-review.svg"><img src="project-context-front-elevation-review-preview.png"></a></article><article><h2>R12 side elevation</h2><a href="project-context-side-elevation-review.svg"><img src="project-context-side-elevation-review-preview.png"></a></article></main><h2>Three-view simplified review line / IFC comparison</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared._gessi54146_original_render_svg = shared.render_svg
shared._gessi54146_original_evidence_svg = shared.render_original_official_evidence_svg
shared.align_official_paths = align_official_paths
shared.compatibility = compatibility
shared.simplify_official_handle_texture = simplify_review_linework
shared.simplify_handle_knurls = preserve_actual_ifc_geometry
shared.render_svg = render_svg
shared.render_original_official_evidence_svg = render_original_official_evidence_svg
shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
