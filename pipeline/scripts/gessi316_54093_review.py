#!/usr/bin/env python3
"""Generate Gessi316 54093 review views from the exact official G000 native DWG."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "gessi316-54093"
shared.ARTICLE_NUMBER = "54093"
shared.DRAWING_PRODUCT_CODE = "54093"
shared.GENERATOR = "pipeline/scripts/gessi316_54093_review.py"
shared.OFFICIAL_CAD_STATUS = "public_official_api_exact_54093_g000_native_dwg_acquired"
shared.OFFICIAL_CAD_ACQUIRED = True
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = True
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/gessi316-54093"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.LINEWORK = shared.OUTPUT_DIR / "official-native-dwg-linework.json"
shared.PDF_VERIFICATION = shared.OUTPUT_DIR / "official-source/official-pdf-verification.json"
shared.REVALIDATION = shared.OUTPUT_DIR / "official-source/official-source-revalidation.json"
shared.SOURCE_DWG = shared.OUTPUT_DIR / "official-source/GPF5409300000G000_3.dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54093 G000 原生 DWG 图纸表达"
shared.SOURCE_LABEL_EN = "drawing representation from the exact Gessi 54093 G000 official native DWG"
shared.SOURCE_DWG_SHA256 = "a9308fc8498d34c8cf2f68fd28aaf90b11e62bfc59424e0a5a3ac7f0d28d5047"
shared.EXPECTED_DESCRIPTION = "High spout basin"
shared.EXPECTED_PATH_COUNTS = {"plan": 7, "front": 914, "side": 21}
shared.VIEW_TOLERANCES_MM = {"plan": 2.0, "front": 2.0, "side": 2.0}
shared.SCOPE = "official Gessi exact 54093 G000 family reference; not a project shop drawing"
shared.IDENTITY_POLICY_KEY = "g001_or_a004_variant_cad_used"
shared.PUBLIC_ACCESS_KEY = "exact_54093_g000_native_dwg_publicly_downloadable"
shared.CANDIDATE_IDENTITY_FLAG = "g001_or_a004_variant_cad_used"


def align_official_paths(view, paths, minimum, maximum):
    center_x = (minimum[0] + maximum[0]) / 2.0
    official_bounds = shared.path_bounds(paths)
    if view == "plan":
        # Native sheet Plan points toward -Y; Side points toward +Y.
        # The semicircle (path 4) joins the vertical stem at its installation
        # axis. Register that semantic axis to Side's 25 mm base radius.
        source_axis_y = min(point[1] for point in paths[4])
        target_axis_y = minimum[1] + 25.0
        return [[[center_x + x, target_axis_y + source_axis_y - y] for x, y in path] for path in paths]
    top_offset = maximum[2] - official_bounds["maximum"][1]
    if view == "front":
        return [[[center_x + x, top_offset + z] for x, z in path] for path in paths]
    return [[[minimum[1] + y, top_offset + z] for y, z in path] for path in paths]


def compatibility(view, official, minimum, maximum):
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


def write_index(manifest: dict) -> None:
    output = shared.OUTPUT_DIR
    cards = "".join(f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>' for item in manifest["views"])
    bonsai = "".join(f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{name}.png"><img src="bonsai-camera-{name}.png"></a></article>' for view, name in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso")))
    (output / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Gessi316 54093 G000 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Gessi316 Meccanica / 54093 G000 high counter basin spout</h1><p>Grey = actual IFC Body; black = simplified proxy; blue = exact official Gessi 54093 G000 native DWG. G001/A004 variants and adjacent products are excluded. The project contains one native sanitary-plan context for this GlobalId; no project elevation is invented.</p><nav><a href="review-contact-sheet.png">Review contact sheet</a><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-native-dwg-linework.json">Native DWG linework</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/official-pdf-verification.json">PDF verification</a><a href="project-context-sanitary-plan.svg">Project sanitary plan</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://areapro.gessi.com/en/product/54093">Official product</a></nav><h2>Project drawing context</h2><main><article><h2>Sanitary plan</h2><a href="project-context-sanitary-plan-review.svg"><img src="project-context-sanitary-plan-review-preview.png"></a></article></main><h2>Three-view native DWG / IFC comparison</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.align_official_paths = align_official_paths
shared.compatibility = compatibility
shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
