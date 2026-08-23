#!/usr/bin/env python3
"""Generate the unresolved electric flue check-valve proxy review."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "electric-flue-check-valve"
shared.ARTICLE_NUMBER = "unresolved"
shared.GENERATOR = "pipeline/scripts/electric_flue_check_valve_review.py"
shared.OFFICIAL_CAD_STATUS = "manufacturer_and_model_unresolved_no_exact_official_CAD_match"
shared.OFFICIAL_CAD_ACQUIRED = False
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
shared.IDENTITY_EVIDENCE_ONLY = False
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/electric-flue-check-valve"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = None


original_render_svg = shared.render_svg


def render_svg(profile, view, raw_edges, proxy, metadata):
    content = original_render_svg(profile, view, raw_edges, proxy, metadata)
    return (
        content
        .replace(
            "Blue = none (exact project configuration not matched)",
            "Blue = none (manufacturer/model unresolved; no official CAD match)",
        )
        .replace("Official identity: unresolved", "Manufacturer/model: unresolved")
        .replace("Official family CAD acquired", "Manufacturer CAD acquired")
        .replace(
            "Official catalogue is identity evidence only.",
            "Project IFC does not identify a manufacturer or model.",
        )
    )


def write_index(manifest: dict) -> None:
    output = shared.OUTPUT_DIR
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (output / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Electric flue check valve review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Electric flue check valve / unresolved project proxy</h1><p>Black line = {shared.SOURCE_LABEL_EN}. No blue line is shown because the project IFC supplies no manufacturer or model and no exact official CAD can be attributed. The Novy 906271 and patent CN209977425U are recorded as rejected type/shape references, not product identity.</p><nav><a href="review-contact-sheet.png">Review contact sheet</a><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-r04-px-elevation.svg">R04 +X elevation</a><a href="project-context-r04-ny-elevation.svg">R04 -Y elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a></nav><h2>Project drawing context</h2><main><article><h2>R04 +X elevation</h2><a href="project-context-r04-px-elevation-review.svg"><img src="project-context-r04-px-elevation-review-preview.png"></a></article><article><h2>R04 -Y elevation</h2><a href="project-context-r04-ny-elevation-review.svg"><img src="project-context-r04-ny-elevation-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from one isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.render_svg = render_svg
shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
