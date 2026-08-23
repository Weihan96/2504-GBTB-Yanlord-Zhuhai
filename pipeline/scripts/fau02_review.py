#!/usr/bin/env python3
"""Generate the FAU02 review with GH2 CAD archived but not misused."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "fau02"
shared.ARTICLE_NUMBER = "GH2 nearest candidate / exact FAU02 article unresolved"
shared.GENERATOR = "pipeline/scripts/fau02_review.py"
shared.OFFICIAL_CAD_STATUS = "official_GH2_acquired_nearest_candidate_exact_project_mismatch"
shared.OFFICIAL_CAD_ACQUIRED = True
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/fau02"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.EXPECTED_DESCRIPTION = "Falper SORGENTE Faucet"


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (shared.OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Falper FAU02 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Falper floor-mounted basin spout / project FAU02</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Falper Cilindro GH2 native 2D/technical DWG, 3D DWG and PDF are archived as the nearest family candidate only. Independent height, stem and flange dimensions fail the exact-project match gate, so no official blue line is shown.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/cad-match-analysis.json">CAD match analysis</a><a href="official-source/Falper-Cilindro-GH2-native-dwg-reference.svg">Archived GH2 native-DWG SVG</a><a href="project-context-furniture-plan.svg">Project plan</a><a href="project-context-front-elevation.svg">Project front</a><a href="project-context-side-elevation.svg">Project side</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://falper.it/rubinetteria-cilindro/">Official page</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture + FFL plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article><article><h2>Front elevation</h2><a href="project-context-front-elevation-review.svg"><img src="project-context-front-elevation-review-preview.png"></a></article><article><h2>Side elevation</h2><a href="project-context-side-elevation-review.svg"><img src="project-context-side-elevation-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
