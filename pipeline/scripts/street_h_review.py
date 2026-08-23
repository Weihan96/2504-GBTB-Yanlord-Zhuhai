#!/usr/bin/env python3
"""Generate the STREET-H review without substituting the parent Street DXF."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "street-h"
shared.ARTICLE_NUMBER = "project STREET-H component / parent family STREET240 + STREET4054"
shared.GENERATOR = "pipeline/scripts/street_h_review.py"
shared.OFFICIAL_CAD_STATUS = "official_complete_parent_top_dxf_acquired_STREET_H_component_match_rejected"
shared.OFFICIAL_CAD_ACQUIRED = True
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/street-h"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "Sink holder"


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png" alt="actual Bonsai IFC Body {view} camera render"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (shared.OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>STREET-H review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>antoniolupi Street sink-holder subcomponent / project STREET-H</h1><p>Black line = {shared.SOURCE_LABEL_EN}. The official native DXF for the complete STREET240 + STREET4054 parent top is archived and mechanically verified, but it is not the same geometry as the repeated 300 × 150 × 100 mm STREET-H subcomponent. No blue parent-top line is substituted into the project representation.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-native-dxf-linework.json">Rejected parent-family DXF linework</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/AL_Street.dxf">Official native DXF</a><a href="official-source/ANTONIOLUPI-official-Street-technical.pdf">Official technical PDF</a><a href="project-context-sanitary-plan.svg">Project sanitary plan</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.antoniolupi.it/en/products/sinks/street">Official product page</a></nav><h2>Project drawing context</h2><main><article><h2>Sanitary plan</h2><a href="project-context-sanitary-plan-review.svg"><img src="project-context-sanitary-plan-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from one isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
