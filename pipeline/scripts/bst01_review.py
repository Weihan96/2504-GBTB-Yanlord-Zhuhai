#!/usr/bin/env python3
"""Generate the Baxter Ninfea / BST01 review without claiming official CAD."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared

shared.PROFILE_KEY = "bst01"
shared.ARTICLE_NUMBER = "Baxter Ninfea bedside table opening side unresolved / BST01"
shared.GENERATOR = "pipeline/scripts/bst01_review.py"
shared.OFFICIAL_CAD_STATUS = "manufacturer_login_required_not_acquired_opening_side_unresolved"
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/bst01"
shared.REGISTER = ROOT / "pipeline/decisions/highpoly-drawing-profile-register.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "Nifea Comodino diam42xh45"

def write_index(manifest):
    cards = "".join(f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>' for item in manifest["views"])
    bonsai = "".join(f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{name}.png"><img src="bonsai-camera-{name}.png"></a></article>' for view, name in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso")))
    (shared.OUTPUT_DIR / "index.html").write_text(f'''<!doctype html><html><meta charset="utf-8"><title>Baxter Ninfea BST01 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Baxter Ninfea / project BST01</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Both official left- and right-opening vectors are archived as evidence, but the project type does not identify the opening side and native CAD requires login; therefore no blue line or variant substitution is shown.</p><nav><a href="review-contact-sheet.png">Review contact sheet</a><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/NINFCLE40D-right-opening.svg">Official R evidence</a><a href="official-source/NINFCLE40S-left-opening.svg">Official L evidence</a><a href="project-context-furniture-plan.svg">Project plan</a><a href="project-context-front-elevation.svg">Project front</a><a href="project-context-side-elevation.svg">Project side</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.baxter.it/en/products/ninfea-tables-and-coffee-tables">Official page</a></nav><h2>Project context</h2><main><article><h2>Furniture Plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article><article><h2>R09 Front</h2><a href="project-context-front-elevation-review.svg"><img src="project-context-front-elevation-review-preview.png"></a></article><article><h2>R09 Side</h2><a href="project-context-side-elevation-review.svg"><img src="project-context-side-elevation-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''', encoding="utf-8")

shared.write_index = write_index

if __name__ == "__main__":
    shared.main()
