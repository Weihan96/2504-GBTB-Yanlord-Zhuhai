#!/usr/bin/env python3
"""Generate the Baxter Beside / BST02 review without claiming official CAD."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "bst02"
shared.ARTICLE_NUMBER = "Baxter Beside 55x58xh35 / BST02"
shared.GENERATOR = "pipeline/scripts/bst02_review.py"
shared.OFFICIAL_CAD_STATUS = "manufacturer_login_required_not_acquired"
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/bst02"
shared.REGISTER = ROOT / "pipeline/decisions/highpoly-drawing-profile-register.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "Beside 55x58xh35"


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
        f'''<!doctype html><html><meta charset="utf-8"><title>Baxter Beside BST02 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Baxter Beside 55x58xh35 / project BST02</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Native 2D/3D/BIM downloads require Baxter login, so no blue CAD line is shown. The official public PDF/SVG are identity and nominal-dimension evidence only. The actual IFC Body is 250 mm high versus the official 350 mm and is not stretched.</p><nav><a href="review-contact-sheet.png">Review contact sheet</a><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/BESICOBS55.svg">Official measurement SVG</a><a href="official-source/Baxter_Beside_TechnicalSheet-page-7.png">Official technical sheet page 7</a><a href="project-context-furniture-plan.svg">Project furniture plan</a><a href="project-context-side-elevation.svg">Project R09 side elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.baxter.it/en/products/beside-tables-and-coffee-tables">Official page</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture Plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article><article><h2>R09 Side Elevation</h2><a href="project-context-side-elevation-review.svg"><img src="project-context-side-elevation-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
