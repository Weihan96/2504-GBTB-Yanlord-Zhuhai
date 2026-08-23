#!/usr/bin/env python3
"""Generate the Baxter Stone / BST03 review without claiming official CAD."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "bst03"
shared.ARTICLE_NUMBER = "Baxter Stone L drawer 45 / BST03"
shared.GENERATOR = "pipeline/scripts/bst03_review.py"
shared.OFFICIAL_CAD_STATUS = "manufacturer_login_required_not_acquired"
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/bst03"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "Stone Bedside Table with Drawer D45"

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
        f'''<!doctype html><html><meta charset="utf-8"><title>Baxter Stone BST03 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Baxter Stone freestanding bedside table with L drawer / project BST03</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Baxter places native 2D/3D/BIM downloads behind login, so no blue CAD line is shown. The archived official vector PDF is identity and nominal-dimension evidence only.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/Baxter_Stone_technical-sheet-page-10-preview.png">Official technical-sheet page 10</a><a href="project-context-furniture-plan.svg">Project furniture plan</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.baxter.it/en/products/stone-beds">Official page</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture Plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
