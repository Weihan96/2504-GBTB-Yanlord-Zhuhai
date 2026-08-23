#!/usr/bin/env python3
"""Generate the Poliform Senzafine / WD02 review without claiming official CAD."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "wd02"
shared.ARTICLE_NUMBER = "Poliform Senzafine glass wardrobe / project WD02"
shared.GENERATOR = "pipeline/scripts/wd02_review.py"
shared.OFFICIAL_CAD_STATUS = "gated_resource_download_no_public_exact_WD02_asset_located_not_acquired"
shared.OFFICIAL_CAD_ACQUIRED = False
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/wd02"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "Poliform Glass Wardrobe"


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
        f'''<!doctype html><html><meta charset="utf-8"><title>Poliform Senzafine WD02 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Poliform Senzafine glass wardrobe / project WD02</h1><p>Black line = {shared.SOURCE_LABEL_EN}. The official page and project registers establish the Senzafine modular wardrobe family, but an exact native CAD asset for this glass-wardrobe configuration was not acquired; no blue line is shown and no catalogue composition is substituted.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/official-product-page-evidence.json">Official-page evidence</a><a href="project-context-furniture-plan.svg">Project furniture plan</a><a href="project-context-front-elevation.svg">R22 +X front elevation</a><a href="project-context-side-elevation.svg">R22 -Y side elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.poliform.it/en/products/senzafine-wardrobe/">Official product page</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article><article><h2>R22 +X front elevation</h2><a href="project-context-front-elevation-review.svg"><img src="project-context-front-elevation-review-preview.png"></a></article><article><h2>R22 -Y side elevation</h2><a href="project-context-side-elevation-review.svg"><img src="project-context-side-elevation-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
