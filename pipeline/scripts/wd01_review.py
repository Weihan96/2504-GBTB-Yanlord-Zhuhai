#!/usr/bin/env python3
"""Generate the Poliform Pivot + Senzafine / WD01 review without claiming official CAD."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "wd01"
shared.ARTICLE_NUMBER = "Poliform Pivot + Senzafine custom wardrobe / project WD01"
shared.GENERATOR = "pipeline/scripts/wd01_review.py"
shared.OFFICIAL_CAD_STATUS = "gated_resource_download_no_public_exact_WD01_configuration_asset_located_not_acquired"
shared.OFFICIAL_CAD_ACQUIRED = False
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/wd01"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "Poliform Pivot Senzafine"


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
        f'''<!doctype html><html><meta charset="utf-8"><title>Poliform Pivot + Senzafine WD01 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Poliform Pivot + Senzafine custom wardrobe / project WD01</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Official Poliform pages and the remote Pivot technical document establish the integrated Pivot/Senzafine system, but the 1202.879 × 673.524 × 2390.023 mm project Body is a custom arrangement and no exact native CAD asset was acquired. No blue line is shown and no catalogue composition is substituted.</p><nav><a href="review-contact-sheet.png">Review contact sheet</a><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/official-product-page-evidence.json">Official-page evidence</a><a href="project-context-furniture-plan.svg">Project furniture plan</a><a href="project-context-front-elevation.svg">R09 +X front elevation</a><a href="project-context-side-elevation.svg">R09 +Y side elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.poliform.it/en/products/pivot/">Official Pivot page</a><a href="https://www.poliform.it/en/products/senzafine-wardrobe/">Official Senzafine page</a><a href="https://www.poliform.it/wp-content/uploads/pdf/278373-pivot-poliform-en-us.pdf">Official Pivot technical PDF</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article><article><h2>R09 +X front elevation</h2><a href="project-context-front-elevation-review.svg"><img src="project-context-front-elevation-review-preview.png"></a></article><article><h2>R09 +Y side elevation</h2><a href="project-context-side-elevation-review.svg"><img src="project-context-side-elevation-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
