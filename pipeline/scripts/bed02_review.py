#!/usr/bin/env python3
"""Generate the Baxter Viktor / BED02 review without claiming official CAD."""

from pathlib import Path

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "bed02"
shared.ARTICLE_NUMBER = "Baxter Viktor / BED02"
shared.GENERATOR = "pipeline/scripts/bed02_review.py"
shared.OFFICIAL_CAD_STATUS = "manufacturer_login_required_not_acquired"
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/bed02"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "VIKTOR 162x234xh106"


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
        f'''<!doctype html><html><meta charset="utf-8"><title>Baxter Viktor BED02 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Baxter Viktor / project BED02</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Official Baxter 2D/3D/BIM downloads require login, so no blue line is shown. The official 1720 mm width, IFC description 1620 mm width and IFC Body 1675.711 mm width remain an explicit review difference.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-furniture-plan.svg">Project furniture plan</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.baxter.it/en/products/viktor-beds">Official page</a></nav><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
