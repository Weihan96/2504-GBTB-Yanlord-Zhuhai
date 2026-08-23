#!/usr/bin/env python3
"""Generate the RODA Bernardo 367 / TAB02 review without claiming official CAD."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "tab02"
shared.ARTICLE_NUMBER = "RODA BERNARDO 367 / TAB02"
shared.GENERATOR = "pipeline/scripts/tab02_review.py"
shared.OFFICIAL_CAD_STATUS = "reserved_area_authentication_required_not_acquired"
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/tab02"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.EXPECTED_DESCRIPTION = "Bernardo 367"


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
        f'''<!doctype html><html><meta charset="utf-8"><title>RODA Bernardo 367 TAB02 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>RODA Bernardo 367 side table / project TAB02</h1><p>Black line = {shared.SOURCE_LABEL_EN}. The official RODA catalogue confirms the exact 500 x 500 x 670 mm model, but native 2D/3D files require reserved-area login and were not acquired; no blue line is shown.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/RODA-official-Bernardo-367-catalogue-extract.pdf">Official catalogue extract</a><a href="project-context-furniture-plan.svg">Project furniture plan</a><a href="project-context-furniture-plan-review.svg">Project plan review</a><a href="project-context-manifest.json">Context manifest</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.rodaonline.com/en/collections/bernardo/">Official page</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture Plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
