#!/usr/bin/env python3
"""Generate the Molteni Sistema 7 / SIS04 review without claiming official CAD."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "sis04"
shared.ARTICLE_NUMBER = "Sistema 7 Wall Unit 4 Doors / SIS04"
shared.GENERATOR = "pipeline/scripts/sis04_review.py"
shared.OFFICIAL_CAD_STATUS = "exact_wall_unit_native_cad_not_acquired"
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/sis04"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.EXPECTED_DESCRIPTION = "Sistema 7 Wall Unit 4 Doors"


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    context_cards = "".join(
        f'<article><h2>Project {label}</h2><a href="{filename}"><img src="{preview}"></a></article>'
        for label, filename, preview in (
            ("Plan", "project-context-furniture-plan.svg", "project-context-furniture-plan-review-preview.png"),
            ("Front", "project-context-r04-front-elevation.svg", "project-context-r04-front-elevation-review-preview.png"),
            ("Side", "project-context-r04-side-elevation.svg", "project-context-r04-side-elevation-review-preview.png"),
        )
    )
    (shared.OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Molteni Sistema 7 SIS04 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Molteni&amp;C Sistema 7 Wall Unit 4 Doors / SIS04</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Official catalogue pages establish identity and the 1960 x 722 x 331 mm standard size class, but exact native Wall Unit CAD was not acquired. The near-name full-height Sistema 7 Doors DWG is explicitly excluded, so no blue line is shown.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-manifest.json">Project context evidence</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.molteni.it/ap/product/sistema-7-wall-unit">Official page</a></nav><h2>Complete project drawing context</h2><main>{context_cards}</main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
