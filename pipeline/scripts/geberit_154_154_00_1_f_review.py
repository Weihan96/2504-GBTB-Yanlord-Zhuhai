#!/usr/bin/env python3
"""Generate the Geberit 154.154.00.1.F flange-component review."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "geberit-154-154-00-1-f"
shared.ARTICLE_NUMBER = "154.154.00.1 / project component .F"
shared.GENERATOR = "pipeline/scripts/geberit_154_154_00_1_f_review.py"
shared.OFFICIAL_CAD_STATUS = "standalone_component_not_published_parent_cadDrawings_empty_native_parent_DWG_urls_404"
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/geberit-154-154-00-1-f"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "Flange for CleanLine shower channel, for screed height at inlet 90–220 mm"


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    context_cards = "".join(
        f'<article><h2>{label}</h2><a href="{svg}"><img src="{preview}"></a></article>'
        for label, svg, preview in (
            ("Sanitary plan", "project-context-sanitary-plan-review.svg", "project-context-sanitary-plan-review-preview.png"),
            ("R12 front elevation", "project-context-r12-front-elevation-review.svg", "project-context-r12-front-elevation-review-preview.png"),
            ("R12 side elevation", "project-context-r12-side-elevation-review.svg", "project-context-r12-side-elevation-review-preview.png"),
        )
    )
    bonsai = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (shared.OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Geberit 154.154.00.1.F review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Geberit CleanLine 154.154.00.1 / project flange component .F</h1><p>Black line = {shared.SOURCE_LABEL_EN}. The official page and EPS identify the complete parent installation set. The project suffix .F is not a separately published manufacturer article and no standalone flange CAD was acquired, so the parent assembly is never shown as blue component linework.</p><nav><a href="review-contact-sheet.png">Review contact sheet</a><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/geberit-154-154-00-1-product-page.html">Archived official parent page</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://catalog.geberit.co.uk/en-GB/product/PRO_199058">Official parent page</a></nav><h2>Project drawing context</h2><main>{context_cards}</main><h2>High-poly / simplified proxy review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
