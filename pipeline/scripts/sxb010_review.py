#!/usr/bin/env python3
"""Generate the sxb010 Venetian-blind review without claiming official CAD."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared


shared.PROFILE_KEY = "sxb010"
shared.ARTICLE_NUMBER = "project sxb010 / Hunter Douglas 25 mm family direction"
shared.GENERATOR = "pipeline/scripts/sxb010_review.py"
shared.OFFICIAL_CAD_STATUS = "not_published_on_verified_manufacturer_surfaces"
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/sxb010"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.EXPECTED_DESCRIPTION = None
shared.VIEW_DEFINITIONS = {
    "plan": {"axes": (0, 2), "label": "PLAN / local XZ (width × depth)"},
    "front": {"axes": (0, 1), "label": "FRONT / local XY (width × height)"},
    "side": {"axes": (2, 1), "label": "SIDE / local ZY (depth × height)"},
}


def write_index(manifest):
    cards = "".join(
        f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>'
        for item in manifest["views"]
    )
    semantic_bonsai = (
        ("Plan", "bonsai-camera-semantic-plan.png", "semantic Plan / local XZ"),
        ("Front", "bonsai-camera-semantic-front-elevation.png", "semantic Front / local XY"),
        ("Side", "bonsai-camera-semantic-side-elevation.png", "semantic Side / local ZY"),
        ("Iso", "bonsai-camera-iso.png", "isometric"),
    )
    bonsai_cards = "".join(
        f'<article><h2>Bonsai {label}</h2><p>{note}</p><a href="{filename}"><img src="{filename}"></a></article>'
        for label, filename, note in semantic_bonsai
    )
    context_cards = "".join(
        f'<article><h2>Project {label}</h2><a href="{filename}"><img src="{preview}"></a></article>'
        for label, filename, preview in (
            ("Plan", "project-context-ffl-plan.svg", "project-context-ffl-plan-review-preview.png"),
            ("Front", "project-context-r07-front-elevation.svg", "project-context-r07-front-elevation-review-preview.png"),
            ("Side", "project-context-r07-side-elevation.svg", "project-context-r07-side-elevation-review-preview.png"),
        )
    )
    (shared.OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>sxb010 Venetian blind review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Project sxb010 / Hunter Douglas 25 mm Venetian-blind direction</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Owner evidence confirms the Hunter Douglas 25 mm family and matte graphite/grey-black direction only. Exact SKU, control type, finish code and shop dimensions remain pending. The verified official download catalogue returns PDF only, so no blue line is shown.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-manifest.json">Project context evidence</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.hunterdouglas.cn/product/venetian-blind/16mm-25mm-venetian-blinds">Official page</a></nav><h2>Three-view review</h2><main>{cards}</main><h2>Complete project drawing context</h2><main>{context_cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
