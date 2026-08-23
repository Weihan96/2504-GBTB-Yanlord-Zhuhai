#!/usr/bin/env python3
"""Generate the Baxter Miami Soft H01 cushion review package."""

from falper_sorgente_linework import ROOT
import gessi316_54294_review as shared
from miamisoft_e09_review import upholstered_proxy_builder


shared.PROFILE_KEY = "miamisoft-h01"
shared.ARTICLE_NUMBER = "Baxter Miami Soft H01 cushion"
shared.GENERATOR = "pipeline/scripts/miamisoft_h01_review.py"
shared.OFFICIAL_CAD_STATUS = "manufacturer_login_required_not_acquired_exact_vector_references_archived"
shared.OFFICIAL_CAD_ACQUIRED = False
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/miamisoft-h01"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "Cushion 80 x 80 cm"
shared.PROXY_BUILDER = upholstered_proxy_builder
original_render_svg = shared.render_svg


def render_svg(profile, view, raw_edges, proxy, metadata):
    return original_render_svg(profile, view, raw_edges, proxy, metadata).replace(
        "Blue = none (exact project configuration not matched)",
        "Blue = none (native CAD not acquired)",
        1,
    )


shared.render_svg = render_svg


def write_index(manifest):
    cards = "".join(f'<article><h2>{item["view"].title()}</h2><a href="{item["view"]}.svg"><img src="{item["view"]}.svg"></a></article>' for item in manifest["views"])
    bonsai_cards = "".join(f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>' for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso")))
    (shared.OUTPUT_DIR / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Baxter Miami Soft H01 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Baxter Miami Soft H01 cushion</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Baxter's official PDF and H01 measurement SVG prove the exact flexible 80 x 80 cm cushion model, but authenticated native 2D/3D/BIM files were not acquired; they are evidence only and no blue CAD line is shown. The installed IFC Body remains leaning and compressed and is not flattened to the catalogue square.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/Baxter_MiamiSoft_technical-sheet-page-8-preview.png">Official sheet page 8</a><a href="official-source/MIAMSOCUSH01.svg">Official H01 measurement SVG</a><a href="project-context-furniture-plan-review.svg">Project plan</a><a href="project-context-r20-side-elevation-review.svg">Project R20 side elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.baxter.it/gb/prodotti/miami-soft-divani-e-poltrone">Official page</a></nav><h2>Complete project drawing context</h2><main><article><h2>Furniture Plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article><article><h2>R20 +Y side elevation</h2><a href="project-context-r20-side-elevation-review.svg"><img src="project-context-r20-side-elevation-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
