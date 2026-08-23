#!/usr/bin/env python3
"""Generate the Baxter Miami Soft E09 right-terminal review package."""

import math

from falper_sorgente_linework import ROOT
from int1_highpoly_type_review import projected_silhouette
import gessi316_54294_review as shared


def upholstered_proxy_builder(vertices, faces, axes, simplify_mm, view):
    """Keep the soft envelope and only long, mechanically stable upholstery creases."""
    silhouette = projected_silhouette(vertices, faces, axes, simplify_mm)
    edge_normals = {}
    for face in faces:
        a, b, c = (vertices[index] for index in face)
        ab = tuple(b[index] - a[index] for index in range(3))
        ac = tuple(c[index] - a[index] for index in range(3))
        normal = (
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        )
        length = math.sqrt(sum(value * value for value in normal))
        if length <= 1e-9:
            continue
        normal = tuple(value / length for value in normal)
        for start, end in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge_normals.setdefault(tuple(sorted((start, end))), []).append(normal)

    crease_limit = math.cos(math.radians(20.0))
    axis_intervals = {"h": {}, "v": {}}
    for (start, end), normals in edge_normals.items():
        sharp = len(normals) == 1 or any(
            abs(sum(a * b for a, b in zip(normals[i], normals[j]))) < crease_limit
            for i in range(len(normals))
            for j in range(i + 1, len(normals))
        )
        if not sharp:
            continue
        p1 = (vertices[start][axes[0]], vertices[start][axes[1]])
        p2 = (vertices[end][axes[0]], vertices[end][axes[1]])
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        if abs(dy) <= 1.5 and abs(dx) >= 80.0:
            coordinate = round(((p1[1] + p2[1]) / 2.0) / 2.0) * 2.0
            axis_intervals["h"].setdefault(coordinate, []).append(sorted((p1[0], p2[0])))
        elif abs(dx) <= 1.5 and abs(dy) >= 80.0:
            coordinate = round(((p1[0] + p2[0]) / 2.0) / 2.0) * 2.0
            axis_intervals["v"].setdefault(coordinate, []).append(sorted((p1[1], p2[1])))

    structural = []
    for direction, groups in axis_intervals.items():
        for coordinate, intervals in groups.items():
            intervals.sort()
            merged = []
            for start, end in intervals:
                if merged and start <= merged[-1][1] + 8.0:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            for start, end in merged:
                if end - start < 100.0:
                    continue
                structural.append(
                    [(start, coordinate), (end, coordinate)]
                    if direction == "h"
                    else [(coordinate, start), (coordinate, end)]
                )
    return silhouette + structural


shared.PROFILE_KEY = "miamisoft-e09"
shared.ARTICLE_NUMBER = "Baxter Miami Soft E09 dx/r"
shared.GENERATOR = "pipeline/scripts/miamisoft_e09_review.py"
shared.OFFICIAL_CAD_STATUS = "manufacturer_login_required_not_acquired_exact_vector_references_archived"
shared.OFFICIAL_CAD_ACQUIRED = False
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/miamisoft-e09"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = "R terminal module 130 x 108 h70/80 cm"
shared.PROXY_BUILDER = upholstered_proxy_builder


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
        f'''<!doctype html><html><meta charset="utf-8"><title>Baxter Miami Soft E09 review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Baxter Miami Soft E09 dx/r right terminal module</h1><p>Black line = {shared.SOURCE_LABEL_EN}. Baxter's official PDF and E09 measurement SVG prove the exact handed variant and nominal dimensions, but authenticated native 2D/3D/BIM files were not acquired; they are evidence only and no blue CAD line is shown.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/Baxter_MiamiSoft_technical-sheet-page-6-preview.png">Official sheet page 6</a><a href="official-source/MIAMSOESE09D.svg">Official E09 measurement SVG</a><a href="project-context-furniture-plan-review.svg">Project plan</a><a href="project-context-r20-side-elevation-review.svg">Project R20 side elevation</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="https://www.baxter.it/gb/prodotti/miami-soft-divani-e-poltrone">Official page</a></nav><h2>Complete project drawing context</h2><main><article><h2>Furniture Plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article><article><h2>R20 +Y side elevation</h2><a href="project-context-r20-side-elevation-review.svg"><img src="project-context-r20-side-elevation-review-preview.png"></a></article></main><h2>Three-view review</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
