#!/usr/bin/env python3
"""Generate the Molteni 505 UP project-composition review."""

import math

from falper_sorgente_linework import ROOT
from int1_highpoly_type_review import projected_silhouette
import gessi316_54294_review as shared


shared.PROFILE_KEY = "505-up-v1-lp-s"
shared.ARTICLE_NUMBER = "505 UP System / project V1.LP.S"
shared.GENERATOR = "pipeline/scripts/505_up_v1_lp_s_review.py"
shared.OFFICIAL_CAD_STATUS = "official_native_family_dwg_archived_exact_project_configuration_not_matched"
shared.OFFICIAL_CAD_ACQUIRED = True
shared.OFFICIAL_CAD_EXACT_PROJECT_CONFIGURATION_MATCH = False
shared.OUTPUT_DIR = ROOT / "output/review/highpoly-types/505-up-v1-lp-s"
shared.REGISTER = shared.OUTPUT_DIR / "profile.json"
shared.ACCESS_RECORD = shared.OUTPUT_DIR / "official-source/source-access-record.json"
shared.EXPECTED_DESCRIPTION = None


def structural_proxy_builder(vertices, faces, axes, simplify_mm, view):
    """Keep the Body silhouette plus true mesh creases, without triangle diagonals."""
    silhouette = projected_silhouette(vertices, faces, axes, simplify_mm)
    if view == "plan":
        return silhouette
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

    crease_limit = math.cos(math.radians(3.0))
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
        if abs(dy) <= 0.8 and abs(dx) >= 24.0:
            coordinate = round(((p1[1] + p2[1]) / 2.0) * 2.0) / 2.0
            axis_intervals["h"].setdefault(coordinate, []).append(sorted((p1[0], p2[0])))
        elif abs(dx) <= 0.8 and abs(dy) >= 24.0:
            coordinate = round(((p1[0] + p2[0]) / 2.0) * 2.0) / 2.0
            axis_intervals["v"].setdefault(coordinate, []).append(sorted((p1[1], p2[1])))

    structural = []
    for direction, groups in axis_intervals.items():
        for coordinate, intervals in groups.items():
            intervals.sort()
            merged = []
            for start, end in intervals:
                if merged and start <= merged[-1][1] + 1.0:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            for start, end in merged:
                if end - start < 30.0:
                    continue
                structural.append(
                    [(start, coordinate), (end, coordinate)]
                    if direction == "h"
                    else [(coordinate, start), (coordinate, end)]
                )
    return silhouette + structural


shared.PROXY_BUILDER = structural_proxy_builder


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
        f'''<!doctype html><html><meta charset="utf-8"><title>Molteni 505 UP review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Molteni&C 505 UP System / project 505 UP V1.LP.S</h1><p>Official native family DWGs are archived below. No published catalogue composition matches the project-specific slatted-left / display-right arrangement, so the candidate drawing line remains {shared.SOURCE_LABEL_EN}; official CAD is not rearranged or falsely shown as blue project geometry.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="profile.json">Profile</a><a href="official-source/source-access-record.json">Source record</a><a href="project-context-furniture-plan-review.svg">Project furniture plan</a><a href="bonsai-review-manifest.json">Bonsai evidence</a><a href="official-source/molteni-505-up-technical-library-full-preview.svg">Official 2021 CAD SVG</a><a href="official-source/molteni-505-up-inspiring-solution-full-preview.svg">Official inspiring solutions CAD SVG</a><a href="https://www.molteni.it/en/ap/product/505-up-system">Official page</a></nav><h2>Project drawing context</h2><main><article><h2>Furniture Plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article></main><h2>Three-view project candidate</h2><main>{cards}</main><h2>Actual Bonsai camera renders from isolated IFC Body</h2><main>{bonsai_cards}</main><code>{manifest["formal_ifc_sha256"]}</code></html>''',
        encoding="utf-8",
    )


shared.write_index = write_index


if __name__ == "__main__":
    shared.main()
