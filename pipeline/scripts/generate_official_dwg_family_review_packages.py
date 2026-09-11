#!/usr/bin/env python3
"""Render four supplemental official-DWG/IFC review packages without IFC writes."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell

from falper_sorgente_linework import ROOT
from int1_highpoly_type_review import (
    VIEWS,
    bounds_3d,
    display_edge_sample,
    mesh_for_one_product,
    projected_raw_edges,
    svg_path,
)


FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
BLUE = "#1677c8"
PRODUCTS = {
    "bed01": {
        "title": "Baxter Casablanca 180 / BED01",
        "scope": "Official independent 2D DWG; mechanically selected native 180×200 Plan/Front/Side linework. The prior ACIS-solid envelope is superseded and excluded from approval.",
        "fit": "Native 2D paths remain 1:1 and receive translation only. The official 2100×2520×900 mm-labelled configuration does not exactly match the project Body; no fit or scaling is applied.",
    },
    "bed02": {
        "title": "Baxter Viktor 160×200 family / BED02",
        "scope": "Official six-spec family CAD; selected 160×200 group is identity/family reference only.",
        "fit": "Project high-poly Body was non-uniformly reduced in width and is not an exact official configuration; no fit or scaling applied.",
    },
    "sis04": {
        "title": "Molteni&C Sistema 7 Wall Unit / SIS04",
        "scope": "Exact official Wall Unit native DWG; selected project configuration is 1962×370×722 mm.",
        "fit": "Project Body is about 1961.975×369×721.991 mm. Overall depth is 370 mm; catalogue 331 mm is the internal depth.",
    },
    "hima01": {
        "title": "Poliform HIMA PVA11 / HIMA01",
        "scope": "Official PVA11 three-element, 1000 mm-high family reference.",
        "fit": "The project Body height agrees, but its folded pose is not the published straight PVA11 arrangement; no pose deformation or scaling applied.",
    },
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def path_bounds(paths: list[list[list[float]]]) -> dict:
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return {
        "minimum": minimum,
        "maximum": maximum,
        "size": [maximum[axis] - minimum[axis] for axis in range(2)],
        "centre": [(minimum[axis] + maximum[axis]) / 2.0 for axis in range(2)],
    }


def align_reference(view: str, paths: list, proxy: list, reflect_x: bool) -> tuple[list, dict]:
    source = [[list(point) for point in path] for path in paths]
    source_bounds = path_bounds(source)
    if reflect_x:
        centre_x = source_bounds["centre"][0]
        source = [
            [[2.0 * centre_x - x, y] for x, y in path]
            for path in source
        ]
        source_bounds = path_bounds(source)
    proxy_bounds = path_bounds(proxy)
    translate_x = proxy_bounds["centre"][0] - source_bounds["centre"][0]
    if view == "plan":
        translate_y = proxy_bounds["centre"][1] - source_bounds["centre"][1]
        anchor = "projection_centre_to_projection_centre"
    else:
        translate_y = proxy_bounds["minimum"][1] - source_bounds["minimum"][1]
        anchor = "horizontal_centre_and_finished_bottom"
    aligned = [
        [[round(x + translate_x, 6), round(y + translate_y, 6)] for x, y in path]
        for path in source
    ]
    return aligned, {
        "mode": "rigid_reflection_and_translation_only" if reflect_x else "translation_only",
        "view_direction_reflection_x": reflect_x,
        "translation_mm": [round(translate_x, 6), round(translate_y, 6)],
        "anchor": anchor,
        "uniform_scale": 1.0,
        "anisotropic_scale_used": False,
        "source_geometry_deformed": False,
        "source_bounds_before_alignment_mm": {
            "minimum": [round(value, 6) for value in source_bounds["minimum"]],
            "maximum": [round(value, 6) for value in source_bounds["maximum"]],
            "size": [round(value, 6) for value in source_bounds["size"]],
        },
        "project_proxy_bounds_mm": {
            "minimum": [round(value, 6) for value in proxy_bounds["minimum"]],
            "maximum": [round(value, 6) for value in proxy_bounds["maximum"]],
            "size": [round(value, 6) for value in proxy_bounds["size"]],
        },
    }


def render_svg(title: str, global_id: str, view: str, raw_edges: list, proxy: list, official: list, source: dict, notes: dict) -> str:
    width, height = 1500, 1000
    plot_x, plot_y, plot_w, plot_h = 55, 185, 1050, 750
    points = [point for edge in raw_edges for point in edge]
    points.extend(point for path in proxy for point in path)
    points.extend(point for path in official for point in path)
    minimum_x = min(point[0] for point in points)
    maximum_x = max(point[0] for point in points)
    minimum_y = min(point[1] for point in points)
    maximum_y = max(point[1] for point in points)
    padding = max(maximum_x - minimum_x, maximum_y - minimum_y) * 0.06 or 1.0
    minimum_x -= padding
    maximum_x += padding
    minimum_y -= padding
    maximum_y += padding
    scale = min(plot_w / (maximum_x - minimum_x), plot_h / (maximum_y - minimum_y))

    def transform(point):
        return (
            plot_x + (point[0] - minimum_x) * scale,
            plot_y + plot_h - (point[1] - minimum_y) * scale,
        )

    raw_path = svg_path([[start, end] for start, end in raw_edges], transform)
    proxy_path = svg_path(proxy, transform, close=True)
    official_path = svg_path(official, transform)
    dashed = source["blue_stroke_style"] == "dashed"
    dash_attribute = ' stroke-dasharray="11 7"' if dashed else ""
    kind_label = {
        "native_dwg_2d_linework": "solid blue = mechanically selected native 2D DWG paths",
        "native_dwg_3d_solid_projected_envelope": "dashed blue = native DWG ACIS-solid projected envelope; not a 2D manufacturer view",
        "native_dwg_dimension_envelope_not_dedicated_2d_view": "dashed blue = dimension envelope from the same native DWG; no dedicated 2D path set",
    }[source["geometry_kind"]]
    source_note = html.escape(notes["scope"])
    fit_note = html.escape(notes["fit"])
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="1500" height="1000" fill="#fbfaf7"/>
<text x="55" y="54" font-family="Arial,sans-serif" font-size="29" font-weight="700" fill="#1f2d3d">{html.escape(title)}</text>
<text x="55" y="92" font-family="Arial,sans-serif" font-size="19" fill="#41566d">{VIEWS[view]["label"]} · project representative {html.escape(global_id)}</text>
<text x="55" y="126" font-family="Arial,sans-serif" font-size="16" fill="#68798a">Grey = actual IFC Body · Black = current geometry-derived review proxy · Blue = official-DWG review reference</text>
<text x="55" y="154" font-family="Arial,sans-serif" font-size="15" fill="{BLUE}">{html.escape(kind_label)}</text>
<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="actual-ifc-body" data-source-kind="project_ifc_body" d="{raw_path}" fill="none" stroke="#87929c" stroke-width="0.5" stroke-opacity="0.3" vector-effect="non-scaling-stroke"/>
<path class="geometry-derived-simplified-proxy" data-source-kind="geometry_derived_simplified_proxy" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="2.4" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-reference native-dwg" data-source-role="review_reference_only" data-geometry-kind="{source["geometry_kind"]}" data-source-scaled="false" d="{official_path}" fill="none" stroke="{BLUE}" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"{dash_attribute} vector-effect="non-scaling-stroke"/>
<text x="1140" y="215" font-family="Arial,sans-serif" font-size="19" font-weight="700" fill="#1f2d3d">Review scope</text>
<foreignObject x="1140" y="240" width="320" height="150"><div xmlns="http://www.w3.org/1999/xhtml" style="font:15px Arial;color:#41566d;line-height:1.4">{source_note}</div></foreignObject>
<text x="1140" y="410" font-family="Arial,sans-serif" font-size="19" font-weight="700" fill="#1f2d3d">Mechanical fit</text>
<foreignObject x="1140" y="435" width="320" height="160"><div xmlns="http://www.w3.org/1999/xhtml" style="font:15px Arial;color:#41566d;line-height:1.4">{fit_note}</div></foreignObject>
<text x="1140" y="640" font-family="Arial,sans-serif" font-size="18" font-weight="700" fill="#1f2d3d">Source semantics</text>
<text x="1140" y="672" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Official reference role only</text>
<text x="1140" y="700" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Primary candidate unchanged</text>
<text x="1140" y="728" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Scale: 1.0 · no fitting</text>
<text x="1140" y="756" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Third-party CAD: false</text>
<text x="1140" y="784" font-family="Arial,sans-serif" font-size="14" fill="#41566d">IFC write: none</text>
</svg>'''


def context_cards(folder: Path) -> str:
    context = load_json(folder / "project-context-manifest.json")
    cards = []
    for record in context.get("views", [])[:3]:
        target = Path(record.get("output") or record.get("path") or "").name
        preview = record.get("review_preview")
        if preview:
            image = Path(preview).name
        else:
            image = target
        cards.append(f'<article><h2>Project {html.escape(str(record["view"]))}</h2><a href="{html.escape(target)}"><img src="{html.escape(image)}"></a></article>')
    return "".join(cards)


def write_index(folder: Path, config: dict, reference: dict, package: dict) -> Path:
    cards = "".join(
        f'<article><h2>{view.title()}</h2><a href="official-dwg-{view}.svg"><img src="official-dwg-{view}.svg"></a></article>'
        for view in ("plan", "front", "side")
    )
    bonsai = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    product_prefix = Path("output/review/highpoly-types") / folder.name
    source_evidence_paths = [
        *reference.get("autocad_source_screenshots", []),
        reference.get("source_dwg"),
        reference.get("source_archive"),
        reference.get("source_archive_inventory"),
    ]
    source_evidence = "".join(
        f'<a href="{html.escape(Path(path).relative_to(product_prefix).as_posix())}">{html.escape(Path(path).name)}</a>'
        for path in source_evidence_paths if path
    )
    superseded_note = (
        '<p><strong>Superseded/error candidate:</strong> the earlier ACIS 3DSOLID projected-envelope package is excluded from approval and retained only in the mechanical record.</p>'
        if reference.get("superseded_error_candidates") else ""
    )
    target = folder / "official-dwg-review-index.html"
    target.write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>{html.escape(config["title"])} official DWG review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>{html.escape(config["title"])}</h1><p>Grey = actual IFC Body; black = current <strong>基于原始高模几何生成的简化图纸表达</strong>; solid blue = mechanically selected native DWG linework; dashed blue = explicitly labelled native-DWG solid/dimension envelope where no dedicated 2D view exists.</p><p>{html.escape(reference["note"])}</p>{superseded_note}<nav><a href="official-dwg-review-contact-sheet.png">Contact sheet</a><a href="official-dwg-review-manifest.json">Review manifest</a><a href="official-dwg-review-reference.json">Mechanical selection record</a><a href="manifest.json">Primary manifest (unchanged source semantics)</a><a href="candidate-representations.json">Primary candidate</a><a href="official-source/source-access-record.json">Source record</a><a href="{html.escape(reference["official_download_url"])}">Official download</a>{source_evidence}</nav><h2>Official-DWG / actual IFC comparison</h2><main>{cards}</main><h2>Actual Bonsai saved-camera renders</h2><main>{bonsai}</main><h2>Available project drawing context</h2><main>{context_cards(folder)}</main><p>Formal IFC SHA-256: <code>{package["formal_ifc_sha256"]}</code>. No derived or authoritative IFC was written.</p></html>''',
        encoding="utf-8",
    )
    primary_index = folder / "index.html"
    if primary_index.is_file():
        content = primary_index.read_text(encoding="utf-8")
        if "official-dwg-review-index.html" not in content:
            link = '<a href="official-dwg-review-index.html">Official DWG SVG review supplement</a>'
            content = content.replace("</nav>", f"{link}</nav>", 1) if "</nav>" in content else content.replace("</body>", f"<p>{link}</p></body>")
            primary_index.write_text(content, encoding="utf-8")
    return target


def write_contact_sheet(folder: Path, config: dict, package: dict) -> tuple[Path, Path]:
    bonsai = {
        view: filename for view, filename in (
            ("plan", "bonsai-camera-plan.png"),
            ("front", "bonsai-camera-front-elevation.png"),
            ("side", "bonsai-camera-side-elevation.png"),
            ("iso", "bonsai-camera-iso.png"),
        )
    }
    context = load_json(folder / "project-context-manifest.json")
    context_html = []
    for record in context.get("views", [])[:3]:
        preview = record.get("review_preview") or record.get("output") or record.get("path")
        context_html.append(f'<article><h2>Project · {html.escape(str(record["view"]))}</h2><img src="{html.escape(Path(preview).name)}"></article>')
    while len(context_html) < 3:
        context_html.append('<article class="empty"><h2>Project context unavailable</h2></article>')
    views = "".join(f'<article><h2>DWG review · {view.title()}</h2><img src="official-dwg-{view}.svg"></article>' for view in ("plan", "front", "side"))
    renders = "".join(f'<article class="bonsai"><h2>Bonsai · {view.title()}</h2><img src="{filename}"></article>' for view, filename in bonsai.items())
    html_path = folder / "official-dwg-review-contact-sheet.html"
    png_path = folder / "official-dwg-review-contact-sheet.png"
    html_path.write_text(
        f'''<!doctype html><html><meta charset="utf-8"><style>*{{box-sizing:border-box}}html,body{{margin:0;width:2400px;height:1800px;overflow:hidden;background:#e9e6df;color:#182430;font-family:Arial,sans-serif}}header{{height:100px;padding:18px 28px;background:#f8f7f3;border-bottom:2px solid #b8bec4}}h1{{margin:0 0 5px;font-size:30px}}header p{{margin:0;font-size:16px}}main{{display:grid;grid-template-columns:repeat(4,1fr);grid-template-rows:repeat(3,556px);gap:10px;padding:10px}}article{{position:relative;overflow:hidden;background:#fff;border:1px solid #c8ced4}}h2{{position:absolute;z-index:2;left:12px;top:10px;margin:0;padding:6px 9px;background:rgba(255,255,255,.92);font-size:18px}}img{{width:100%;height:100%;object-fit:contain}}.bonsai{{background:#303030}}.summary{{padding:75px 25px;font-size:20px;line-height:1.45}}.empty{{background:#f5f3ee}}</style><header><h1>{html.escape(config["title"])}</h1><p>Grey actual IFC Body · black geometry-derived review proxy · solid blue native 2D DWG · dashed blue explicitly labelled DWG envelope · formal IFC unchanged</p></header><main>{views}<article class="summary"><h2>Source / scope</h2><p>{html.escape(config["scope"])}</p><p>{html.escape(config["fit"])}</p><p>No derived IFC · no authoritative IFC write.</p></article>{renders}{''.join(context_html)}<article class="summary"><h2>Evidence</h2><p>DWG SHA-256<br>{html.escape(package["source_dwg_sha256"])}</p><p>Formal IFC<br>{html.escape(package["formal_ifc_sha256"])}</p></article></main></html>''',
        encoding="utf-8",
    )
    with tempfile.TemporaryDirectory(prefix="official-dwg-contact-sheet-") as profile:
        png_path.unlink(missing_ok=True)
        process = subprocess.Popen([
            str(CHROME), "--headless=new", "--disable-gpu", "--disable-background-networking",
            "--no-first-run", "--hide-scrollbars", "--force-device-scale-factor=1",
            "--window-size=2400,1800", f"--user-data-dir={profile}",
            f"--screenshot={png_path}", html_path.as_uri(),
        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(300):
            if png_path.is_file() and png_path.stat().st_size:
                break
            if process.poll() is not None:
                break
            time.sleep(0.1)
        if process.poll() is None:
            process.terminate()
            process.wait(timeout=3)
    if not png_path.is_file() or png_path.stat().st_size == 0:
        raise RuntimeError(f"contact sheet render failed: {folder.name}")
    return html_path, png_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--slug", choices=PRODUCTS)
    args = parser.parse_args()
    if not CHROME.is_file():
        raise RuntimeError("Google Chrome is required for contact-sheet evidence")
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before review generation")
    model = ifcopenshell.open(FORMAL_IFC)
    selected_products = PRODUCTS.items() if args.slug is None else [(args.slug, PRODUCTS[args.slug])]
    for slug, config in selected_products:
        folder = ROOT / "output/review/highpoly-types" / slug
        reference_path = folder / "official-dwg-review-reference.json"
        reference = load_json(reference_path)
        if sha256(ROOT / reference["source_dwg"]) != reference["source_dwg_sha256"]:
            raise RuntimeError(f"{slug} source DWG hash mismatch")
        primary_manifest = load_json(folder / "manifest.json")
        candidate = load_json(folder / "candidate-representations.json")
        if candidate.get("source_kind") != "geometry_derived_simplified_proxy" or candidate.get("formal_ifc_write_allowed") is not False:
            raise RuntimeError(f"{slug} primary candidate source/write gate drifted")
        product, vertices, faces = mesh_for_one_product(model, primary_manifest["representative_global_id"])
        minimum, maximum = bounds_3d(vertices)
        package_views = []
        for view, definition in VIEWS.items():
            proxy = candidate["views"][view]["proxy_paths_mm"]
            source_view = reference["views"][view]
            official, alignment = align_reference(
                view,
                source_view["paths_mm"],
                proxy,
                source_view.get("reflect_x_for_project_view_direction", False),
            )
            all_edges = projected_raw_edges(vertices, faces, definition["axes"])
            edges = display_edge_sample(all_edges, maximum=850)
            target = folder / f"official-dwg-{view}.svg"
            target.write_text(render_svg(config["title"], product.GlobalId, view, edges, proxy, official, source_view, config), encoding="utf-8")
            package_views.append({
                "view": view,
                "svg": relative(target),
                "svg_sha256": sha256(target),
                "official_reference_geometry_kind": source_view["geometry_kind"],
                "official_reference_path_count": source_view["path_count"],
                "blue_stroke_style": source_view["blue_stroke_style"],
                "alignment": alignment,
                "project_minus_official_bounds_size_mm": [
                    round(alignment["project_proxy_bounds_mm"]["size"][axis] - alignment["source_bounds_before_alignment_mm"]["size"][axis], 6)
                    for axis in range(2)
                ],
                "primary_candidate_source_kind": candidate["source_kind"],
                "primary_candidate_replaced": False,
            })
        package = {
            "schema_version": 1,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "generator": "pipeline/scripts/generate_official_dwg_family_review_packages.py",
            "product_slug": slug,
            "display_name": config["title"],
            "representative_global_id": product.GlobalId,
            "formal_ifc_sha256": FORMAL_SHA256,
            "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
            "source_dwg": reference["source_dwg"],
            "source_dwg_absolute_path": reference["source_dwg_absolute_path"],
            "source_dwg_sha256": reference["source_dwg_sha256"],
            "official_download_url": reference["official_download_url"],
            "source_archive": reference.get("source_archive"),
            "source_archive_sha256": reference.get("source_archive_sha256"),
            "source_archive_member": reference.get("source_archive_member"),
            "source_archive_member_sha256": reference.get("source_archive_member_sha256"),
            "source_archive_inventory": reference.get("source_archive_inventory"),
            "source_archive_inventory_sha256": reference.get("source_archive_inventory_sha256"),
            "local_byte_identical_archive_copies": reference.get("local_byte_identical_archive_copies", []),
            "selected_geometry_layers": reference.get("selected_geometry_layers"),
            "excluded_layers": reference.get("excluded_layers"),
            "autocad_source_screenshots": reference.get("autocad_source_screenshots", []),
            "superseded_error_candidates": reference.get("superseded_error_candidates", []),
            "configuration": reference["configuration"],
            "reference_role": reference["reference_role"],
            "scope": config["scope"],
            "fit_note": config["fit"],
            "official_dwg_review_reference": relative(reference_path),
            "official_dwg_review_reference_sha256": sha256(reference_path),
            "primary_candidate_source_kind": candidate["source_kind"],
            "primary_candidate_source_label_zh": candidate["source_label_zh"],
            "primary_candidate_replaced": False,
            "official_reference_used_as_primary_candidate": False,
            "third_party_cad_used": False,
            "derived_ifc_write_performed": False,
            "formal_authoritative_ifc_write_performed": False,
            "bounds_mm": {
                "minimum": [round(value, 6) for value in minimum],
                "maximum": [round(value, 6) for value in maximum],
                "size": [round(maximum[index] - minimum[index], 6) for index in range(3)],
            },
            "views": package_views,
            "pass": True,
        }
        package_path = folder / "official-dwg-review-manifest.json"
        write_json(package_path, package)
        index_path = write_index(folder, config, reference, package)
        contact_html, contact_png = write_contact_sheet(folder, config, package)
        package["index"] = relative(index_path)
        package["index_sha256"] = sha256(index_path)
        package["contact_sheet_html"] = relative(contact_html)
        package["contact_sheet_html_sha256"] = sha256(contact_html)
        package["contact_sheet_png"] = relative(contact_png)
        package["contact_sheet_png_sha256"] = sha256(contact_png)
        package["formal_ifc_bytes_unchanged"] = sha256(FORMAL_IFC) == FORMAL_SHA256
        write_json(package_path, package)
        print(relative(package_path))
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during review generation")


if __name__ == "__main__":
    main()
