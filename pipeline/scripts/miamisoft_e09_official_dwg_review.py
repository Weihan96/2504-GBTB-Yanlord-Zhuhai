#!/usr/bin/env python3
"""Build the pending Baxter Miami Soft E09 native-DWG review package.

This script is intentionally review-only: it never writes IFC.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import subprocess
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import ezdxf
import ifcopenshell
from ezdxf import disassemble

from falper_sorgente_linework import ROOT
from int1_highpoly_type_review import VIEWS, display_edge_sample, mesh_for_one_product, projected_raw_edges, svg_path


FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
PRODUCT = ROOT / "output/review/highpoly-types/miamisoft-e09"
SOURCE = PRODUCT / "official-source"
DOWNLOAD = SOURCE / "official-download"
ARCHIVE = DOWNLOAD / "Baxter_MiamiSoft_Sofa_2D_3D.zip"
DWG = DOWNLOAD / "Miami Soft_Divano_Abaco.dwg"
DXF_MEMBER = "MIAMI SOFT/2D/Miami Soft_Divano_Abaco.dwg"
OFFICIAL_URL = "https://dam.baxter.it/asset/b9c9485c-9716-422d-9c1f-7851c7bc9530/Baxter_MiamiSoft_Sofa_2D_3D.zip"
PRODUCT_URL = "https://www.baxter.it/gb/prodotti/miami-soft-divani-e-poltrone"
API_URL = "https://www.baxter.it/api/catalog/get_pdp_drawings?baxterId=MIAMSO&outdoor=false&indoor=true&open_air=false&settore=divani-e-poltrone"
ZIP_SHA256 = "f422e836ed7d666616340d696af93002fac6c9deff5c9e4a714ea9b8e1cbafd8"
DWG_SHA256 = "efc3315feed0a2dac8938aebbd7336856dabba6721395e4520783c39ba0cc51a"
GLOBAL_ID = "3osoWufdD1mhDDAM6lcix4"
ALLOWED_LAYERS = {"_ARREDO", "_CUSCINO SEDUTA", "_RULLO", "_PIEDINI"}
VIEW_BOXES = {
    "plan": [13700.0, 25000.0, 15900.0, 27400.0],
    "front": [13700.0, 27800.0, 15900.0, 29200.0],
    "side": [11500.0, 27800.0, 13700.0, 29200.0],
}
EXPECTED_PATH_COUNTS = {"plan": 64, "front": 130, "side": 85}
PINNED_PROXY_REVISION = "139d9ec178889e4f6bc0b81095b9bc79f33b3b79"
PINNED_PROXY_CANDIDATE = "output/review/highpoly-types/miamisoft-e09/candidate-representations.json"
PINNED_PROXY_CANDIDATE_SHA256 = "bfd6da993c754be473ee0296bd356ba7597c2918fac1e4372445268a3b90f305"
PINNED_PLAN_PROXY_PATH_COUNT = 25
BLUE = "#1677c8"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def pinned_geometry_derived_plan_proxy() -> list:
    result = subprocess.run(
        ["git", "show", f"{PINNED_PROXY_REVISION}:{PINNED_PROXY_CANDIDATE}"],
        cwd=ROOT,
        capture_output=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"cannot read pinned pre-DWG candidate: {result.stderr.decode(errors='replace')}")
    if hashlib.sha256(result.stdout).hexdigest() != PINNED_PROXY_CANDIDATE_SHA256:
        raise RuntimeError("pinned pre-DWG candidate hash mismatch")
    candidate = json.loads(result.stdout)
    if candidate.get("source_kind") != "geometry_derived_simplified_proxy":
        raise RuntimeError("pinned Plan comparison is not the geometry-derived candidate")
    paths = candidate["views"]["plan"]["proxy_paths_mm"]
    if len(paths) != PINNED_PLAN_PROXY_PATH_COUNT:
        raise RuntimeError("pinned geometry-derived Plan proxy path count drifted")
    return paths


def rounded_paths(paths: list) -> list:
    return [[[round(float(x), 6), round(float(y), 6)] for x, y in path] for path in paths]


def path_bounds(paths: list) -> dict:
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return {
        "minimum": [round(value, 6) for value in minimum],
        "maximum": [round(value, 6) for value in maximum],
        "size": [round(maximum[axis] - minimum[axis], 6) for axis in range(2)],
        "centre": [round((minimum[axis] + maximum[axis]) / 2.0, 6) for axis in range(2)],
    }


def convert_dwg() -> Path:
    temporary = Path(tempfile.mkdtemp(prefix="miamisoft-e09-native-dwg-"))
    result = subprocess.run(["dwg2dxf", str(DWG)], cwd=temporary, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"dwg2dxf failed: {result.stderr[-2000:]}")
    candidates = sorted(temporary.glob("*.dxf"))
    if len(candidates) != 1:
        raise RuntimeError(f"expected one DXF, found {candidates}")
    return candidates[0]


def extract_e09_paths(dxf: Path) -> tuple[dict, dict]:
    document = ezdxf.readfile(dxf)
    modelspace = document.modelspace()
    selected = {view: [] for view in VIEW_BOXES}
    records = {view: [] for view in VIEW_BOXES}
    for entity in modelspace:
        layer = entity.dxf.layer
        if layer not in ALLOWED_LAYERS:
            continue
        try:
            primitives = list(disassemble.to_primitives(disassemble.recursive_decompose([entity])))
        except Exception:
            continue
        for primitive in primitives:
            try:
                vertices = list(primitive.vertices(distance=0.5))
            except TypeError:
                vertices = list(primitive.vertices())
            path = [(float(vertex.x), float(vertex.y)) for vertex in vertices]
            if len(path) < 2:
                continue
            centre = [sum(point[axis] for point in path) / len(path) for axis in range(2)]
            for view, box in VIEW_BOXES.items():
                if box[0] <= centre[0] <= box[2] and box[1] <= centre[1] <= box[3]:
                    selected[view].append(path)
                    records[view].append({
                        "source_entity_type": entity.dxftype(),
                        "source_entity_handle": entity.dxf.get("handle"),
                        "source_layer": layer,
                    })
    for view, expected in EXPECTED_PATH_COUNTS.items():
        if len(selected[view]) != expected:
            raise RuntimeError(f"{view}: expected {expected} E09 paths, found {len(selected[view])}")
    return {view: rounded_paths(paths) for view, paths in selected.items()}, records


def align_reference(view: str, paths: list, proxy: list) -> tuple[list, dict]:
    source = [[list(point) for point in path] for path in paths]
    before = path_bounds(source)
    reflect_x = view == "side"
    if reflect_x:
        centre_x = before["centre"][0]
        source = [[[round(2.0 * centre_x - x, 6), y] for x, y in path] for path in source]
    after_reflection = path_bounds(source)
    target = path_bounds(proxy)
    translate_x = target["centre"][0] - after_reflection["centre"][0]
    if view == "plan":
        translate_y = target["centre"][1] - after_reflection["centre"][1]
        anchor = "projection_centre_to_projection_centre"
    else:
        translate_y = target["minimum"][1] - after_reflection["minimum"][1]
        anchor = "horizontal_centre_and_finished_bottom"
    aligned = rounded_paths([
        [(x + translate_x, y + translate_y) for x, y in path]
        for path in source
    ])
    aligned_bounds = path_bounds(aligned)
    residual = [
        round(aligned_bounds["size"][axis] - target["size"][axis], 6)
        for axis in range(2)
    ]
    return aligned, {
        "mode": "view_direction_reflection_and_translation_only" if reflect_x else "translation_only",
        "view_direction_reflection_x": reflect_x,
        "reflection_reason": "Official left-looking side is reflected to the project right-side review direction; this does not alter size." if reflect_x else None,
        "translation_mm": [round(translate_x, 6), round(translate_y, 6)],
        "anchor": anchor,
        "uniform_scale": 1.0,
        "anisotropic_scale_used": False,
        "source_geometry_deformed": False,
        "source_bounds_before_alignment_mm": before,
        "source_bounds_after_view_direction_mm": after_reflection,
        "aligned_official_bounds_mm": aligned_bounds,
        "project_proxy_bounds_mm": target,
        "official_minus_proxy_envelope_mm": residual,
    }


def render_svg(view: str, raw_edges: list, proxy: list, official: list, alignment: dict) -> str:
    width, height = 1500, 1000
    plot_x, plot_y, plot_w, plot_h = 55, 185, 1050, 750
    points = [point for edge in raw_edges for point in edge]
    points.extend(point for path in proxy for point in path)
    points.extend(point for path in official for point in path)
    minimum_x = min(point[0] for point in points)
    maximum_x = max(point[0] for point in points)
    minimum_y = min(point[1] for point in points)
    maximum_y = max(point[1] for point in points)
    padding = max(maximum_x - minimum_x, maximum_y - minimum_y) * 0.06
    minimum_x -= padding
    maximum_x += padding
    minimum_y -= padding
    maximum_y += padding
    scale = min(plot_w / (maximum_x - minimum_x), plot_h / (maximum_y - minimum_y))

    def transform(point):
        return plot_x + (point[0] - minimum_x) * scale, plot_y + plot_h - (point[1] - minimum_y) * scale

    raw_path = svg_path([[start, end] for start, end in raw_edges], transform)
    proxy_path = svg_path(proxy, transform, close=True)
    official_path = svg_path(official, transform)
    official_size = alignment["aligned_official_bounds_mm"]["size"]
    proxy_size = alignment["project_proxy_bounds_mm"]["size"]
    residual = alignment["official_minus_proxy_envelope_mm"]
    side_note = "Side uses view-direction X reflection only; size remains 1:1." if view == "side" else "No reflection; official paths use translation only."
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="1500" height="1000" fill="#fbfaf7"/>
<text x="55" y="54" font-family="Arial,sans-serif" font-size="29" font-weight="700" fill="#1f2d3d">Baxter Miami Soft E09 dx/r · 130×108 h70/80 cm</text>
<text x="55" y="92" font-family="Arial,sans-serif" font-size="19" fill="#41566d">{VIEWS[view]["label"]} · project representative {GLOBAL_ID}</text>
<text x="55" y="126" font-family="Arial,sans-serif" font-size="16" fill="#68798a">Grey = actual IFC Body · Black = prior high-poly-derived comparison · Blue = exact official E09 native 2D DWG candidate</text>
<text x="55" y="154" font-family="Arial,sans-serif" font-size="15" fill="{BLUE}">Solid blue · 1:1 · no scaling · pending individual approval · no IFC write</text>
<rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" rx="10" fill="#fff" stroke="#cad2d9" stroke-width="2"/>
<path class="actual-ifc-body" data-source-kind="project_ifc_body" d="{raw_path}" fill="none" stroke="#87929c" stroke-width="0.5" stroke-opacity="0.28" vector-effect="non-scaling-stroke"/>
<path class="geometry-derived-simplified-proxy comparison-only" data-source-kind="geometry_derived_simplified_proxy" d="{proxy_path}" fill="none" stroke="#111820" stroke-width="2.1" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<path class="official-candidate native-dwg" data-source-kind="native_dwg" data-source-role="pending_candidate" data-source-scaled="false" d="{official_path}" fill="none" stroke="{BLUE}" stroke-width="2.3" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>
<text x="1140" y="220" font-family="Arial,sans-serif" font-size="19" font-weight="700" fill="#1f2d3d">Official selection</text>
<text x="1140" y="255" font-family="Arial,sans-serif" font-size="15" fill="#41566d">E09 · TERMINAL MODULE R</text>
<text x="1140" y="283" font-family="Arial,sans-serif" font-size="15" fill="#41566d">width 130 · depth 108</text>
<text x="1140" y="311" font-family="Arial,sans-serif" font-size="15" fill="#41566d">height 70/80 · seat h. 40</text>
<text x="1140" y="360" font-family="Arial,sans-serif" font-size="19" font-weight="700" fill="#1f2d3d">Mechanical fit</text>
<text x="1140" y="395" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Official envelope: {official_size[0]:.1f} × {official_size[1]:.1f} mm</text>
<text x="1140" y="423" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Black proxy: {proxy_size[0]:.1f} × {proxy_size[1]:.1f} mm</text>
<text x="1140" y="451" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Residual: {residual[0]:+.1f} × {residual[1]:+.1f} mm</text>
<foreignObject x="1140" y="490" width="315" height="95"><div xmlns="http://www.w3.org/1999/xhtml" style="font:14px Arial;color:#41566d;line-height:1.4">{html.escape(side_note)}</div></foreignObject>
<text x="1140" y="625" font-family="Arial,sans-serif" font-size="19" font-weight="700" fill="#1f2d3d">Write gate</text>
<text x="1140" y="660" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Review: pending</text>
<text x="1140" y="688" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Derived IFC allowed: false</text>
<text x="1140" y="716" font-family="Arial,sans-serif" font-size="14" fill="#41566d">Formal IFC write: none</text>
</svg>'''


def write_index() -> None:
    cards = "".join(
        f'<article><h2>{view.title()}</h2><a href="{view}.svg"><img src="{view}.svg"></a></article>'
        for view in ("plan", "front", "side")
    )
    bonsai = "".join(
        f'<article><h2>Bonsai {view.title()}</h2><a href="bonsai-camera-{filename}.png"><img src="bonsai-camera-{filename}.png"></a></article>'
        for view, filename in (("plan", "plan"), ("front", "front-elevation"), ("side", "side-elevation"), ("iso", "iso"))
    )
    (PRODUCT / "index.html").write_text(
        f'''<!doctype html><html><meta charset="utf-8"><title>Baxter Miami Soft E09 official DWG review</title><style>body{{font:16px Arial;margin:30px;background:#f5f3ef;color:#1f2d3d}}main{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:18px}}article{{background:#fff;padding:14px;border-radius:10px}}img{{width:100%}}a{{margin-right:18px}}</style><h1>Baxter Miami Soft E09 dx/r · 130×108 h70/80 cm</h1><p>Grey = actual IFC Body; black = prior <strong>基于原始高模几何生成的简化图纸表达</strong> for comparison; solid blue = mechanically selected E09 native Plan/Front/Side from Baxter's official DWG. Blue stays 1:1 and receives translation only; Side additionally uses a documented view-direction reflection. Review remains pending and no IFC was written.</p><nav><a href="manifest.json">Manifest</a><a href="candidate-representations.json">Candidate</a><a href="official-dwg-review-reference.json">Mechanical record</a><a href="official-source/source-access-record.json">Source record</a><a href="official-source/autocad-miami-soft-e09-three-view.png">Official DWG E09 screenshot</a><a href="official-source/autocad-miami-soft-full-modelspace.png">Official DWG full modelspace</a><a href="official-source/official-download/archive-inventory.json">ZIP inventory</a><a href="{OFFICIAL_URL}">Official ZIP</a><a href="project-context-furniture-plan-review.svg">Project plan</a><a href="project-context-r20-side-elevation-review.svg">Project elevation</a><a href="bonsai-review-manifest.json">Bonsai saved-camera evidence</a><a href="{PRODUCT_URL}">Official product page</a></nav><h2>Exact E09 native-DWG candidate</h2><main>{cards}</main><h2>Actual Bonsai saved-camera renders of the current IFC Body</h2><main>{bonsai}</main><h2>Project drawing context</h2><main><article><h2>Furniture Plan</h2><a href="project-context-furniture-plan-review.svg"><img src="project-context-furniture-plan-review-preview.png"></a></article><article><h2>R20 side elevation</h2><a href="project-context-r20-side-elevation-review.svg"><img src="project-context-r20-side-elevation-review-preview.png"></a></article></main><h2>Official DWG evidence</h2><main><article><h2>E09 three-view</h2><a href="official-source/autocad-miami-soft-e09-three-view.png"><img src="official-source/autocad-miami-soft-e09-three-view.png"></a></article><article><h2>Full family modelspace</h2><a href="official-source/autocad-miami-soft-full-modelspace.png"><img src="official-source/autocad-miami-soft-full-modelspace.png"></a></article></main><p>Formal IFC SHA-256: <code>{FORMAL_SHA256}</code>. No derived or authoritative IFC write.</p></html>''',
        encoding="utf-8",
    )


def _svg_path_data(svg_text: str, class_name: str) -> tuple[re.Match, str]:
    pattern = re.compile(rf'(<path class="{re.escape(class_name)}"[^>]*?\sd=")([^"]*)(")')
    match = pattern.search(svg_text)
    if match is None:
        raise RuntimeError(f"SVG layer not found: {class_name}")
    return match, match.group(2)


def _fit_existing_svg_transform(source_paths: list, rendered_path_data: str):
    source_points = [point for path in source_paths for point in path]
    rendered_points = [
        (float(x), float(y))
        for x, y in re.findall(r'[ML]\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)', rendered_path_data)
    ]
    if len(source_points) != len(rendered_points):
        raise RuntimeError("official SVG/source point correspondence drifted")

    def fit(source_values, rendered_values):
        source_mean = sum(source_values) / len(source_values)
        rendered_mean = sum(rendered_values) / len(rendered_values)
        variance = sum((value - source_mean) ** 2 for value in source_values)
        slope = sum(
            (source - source_mean) * (rendered - rendered_mean)
            for source, rendered in zip(source_values, rendered_values)
        ) / variance
        intercept = rendered_mean - slope * source_mean
        return slope, intercept

    scale_x, translate_x = fit(
        [point[0] for point in source_points],
        [point[0] for point in rendered_points],
    )
    scale_y, translate_y = fit(
        [point[1] for point in source_points],
        [point[1] for point in rendered_points],
    )
    if abs(scale_x + scale_y) > 2e-6:
        raise RuntimeError("existing official SVG transform is not uniform")

    def transform(point):
        return translate_x + scale_x * point[0], translate_y + scale_y * point[1]

    regenerated = svg_path(source_paths, transform)
    if regenerated != rendered_path_data:
        raise RuntimeError("fitted transform does not byte-reproduce the existing official blue path")
    return transform


def _rebuild_existing_plan_transform(candidate: dict, rendered_official_path_data: str):
    model = ifcopenshell.open(str(FORMAL_IFC))
    _, vertices, faces = mesh_for_one_product(model, GLOBAL_ID)
    raw_edges = display_edge_sample(
        projected_raw_edges(vertices, faces, VIEWS["plan"]["axes"]),
        maximum=1000,
    )
    official_paths = candidate["views"]["plan"]["official_cad_paths_mm"]
    current_black_paths = candidate["views"]["plan"]["comparison_geometry_derived_proxy_paths_mm"]
    points = [point for edge in raw_edges for point in edge]
    points.extend(point for path in current_black_paths for point in path)
    points.extend(point for path in official_paths for point in path)
    minimum_x = min(point[0] for point in points)
    maximum_x = max(point[0] for point in points)
    minimum_y = min(point[1] for point in points)
    maximum_y = max(point[1] for point in points)
    padding = max(maximum_x - minimum_x, maximum_y - minimum_y) * 0.06
    minimum_x -= padding
    maximum_x += padding
    minimum_y -= padding
    maximum_y += padding
    scale = min(1050 / (maximum_x - minimum_x), 750 / (maximum_y - minimum_y))

    def transform(point):
        return 55 + (point[0] - minimum_x) * scale, 185 + 750 - (point[1] - minimum_y) * scale

    if svg_path(official_paths, transform) != rendered_official_path_data:
        raise RuntimeError("exact review viewport does not byte-reproduce the existing official blue path")
    return transform


def repair_plan_proxy_overlay_only() -> None:
    """Replace only the Plan black path; preserve its grey and blue path bytes."""
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash changed before Plan overlay repair")
    plan_svg = PRODUCT / "plan.svg"
    front_svg = PRODUCT / "front.svg"
    side_svg = PRODUCT / "side.svg"
    before_file_hashes = {
        "plan": sha256(plan_svg),
        "front": sha256(front_svg),
        "side": sha256(side_svg),
    }
    svg_text = plan_svg.read_text(encoding="utf-8")
    _, actual_path_data = _svg_path_data(svg_text, "actual-ifc-body")
    _, official_path_data = _svg_path_data(svg_text, "official-candidate native-dwg")
    black_match, _ = _svg_path_data(svg_text, "geometry-derived-simplified-proxy comparison-only")

    candidate_path = PRODUCT / "candidate-representations.json"
    candidate = load_json(candidate_path)
    official_paths = candidate["views"]["plan"]["official_cad_paths_mm"]
    if len(official_paths) != EXPECTED_PATH_COUNTS["plan"]:
        raise RuntimeError("official Plan candidate path count drifted")
    official_semantics_before = hashlib.sha256(json.dumps({
        "paths": official_paths,
        "alignment": candidate["views"]["plan"]["alignment"],
    }, separators=(",", ":"), sort_keys=True).encode()).hexdigest()
    transform = _rebuild_existing_plan_transform(candidate, official_path_data)
    proxy_paths = pinned_geometry_derived_plan_proxy()
    black_path_data = svg_path(proxy_paths, transform, close=True)
    new_svg_text = svg_text[:black_match.start(2)] + black_path_data + svg_text[black_match.end(2):]
    plan_svg.write_text(new_svg_text, encoding="utf-8")

    _, actual_after = _svg_path_data(new_svg_text, "actual-ifc-body")
    _, official_after = _svg_path_data(new_svg_text, "official-candidate native-dwg")
    if actual_after != actual_path_data or official_after != official_path_data:
        raise RuntimeError("Plan repair changed the grey or blue path bytes")
    if sha256(front_svg) != before_file_hashes["front"] or sha256(side_svg) != before_file_hashes["side"]:
        raise RuntimeError("Plan-only repair changed Front or Side")

    candidate["views"]["plan"]["comparison_geometry_derived_proxy_paths_mm"] = proxy_paths
    candidate["views"]["plan"]["comparison_geometry_derived_proxy_provenance"] = {
        "source_kind": "geometry_derived_simplified_proxy",
        "source_revision": PINNED_PROXY_REVISION,
        "source_candidate": PINNED_PROXY_CANDIDATE,
        "source_candidate_sha256": PINNED_PROXY_CANDIDATE_SHA256,
        "path_count": PINNED_PLAN_PROXY_PATH_COUNT,
        "recomputed": False,
        "official_blue_paths_reused_as_black": False,
    }
    official_semantics_after = hashlib.sha256(json.dumps({
        "paths": candidate["views"]["plan"]["official_cad_paths_mm"],
        "alignment": candidate["views"]["plan"]["alignment"],
    }, separators=(",", ":"), sort_keys=True).encode()).hexdigest()
    if official_semantics_after != official_semantics_before:
        raise RuntimeError("official Plan paths or 1:1 alignment semantics changed")

    reference_path = PRODUCT / "official-dwg-review-reference.json"
    reference = load_json(reference_path)
    reference["plan_geometry_derived_proxy_overlay"] = {
        "source_kind": "geometry_derived_simplified_proxy",
        "source_revision": PINNED_PROXY_REVISION,
        "source_candidate": PINNED_PROXY_CANDIDATE,
        "source_candidate_sha256": PINNED_PROXY_CANDIDATE_SHA256,
        "path_count": PINNED_PLAN_PROXY_PATH_COUNT,
        "role": "black comparison layer only",
        "recomputed": False,
        "official_blue_paths_reused_as_black": False,
        "actual_body_path_bytes_unchanged": True,
        "official_blue_path_bytes_unchanged": True,
    }
    write_json(reference_path, reference)
    candidate["official_dwg_review_reference_sha256"] = sha256(reference_path)
    write_json(candidate_path, candidate)

    manifest_path = PRODUCT / "manifest.json"
    manifest = load_json(manifest_path)
    manifest["generated_at"] = datetime.now(timezone.utc).isoformat()
    manifest["candidate_representations_sha256"] = sha256(candidate_path)
    manifest["official_reference"]["mechanical_selection_record_sha256"] = sha256(reference_path)
    plan_record = next(record for record in manifest["views"] if record["view"] == "plan")
    plan_record["svg_sha256"] = sha256(plan_svg)
    plan_record["comparison_proxy_path_count"] = PINNED_PLAN_PROXY_PATH_COUNT
    plan_record["layer_path_counts"] = {
        "actual_ifc_body_displayed_edges": plan_record["displayed_raw_edge_count"],
        "geometry_derived_simplified_proxy": PINNED_PLAN_PROXY_PATH_COUNT,
        "official_native_dwg": EXPECTED_PATH_COUNTS["plan"],
    }
    manifest["review_status"] = "visual_review_pending"
    manifest["approved_for_drawing_ifc"] = False
    manifest["derived_ifc_write_allowed"] = False
    manifest["formal_ifc_write"] = "not performed"
    manifest["formal_ifc_sha256_after_generation"] = sha256(FORMAL_IFC)
    manifest["formal_ifc_bytes_unchanged"] = True
    write_json(manifest_path, manifest)

    approval_path = ROOT / "pipeline/decisions/miamisoft-e09-drawing-approval.json"
    approval = load_json(approval_path)
    approval["candidate_manifest_sha256"] = sha256(manifest_path)
    approval["status"] = "pending"
    approval["approved_views"] = []
    approval["derived_ifc_write_allowed"] = False
    approval["formal_authoritative_ifc_write_allowed"] = False
    write_json(approval_path, approval)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash changed during Plan overlay repair")
    print(json.dumps({
        "mode": "plan_proxy_overlay_only",
        "plan_svg_before_sha256": before_file_hashes["plan"],
        "plan_svg_after_sha256": sha256(plan_svg),
        "front_svg_sha256_unchanged": sha256(front_svg),
        "side_svg_sha256_unchanged": sha256(side_svg),
        "plan_layer_path_counts": {
            "actual_ifc_body_displayed_edges": plan_record["displayed_raw_edge_count"],
            "geometry_derived_simplified_proxy": len(proxy_paths),
            "official_native_dwg": len(official_paths),
        },
        "actual_body_path_data_sha256": hashlib.sha256(actual_after.encode()).hexdigest(),
        "official_blue_path_data_sha256": hashlib.sha256(official_after.encode()).hexdigest(),
        "official_paths_and_alignment_semantics_sha256": official_semantics_after,
        "formal_ifc_sha256": sha256(FORMAL_IFC),
        "derived_ifc_written": False,
    }, indent=2))


def main() -> None:
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash changed before review generation")
    if sha256(ARCHIVE) != ZIP_SHA256 or sha256(DWG) != DWG_SHA256:
        raise RuntimeError("official Baxter source hash mismatch")

    inventory = []
    with zipfile.ZipFile(ARCHIVE) as archive:
        for item in archive.infolist():
            inventory.append({
                "path": item.filename,
                "uncompressed_size_bytes": item.file_size,
                "compressed_size_bytes": item.compress_size,
                "crc32_hex": f"{item.CRC:08x}",
                "is_directory": item.is_dir(),
            })
    archive_inventory = DOWNLOAD / "archive-inventory.json"
    write_json(archive_inventory, {
        "schema_version": 1,
        "official_download_url": OFFICIAL_URL,
        "archive": relative(ARCHIVE),
        "archive_sha256": ZIP_SHA256,
        "member_count": len(inventory),
        "selected_native_2d_member": DXF_MEMBER,
        "selected_native_2d_local_path": relative(DWG),
        "selected_native_2d_sha256": DWG_SHA256,
        "members": inventory,
    })

    dxf = convert_dwg()
    source_paths, entity_records = extract_e09_paths(dxf)
    old_candidate = load_json(PRODUCT / "candidate-representations.json")
    proxy_paths = {
        view: old_candidate["views"][view].get(
            "comparison_geometry_derived_proxy_paths_mm",
            old_candidate["views"][view].get("proxy_paths_mm"),
        )
        for view in VIEWS
    }
    if any(paths is None for paths in proxy_paths.values()):
        raise RuntimeError("the prior geometry-derived comparison paths are missing")
    proxy_paths["plan"] = pinned_geometry_derived_plan_proxy()

    model = ifcopenshell.open(str(FORMAL_IFC))
    product, vertices, faces = mesh_for_one_product(model, GLOBAL_ID)
    aligned_paths = {}
    alignment_records = {}
    view_records = []
    for view, view_data in VIEWS.items():
        aligned, alignment = align_reference(view, source_paths[view], proxy_paths[view])
        raw = display_edge_sample(projected_raw_edges(vertices, faces, view_data["axes"]), maximum=1000)
        svg = PRODUCT / f"{view}.svg"
        svg.write_text(render_svg(view, raw, proxy_paths[view], aligned, alignment), encoding="utf-8")
        aligned_paths[view] = aligned
        alignment_records[view] = alignment
        view_records.append({
            "view": view,
            "svg": relative(svg),
            "svg_sha256": sha256(svg),
            "projection_axes": list(view_data["axes"]),
            "raw_edge_count": len(projected_raw_edges(vertices, faces, view_data["axes"])),
            "displayed_raw_edge_count": len(raw),
            "candidate_path_count": len(aligned),
            "comparison_proxy_path_count": len(proxy_paths[view]),
            "drawing_line_source_kind": "native_dwg",
            "drawing_line_source_label_zh": "Baxter 官方 Miami Soft E09 原生二维 DWG 候选",
            "official_cad_path_count": len(aligned),
            "blue_line_present": True,
            "uniform_scale": 1.0,
            "view_direction_reflection_x": alignment["view_direction_reflection_x"],
        })

    reference_path = PRODUCT / "official-dwg-review-reference.json"
    reference = {
        "schema_version": 1,
        "profile_key": "miamisoft-e09",
        "manufacturer": "Baxter",
        "family": "Miami Soft",
        "model_code": "E09",
        "variant": "dx/r - right terminal module - 130 x 108 h70/80 cm",
        "source_kind": "native_dwg",
        "source_geometry_kind": "native_dwg_2d_linework",
        "official_product_page": PRODUCT_URL,
        "official_api_url": API_URL,
        "official_download_url": OFFICIAL_URL,
        "restricted_asset_name": "Baxter_MiamiSoft_Sofa_2D_3D",
        "known_catalog_asset_alias": "Baxter_MiamiSoft_divano_2D_3D",
        "restricted_asset_database_id": "B9C9485C-9716-422D-9C1F7851C7BC9530",
        "source_archive": relative(ARCHIVE),
        "source_archive_sha256": ZIP_SHA256,
        "source_archive_inventory": relative(archive_inventory),
        "source_archive_inventory_sha256": sha256(archive_inventory),
        "source_dwg": relative(DWG),
        "source_dwg_sha256": DWG_SHA256,
        "source_dwg_format": "AutoCAD 2018/2019/2020 DWG",
        "autocad_source_screenshots": [
            relative(SOURCE / "autocad-miami-soft-e09-three-view.png"),
            relative(SOURCE / "autocad-miami-soft-full-modelspace.png"),
        ],
        "autocad_source_screenshot_sha256": {
            "e09_three_view": sha256(SOURCE / "autocad-miami-soft-e09-three-view.png"),
            "full_modelspace": sha256(SOURCE / "autocad-miami-soft-full-modelspace.png"),
        },
        "mechanical_extraction": {
            "converter": "LibreDWG dwg2dxf",
            "reader": f"ezdxf {ezdxf.__version__}",
            "curve_flattening_distance_mm": 0.5,
            "selection_method": "entity-centre inside exact E09 view rectangles after AutoCAD visual localization",
            "geometry_layers_included": sorted(ALLOWED_LAYERS),
            "annotation_layers_excluded": ["_QUOTE", "TESTO", "LINEA"],
            "view_selection_boxes_dwg_mm": VIEW_BOXES,
            "view_path_counts": EXPECTED_PATH_COUNTS,
            "view_bounds_dwg_mm": {view: path_bounds(paths) for view, paths in source_paths.items()},
            "per_path_source_entity_records": entity_records,
        },
        "alignment": alignment_records,
        "official_identity_text_mechanically_located": {
            "insertion_point_dwg_mm": [14169.291, 24660.779],
            "text": ["E09", "TERMINAL MODULE R", "width 130", "depth 108", "height 70/80", "seat h. 40"],
        },
        "exact_project_configuration_match": True,
        "third_party_cad_used": False,
        "review_status": "visual_review_pending",
        "derived_ifc_write_allowed": False,
        "formal_authoritative_ifc_write_allowed": False,
        "pass": True,
    }
    write_json(reference_path, reference)

    candidate = {
        "schema_version": 1,
        "profile_key": "miamisoft-e09",
        "representative_global_id": GLOBAL_ID,
        "ifc_type_name": "MiamiSoft E09",
        "article_number": "Baxter Miami Soft E09 dx/r",
        "units": "mm",
        "source_kind": "native_dwg",
        "source_label_zh": "Baxter 官方 Miami Soft E09 原生二维 DWG 候选",
        "source_label_en": "exact E09 native 2D linework from Baxter's official Miami Soft DWG",
        "official_cad_used": True,
        "third_party_cad_used": False,
        "official_dwg_review_reference": relative(reference_path),
        "official_dwg_review_reference_sha256": sha256(reference_path),
        "derived_ifc_write_allowed": False,
        "formal_ifc_write_allowed": False,
        "review_status": "visual_review_pending",
        "views": {
            view: {
                "projection_axes": list(VIEWS[view]["axes"]),
                "candidate_paths_mm": aligned_paths[view],
                "official_cad_paths_mm": aligned_paths[view],
                "comparison_geometry_derived_proxy_paths_mm": proxy_paths[view],
                "alignment": alignment_records[view],
                "candidate_path_count": len(aligned_paths[view]),
                "source_kind": "native_dwg",
                "source_scaled": False,
            }
            for view in VIEWS
        },
    }
    candidate_path = PRODUCT / "candidate-representations.json"
    write_json(candidate_path, candidate)

    access_path = SOURCE / "source-access-record.json"
    access = load_json(access_path)
    access["research_date"] = datetime.now(timezone.utc).date().isoformat()
    access["official_product_cad"] = {
        "published_access_surface": "Baxter product drawings API and official DAM",
        "restricted_asset_name": "Baxter_MiamiSoft_Sofa_2D_3D",
        "known_catalog_asset_alias": "Baxter_MiamiSoft_divano_2D_3D",
        "restricted_asset_database_id": "B9C9485C-9716-422D-9C1F7851C7BC9530",
        "official_api_url": API_URL,
        "official_download_url": OFFICIAL_URL,
        "authentication_required_on_user_download_page": True,
        "official_asset_url_returned_by_product_drawings_api": True,
        "acquired": True,
        "exact_project_configuration_match": True,
        "selected_configuration": "E09 dx/r right terminal module 130 x 108 h70/80 cm",
        "family_archive": relative(ARCHIVE),
        "family_archive_sha256": ZIP_SHA256,
        "archive_inventory": relative(archive_inventory),
        "archive_inventory_sha256": sha256(archive_inventory),
        "local_cad_files": [{
            "path": relative(DWG),
            "sha256": DWG_SHA256,
            "archive_member": DXF_MEMBER,
            "format": "AutoCAD 2018/2019/2020 DWG",
        }],
        "selected_native_views": {view: {"path_count": len(paths), "bounds_dwg_mm": path_bounds(paths)} for view, paths in source_paths.items()},
        "third_party_cad_used": False,
        "official_vector_pdf_archived": True,
        "official_measurement_svg_archived": True,
        "official_vector_references_used_as_cad_geometry": False,
    }
    access["drawing_geometry_source"] = {
        "source_kind": "native_dwg",
        "source_label_zh": "Baxter 官方 Miami Soft E09 原生二维 DWG 候选",
        "source_label_en": "exact E09 native 2D linework from Baxter's official Miami Soft DWG",
        "official_cad_used": True,
        "third_party_cad_used": False,
        "review_status": "visual_review_pending",
        "official_dwg_review_reference": relative(reference_path),
    }
    access["scope"] = "Exact Baxter Miami Soft E09 dx/r native Plan/Front/Side selected from the official family DWG; 1:1, translation only, with Side view-direction reflection; pending review and not yet written to IFC."
    access["pass"] = True
    write_json(access_path, access)

    old_manifest = load_json(PRODUCT / "manifest.json")
    manifest = {
        **old_manifest,
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/miamisoft_e09_official_dwg_review.py",
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(FORMAL_IFC),
        "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
        "source_kind": "native_dwg",
        "source_label_zh": "Baxter 官方 Miami Soft E09 原生二维 DWG 候选",
        "drawing_source": {
            "source_kind": "native_dwg",
            "source_label_zh": "Baxter 官方 Miami Soft E09 原生二维 DWG 候选",
            "source_label_en": "exact E09 native 2D linework from Baxter's official Miami Soft DWG",
            "official_cad_used": True,
            "third_party_cad_used": False,
            "official_download_url": OFFICIAL_URL,
            "official_archive": relative(ARCHIVE),
            "official_archive_sha256": ZIP_SHA256,
            "official_dwg": relative(DWG),
            "official_dwg_sha256": DWG_SHA256,
            "official_source_access_record": relative(access_path),
            "official_source_access_record_sha256": sha256(access_path),
            "official_product_cad_status": "official_family_archive_acquired_exact_E09_native_2d_selected_pending_approval",
        },
        "official_reference": {
            "source_id": "BAXTER-MIAMI-SOFT-E09-NATIVE-DWG-001",
            "source_kind": "native_dwg",
            "manufacturer": "Baxter",
            "family": "Miami Soft",
            "model_code": "E09",
            "variant": "dx/r right terminal module, 130 x 108 h70/80 cm",
            "scope": access["scope"],
            "product_page": PRODUCT_URL,
            "official_download_url": OFFICIAL_URL,
            "source_archive": relative(ARCHIVE),
            "source_archive_sha256": ZIP_SHA256,
            "source_dwg": relative(DWG),
            "source_dwg_sha256": DWG_SHA256,
            "source_access_record": relative(access_path),
            "mechanical_selection_record": relative(reference_path),
            "mechanical_selection_record_sha256": sha256(reference_path),
            "official_product_cad_status": "acquired_exact_E09_native_2d_selected_pending_approval",
            "official_product_cad_used": True,
        },
        "official_identity_evidence_only": False,
        "official_cad_acquired": True,
        "official_cad_exact_project_configuration_match": True,
        "official_cad_used": True,
        "third_party_cad_used": False,
        "blue_line_present": True,
        "official_source_access_record": relative(access_path),
        "official_source_access_record_sha256": sha256(access_path),
        "candidate_representations": relative(candidate_path),
        "candidate_representations_sha256": sha256(candidate_path),
        "review_status": "visual_review_pending",
        "approved_for_drawing_ifc": False,
        "derived_ifc_write_allowed": False,
        "formal_ifc_write": "not performed",
        "views": view_records,
        "pass": True,
    }
    manifest_path = PRODUCT / "manifest.json"
    write_json(manifest_path, manifest)

    approval_path = ROOT / "pipeline/decisions/miamisoft-e09-drawing-approval.json"
    write_json(approval_path, {
        "schema_version": 1,
        "profile_key": "miamisoft-e09",
        "ifc_type_name": "MiamiSoft E09",
        "candidate_manifest_sha256": sha256(manifest_path),
        "status": "pending",
        "reviewer": None,
        "review_date": None,
        "approved_views": [],
        "derived_ifc_write_allowed": False,
        "formal_authoritative_ifc_write_allowed": False,
        "scope": access["scope"],
        "approval_evidence": None,
    })

    write_index()
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash changed during review generation")
    print(json.dumps({
        "profile_key": "miamisoft-e09",
        "source_archive_sha256": ZIP_SHA256,
        "source_dwg_sha256": DWG_SHA256,
        "candidate_sha256": sha256(candidate_path),
        "manifest_sha256": sha256(manifest_path),
        "view_path_counts": EXPECTED_PATH_COUNTS,
        "formal_ifc_sha256": sha256(FORMAL_IFC),
        "derived_ifc_written": False,
    }, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan-proxy-repair-only", action="store_true")
    arguments = parser.parse_args()
    repair_plan_proxy_overlay_only() if arguments.plan_proxy_repair_only else main()
