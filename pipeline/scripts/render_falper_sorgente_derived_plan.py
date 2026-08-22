#!/usr/bin/env python3
"""Render the approved Falper WFB plan from the corrected derived IFC."""

from __future__ import annotations

import argparse
import html
from pathlib import Path

import ifcopenshell
import ifcopenshell.util.element

from falper_sorgente_linework import EXPECTED, ROOT, SCOPE, load_json, sha256


GLOBAL_ID = "350tdaubr8QP3Cu2YMQZIN"
BLUE = "#1677c8"
DEFAULT_IFC = (
    ROOT
    / "output/review/highpoly-types/falper-sorgente/"
    "Falper-Sorgente-WFB-derived-drawing-corrected.ifc"
)
DEFAULT_MANIFEST = (
    ROOT
    / "output/review/highpoly-types/falper-sorgente/derived-ifc-corrected-manifest.json"
)
DEFAULT_OUTPUT = (
    ROOT / "output/review/highpoly-types/falper-sorgente/plan-approved-corrected.svg"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_IFC)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    derived_path = args.input.resolve()
    manifest = load_json(args.manifest.resolve())
    if (
        manifest.get("pass") is not True
        or manifest.get("approval", {}).get("status") != "approved"
        or manifest.get("representation_path_counts", {}).get("plan") != 5
        or manifest.get("proxy_geometry_included") is not False
        or sha256(derived_path) != manifest.get("derived_ifc_sha256")
    ):
        raise RuntimeError("corrected derived IFC manifest gate failed")

    model = ifcopenshell.open(derived_path)
    product = model.by_guid(GLOBAL_ID)
    representation = next(
        (
            item
            for item in product.Representation.Representations
            if item.RepresentationIdentifier == "FalperWFBPlan"
        ),
        None,
    )
    if (
        representation is None
        or representation.RepresentationType != "GeometricCurveSet"
        or representation.ContextOfItems.TargetView != "PLAN_VIEW"
    ):
        raise RuntimeError("corrected FalperWFBPlan representation gate failed")
    paths = [
        [
            (float(point.Coordinates[0]), float(point.Coordinates[1]))
            for point in polyline.Points
        ]
        for curve_set in representation.Items
        for polyline in curve_set.Elements
    ]
    if len(paths) != 5 or any(path[0] != path[-1] for path in paths):
        raise RuntimeError("official native DWG plan path identity gate failed")

    source = ifcopenshell.util.element.get_pset(
        product, "Pset_FalperSorgenteDrawingSource"
    ) or {}
    if (
        source.get("SourceKind") != "native_dwg"
        or source.get("SourceDwgSha256") != EXPECTED["wfb_2d"]
        or source.get("SourcePdfSha256") != EXPECTED["pdf"]
        or source.get("EvidenceScope") != SCOPE
        or source.get("RepresentationGeometrySource") != "official_native_dwg_paths_mm"
        or source.get("ProxyGeometryIncluded") != "false"
        or source.get("NativeDwgPlanPathCount") != "5"
    ):
        raise RuntimeError("corrected plan source metadata gate failed")

    all_points = [point for path in paths for point in path]
    min_x = min(point[0] for point in all_points)
    max_x = max(point[0] for point in all_points)
    min_y = min(point[1] for point in all_points)
    max_y = max(point[1] for point in all_points)
    extent_x, extent_y = max_x - min_x, max_y - min_y
    panel_x, panel_y, panel_width, panel_height = 60.0, 155.0, 1080.0, 740.0
    drawing_size = 620.0
    scale = min(drawing_size / extent_x, drawing_size / extent_y)
    origin_x = panel_x + panel_width / 2 - (min_x + max_x) * scale / 2
    origin_y = panel_y + panel_height / 2 + (min_y + max_y) * scale / 2

    def svg_path(points: list[tuple[float, float]]) -> str:
        return " ".join(
            f"{'M' if index == 0 else 'L'} {origin_x + x * scale:.3f} {origin_y - y * scale:.3f}"
            for index, (x, y) in enumerate(points)
        )

    linework = "\n".join(f'      <path d="{svg_path(path)}"/>' for path in paths)
    reviewer = html.escape(str(source["Reviewer"]))
    review_date = html.escape(str(source["ReviewDate"]))
    dwg_hash = html.escape(source["SourceDwgSha256"])
    derived_hash = html.escape(manifest["derived_ifc_sha256"])
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="1080" viewBox="0 0 1200 1080">
  <rect width="1200" height="1080" fill="#f8f7f3"/>
  <text x="60" y="58" font-family="Arial, sans-serif" font-size="31" font-weight="700" fill="#243447">Falper Sorgente WFB · approved plan</text>
  <text x="60" y="96" font-family="Arial, sans-serif" font-size="18" fill="#4e6478">Derived IFC representation · isolated BS01 representative {GLOBAL_ID}</text>
  <rect x="922" y="36" width="218" height="64" rx="8" fill="#e9f4fc" stroke="{BLUE}" stroke-width="2"/>
  <text x="1031" y="63" text-anchor="middle" font-family="Arial, sans-serif" font-size="16" font-weight="700" fill="{BLUE}">APPROVED</text>
  <text x="1031" y="85" text-anchor="middle" font-family="Arial, sans-serif" font-size="13" fill="#42647f">DERIVED IFC</text>
  <rect x="{panel_x}" y="{panel_y}" width="{panel_width}" height="{panel_height}" rx="10" fill="#ffffff" stroke="#cbd4dc" stroke-width="2"/>
  <text x="88" y="199" font-family="Arial, sans-serif" font-size="21" font-weight="700" fill="#243447">PLAN / XY</text>
  <text x="1112" y="199" text-anchor="end" font-family="Arial, sans-serif" font-size="14" fill="#60758a">FalperWFBPlan · IFC GeometricCurveSet · 5 native DWG paths</text>
  <g class="official-wfb-native-dwg" data-source-kind="native_dwg" fill="none" stroke="{BLUE}" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke">
{linework}
  </g>
  <line x1="220" y1="854" x2="980" y2="854" stroke="#d7dee5" stroke-width="1"/>
  <text x="600" y="879" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" fill="#60758a">source extents {extent_x:.0f} × {extent_y:.0f} mm · official blue linework</text>
  <rect x="60" y="925" width="1080" height="113" rx="8" fill="#ffffff" stroke="#cbd4dc" stroke-width="2"/>
  <text x="84" y="956" font-family="Arial, sans-serif" font-size="15" font-weight="700" fill="#243447">SOURCE / APPROVAL</text>
  <text x="84" y="982" font-family="Arial, sans-serif" font-size="14" fill="#4e6478">Falper official Sorgente WFB native DWG · PDF vector cross-check passed · family reference, not project shop drawing</text>
  <text x="84" y="1007" font-family="Arial, sans-serif" font-size="13" fill="#60758a">DWG SHA-256 {dwg_hash}</text>
  <text x="84" y="1030" font-family="Arial, sans-serif" font-size="13" fill="#60758a">Derived IFC SHA-256 {derived_hash}</text>
  <text x="1116" y="1007" text-anchor="end" font-family="Arial, sans-serif" font-size="14" fill="#4e6478">reviewer {reviewer}</text>
  <text x="1116" y="1030" text-anchor="end" font-family="Arial, sans-serif" font-size="14" fill="#4e6478">review date {review_date}</text>
</svg>
"""
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
