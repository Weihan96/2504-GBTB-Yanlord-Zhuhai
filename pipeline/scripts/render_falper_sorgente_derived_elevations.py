#!/usr/bin/env python3
"""Render approved Falper WFB front and side elevations from the derived IFC."""

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
    ROOT / "output/review/highpoly-types/falper-sorgente/elevations-approved-corrected.svg"
)


def extract_paths(representation, axis: int) -> list[list[tuple[float, float]]]:
    if representation.RepresentationType != "GeometricCurveSet":
        raise RuntimeError("approved elevation is not a GeometricCurveSet")
    paths = []
    for curve_set in representation.Items:
        for polyline in curve_set.Elements:
            coordinates = [
                (float(point.Coordinates[axis]), float(point.Coordinates[2]))
                for point in polyline.Points
            ]
            if len(coordinates) >= 2:
                paths.append(coordinates)
    return paths


def panel_svg(
    label: str,
    representation_identifier: str,
    paths: list[list[tuple[float, float]]],
    panel_x: float,
) -> str:
    panel_y, panel_width, panel_height = 170.0, 710.0, 720.0
    all_points = [point for path in paths for point in path]
    min_x = min(point[0] for point in all_points)
    max_x = max(point[0] for point in all_points)
    min_z = min(point[1] for point in all_points)
    max_z = max(point[1] for point in all_points)
    width_mm = max_x - min_x
    height_mm = max_z - min_z
    drawing_width, drawing_height = 520.0, 610.0
    scale = min(drawing_width / width_mm, drawing_height / height_mm)
    origin_x = panel_x + panel_width / 2 - (min_x + max_x) * scale / 2
    origin_y = panel_y + 50 + (max_z * scale)

    def mapped_path(points: list[tuple[float, float]]) -> str:
        commands = []
        for index, (x, z) in enumerate(points):
            prefix = "M" if index == 0 else "L"
            commands.append(f"{prefix} {origin_x + x * scale:.3f} {origin_y - z * scale:.3f}")
        return " ".join(commands)

    linework = "\n".join(
        f'    <path d="{mapped_path(path)}"/>' for path in paths
    )
    return f"""
  <g class="elevation-panel" data-representation="{representation_identifier}">
    <rect x="{panel_x}" y="{panel_y}" width="{panel_width}" height="{panel_height}" rx="10" fill="#ffffff" stroke="#cbd4dc" stroke-width="2"/>
    <text x="{panel_x + 28}" y="{panel_y + 40}" font-family="Arial, sans-serif" font-size="21" font-weight="700" fill="#243447">{label}</text>
    <text x="{panel_x + panel_width - 28}" y="{panel_y + 40}" text-anchor="end" font-family="Arial, sans-serif" font-size="14" fill="#60758a">{representation_identifier} · IFC GeometricCurveSet</text>
    <g class="official-wfb-native-dwg" data-source-kind="native_dwg" fill="none" stroke="{BLUE}" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke">
{linework}
    </g>
    <line x1="{panel_x + 88}" y1="{panel_y + panel_height - 46}" x2="{panel_x + panel_width - 88}" y2="{panel_y + panel_height - 46}" stroke="#d7dee5" stroke-width="1"/>
    <text x="{panel_x + panel_width / 2}" y="{panel_y + panel_height - 20}" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" fill="#60758a">source extents {width_mm:.0f} × {height_mm:.0f} mm · official blue linework</text>
  </g>"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_IFC)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    derived_path = args.input.resolve()
    manifest = load_json(args.manifest.resolve())
    if manifest.get("pass") is not True or manifest.get("approval", {}).get("status") != "approved":
        raise RuntimeError("derived IFC manifest is not approved")
    if sha256(derived_path) != manifest.get("derived_ifc_sha256"):
        raise RuntimeError("derived IFC hash does not match its manifest")

    model = ifcopenshell.open(derived_path)
    product = model.by_guid(GLOBAL_ID)
    if product is None:
        raise RuntimeError("approved Falper representative is missing")
    representations = {
        representation.RepresentationIdentifier: representation
        for representation in product.Representation.Representations
    }
    required = {
        "FalperWFBFront": ("ELEVATION_VIEW", 0),
        "FalperWFBSide": ("ELEVATION_VIEW", 1),
    }
    for identifier, (target_view, _) in required.items():
        representation = representations.get(identifier)
        if representation is None or representation.ContextOfItems.TargetView != target_view:
            raise RuntimeError(f"approved IFC representation gate failed: {identifier}")

    source = ifcopenshell.util.element.get_pset(
        product, "Pset_FalperSorgenteDrawingSource"
    ) or {}
    if (
        source.get("SourceKind") != "native_dwg"
        or source.get("SourceDwgSha256") != EXPECTED["wfb_2d"]
        or source.get("SourcePdfSha256") != EXPECTED["pdf"]
        or source.get("EvidenceScope") != SCOPE
        or source.get("ReviewStatus") != "APPROVED"
        or source.get("RepresentationGeometrySource") != "official_native_dwg_paths_mm"
        or source.get("ProxyGeometryIncluded") != "false"
    ):
        raise RuntimeError("approved IFC source metadata gate failed")

    front = extract_paths(representations["FalperWFBFront"], 0)
    side = extract_paths(representations["FalperWFBSide"], 1)
    if len(front) != 4 or len(side) != 4:
        raise RuntimeError("approved native DWG elevation path count gate failed")
    reviewer = html.escape(str(source["Reviewer"]))
    review_date = html.escape(str(source["ReviewDate"]))
    derived_hash = html.escape(manifest["derived_ifc_sha256"])
    dwg_hash = html.escape(source["SourceDwgSha256"])
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="1080" viewBox="0 0 1600 1080">
  <rect width="1600" height="1080" fill="#f8f7f3"/>
  <text x="60" y="58" font-family="Arial, sans-serif" font-size="31" font-weight="700" fill="#243447">Falper Sorgente WFB · approved elevations</text>
  <text x="60" y="96" font-family="Arial, sans-serif" font-size="18" fill="#4e6478">Derived IFC representation · isolated BS01 representative {GLOBAL_ID}</text>
  <rect x="1322" y="36" width="218" height="64" rx="8" fill="#e9f4fc" stroke="{BLUE}" stroke-width="2"/>
  <text x="1431" y="63" text-anchor="middle" font-family="Arial, sans-serif" font-size="16" font-weight="700" fill="{BLUE}">APPROVED</text>
  <text x="1431" y="85" text-anchor="middle" font-family="Arial, sans-serif" font-size="13" fill="#42647f">DERIVED IFC</text>
{panel_svg("FRONT ELEVATION / XZ", "FalperWFBFront", front, 60.0)}
{panel_svg("SIDE ELEVATION / YZ", "FalperWFBSide", side, 830.0)}
  <rect x="60" y="920" width="1480" height="118" rx="8" fill="#ffffff" stroke="#cbd4dc" stroke-width="2"/>
  <text x="84" y="952" font-family="Arial, sans-serif" font-size="15" font-weight="700" fill="#243447">SOURCE / APPROVAL</text>
  <text x="84" y="979" font-family="Arial, sans-serif" font-size="14" fill="#4e6478">Falper official Sorgente WFB native DWG · PDF vector cross-check passed · family reference, not project shop drawing</text>
  <text x="84" y="1004" font-family="Arial, sans-serif" font-size="13" fill="#60758a">DWG SHA-256 {dwg_hash}</text>
  <text x="84" y="1027" font-family="Arial, sans-serif" font-size="13" fill="#60758a">Derived IFC SHA-256 {derived_hash}</text>
  <text x="1516" y="979" text-anchor="end" font-family="Arial, sans-serif" font-size="14" fill="#4e6478">reviewer {reviewer}</text>
  <text x="1516" y="1004" text-anchor="end" font-family="Arial, sans-serif" font-size="14" fill="#4e6478">review date {review_date}</text>
</svg>
"""
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
