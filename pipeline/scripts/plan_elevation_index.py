#!/usr/bin/env python3
"""Replace per-Drawing plan markers with the official 12-anchor elevation index.

Bonsai emits one marker for every IFC Drawing.  This project retains the 36
official elevation Drawings plus project-only auxiliary views, but the plan
index follows the developer DWG: twelve anchor points carrying 36 official
view directions.  Auxiliary EL-P01/P02 views are deliberately absent.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REGISTER = ROOT / "pipeline/decisions/int1-elevation-view-register.csv"
START = "<!-- official-elevation-index:start -->"
END = "<!-- official-elevation-index:end -->"
RAW_MARKER = re.compile(
    r'^\s*<(?:use\b[^>]*(?:xlink:href|href)="#elevation-(?:arrow|tag)"[^>]*/>|'
    r'text\b[^>]*class="ELEVATION"[^>]*>.*?</text>)\s*$',
    re.MULTILINE,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_index(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 36:
        raise RuntimeError(f"official elevation register must contain 36 views: {path}")
    anchors = {row["anchor_id"] for row in rows}
    views = {row["view_id"] for row in rows}
    if anchors != {f"A{number}" for number in range(1, 13)} or len(views) != 36:
        raise RuntimeError("official elevation register identity is incomplete")
    if any(row["sheet_id"].startswith("EL-P") for row in rows):
        raise RuntimeError("project auxiliary views must not enter the official plan index")
    return rows


def svg_point(row: dict[str, str]) -> tuple[float, float]:
    # The native 1:50 plan Drawing uses a 400 mm square viewBox.  IFC +X maps
    # right and IFC +Y maps up; the DWG-derived coordinates are millimetres.
    return 200.0 + float(row["ifc_x_mm"]) / 50.0, 200.0 - float(row["ifc_y_mm"]) / 50.0


def direction_graphic(direction: str, view_id: str) -> str:
    vectors = {
        "+Y": (0.0, -1.0),
        "+X": (1.0, 0.0),
        "-Y": (0.0, 1.0),
        "-X": (-1.0, 0.0),
    }
    dx, dy = vectors[direction]
    start_x, start_y = dx * 3.0, dy * 3.0
    end_x, end_y = dx * 7.0, dy * 7.0
    label_x, label_y = dx * 9.4, dy * 9.4 + 0.55
    if direction in {"+X", "-X"}:
        label_y += 0.15
    return (
        f'<g class="official-elevation-direction" data-direction="{direction}" data-view-id="{view_id}">'
        f'<line x1="{start_x:.2f}" y1="{start_y:.2f}" x2="{end_x:.2f}" y2="{end_y:.2f}"/>'
        f'<circle cx="{end_x:.2f}" cy="{end_y:.2f}" r="0.75"/>'
        f'<text x="{label_x:.2f}" y="{label_y:.2f}">{view_id}</text></g>'
    )


def render_index(rows: list[dict[str, str]]) -> str:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["anchor_id"]].append(row)
    parts = [
        START,
        '<g id="official-elevation-index" data-source="official-DWG-KP-01" '
        'data-anchor-count="12" data-view-count="36">',
        '<style>.official-elevation-anchor{font-family:-apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif;}'
        '.official-elevation-anchor>circle{fill:#fffdf8;stroke:#8A6746;stroke-width:.34;}'
        '.official-elevation-anchor>.anchor-label{font-size:1.45px;font-weight:700;fill:#6F5035;text-anchor:middle;dominant-baseline:central;}'
        '.official-elevation-direction line{stroke:#B88A5A;stroke-width:.34;}'
        '.official-elevation-direction circle{fill:#B88A5A;stroke:none;}'
        '.official-elevation-direction text{font-size:1.55px;font-weight:650;fill:#6F5035;stroke:#fffdf8;stroke-width:.42;paint-order:stroke;'
        'text-anchor:middle;dominant-baseline:central;}</style>',
    ]
    for anchor_id in sorted(grouped, key=lambda value: int(value[1:])):
        anchor_rows = sorted(grouped[anchor_id], key=lambda row: int(row["view_id"]))
        x, y = svg_point(anchor_rows[0])
        view_ids = ",".join(row["view_id"] for row in anchor_rows)
        directions = "".join(
            direction_graphic(row["direction"], row["view_id"])
            for row in anchor_rows
        )
        parts.append(
            f'<g class="official-elevation-anchor" data-anchor-id="{anchor_id}" '
            f'data-view-ids="{view_ids}" data-ifc-x-mm="{anchor_rows[0]["ifc_x_mm"]}" '
            f'data-ifc-y-mm="{anchor_rows[0]["ifc_y_mm"]}" transform="translate({x:.3f} {y:.3f})">'
            f'<circle r="2.65"/><text class="anchor-label" x="0" y="0">{anchor_id}</text>{directions}</g>'
        )
    parts.extend(["</g>", END])
    return "\n".join(parts)


def apply_official_elevation_index(
    svg_path: Path,
    register_path: Path = DEFAULT_REGISTER,
    *,
    write: bool = True,
) -> dict[str, Any]:
    rows = read_index(register_path)
    source = svg_path.read_text(encoding="utf-8")
    had_managed_index = START in source or END in source
    if had_managed_index:
        if source.count(START) != 1 or source.count(END) != 1:
            raise RuntimeError(f"broken managed elevation-index block: {svg_path}")
        source = re.sub(
            r"\s*" + re.escape(START) + r".*?" + re.escape(END),
            "",
            source,
            flags=re.DOTALL,
        )
    raw_uses = len(re.findall(r'(?:xlink:href|href)="#elevation-(?:arrow|tag)"', source))
    raw_text = len(re.findall(r'<text\b[^>]*class="ELEVATION"', source))
    if raw_uses not in {0, 88} or raw_text not in {0, 88}:
        raise RuntimeError(
            f"unexpected native elevation-marker composition in {svg_path}: "
            f"uses={raw_uses}, labels={raw_text}"
        )
    cleaned = RAW_MARKER.sub("", source)
    if "</svg>" not in cleaned:
        raise RuntimeError(f"SVG closing tag missing: {svg_path}")
    cleaned = cleaned.replace("</svg>", "", 1).rstrip()
    output = cleaned + "\n" + render_index(rows) + "\n</svg>\n"
    if not write and output != svg_path.read_text(encoding="utf-8"):
        raise RuntimeError(f"official elevation index is stale: {svg_path}")
    if write:
        svg_path.write_text(output, encoding="utf-8")
    return {
        "svg": str(svg_path),
        # A managed file has already had the 44 native Drawing markers
        # removed, so an idempotent rerun sees zero raw references.  Preserve
        # the actual suppression count instead of degrading its provenance.
        "suppressed_ifc_drawing_marker_count": (
            raw_uses // 2 if raw_uses else 44 if had_managed_index else 0
        ),
        "anchor_count": 12,
        "official_view_count": 36,
        "auxiliary_marker_count": 0,
        "svg_sha256": sha256(svg_path),
        "pass": output.count('class="official-elevation-anchor"') == 12,
    }


def update_source_manifest(
    svg_path: Path, report: dict[str, Any], *, write: bool = True
) -> None:
    manifest_path = svg_path.with_name(svg_path.stem + "-source.json")
    if not manifest_path.is_file():
        return
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_index = {
        "source": "pipeline/decisions/int1-elevation-view-register.csv",
        "anchor_count": 12,
        "official_view_count": 36,
        "auxiliary_marker_count": 0,
        "suppressed_ifc_drawing_marker_count": report["suppressed_ifc_drawing_marker_count"],
    }
    if not write:
        if (
            payload.get("source_svg_sha256") != report["svg_sha256"]
            or payload.get("source_svg_bytes") != svg_path.stat().st_size
            or payload.get("elevation_index") != expected_index
        ):
            raise RuntimeError(f"plan source manifest is stale: {manifest_path}")
        return
    payload["source_svg_sha256"] = report["svg_sha256"]
    payload["source_svg_bytes"] = svg_path.stat().st_size
    payload["elevation_index"] = expected_index
    manifest_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("svg", nargs="*", type=Path)
    parser.add_argument("--register", type=Path, default=DEFAULT_REGISTER)
    parser.add_argument("--scan", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify the managed index and source manifests without writing them.",
    )
    args = parser.parse_args()
    paths = list(args.svg)
    if args.scan:
        paths.extend(
            path for path in args.scan.rglob("*.svg")
            if '#elevation-arrow' in path.read_text(encoding="utf-8", errors="ignore")
            or START in path.read_text(encoding="utf-8", errors="ignore")
        )
    paths = sorted({path.resolve() for path in paths})
    if not paths:
        raise RuntimeError("no SVGs selected")
    reports = []
    for path in paths:
        report = apply_official_elevation_index(
            path, args.register.resolve(), write=not args.check
        )
        update_source_manifest(path, report, write=not args.check)
        reports.append(report)
    payload = {
        "mode": (
            "official_12_anchor_plan_elevation_index_check"
            if args.check
            else "official_12_anchor_plan_elevation_index"
        ),
        "drawing_count": len(reports),
        "drawings": reports,
        "pass": all(report["pass"] for report in reports),
    }
    if args.report and not args.check:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
