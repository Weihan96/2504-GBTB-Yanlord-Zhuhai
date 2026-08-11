#!/usr/bin/env python3
"""Compile official-indexed IFC elevation renders into CAD-like review sheets.

The script is read-only with respect to IFC.  It consumes the camera/render
manifest produced inside the open Bonsai session, checks every raster hash, and
assembles nine SVG + PNG sheets.  Raw non-integer projected dimensions are
preserved and highlighted; the compiler never rounds or writes IFC geometry.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import os
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPECTED_SHEET_COUNTS = {
    "EL-01": 3,
    "EL-02": 2,
    "EL-03": 4,
    "EL-04": 4,
    "EL-05": 4,
    "EL-06": 5,
    "EL-07": 4,
    "EL-08": 6,
    "EL-09": 4,
}
PAGE_WIDTH_MM = 500.0
PAGE_HEIGHT_MM = 400.0
DRAWING_X = 7.0
DRAWING_Y = 7.0
DRAWING_WIDTH = 397.0
DRAWING_HEIGHT = 386.0
PANEL_X = 409.0
PANEL_WIDTH = 84.0
FLOAT_TOLERANCE_MM = 0.01
HIGH_TOLERANCE_MM = 0.1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def resolve(root: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def nearest_integer(value: float) -> tuple[int, float]:
    target = int(round(value))
    return target, value - target


def severity(value: float) -> str:
    _, residual = nearest_integer(value)
    absolute = abs(residual)
    if absolute <= FLOAT_TOLERANCE_MM:
        return "integer"
    if absolute <= HIGH_TOLERANCE_MM:
        return "amber"
    return "magenta"


def format_mm(value: float) -> str:
    if severity(value) == "integer":
        return str(int(round(value)))
    return f"{value:.3f}".rstrip("0").rstrip(".")


def noninteger_values(record: dict[str, Any]) -> list[dict[str, Any]]:
    projected = record["projected"]
    values = {
        "宽": float(projected["u_max_mm"]) - float(projected["u_min_mm"]),
        "高": float(projected["z_max_mm"]) - float(projected["z_min_mm"]),
        "底": float(projected["z_min_mm"]),
        "顶": float(projected["z_max_mm"]),
    }
    result = []
    for label, value in values.items():
        value_severity = severity(value)
        if value_severity == "integer":
            continue
        target, residual = nearest_integer(value)
        result.append(
            {
                "label": label,
                "value_mm": value,
                "nearest_integer_mm": target,
                "residual_mm": residual,
                "severity": value_severity,
            }
        )
    return result


def find_chrome(explicit: Path | None) -> Path:
    candidates = [
        explicit,
        Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"),
        Path("/Applications/Chromium.app/Contents/MacOS/Chromium"),
    ]
    for candidate in candidates:
        if candidate and candidate.is_file():
            return candidate
    raise RuntimeError("Chrome/Chromium is required to render elevation proof PNGs")


def render_png(svg_path: Path, png_path: Path, chrome: Path) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(chrome),
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            "--allow-file-access-from-files",
            "--force-device-scale-factor=1",
            f"--screenshot={png_path}",
            "--window-size=2000,1600",
            svg_path.resolve().as_uri(),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if not png_path.is_file() or png_path.stat().st_size == 0:
        raise RuntimeError(f"proof PNG was not created: {png_path}")
    if png_path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError(f"proof output is not PNG: {png_path}")


def grid_for(count: int) -> tuple[int, int]:
    if count <= 2:
        return 1, count
    if count <= 4:
        return 2, 2
    return 2, 3


def map_point(
    value_u: float,
    value_z: float,
    frame: dict[str, Any],
    x: float,
    y: float,
    width: float,
    height: float,
) -> tuple[float, float]:
    u_min = float(frame["u_min_mm"])
    u_max = float(frame["u_max_mm"])
    z_min = float(frame["z_min_mm"])
    z_max = float(frame["z_max_mm"])
    span_u = max(u_max - u_min, 1.0)
    span_z = max(z_max - z_min, 1.0)
    return (
        x + (value_u - u_min) / span_u * width,
        y + (z_max - value_z) / span_z * height,
    )


def object_highlights(
    view: dict[str, Any],
    image_x: float,
    image_y: float,
    image_width: float,
    image_height: float,
) -> tuple[str, list[dict[str, Any]]]:
    markup: list[str] = []
    issues: list[dict[str, Any]] = []
    frame = view["frame"]
    for record in view.get("objects", []):
        values = noninteger_values(record)
        if not values:
            continue
        record_severity = "magenta" if any(item["severity"] == "magenta" for item in values) else "amber"
        projected = record["projected"]
        x0, y1 = map_point(
            float(projected["u_min_mm"]),
            float(projected["z_min_mm"]),
            frame,
            image_x,
            image_y,
            image_width,
            image_height,
        )
        x1, y0 = map_point(
            float(projected["u_max_mm"]),
            float(projected["z_max_mm"]),
            frame,
            image_x,
            image_y,
            image_width,
            image_height,
        )
        if x1 < image_x or x0 > image_x + image_width or y1 < image_y or y0 > image_y + image_height:
            continue
        x0 = max(image_x, min(x0, image_x + image_width))
        x1 = max(image_x, min(x1, image_x + image_width))
        y0 = max(image_y, min(y0, image_y + image_height))
        y1 = max(image_y, min(y1, image_y + image_height))
        if x1 - x0 < 0.25 or y1 - y0 < 0.25:
            continue
        markup.append(
            f'<rect class="noninteger-{record_severity}" x="{x0:.3f}" y="{y0:.3f}" '
            f'width="{x1-x0:.3f}" height="{y1-y0:.3f}" data-global-id="{esc(record["global_id"])}"/>'
        )
        maximum = max(values, key=lambda item: abs(float(item["residual_mm"])))
        issues.append(
            {
                "view_id": view["view_id"],
                "global_id": record["global_id"],
                "ifc_class": record["ifc_class"],
                "name": record.get("name", ""),
                "values": values,
                "maximum_abs_residual_mm": abs(float(maximum["residual_mm"])),
                "severity": record_severity,
            }
        )
    return "".join(markup), issues


def dimension_markup(
    x: float,
    y: float,
    width: float,
    height: float,
    frame: dict[str, Any],
) -> str:
    raw_width = float(frame["u_max_mm"]) - float(frame["u_min_mm"])
    raw_height = float(frame["z_max_mm"]) - float(frame["z_min_mm"])
    width_class = f'dim-{severity(raw_width)}'
    height_class = f'dim-{severity(raw_height)}'
    baseline_y = y + height + 4.0
    vertical_x = x - 3.2
    return (
        f'<line class="dim-line" x1="{x:.3f}" y1="{baseline_y:.3f}" x2="{x+width:.3f}" y2="{baseline_y:.3f}"/>'
        f'<line class="dim-tick" x1="{x:.3f}" y1="{baseline_y-1.2:.3f}" x2="{x:.3f}" y2="{baseline_y+1.2:.3f}"/>'
        f'<line class="dim-tick" x1="{x+width:.3f}" y1="{baseline_y-1.2:.3f}" x2="{x+width:.3f}" y2="{baseline_y+1.2:.3f}"/>'
        f'<text class="{width_class}" x="{x+width/2:.3f}" y="{baseline_y-0.7:.3f}" text-anchor="middle">{format_mm(raw_width)}</text>'
        f'<line class="dim-line" x1="{vertical_x:.3f}" y1="{y:.3f}" x2="{vertical_x:.3f}" y2="{y+height:.3f}"/>'
        f'<line class="dim-tick" x1="{vertical_x-1.2:.3f}" y1="{y:.3f}" x2="{vertical_x+1.2:.3f}" y2="{y:.3f}"/>'
        f'<line class="dim-tick" x1="{vertical_x-1.2:.3f}" y1="{y+height:.3f}" x2="{vertical_x+1.2:.3f}" y2="{y+height:.3f}"/>'
        f'<text class="{height_class}" x="{vertical_x-1.1:.3f}" y="{y+height/2:.3f}" text-anchor="middle" transform="rotate(-90 {vertical_x-1.1:.3f} {y+height/2:.3f})">{format_mm(raw_height)}</text>'
    )


def view_markup(
    root: Path,
    sheet_svg: Path,
    view: dict[str, Any],
    register: dict[str, str],
    x: float,
    y: float,
    width: float,
    height: float,
) -> tuple[str, list[dict[str, Any]]]:
    header_height = 9.0
    dimension_margin = 9.0
    image_x = x + 7.0
    image_y = y + header_height + 1.0
    image_width = width - 10.0
    image_height = height - header_height - dimension_margin - 2.0
    png_path = resolve(root, view["png"]["path"])
    href = os.path.relpath(png_path, sheet_svg.parent).replace(os.sep, "/")
    highlights, issues = object_highlights(view, image_x, image_y, image_width, image_height)
    ceiling_levels = sorted(
        {
            int(round(float(record["projected"]["z_min_mm"])))
            for record in view.get("objects", [])
            if "ceiling" in str(record.get("name", "")).lower()
            and float(record["projected"]["z_min_mm"]) >= 2000
        }
    )
    ceiling_label = "/".join(str(value) for value in ceiling_levels[:4]) or "未建模"
    markup = f'''
<g class="elevation-view" data-view-id="{esc(view['view_id'])}" data-sheet-id="{esc(view['sheet_id'])}">
  <rect class="view-border" x="{x:.3f}" y="{y:.3f}" width="{width:.3f}" height="{height:.3f}"/>
  <text class="view-title" x="{x+4:.3f}" y="{y+5.7:.3f}">{esc(view['view_id'])} · {esc(register['direction'])} · {esc(register['space_reference'])}</text>
  <text class="view-meta" x="{x+width-4:.3f}" y="{y+5.7:.3f}" text-anchor="end">天花 {esc(ceiling_label)} mm</text>
  <rect class="image-bg" x="{image_x:.3f}" y="{image_y:.3f}" width="{image_width:.3f}" height="{image_height:.3f}"/>
  <image x="{image_x:.3f}" y="{image_y:.3f}" width="{image_width:.3f}" height="{image_height:.3f}" preserveAspectRatio="none" href="{esc(href)}"/>
  {highlights}
</g>
'''
    return markup, issues


def make_sheet_svg(
    root: Path,
    svg_path: Path,
    sheet_id: str,
    rows: list[dict[str, str]],
    views: list[dict[str, Any]],
    ifc_hash: str,
    register_hash: str,
) -> tuple[str, list[dict[str, Any]]]:
    columns, grid_rows = grid_for(len(views))
    gap = 5.0
    cell_width = (DRAWING_WIDTH - gap * (columns - 1)) / columns
    cell_height = (DRAWING_HEIGHT - gap * (grid_rows - 1)) / grid_rows
    register_by_id = {row["view_id"]: row for row in rows}
    markup: list[str] = []
    issues: list[dict[str, Any]] = []
    for index, view in enumerate(views):
        column = index % columns
        row_index = index // columns
        x = DRAWING_X + column * (cell_width + gap)
        y = DRAWING_Y + row_index * (cell_height + gap)
        view_svg, view_issues = view_markup(
            root,
            svg_path,
            view,
            register_by_id[view["view_id"]],
            x,
            y,
            cell_width,
            cell_height,
        )
        markup.append(view_svg)
        issues.extend(view_issues)
    title = rows[0]["official_title"]
    magenta_count = sum(item["severity"] == "magenta" for item in issues)
    amber_count = sum(item["severity"] == "amber" for item in issues)
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="500mm" height="400mm" viewBox="0 0 500 400"
 data-sheet-id="{esc(sheet_id)}" data-source-ifc-sha256="{ifc_hash}" data-view-register-sha256="{register_hash}">
<style>
 @page {{ size:500mm 400mm; margin:0; }}
 * {{ box-sizing:border-box; }}
 text {{ font-family:Arial,"PingFang SC",sans-serif; fill:#111827; }}
 .view-border {{ fill:#fff; stroke:#111827; stroke-width:.42; }}
 .image-bg {{ fill:#fafafa; stroke:#94a3b8; stroke-width:.18; }}
 .view-title {{ font-size:3.2px; font-weight:700; }}
 .view-meta {{ font-size:2.3px; fill:#475569; }}
 .dim-line,.dim-tick {{ stroke:#374151; stroke-width:.22; fill:none; }}
 .dim-integer {{ font-size:2.25px; fill:#111827; }}
 .dim-amber {{ font-size:2.25px; font-weight:700; fill:#b45309; }}
 .dim-magenta {{ font-size:2.25px; font-weight:700; fill:#c026d3; }}
 .noninteger-amber {{ fill:#f59e0b; fill-opacity:.08; stroke:#d97706; stroke-width:.45; stroke-dasharray:2 1; }}
 .noninteger-magenta {{ fill:#d946ef; fill-opacity:.09; stroke:#c026d3; stroke-width:.6; stroke-dasharray:2 1; }}
 .panel {{ fill:#f8fafc; stroke:#111827; stroke-width:.45; }}
 .sheet {{ font-size:7px; font-weight:800; fill:#0f4c81; }}
 .title {{ font-size:4.3px; font-weight:700; }}
 .subtitle {{ font-size:2.7px; font-weight:700; fill:#991b1b; }}
 .note {{ font-size:2.45px; fill:#334155; }}
 .legend-amber {{ fill:#d97706; }} .legend-magenta {{ fill:#c026d3; }}
</style>
<rect width="500" height="400" fill="white"/>
{''.join(markup)}
<rect class="panel" x="{PANEL_X}" y="7" width="{PANEL_WIDTH}" height="386"/>
<text class="sheet" x="414" y="20">{esc(sheet_id)}</text>
<text class="title" x="414" y="30">{esc(title)}</text>
<text class="subtitle" x="414" y="39">IFC 立面协调候选 / 非施工发布</text>
<text class="note" x="414" y="51">官方平面索引：1:50@A2</text>
<text class="note" x="414" y="57">视图：{len(views)} / 方向按 DWG 箭头</text>
<text class="note" x="414" y="63">天花：当前 DCL IFC 下表面</text>
<text class="note" x="414" y="69">材质：当前 IFC 未绑定，不推断</text>
<rect class="legend-magenta" x="414" y="82" width="4" height="4"/><text class="note" x="421" y="85.4">残差 &gt; 0.1 mm：{magenta_count}</text>
<rect class="legend-amber" x="414" y="91" width="4" height="4"/><text class="note" x="421" y="94.4">0.01–0.1 mm：{amber_count}</text>
<text class="note" x="414" y="106">高亮保留原始值；不自动归整 IFC。</text>
<text class="note" x="414" y="115">完成面缺失区不得用 Space 顶面补齐。</text>
<text class="note" x="414" y="124">DEMOLISH 墙已排除；真实深度遮挡。</text>
<text class="note" x="414" y="344">IFC SHA</text>
<text class="note" x="414" y="351">{ifc_hash[:24]}…</text>
<text class="note" x="414" y="362">View register SHA</text>
<text class="note" x="414" y="369">{register_hash[:24]}…</text>
<text class="subtitle" x="414" y="382">AUTOMATIC IFC WRITE = FALSE</text>
</svg>
''', issues


def validate_inputs(
    root: Path,
    ifc_path: Path,
    register_path: Path,
    manifest_path: Path,
) -> tuple[list[dict[str, str]], dict[str, Any], dict[str, dict[str, Any]], str]:
    ifc_hash = sha256(ifc_path)
    rows = read_csv(register_path)
    if len(rows) != 36 or len({row["view_id"] for row in rows}) != 36:
        raise RuntimeError("elevation view register must contain 36 unique views")
    if Counter(row["sheet_id"] for row in rows) != Counter(EXPECTED_SHEET_COUNTS):
        raise RuntimeError("elevation view register sheet counts do not match EL-01..EL-09")
    if {row["direction"] for row in rows} != {"+Y", "+X", "-Y", "-X"}:
        raise RuntimeError("elevation view register must contain four official directions")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("source_ifc_sha256") != ifc_hash:
        raise RuntimeError("elevation render manifest is stale against the formal IFC")
    if manifest.get("view_register_sha256") != sha256(register_path):
        raise RuntimeError("elevation render manifest is stale against the view register")
    views = manifest.get("views", [])
    if len(views) != 36 or len({str(view["view_id"]) for view in views}) != 36:
        raise RuntimeError("elevation render manifest must contain 36 unique views")
    by_id = {str(view["view_id"]): view for view in views}
    if set(by_id) != {row["view_id"] for row in rows}:
        raise RuntimeError("rendered views do not match the official view register")
    for view in views:
        if int(view.get("demolish_visible_count", -1)) != 0:
            raise RuntimeError(f"DEMOLISH wall leaked into elevation {view['view_id']}")
        png = resolve(root, view["png"]["path"])
        if not png.is_file() or sha256(png) != view["png"]["sha256"]:
            raise RuntimeError(f"stale or missing raw elevation PNG: {view['view_id']}")
        if png.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
            raise RuntimeError(f"raw elevation is not PNG: {view['view_id']}")
        for key in ("u_min_mm", "u_max_mm", "z_min_mm", "z_max_mm"):
            if not math.isfinite(float(view["frame"][key])):
                raise RuntimeError(f"invalid frame {key}: {view['view_id']}")
    return rows, manifest, by_id, ifc_hash


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--input-ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument("--views", type=Path, default=Path("pipeline/decisions/int1-elevation-view-register.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("build/int1/elevation-render-manifest.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("drawings/elevations"))
    parser.add_argument("--proof-dir", type=Path, default=Path("output/images/elevations"))
    parser.add_argument("--report", type=Path, default=Path("build/int1/elevation-sheet-candidate.json"))
    parser.add_argument("--chrome", type=Path)
    parser.add_argument("--skip-render-png", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    ifc_path = resolve(root, args.input_ifc)
    register_path = resolve(root, args.views)
    manifest_path = resolve(root, args.manifest)
    output_dir = resolve(root, args.output_dir)
    proof_dir = resolve(root, args.proof_dir)
    report_path = resolve(root, args.report)
    rows, manifest, views_by_id, ifc_hash = validate_inputs(root, ifc_path, register_path, manifest_path)
    register_hash = sha256(register_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    proof_dir.mkdir(parents=True, exist_ok=True)
    chrome = None if args.skip_render_png else find_chrome(args.chrome)
    sheet_reports = []
    all_issues: list[dict[str, Any]] = []
    for sheet_id in EXPECTED_SHEET_COUNTS:
        sheet_rows = [row for row in rows if row["sheet_id"] == sheet_id]
        sheet_views = [views_by_id[row["view_id"]] for row in sheet_rows]
        svg_path = output_dir / f"{sheet_id}-ifc-elevation-candidate.svg"
        png_path = proof_dir / f"{sheet_id}-ifc-elevation-candidate.png"
        svg, issues = make_sheet_svg(root, svg_path, sheet_id, sheet_rows, sheet_views, ifc_hash, register_hash)
        svg_path.write_text(svg, encoding="utf-8")
        if chrome:
            render_png(svg_path, png_path, chrome)
        all_issues.extend(issues)
        sheet_reports.append(
            {
                "sheet_id": sheet_id,
                "official_title": sheet_rows[0]["official_title"],
                "view_ids": [row["view_id"] for row in sheet_rows],
                "svg": {"path": str(svg_path.relative_to(root)), "sha256": sha256(svg_path)},
                "proof_png": None if not png_path.is_file() else {
                    "path": str(png_path.relative_to(root)),
                    "sha256": sha256(png_path),
                },
                "noninteger_object_count": len(issues),
                "magenta_object_count": sum(item["severity"] == "magenta" for item in issues),
                "amber_object_count": sum(item["severity"] == "amber" for item in issues),
            }
        )
    report = {
        "mode": "read_only_official_indexed_ifc_elevation_candidate",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_ifc": {"path": str(ifc_path), "sha256": ifc_hash},
        "view_register": {"path": str(register_path), "sha256": register_hash},
        "render_manifest": {"path": str(manifest_path), "sha256": sha256(manifest_path)},
        "rules": {
            "float_tolerance_mm": FLOAT_TOLERANCE_MM,
            "magenta_threshold_mm": HIGH_TOLERANCE_MM,
            "integer_geometry_is_not_modified": True,
            "ceiling_basis": "current DCL IFC downward finish faces along each view",
            "space_top_used_as_ceiling": False,
            "material_inference_allowed": False,
        },
        "summary": {
            "sheet_count": len(sheet_reports),
            "view_count": len(rows),
            "noninteger_object_count": len(all_issues),
            "magenta_object_count": sum(item["severity"] == "magenta" for item in all_issues),
            "amber_object_count": sum(item["severity"] == "amber" for item in all_issues),
        },
        "sheets": sheet_reports,
        "noninteger_issues": sorted(
            all_issues,
            key=lambda item: (item["view_id"], -float(item["maximum_abs_residual_mm"]), item["global_id"]),
        ),
        "gates": {
            "official_view_count_36": len(rows) == 36,
            "official_sheet_count_9": len(sheet_reports) == 9,
            "render_manifest_matches_ifc_and_register": True,
            "all_raw_png_hashes_match": True,
            "demolish_wall_visible_count": 0,
            "all_svg_outputs_present": all((root / sheet["svg"]["path"]).is_file() for sheet in sheet_reports),
            "all_proof_png_outputs_present": args.skip_render_png or all(sheet["proof_png"] for sheet in sheet_reports),
            "automatic_ifc_write_allowed": False,
            "construction_release_ready": False,
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(report_path), **report["summary"], "gates": report["gates"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
