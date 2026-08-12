#!/usr/bin/env python3
"""Compose public-space Bonsai Drawing SVGs and render phone-review PNGs."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REGISTER = PROJECT_ROOT / "pipeline/decisions/int1-public-elevation-register.csv"
NATIVE_DIR = PROJECT_ROOT / "drawings/elevations/native"
SHEET_DIR = PROJECT_ROOT / "drawings/elevations"
IMAGE_DIR = PROJECT_ROOT / "output/images/elevations"
WIDTH = 2400
HEIGHT = 1500
MARGIN = 54
HEADER = 154
GAP = 34
SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def chrome() -> Path:
    candidates = [
        Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"),
        Path("/Applications/Chromium.app/Contents/MacOS/Chromium"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    found = shutil.which("google-chrome") or shutil.which("chromium")
    if found:
        return Path(found)
    raise RuntimeError("Chrome/Chromium is required to render review PNGs")


def root(title: str, subtitle: str) -> ET.Element:
    svg = ET.Element(
        f"{{{SVG_NS}}}svg",
        {
            "width": str(WIDTH),
            "height": str(HEIGHT),
            "viewBox": f"0 0 {WIDTH} {HEIGHT}",
            "data-source": "native Bonsai Drawing SVG composition",
        },
    )
    style = ET.SubElement(svg, f"{{{SVG_NS}}}style")
    style.text = (
        "text{font-family:-apple-system,BlinkMacSystemFont,'PingFang SC',sans-serif;}"
        ".int1-public-title{font-size:40px;font-weight:750;fill:#102a43;}"
        ".int1-public-meta{font-size:20px;fill:#486581;}"
        ".int1-public-label{font-size:20px;font-weight:700;fill:#102a43;}"
        ".int1-public-note{font-size:16px;fill:#7c2d12;}"
        ".int1-public-frame{fill:#fff;stroke:#bcccdc;stroke-width:2;}"
        ".int1-public-break{fill:none;stroke:#d97706;stroke-width:4;}"
    )
    ET.SubElement(svg, f"{{{SVG_NS}}}rect", {"width": str(WIDTH), "height": str(HEIGHT), "fill": "#f7f9fc"})
    heading = ET.SubElement(svg, f"{{{SVG_NS}}}text", {"x": str(MARGIN), "y": "58", "class": "int1-public-title"})
    heading.text = title
    meta = ET.SubElement(svg, f"{{{SVG_NS}}}text", {"x": str(MARGIN), "y": "96", "class": "int1-public-meta"})
    meta.text = subtitle
    warning = ET.SubElement(svg, f"{{{SVG_NS}}}text", {"x": str(MARGIN), "y": "130", "class": "int1-public-note"})
    warning.text = "项目编译审核图，不冒充开发商原图；红/橙框为非整数世界几何。"
    return svg


def add_native(svg: ET.Element, row: dict[str, str], x: float, y: float, width: float, height: float) -> None:
    ET.SubElement(
        svg,
        f"{{{SVG_NS}}}rect",
        {"x": f"{x:.2f}", "y": f"{y:.2f}", "width": f"{width:.2f}", "height": f"{height:.2f}", "rx": "9", "class": "int1-public-frame"},
    )
    label = ET.SubElement(svg, f"{{{SVG_NS}}}text", {"x": f"{x+16:.2f}", "y": f"{y+31:.2f}", "class": "int1-public-label"})
    label.text = f'{row["view_id"]} · {row["space_reference"]} · 视向 {row["direction"]}'
    source_path = NATIVE_DIR / f'{row["drawing_name"]}.svg'
    normalize_svg(source_path)
    source = ET.parse(source_path).getroot()
    view_box = source.attrib.get("viewBox")
    if not view_box:
        raise RuntimeError(f'missing viewBox: {row["drawing_name"]}')
    nested = ET.SubElement(
        svg,
        f"{{{SVG_NS}}}svg",
        {"x": f"{x+12:.2f}", "y": f"{y+46:.2f}", "width": f"{width-24:.2f}", "height": f"{height-58:.2f}", "viewBox": view_box, "preserveAspectRatio": "xMidYMid meet"},
    )
    for child in source:
        nested.append(deepcopy(child))


def write_and_render(svg: ET.Element, stem: str) -> dict[str, object]:
    svg_path = SHEET_DIR / f"{stem}.svg"
    png_path = IMAGE_DIR / f"{stem}.png"
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(svg).write(svg_path, encoding="utf-8", xml_declaration=True)
    normalize_svg(svg_path)
    subprocess.run(
        [
            str(chrome()),
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            "--force-device-scale-factor=1",
            f"--screenshot={png_path}",
            f"--window-size={WIDTH},{HEIGHT}",
            svg_path.resolve().as_uri(),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if png_path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError(f"invalid PNG: {png_path}")
    parsed = ET.parse(svg_path).getroot()
    image_count = sum(element.tag.endswith("image") for element in parsed.iter())
    path_count = sum(element.tag.endswith("path") for element in parsed.iter())
    if image_count or path_count < 100:
        raise RuntimeError(
            f"sheet QA failed for {stem}: images={image_count}, paths={path_count}"
        )
    print(svg_path.relative_to(PROJECT_ROOT))
    print(png_path.relative_to(PROJECT_ROOT))
    return {
        "stem": stem,
        "svg": str(svg_path.relative_to(PROJECT_ROOT)),
        "svg_sha256": sha256(svg_path),
        "png": str(png_path.relative_to(PROJECT_ROOT)),
        "png_sha256": sha256(png_path),
        "image_count": image_count,
        "path_count": path_count,
        "pixel_size": [WIDTH, HEIGHT],
    }


def main() -> None:
    rows = list(csv.DictReader(REGISTER.open(encoding="utf-8-sig")))
    unfolded = [row for row in rows if row["sheet_id"] == "EL-P01"]
    long_sections = [row for row in rows if row["sheet_id"] == "EL-P02"]
    if len(unfolded) != 6 or len(long_sections) != 2:
        raise RuntimeError("public elevation register must contain 6 unfolded + 2 long views")

    p01 = root(
        "EL-P01  客厅—餐厅—西厨—中厨折线展开立面",
        "6 段 Bonsai 原生 ELEVATION_VIEW · 1:30 · 各段之间以断开符号明确分隔",
    )
    cell_width = (WIDTH - 2 * MARGIN - 2 * GAP) / 3
    cell_height = (HEIGHT - HEADER - MARGIN - GAP) / 2
    for index, row in enumerate(unfolded):
        column, band = index % 3, index // 3
        x = MARGIN + column * (cell_width + GAP)
        y = HEADER + band * (cell_height + GAP)
        add_native(p01, row, x, y, cell_width, cell_height)
        if column < 2:
            bx = x + cell_width + GAP / 2
            by = y + cell_height / 2
            ET.SubElement(p01, f"{{{SVG_NS}}}path", {"d": f"M {bx-9} {by-28} l 18 18 l -18 18 l 18 18", "class": "int1-public-break"})
    reports = [write_and_render(p01, "EL-P01-bonsai-public-unfolded")]

    p02 = root(
        "EL-P02  客厅—餐厅—西厨—中厨互补贯穿长立面",
        "2 条 Bonsai 原生正交长剖 · 1:30 · 相反视向互补，保留公共空间整体关系",
    )
    long_height = (HEIGHT - HEADER - MARGIN - GAP) / 2
    for index, row in enumerate(long_sections):
        add_native(p02, row, MARGIN, HEADER + index * (long_height + GAP), WIDTH - 2 * MARGIN, long_height)
    reports.append(write_and_render(p02, "EL-P02-bonsai-public-long-section"))
    report_path = PROJECT_ROOT / "build/int1/public-elevation-preview-report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(
            {
                "sheet_count": len(reports),
                "native_view_count": len(rows),
                "sheets": reports,
                "pass": len(reports) == 2
                and all(not report["image_count"] for report in reports),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
