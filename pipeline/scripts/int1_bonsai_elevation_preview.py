#!/usr/bin/env python3
"""Compose phone-review SVG sheets from native Bonsai elevation SVGs."""

from __future__ import annotations

import csv
import math
import xml.etree.ElementTree as ET
from collections import defaultdict
from copy import deepcopy
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REGISTER = PROJECT_ROOT / "pipeline/decisions/int1-elevation-view-register.csv"
NATIVE_DIR = PROJECT_ROOT / "drawings/elevations/native"
SHEET_DIR = PROJECT_ROOT / "drawings/elevations"
WIDTH = 1600
HEIGHT = 1200
MARGIN = 44
HEADER = 92
GAP = 28
SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def native_name(row: dict[str, str]) -> str:
    suffix = {"+Y": "PY", "+X": "PX", "-Y": "NY", "-X": "NX"}[
        row["direction"]
    ]
    room = row["space_reference"].split(" ", 1)[0]
    return f'{row["sheet_id"]}-{row["view_id"]}-{room}-{suffix}'


def dimensions(count: int) -> tuple[int, int]:
    columns = 2 if count <= 4 else 3
    return columns, math.ceil(count / columns)


def main() -> None:
    rows = list(csv.DictReader(REGISTER.open(encoding="utf-8-sig")))
    sheets: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        sheets[row["sheet_id"]].append(row)

    for sheet_id, views in sorted(sheets.items()):
        views.sort(key=lambda row: int(row["view_id"]))
        columns, row_count = dimensions(len(views))
        cell_width = (WIDTH - 2 * MARGIN - (columns - 1) * GAP) / columns
        cell_height = (HEIGHT - HEADER - MARGIN - (row_count - 1) * GAP) / row_count
        root = ET.Element(
            f"{{{SVG_NS}}}svg",
            {
                "width": str(WIDTH),
                "height": str(HEIGHT),
                "viewBox": f"0 0 {WIDTH} {HEIGHT}",
                "data-source": "native Bonsai Drawing SVG composition",
            },
        )
        style = ET.SubElement(root, f"{{{SVG_NS}}}style")
        style.text = (
            "text{font-family:-apple-system,BlinkMacSystemFont,'PingFang SC',sans-serif;}"
            ".title{font-size:32px;font-weight:700;fill:#102a43;}"
            ".meta{font-size:16px;fill:#486581;}"
            ".view-title{font-size:18px;font-weight:650;fill:#102a43;}"
            ".review{font-size:14px;fill:#b42318;}"
            ".frame{fill:#fff;stroke:#bcccdc;stroke-width:1.5;}"
        )
        ET.SubElement(root, f"{{{SVG_NS}}}rect", {"width": str(WIDTH), "height": str(HEIGHT), "fill": "#f7f9fc"})
        title = ET.SubElement(root, f"{{{SVG_NS}}}text", {"x": str(MARGIN), "y": "43", "class": "title"})
        title.text = f"{sheet_id}  Bonsai 原生室内立面审核图"
        meta = ET.SubElement(root, f"{{{SVG_NS}}}text", {"x": str(MARGIN), "y": "72", "class": "meta"})
        meta.text = "IFC ELEVATION_VIEW · 1:50 · 无 underlay 纹理 · 红/橙框为非整数世界几何"

        for index, row in enumerate(views):
            column = index % columns
            grid_row = index // columns
            x = MARGIN + column * (cell_width + GAP)
            y = HEADER + grid_row * (cell_height + GAP)
            ET.SubElement(
                root,
                f"{{{SVG_NS}}}rect",
                {
                    "x": f"{x:.2f}",
                    "y": f"{y:.2f}",
                    "width": f"{cell_width:.2f}",
                    "height": f"{cell_height:.2f}",
                    "rx": "8",
                    "class": "frame",
                },
            )
            name = native_name(row)
            source = ET.parse(NATIVE_DIR / f"{name}.svg").getroot()
            view_box = source.attrib.get("viewBox")
            if not view_box:
                raise RuntimeError(f"missing viewBox: {name}")
            label = ET.SubElement(
                root,
                f"{{{SVG_NS}}}text",
                {"x": f"{x+16:.2f}", "y": f"{y+27:.2f}", "class": "view-title"},
            )
            label.text = (
                f'{row["view_id"]} · {row["space_reference"]} · 视向 {row["direction"]}'
            )
            if row["review_status"].startswith("review_required"):
                review = ET.SubElement(
                    root,
                    f"{{{SVG_NS}}}text",
                    {"x": f"{x+16:.2f}", "y": f"{y+49:.2f}", "class": "review"},
                )
                review.text = "官方标题/当前 Space 映射待复核"
            nested = ET.SubElement(
                root,
                f"{{{SVG_NS}}}svg",
                {
                    "x": f"{x+12:.2f}",
                    "y": f"{y+56:.2f}",
                    "width": f"{cell_width-24:.2f}",
                    "height": f"{cell_height-68:.2f}",
                    "viewBox": view_box,
                    "preserveAspectRatio": "xMidYMid meet",
                },
            )
            for child in source:
                nested.append(deepcopy(child))

        output = SHEET_DIR / f"{sheet_id}-bonsai-native.svg"
        ET.ElementTree(root).write(output, encoding="utf-8", xml_declaration=True)
        lines = output.read_text(encoding="utf-8").splitlines()
        output.write_text(
            "\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8"
        )
        print(output.relative_to(PROJECT_ROOT))


if __name__ == "__main__":
    main()
