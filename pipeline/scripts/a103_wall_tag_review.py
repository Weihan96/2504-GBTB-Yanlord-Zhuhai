#!/usr/bin/env python3
"""Render a phone-friendly A-103 wall-tag review image without writing IFC."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import struct
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
SOURCE_SVG = ROOT / "drawings/Wall Plan.svg"
CANDIDATE_REPORT = ROOT / "build/a103/a103-wall-tag-candidate.json"
OUTPUT_SVG = ROOT / "drawings/A103-wall-tag-review.svg"
OUTPUT_PNG = ROOT / "output/images/A-103-wall-tag-review.png"
RENDER_REPORT = ROOT / "build/a103/a103-wall-tag-review.json"
CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def world_to_svg(x_mm: float, y_mm: float) -> tuple[float, float]:
    return ((x_mm + 10000.0) / 50.0, (10000.0 - y_mm) / 50.0)


def intersects(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> bool:
    padding = 0.45
    return not (
        a[2] + padding <= b[0]
        or b[2] + padding <= a[0]
        or a[3] + padding <= b[1]
        or b[3] + padding <= a[1]
    )


def layout(records: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], int]:
    width, height = 8.8, 3.2
    boxes: list[tuple[float, float, float, float]] = []
    placed: list[dict[str, Any]] = []
    unresolved = 0
    for record in records:
        anchor_x, anchor_y = world_to_svg(record["centre_x_mm"], record["centre_y_mm"])
        selected = None
        for radius in (3.0, 5.5, 8.0, 11.0, 14.5, 18.0, 22.0, 27.0):
            for angle in (0, 180, 90, 270, 45, 135, 315, 225):
                radians = math.radians(angle)
                x = anchor_x + math.cos(radians) * radius
                y = anchor_y - math.sin(radians) * radius
                box = (x - width / 2, y - height / 2, x + width / 2, y + height / 2)
                if min(box) < 1.5 or box[2] > 398.5 or box[3] > 398.5:
                    continue
                if any(intersects(box, other) for other in boxes):
                    continue
                selected = (x, y, box)
                break
            if selected:
                break
        if selected is None:
            unresolved += 1
            selected = (anchor_x, anchor_y, (anchor_x, anchor_y, anchor_x, anchor_y))
        x, y, box = selected
        boxes.append(box)
        placed.append({**record, "anchor_x": anchor_x, "anchor_y": anchor_y, "label_x": x, "label_y": y})
    return placed, unresolved


def render_overlay(records: list[dict[str, Any]]) -> str:
    parts = [
        "<style>",
        ".tag-review-line{stroke:#334155;stroke-width:.28;opacity:.8}",
        ".tag-review-anchor{fill:#0f172a;stroke:#fff;stroke-width:.25}",
        ".tag-review-existing{fill:#e0f2fe;stroke:#0369a1;stroke-width:.3}",
        ".tag-review-new{fill:#dcfce7;stroke:#15803d;stroke-width:.3}",
        ".tag-review-text{font:700 1.85px ui-monospace,SFMono-Regular,monospace;text-anchor:middle;dominant-baseline:central;fill:#0f172a}",
        ".tag-review-title{font:700 4px -apple-system,BlinkMacSystemFont,'PingFang SC',sans-serif;fill:#0f172a}",
        ".tag-review-meta{font:2.3px -apple-system,BlinkMacSystemFont,'PingFang SC',sans-serif;fill:#334155}",
        "</style>",
        '<g id="A103_WALL_TAG_REVIEW">',
        '<rect x="4" y="4" width="205" height="17" rx="2" fill="#ffffff" fill-opacity=".94" stroke="#94a3b8" stroke-width=".35"/>',
        '<text class="tag-review-title" x="8" y="10">A-103 墙编号一次性审核候选</text>',
        '<text class="tag-review-meta" x="8" y="16">现状 EW01–EW84｜新建 NW01–NW04｜北→南、同排西→东｜未写正式 IFC</text>',
    ]
    for record in records:
        css = "tag-review-new" if record["status"] == "NEW" else "tag-review-existing"
        parts.extend(
            [
                f'<line class="tag-review-line" x1="{record["anchor_x"]:.3f}" y1="{record["anchor_y"]:.3f}" x2="{record["label_x"]:.3f}" y2="{record["label_y"]:.3f}"/>',
                f'<circle class="tag-review-anchor" cx="{record["anchor_x"]:.3f}" cy="{record["anchor_y"]:.3f}" r=".72"/>',
                f'<rect class="{css}" x="{record["label_x"] - 4.4:.3f}" y="{record["label_y"] - 1.6:.3f}" width="8.8" height="3.2" rx=".7"/>',
                f'<text class="tag-review-text" x="{record["label_x"]:.3f}" y="{record["label_y"]:.3f}">{html.escape(record["candidate_tag"])}</text>',
            ]
        )
    parts.append("</g>")
    return "".join(parts)


def png_dimensions(path: Path) -> tuple[int, int]:
    header = path.read_bytes()[:24]
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError("review output is not a valid PNG")
    return struct.unpack(">II", header[16:24])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=FORMAL_IFC)
    parser.add_argument("--source-svg", type=Path, default=SOURCE_SVG)
    parser.add_argument("--candidate-report", type=Path, default=CANDIDATE_REPORT)
    parser.add_argument("--output-svg", type=Path, default=OUTPUT_SVG)
    parser.add_argument("--output-png", type=Path, default=OUTPUT_PNG)
    parser.add_argument("--report", type=Path, default=RENDER_REPORT)
    parser.add_argument("--chrome", type=Path, default=CHROME)
    args = parser.parse_args()

    candidate = json.loads(args.candidate_report.read_text(encoding="utf-8"))
    if not candidate.get("pass") or candidate["source_ifc_sha256"] != sha256(args.ifc):
        raise RuntimeError("wall-tag candidate is failed or stale")
    records, unresolved = layout(candidate["records"])
    if unresolved:
        raise RuntimeError(f"could not place {unresolved} wall-tag labels without collision")
    source = args.source_svg.read_text(encoding="utf-8")
    if source.count("</svg>") != 1:
        raise RuntimeError("Wall Plan source is not a single SVG document")
    output = source.replace("</svg>", render_overlay(records) + "</svg>")
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(output, encoding="utf-8")
    subprocess.run(
        [
            str(args.chrome),
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            "--force-device-scale-factor=2",
            f"--screenshot={args.output_png.resolve()}",
            "--window-size=1512,1512",
            args.output_svg.resolve().as_uri(),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    width, height = png_dimensions(args.output_png)
    report = {
        "source_ifc_sha256": candidate["source_ifc_sha256"],
        "candidate_report": str(args.candidate_report.resolve()),
        "output_svg": str(args.output_svg.resolve()),
        "output_png": str(args.output_png.resolve()),
        "tag_count": len(records),
        "label_collision_count": unresolved,
        "png_width_px": width,
        "png_height_px": height,
        "formal_ifc_write_allowed": False,
        "pass": len(records) == 88 and unresolved == 0 and width >= 2000 and height >= 2000,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False))
    return 0 if report["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
