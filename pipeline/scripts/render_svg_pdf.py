#!/usr/bin/env python3
"""Render one SVG sheet to PDF and verify its page geometry.

Chrome preserves the vector SVG.  Poppler then supplies a mechanical page
check and a PNG proof for visual QA.  The proof belongs under ignored ``tmp/``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


MM_TO_POINTS = 72.0 / 25.4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_chrome(explicit: Path | None = None) -> Path:
    candidates = [
        explicit,
        Path(os.environ["CHROME_BIN"]) if os.environ.get("CHROME_BIN") else None,
        Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"),
        Path("/Applications/Chromium.app/Contents/MacOS/Chromium"),
    ]
    for name in ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser"):
        binary = shutil.which(name)
        if binary:
            candidates.append(Path(binary))
    for candidate in candidates:
        if candidate and candidate.is_file():
            return candidate
    raise RuntimeError("Chrome/Chromium not found; pass --chrome or set CHROME_BIN")


def parse_pdfinfo(text: str) -> dict[str, float | int]:
    pages_match = re.search(r"^Pages:\s+(\d+)\s*$", text, re.MULTILINE)
    size_match = re.search(r"^Page size:\s+([0-9.]+)\s+x\s+([0-9.]+)\s+pts", text, re.MULTILINE)
    if not pages_match or not size_match:
        raise RuntimeError("unable to parse pdfinfo page metadata")
    return {
        "pages": int(pages_match.group(1)),
        "width_points": float(size_match.group(1)),
        "height_points": float(size_match.group(2)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-svg", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--proof-png", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--chrome", type=Path)
    parser.add_argument("--width-mm", type=float, default=500.0)
    parser.add_argument("--height-mm", type=float, default=400.0)
    parser.add_argument("--page-tolerance-points", type=float, default=1.0)
    parser.add_argument("--proof-dpi", type=int, default=100)
    args = parser.parse_args()

    chrome = find_chrome(args.chrome)
    pdfinfo = shutil.which("pdfinfo")
    pdftoppm = shutil.which("pdftoppm")
    if not pdfinfo or not pdftoppm:
        raise RuntimeError("Poppler pdfinfo and pdftoppm are required")

    input_svg = args.input_svg.resolve()
    output_pdf = args.output_pdf.resolve()
    proof_png = args.proof_png.resolve()
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    proof_png.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            str(chrome),
            "--headless=new",
            "--disable-gpu",
            "--no-pdf-header-footer",
            f"--print-to-pdf={output_pdf}",
            input_svg.as_uri(),
        ],
        check=True,
    )
    info_text = subprocess.run(
        [pdfinfo, str(output_pdf)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    page = parse_pdfinfo(info_text)
    expected_width = args.width_mm * MM_TO_POINTS
    expected_height = args.height_mm * MM_TO_POINTS
    page_pass = (
        page["pages"] == 1
        and abs(float(page["width_points"]) - expected_width) <= args.page_tolerance_points
        and abs(float(page["height_points"]) - expected_height) <= args.page_tolerance_points
    )
    if not page_pass:
        raise RuntimeError(
            f"PDF page gate failed: {page}; expected 1 page at {args.width_mm}x{args.height_mm} mm"
        )

    proof_stem = proof_png.with_suffix("")
    subprocess.run(
        [
            pdftoppm,
            "-png",
            "-r",
            str(args.proof_dpi),
            "-singlefile",
            str(output_pdf),
            str(proof_stem),
        ],
        check=True,
    )
    if not proof_png.is_file() or proof_png.stat().st_size == 0:
        raise RuntimeError("PDF proof PNG was not created")

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_svg": str(input_svg),
        "source_svg_sha256": sha256(input_svg),
        "output_pdf": str(output_pdf),
        "output_pdf_sha256": sha256(output_pdf),
        "proof_png": str(proof_png),
        "proof_png_sha256": sha256(proof_png),
        "chrome": str(chrome),
        "page": page,
        "expected": {
            "pages": 1,
            "width_mm": args.width_mm,
            "height_mm": args.height_mm,
            "tolerance_points": args.page_tolerance_points,
        },
        "pass": True,
    }
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"PDF: {args.output_pdf}")
    print(f"proof PNG: {args.proof_png}")
    print(f"page: {page['pages']} × {page['width_points']} × {page['height_points']} pt")


if __name__ == "__main__":
    main()
