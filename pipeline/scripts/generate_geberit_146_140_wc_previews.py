#!/usr/bin/env python3
"""Render Geberit 146.140 Bonsai scene SVG previews and one PDF review pack."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[2]
PRODUCT = ROOT / "output/review/highpoly-types/geberit-146-140"
SVG_DIR = PRODUCT / "bonsai-drawings/wc"
PDF = PRODUCT / "Geberit-146-140-WC-project-drawings.pdf"
MANIFEST = PRODUCT / "Geberit-146-140-WC-drawing-manifest.json"
FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
VIEWS = {
    "plan": "GEBERIT-146-140-WC-PLAN",
    "front": "GEBERIT-146-140-WC-FRONT",
    "side": "GEBERIT-146-140-WC-SIDE",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def render_svg(source: Path, target: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="geberit-146140-preview-") as temporary:
        result = subprocess.run(
            ["qlmanage", "-t", "-s", "2200", "-o", temporary, str(source)],
            capture_output=True,
            text=True,
        )
        generated = Path(temporary) / f"{source.name}.png"
        if result.returncode != 0 or not generated.is_file():
            raise RuntimeError(f"Quick Look SVG render failed: {result.stderr.strip()}")
        shutil.copy2(generated, target)


def create_pdf(records: list[dict]) -> None:
    page_width, page_height = landscape(A4)
    output = canvas.Canvas(str(PDF), pagesize=(page_width, page_height))
    for record in records:
        output.setTitle("Geberit AquaClean Sela 146.140 WC scene drawings")
        output.setFont("Helvetica-Bold", 15)
        output.drawString(36, page_height - 32, f"Geberit 146.140 — {record['view'].title()} scene SVG")
        output.setFont("Helvetica", 8)
        output.drawString(36, page_height - 46, "Bonsai Create Drawing / representative 1rhZG98PPCSxaLeMFLTYb9 / scale 1:25")
        with Image.open(record["preview"]) as image:
            image_width, image_height = image.size
        max_width, max_height = page_width - 72, page_height - 82
        scale = min(max_width / image_width, max_height / image_height)
        draw_width, draw_height = image_width * scale, image_height * scale
        output.drawImage(
            ImageReader(str(record["preview"])),
            (page_width - draw_width) / 2,
            22,
            draw_width,
            draw_height,
            preserveAspectRatio=True,
            mask="auto",
        )
        output.showPage()
    output.save()


def main() -> None:
    if sha256(FORMAL) != FORMAL_SHA:
        raise RuntimeError("formal IFC changed before preview generation")
    records = []
    for view, stem in VIEWS.items():
        svg = SVG_DIR / f"{stem}.svg"
        preview = PRODUCT / f"geberit-146-140-wc-{view}.png"
        render_svg(svg, preview)
        records.append(
            {
                "view": view,
                "svg": str(svg.relative_to(ROOT)),
                "svg_bytes": svg.stat().st_size,
                "svg_sha256": sha256(svg),
                "preview": preview,
                "preview_path": str(preview.relative_to(ROOT)),
                "preview_bytes": preview.stat().st_size,
                "preview_sha256": sha256(preview),
            }
        )
    create_pdf(records)
    manifest = {
        "schema_version": 1,
        "product": "Geberit AquaClean Sela 146.140",
        "article_number": "146.140.11.1",
        "representative_global_id": "1rhZG98PPCSxaLeMFLTYb9",
        "scene_scope": "representative WC instance only; sibling instance excluded",
        "generator": "Bonsai 0.8.4 bpy.ops.bim.create_drawing with raster previews rendered from the persisted SVG outputs",
        "views": [
            {key: value for key, value in record.items() if key != "preview"}
            for record in records
        ],
        "pdf": str(PDF.relative_to(ROOT)),
        "pdf_bytes": PDF.stat().st_size,
        "pdf_sha256": sha256(PDF),
        "formal_ifc_sha256": sha256(FORMAL),
        "formal_ifc_bytes_unchanged": True,
        "pass": True,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
