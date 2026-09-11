#!/usr/bin/env python3
"""Render TRAP01 cabinet-internal detail SVGs to PNG/PDF and finalize evidence."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/cabinet-internal-detail"
EVIDENCE = OUTPUT_DIR / "TRAP01-cabinet-internal-detail-create-drawing-evidence.json"
MANIFEST = PRODUCT_DIR / "TRAP01-cabinet-internal-detail-manifest.json"
COMBINED_PDF = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-cabinet-internal-detail.pdf"
PREVIEW_PREFIX = "trap01-cabinet-internal-detail"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-derived-drawing.ifc"
CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
PDFTOPPM = Path("/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/pdftoppm")
VIEWS = {
    "plan": "TRAP01-CABINET-INTERNAL-DETAIL-PLAN",
    "front": "TRAP01-CABINET-INTERNAL-DETAIL-FRONT",
    "side": "TRAP01-CABINET-INTERNAL-DETAIL-SIDE",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def record(path: Path) -> dict:
    return {"path": relative(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def svg_to_pdf(svg: Path, pdf: Path):
    with tempfile.TemporaryDirectory(prefix="trap01-detail-chrome-") as profile:
        process = subprocess.Popen(
            [
                str(CHROME), "--headless=new", "--disable-gpu", "--disable-background-networking",
                "--disable-component-update", "--no-first-run", "--no-default-browser-check",
                "--no-pdf-header-footer", f"--user-data-dir={profile}",
                f"--print-to-pdf={pdf}", svg.as_uri(),
            ],
            cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        for _ in range(450):
            if pdf.is_file() and pdf.stat().st_size:
                break
            if process.poll() is not None:
                break
            time.sleep(0.1)
        if process.poll() is None:
            process.terminate()
            process.wait(timeout=3)
    if not pdf.is_file() or not pdf.stat().st_size:
        raise RuntimeError(f"Chrome SVG-to-PDF failed: {svg}")


def inspect_preview(path: Path):
    image = Image.open(path).convert("RGB")
    dark = grey = blue = 0
    for red, green, value_blue in image.getdata():
        dark += int(max(red, green, value_blue) < 90)
        grey += int(abs(red - green) <= 20 and abs(green - value_blue) <= 20 and 90 <= red <= 225)
        blue += int(value_blue >= 110 and value_blue > red * 1.25 and value_blue > green * 1.03)
    return {
        "width": image.width,
        "height": image.height,
        "black_pixel_count": dark,
        "grey_context_pixel_count": grey,
        "blue_pixel_count": blue,
        "current_configuration_visible": dark >= 50,
        "retained_context_visible": grey >= 100,
        "blue_dashed_reference_absent": blue == 0,
        "pass": dark >= 50 and grey >= 100 and blue == 0,
    }


def main():
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch during TRAP01 detail packaging")
    evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
    if evidence.get("status") != "persisted_reloaded_verified":
        raise RuntimeError("TRAP01 detail persisted/reloaded evidence gate failed")
    by_view = {item["view"]: item for item in evidence["outputs"]["views"]}
    records, page_pdfs = [], []
    for view, name in VIEWS.items():
        source = by_view[view]
        svg = OUTPUT_DIR / f"{name}.svg"
        root = ET.parse(svg).getroot()
        if (
            source["create_result"] != ["FINISHED"]
            or source["svg"]["no_target_or_annotation_duplicate"] is not True
            or source["svg"]["all_occluders_absent"] is not True
            or root.attrib.get("data-blue-dashed-reference-displayed") != "false"
        ):
            raise RuntimeError(f"TRAP01 detail {view} scene SVG gate failed")
        pdf = OUTPUT_DIR / f"{name}.pdf"
        svg_to_pdf(svg, pdf)
        page_pdfs.append(pdf)
        prefix = PRODUCT_DIR / f"{PREVIEW_PREFIX}-{view}"
        subprocess.run(
            [str(PDFTOPPM), "-png", "-singlefile", "-r", "220", str(pdf), str(prefix)],
            cwd=ROOT, check=True, capture_output=True,
        )
        preview = prefix.with_suffix(".png")
        visual = inspect_preview(preview)
        if visual["pass"] is not True:
            raise RuntimeError(f"TRAP01 detail {view} visual gate failed: {visual}")
        records.append({
            "view": view,
            "drawing_global_id": source["drawing_global_id"],
            "annotation_global_id": source["annotation_global_id"],
            "svg": record(svg),
            "page_pdf": record(pdf),
            "preview": record(preview),
            "visual": visual,
            "svg_validation": source["svg"],
            "camera": source["camera"],
        })
    writer = PdfWriter()
    for page in page_pdfs:
        reader = PdfReader(str(page))
        if len(reader.pages) != 1:
            raise RuntimeError(f"expected one detail PDF page: {page}")
        writer.add_page(reader.pages[0])
    with COMBINED_PDF.open("wb") as stream:
        writer.write(stream)
    if len(PdfReader(str(COMBINED_PDF)).pages) != 3:
        raise RuntimeError("TRAP01 detail combined PDF must have three pages")
    manifest = {
        "schema_version": 1,
        "mode": "actual_bonsai_create_drawing_cabinet_internal_orthographic_detail",
        "formal_ifc_sha256": sha256(FORMAL_IFC),
        "formal_ifc_bytes_unchanged": sha256(FORMAL_IFC) == FORMAL_SHA256,
        "derived_ifc": record(DERIVED_IFC),
        "create_drawing_evidence": record(EVIDENCE),
        "current_installation_only": True,
        "blue_dashed_reference_displayed": False,
        "filtered_occluder_global_ids": [item["global_id"] for item in evidence["filters"]["excluded_occluders"]],
        "retained_context_element_count": evidence["filters"]["include_count"],
        "no_duplicate_target_or_annotation": all(item["svg_validation"]["no_target_or_annotation_duplicate"] for item in records),
        "all_three_current_configurations_visible": all(item["visual"]["current_configuration_visible"] for item in records),
        "views": records,
        "combined_pdf": record(COMBINED_PDF),
        "pass": all(item["visual"]["pass"] for item in records),
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    evidence["visual"] = {
        "preview_manifest": str(MANIFEST),
        "rendered_previews": [item["preview"] for item in records],
        "all_three_current_configurations_visible": manifest["all_three_current_configurations_visible"],
        "blue_dashed_reference_displayed": False,
    }
    evidence["outputs"]["rendered_views"] = records
    evidence["outputs"]["combined_pdf"] = record(COMBINED_PDF)
    evidence["verdict"] = "pass"
    evidence["pass"] = True
    evidence["status"] = "persisted_reloaded_rendered_verified"
    evidence["tests"]["all_three_rendered_previews_show_current_configuration"] = manifest["all_three_current_configurations_visible"]
    EVIDENCE.write_text(json.dumps(evidence, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    # Refresh evidence hash in the final manifest after its final write.
    manifest["create_drawing_evidence"] = record(EVIDENCE)
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": relative(MANIFEST), "combined_pdf": relative(COMBINED_PDF), "pass": manifest["pass"]}, indent=2))


if __name__ == "__main__":
    main()
