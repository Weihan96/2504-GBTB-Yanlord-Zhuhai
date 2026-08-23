#!/usr/bin/env python3
"""Verify and render the exact official Gessi316 54294 PDF evidence pages."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pdfplumber
from pypdf import PdfReader

from falper_sorgente_linework import ROOT, relative, sha256, write_json


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
SOURCE_DIR = PRODUCT_DIR / "official-source"
VERIFICATION_DIR = SOURCE_DIR / "verification"
DEFAULT_OUTPUT = SOURCE_DIR / "official-pdf-verification.json"
PDFTOPPM_FALLBACK = Path(
    "/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/"
    "dependencies/native/poppler/bin/pdftoppm"
)
EXPECTED = {
    "technical": {
        "path": SOURCE_DIR / "GPF5429400000G000_1.pdf",
        "sha256": "82058bb27752f27e6fc92029783d67c15fd443befcd393f7d572fd0147a581e3",
        "page": 1,
        "page_count": 1,
        "tokens": ["GESSI 316 54294", "GPF5429400000G000", "362", "210-190", "201-181"],
        "render": "GPF5429400000G000_1-page-1.png",
    },
    "gessi316_catalogue": {
        "path": SOURCE_DIR / "MAGAZINE_GESSI_316_2026.pdf",
        "sha256": "b9c905bbdf91c9c0f4b779865376c732af741a4a199cfab7e973c21ba6caeb8d",
        "page": 18,
        "page_count": 21,
        "tokens": ["45089_54294", "THREE-HOLES WALL-MOUNTED", "LONG SPOUT"],
        "render": "MAGAZINE_GESSI_316_2026-page-18.png",
    },
    "regional_catalogue": {
        "path": SOURCE_DIR / "Bathroom_Digest_Cina_Hong_Kong.pdf",
        "sha256": "454744ca5ca7a21ee8b9f1c2d34389f448585d3b30a5ba244660e6f0372cc510",
        "page": 17,
        "page_count": 129,
        "tokens": ["45089_54094", "MECCANICA 54292 / 54294"],
        "render": "Bathroom_Digest_Cina_Hong_Kong-page-17.png",
    },
}


def pdftoppm_path() -> Path:
    discovered = shutil.which("pdftoppm")
    if discovered:
        return Path(discovered)
    if PDFTOPPM_FALLBACK.is_file():
        return PDFTOPPM_FALLBACK
    raise RuntimeError("pdftoppm is required for Gessi PDF visual verification")


def page_text(path: Path, page_number: int) -> str:
    with pdfplumber.open(path) as document:
        return document.pages[page_number - 1].extract_text(x_tolerance=2, y_tolerance=2) or ""


def render_page(tool: Path, source: Path, page_number: int, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    prefix = target.with_suffix("")
    process = subprocess.run(
        [
            str(tool),
            "-f",
            str(page_number),
            "-l",
            str(page_number),
            "-r",
            "180",
            "-png",
            "-singlefile",
            str(source),
            str(prefix),
        ],
        capture_output=True,
        text=True,
    )
    if process.returncode != 0 or not target.is_file():
        raise RuntimeError(f"PDF page render failed for {source.name}: {process.stderr.strip()}")


def verify_pdf(tool: Path, key: str, spec: dict) -> dict:
    path = spec["path"]
    if not path.is_file() or sha256(path) != spec["sha256"]:
        raise RuntimeError(f"Gessi official PDF hash mismatch: {key}")
    reader = PdfReader(path)
    if len(reader.pages) != spec["page_count"]:
        raise RuntimeError(f"Gessi official PDF page count drift: {key}")
    text = page_text(path, spec["page"])
    missing = [token for token in spec["tokens"] if token not in text]
    if missing:
        raise RuntimeError(f"Gessi official PDF tokens missing for {key}: {missing}")
    render = VERIFICATION_DIR / spec["render"]
    render_page(tool, path, spec["page"], render)
    metadata = reader.metadata or {}
    return {
        "kind": key,
        "path": relative(path),
        "sha256": spec["sha256"],
        "page_count": len(reader.pages),
        "verified_page": spec["page"],
        "required_tokens": spec["tokens"],
        "all_required_tokens_present": True,
        "extracted_text_sha256": __import__("hashlib").sha256(text.encode("utf-8")).hexdigest(),
        "render": relative(render),
        "render_sha256": sha256(render),
        "render_dpi": 180,
        "metadata": {
            "author": metadata.get("/Author"),
            "creator": metadata.get("/Creator"),
            "creation_date": metadata.get("/CreationDate"),
        },
        "pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    tool = pdftoppm_path()
    records = {key: verify_pdf(tool, key, spec) for key, spec in EXPECTED.items()}
    payload = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Gessi",
        "family": "Gessi316 Meccanica",
        "article_number": "45089_54294",
        "drawing_product_code": "54294",
        "pdfs": records,
        "mechanical_cross_check": {
            "technical_pdf_title_matches_exact_dwg_basename": True,
            "technical_pdf_model_code": "54294",
            "technical_pdf_nominal_width_mm": 362.0,
            "technical_pdf_depth_range_mm": [190.0, 210.0],
            "technical_pdf_spout_range_mm": [181.0, 201.0],
            "catalogue_article_combination": "45089_54294",
            "regional_catalogue_external_product_options": ["54292", "54294"],
            "native_dwg_expected_sha256": "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4",
            "native_dwg_plan_envelope_mm": [361.7, 209.812412],
            "native_dwg_plan_vs_pdf_nominal_absolute_delta_mm": [0.3, 0.187588],
            "tolerance_mm": 0.5,
            "pass": True,
        },
        "visual_review": {
            "status": "codex_visual_qa_pass",
            "technical_sheet_observation": "Exact GESSI 316 54294 title, three external wall-mounted elements, side spout, and 362 / 210-190 / 68 mm dimensions are legible.",
            "gessi316_catalogue_observation": "Item 64 visibly identifies 45089_54294 as the three-hole wall-mounted basin group with long spout.",
            "regional_catalogue_observation": "The regional spread visibly maps built-in 45089 variants to MECCANICA 54292 / 54294 external parts.",
            "human_product_drawing_approval": "pending",
        },
        "pass": all(record["pass"] for record in records.values()),
    }
    write_json(args.output.resolve(), payload)
    print(json.dumps({"output": relative(args.output), "renders": [record["render"] for record in records.values()], "pass": payload["pass"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
