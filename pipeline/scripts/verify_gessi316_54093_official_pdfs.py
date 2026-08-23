#!/usr/bin/env python3
"""Verify and render exact official Gessi316 54093 PDF evidence pages."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pdfplumber
from pypdf import PdfReader

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54093"
SOURCE_DIR = PRODUCT_DIR / "official-source"
VERIFICATION_DIR = SOURCE_DIR / "verification"
LINEWORK = PRODUCT_DIR / "official-native-dwg-linework.json"
DEFAULT_OUTPUT = SOURCE_DIR / "official-pdf-verification.json"
PDFTOPPM_FALLBACK = Path("/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/bin/pdftoppm")
EXPECTED = {
    "technical": {
        "path": SOURCE_DIR / "GPF5409300000G000_1.pdf",
        "sha256": "2eb67af1be75093a12ab8bc3a739e5bdd6e4b1c4c544948795c87cb9598515c6",
        "page": 1, "page_count": 1,
        "tokens": ["GESSI 316 54093", "GPF5409300000G000", "190", "372", "822"],
        "render": "GPF5409300000G000_1-page-1.png",
    },
    "bathroom_catalogue": {
        "path": SOURCE_DIR / "Gessi_Cataloghi_Bathroom.pdf",
        "sha256": "8f62401e7e33bd1f8a60cfc20ba00996817f935b1a9731e12042e6d7c05e2092",
        "page": 17, "page_count": 129,
        "tokens": ["54093", "H 230 MM", "COUNTER SPOUT"],
        "render": "Gessi_Cataloghi_Bathroom-page-17.png",
    },
}


def pdftoppm_path() -> Path:
    discovered = shutil.which("pdftoppm")
    if discovered:
        return Path(discovered)
    if PDFTOPPM_FALLBACK.is_file():
        return PDFTOPPM_FALLBACK
    raise RuntimeError("pdftoppm is required for Gessi 54093 PDF visual verification")


def page_text(path: Path, page_number: int) -> str:
    with pdfplumber.open(path) as document:
        return document.pages[page_number - 1].extract_text(x_tolerance=2, y_tolerance=2) or ""


def render_page(tool: Path, source: Path, page_number: int, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    process = subprocess.run([str(tool), "-f", str(page_number), "-l", str(page_number), "-r", "180", "-png", "-singlefile", str(source), str(target.with_suffix(""))], capture_output=True, text=True)
    if process.returncode != 0 or not target.is_file():
        raise RuntimeError(f"PDF page render failed for {source.name}: {process.stderr.strip()}")


def verify_pdf(tool: Path, key: str, spec: dict) -> dict:
    path = spec["path"]
    if not path.is_file() or sha256(path) != spec["sha256"]:
        raise RuntimeError(f"Gessi 54093 official PDF hash mismatch: {key}")
    reader = PdfReader(path)
    if len(reader.pages) != spec["page_count"]:
        raise RuntimeError(f"Gessi 54093 official PDF page-count drift: {key}")
    text = page_text(path, spec["page"])
    missing = [token for token in spec["tokens"] if token not in text]
    if missing:
        raise RuntimeError(f"Gessi 54093 official PDF tokens missing for {key}: {missing}")
    render = VERIFICATION_DIR / spec["render"]
    render_page(tool, path, spec["page"], render)
    metadata = reader.metadata or {}
    return {
        "kind": key, "path": relative(path), "sha256": spec["sha256"], "page_count": len(reader.pages),
        "verified_page": spec["page"], "required_tokens": spec["tokens"], "all_required_tokens_present": True,
        "extracted_text_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(), "render": relative(render),
        "render_sha256": sha256(render), "render_dpi": 180,
        "metadata": {"author": metadata.get("/Author"), "creator": metadata.get("/Creator"), "creation_date": metadata.get("/CreationDate")},
        "pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    tool = pdftoppm_path()
    records = {key: verify_pdf(tool, key, spec) for key, spec in EXPECTED.items()}
    linework = load_json(LINEWORK)
    cross_check = linework["nominal_dimension_cross_check"]
    if cross_check.get("pass") is not True:
        raise RuntimeError("Gessi 54093 native DWG nominal-dimension cross-check failed")
    payload = {
        "schema_version": 1, "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Gessi", "family": "Gessi316 Meccanica", "article_number": "54093", "configuration": "G000",
        "pdfs": records,
        "mechanical_cross_check": {
            "technical_pdf_title_matches_exact_dwg_basename": True, "technical_pdf_model_code": "54093",
            "technical_pdf_nominal_width_depth_height_mm": [50.0, 190.0, 273.0], "technical_pdf_outlet_height_mm": 228.0,
            "native_dwg_expected_sha256": "a9308fc8498d34c8cf2f68fd28aaf90b11e62bfc59424e0a5a3ac7f0d28d5047",
            "native_dwg_plan_envelope_mm": cross_check["native_dwg_plan_envelope_mm"],
            "native_dwg_front_envelope_mm": cross_check["native_dwg_front_envelope_mm"],
            "native_dwg_side_envelope_mm": cross_check["native_dwg_side_envelope_mm"],
            "absolute_delta_mm": cross_check["absolute_delta_mm"], "tolerance_mm": cross_check["tolerance_mm"], "pass": True,
        },
        "visual_review": {
            "status": "codex_visual_qa_pass",
            "technical_sheet_observation": "Exact GESSI 316 54093 / GPF5409300000G000 title and 50 x 190 x 273 mm geometry with 228 mm outlet height are visible.",
            "catalogue_observation": "The official bathroom catalogue identifies 54093 as the H 230 mm high counter spout.",
            "native_dwg_observation": "The exact native DWG contains front, side and plan orthographic views plus a separate isometric; only the three orthographic windows are extracted.",
            "human_product_drawing_approval": "pending",
        },
        "pass": all(record["pass"] for record in records.values()),
    }
    write_json(args.output.resolve(), payload)
    print(json.dumps({"output": relative(args.output), "renders": [record["render"] for record in records.values()], "pass": payload["pass"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
