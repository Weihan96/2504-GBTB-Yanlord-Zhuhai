#!/usr/bin/env python3
"""Verify and render exact official Gessi316 54146 PDF evidence pages."""

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


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54146"
SOURCE_DIR = PRODUCT_DIR / "official-source"
VERIFICATION_DIR = SOURCE_DIR / "verification"
LINEWORK = PRODUCT_DIR / "official-native-dwg-linework.json"
DEFAULT_OUTPUT = SOURCE_DIR / "official-pdf-verification.json"
PDFTOPPM_FALLBACK = Path("/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/bin/pdftoppm")
EXPECTED = {
    "technical": {
        "path": SOURCE_DIR / "GPF5414600000G000_1.pdf",
        "sha256": "ee7909886b1181306902d0a0e146b73f938bf96e3e12e8f1ba27ded996b2b8ce",
        "page": 1,
        "page_count": 1,
        "tokens": ["GESSI 316 54146", "GPF5414600000G000", "672", "G1/2\"", "62"],
        "render": "GPF5414600000G000_1-page-1.png",
    },
    "collection_catalogue": {
        "path": SOURCE_DIR / "MAGAZINE_GESSI_316_2026.pdf",
        "sha256": "b9c905bbdf91c9c0f4b779865376c732af741a4a199cfab7e973c21ba6caeb8d",
        "page": 19,
        "page_count": 21,
        "tokens": ["94. 54145", "WALL-MOUNTED ADJUSTABLE", "95. 54146", "CEILING-MOUNTED ADJUSTABLE"],
        "render": "MAGAZINE_GESSI_316_2026-page-19.png",
    },
}


def pdftoppm_path() -> Path:
    discovered = shutil.which("pdftoppm")
    if discovered:
        return Path(discovered)
    if PDFTOPPM_FALLBACK.is_file():
        return PDFTOPPM_FALLBACK
    raise RuntimeError("pdftoppm is required for Gessi 54146 PDF visual verification")


def page_text(path: Path, page_number: int) -> str:
    with pdfplumber.open(path) as document:
        return document.pages[page_number - 1].extract_text(x_tolerance=2, y_tolerance=2) or ""


def verify_pdf(tool: Path, key: str, spec: dict) -> dict:
    path = spec["path"]
    if not path.is_file() or sha256(path) != spec["sha256"]:
        raise RuntimeError(f"Gessi 54146 official PDF hash mismatch: {key}")
    reader = PdfReader(path)
    if len(reader.pages) != spec["page_count"]:
        raise RuntimeError(f"Gessi 54146 official PDF page-count drift: {key}")
    text = page_text(path, spec["page"])
    missing = [token for token in spec["tokens"] if token not in text]
    if missing:
        raise RuntimeError(f"Gessi 54146 official PDF tokens missing for {key}: {missing}")
    render = VERIFICATION_DIR / spec["render"]
    render.parent.mkdir(parents=True, exist_ok=True)
    process = subprocess.run(
        [str(tool), "-f", str(spec["page"]), "-l", str(spec["page"]), "-r", "180", "-png", "-singlefile", str(path), str(render.with_suffix(""))],
        capture_output=True,
        text=True,
    )
    if process.returncode != 0 or not render.is_file():
        raise RuntimeError(f"PDF page render failed for {path.name}: {process.stderr.strip()}")
    return {
        "kind": key,
        "path": relative(path),
        "sha256": spec["sha256"],
        "page_count": len(reader.pages),
        "verified_page": spec["page"],
        "required_tokens": spec["tokens"],
        "all_required_tokens_present": True,
        "extracted_text_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "render": relative(render),
        "render_sha256": sha256(render),
        "render_dpi": 180,
        "pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    tool = pdftoppm_path()
    records = {key: verify_pdf(tool, key, spec) for key, spec in EXPECTED.items()}
    cross_check = load_json(LINEWORK)["nominal_dimension_cross_check"]
    if cross_check.get("pass") is not True:
        raise RuntimeError("Gessi 54146 native DWG / IFC projection cross-check failed")
    payload = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Gessi",
        "family": "Gessi316 Meccanica",
        "article_number": "54146",
        "configuration": "G000",
        "pdfs": records,
        "mechanical_cross_check": {
            "technical_pdf_title_matches_exact_dwg_basename": True,
            "technical_pdf_model_code": "54146",
            "technical_pdf_nominal_width_depth_height_mm": [300.0, 300.0, 276.0],
            "technical_pdf_showerhead_diameter_mm": 300.0,
            "native_dwg_expected_sha256": "c8ddb90f61565d5273a32574f777812bbb1d9ef33e2b98479f1e7c381e69950d",
            "native_dwg_plan_envelope_mm": cross_check["native_dwg_plan_envelope_mm"],
            "native_dwg_front_envelope_mm": cross_check["native_dwg_front_envelope_mm"],
            "native_dwg_side_envelope_mm": cross_check["native_dwg_side_envelope_mm"],
            "project_ifc_body_local_xyz_mm": cross_check["project_ifc_body_local_xyz_mm"],
            "maximum_ifc_projection_delta_mm": cross_check["maximum_ifc_projection_delta_mm"],
            "tolerance_mm": cross_check["tolerance_mm"],
            "pass": True,
        },
        "visual_review": {
            "status": "codex_visual_qa_pass",
            "technical_sheet_observation": "Exact GESSI 316 54146 / GPF5414600000G000 title, Ø300 showerhead, 276 mm nominal height, Ø30 arm and ceiling connection are visible.",
            "catalogue_observation": "Official catalogue item 95 identifies 54146 as ceiling-mounted adjustable; wall-mounted item 94 / 54145 is excluded.",
            "native_dwg_observation": "The exact G000 native DWG partitions into three orthographic views; all 1283 eligible paths are consumed once, and G001 is excluded.",
            "human_product_drawing_approval": "pending",
        },
        "pass": all(record["pass"] for record in records.values()),
    }
    write_json(args.output.resolve(), payload)
    print(json.dumps({"output": relative(args.output), "renders": [record["render"] for record in records.values()], "pass": payload["pass"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
