#!/usr/bin/env python3
"""Render the corrected project-axis-aware TRAP01 detail Drawings."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
import generate_trap01_internal_detail_drawings as generator  # noqa: E402


PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
generator.OUTPUT_DIR = PRODUCT_DIR / "bonsai-drawings/cabinet-internal-detail-v2"
generator.EVIDENCE = generator.OUTPUT_DIR / "TRAP01-cabinet-internal-detail-v2-create-drawing-evidence.json"
generator.MANIFEST = PRODUCT_DIR / "TRAP01-cabinet-internal-detail-v2-manifest.json"
generator.COMBINED_PDF = PRODUCT_DIR / "Geberit-151.116.11.1-TRAP01-cabinet-internal-detail-v2.pdf"
generator.PREVIEW_PREFIX = "trap01-cabinet-internal-detail-v2"
generator.VIEWS = {
    "plan": "TRAP01-CABINET-INTERNAL-DETAIL-V2-PLAN",
    "front": "TRAP01-CABINET-INTERNAL-DETAIL-V2-FRONT",
    "side": "TRAP01-CABINET-INTERNAL-DETAIL-V2-SIDE",
}


if __name__ == "__main__":
    generator.main()
