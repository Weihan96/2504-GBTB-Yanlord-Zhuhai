#!/usr/bin/env python3
"""Extract WFA/WFB linework by decoding the archived native DWG files."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import (
    AUTOCAD_2D_URL,
    EXPECTED,
    PRODUCT_URL,
    ROOT,
    SCOPE,
    TECHNICAL_DWG_URL,
    extract_native_dwg_variant,
    relative,
    sha256,
    write_json,
)


DEFAULT_WFA = ROOT / "drawings/evidence/FALPER-official-Sorgente-WFA-2D.dwg"
DEFAULT_WFB = ROOT / "drawings/evidence/FALPER-official-Sorgente-WFB-2D.dwg"
DEFAULT_ARCHIVE = ROOT / "drawings/evidence/FALPER-official-Lavabi-Freestanding-Autocad-2D.zip"
DEFAULT_OUTPUT = ROOT / "pipeline/decisions/falper-sorgente-official-dwg-linework.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wfa", type=Path, default=DEFAULT_WFA)
    parser.add_argument("--wfb", type=Path, default=DEFAULT_WFB)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    archive_hash = sha256(args.archive)
    variants = {
        "WFA": extract_native_dwg_variant(args.wfa.resolve(), "WFA"),
        "WFB": extract_native_dwg_variant(args.wfb.resolve(), "WFB"),
    }
    if variants["WFB"]["source_dwg_sha256"] == EXPECTED["wfa_3d"]:
        raise RuntimeError("WFA 3D must never be used as WFB drawing linework")
    payload = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source_kind": "native_dwg",
        "scope": SCOPE,
        "units": "mm",
        "product_url": PRODUCT_URL,
        "autocad_2d_zip_url": AUTOCAD_2D_URL,
        "technical_dwg_zip_url": TECHNICAL_DWG_URL,
        "archived_autocad_2d_zip": relative(args.archive),
        "archived_autocad_2d_zip_sha256": archive_hash,
        "variants": variants,
        "identity_guards": {
            "wfa_2d_sha256": EXPECTED["wfa_2d"],
            "wfb_2d_sha256": EXPECTED["wfb_2d"],
            "wfa_3d_sha256_identity_only_not_geometry_source": EXPECTED["wfa_3d"],
            "project_review_model": "WFB",
            "wfa_3d_used_as_wfb": False,
        },
    }
    write_json(args.output.resolve(), payload)
    print(
        f"native DWG linework: {relative(args.output)}; "
        f"WFA={variants['WFA']['source_dwg_sha256']}; "
        f"WFB={variants['WFB']['source_dwg_sha256']}"
    )


if __name__ == "__main__":
    main()
