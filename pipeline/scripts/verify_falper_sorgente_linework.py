#!/usr/bin/env python3
"""Mechanically cross-check WFA/WFB native DWG linework against official PDF."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import (
    EXPECTED,
    ROOT,
    SCOPE,
    load_json,
    pathset_diff_mm,
    relative,
    sha256,
    write_json,
)


DEFAULT_DWG = ROOT / "pipeline/decisions/falper-sorgente-official-dwg-linework.json"
DEFAULT_PDF = ROOT / "pipeline/decisions/falper-sorgente-official-pdf-linework.json"
DEFAULT_OUTPUT = ROOT / "pipeline/decisions/falper-sorgente-dwg-pdf-verification.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dwg", type=Path, default=DEFAULT_DWG)
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--plan-tolerance-mm", type=float, default=0.75)
    parser.add_argument("--elevation-tolerance-mm", type=float, default=6.00)
    args = parser.parse_args()
    dwg = load_json(args.dwg.resolve())
    pdf = load_json(args.pdf.resolve())
    if dwg["source_kind"] != "native_dwg":
        raise RuntimeError("blue-line source must be native_dwg")
    if dwg["scope"] != SCOPE or pdf["scope"] != SCOPE:
        raise RuntimeError("Falper family-reference scope drifted")
    if pdf["source_pdf_sha256"] != EXPECTED["pdf"]:
        raise RuntimeError("official PDF hash drifted")
    records = []
    for code in ("WFA", "WFB"):
        expected_hash = EXPECTED[f"{code.lower()}_2d"]
        variant = dwg["variants"][code]
        if variant["source_dwg_sha256"] != expected_hash:
            raise RuntimeError(f"{code} DWG identity drifted")
        if variant["model_code"] != code or pdf["variants"][code]["model_code"] != code:
            raise RuntimeError(f"{code} identity mismatch between DWG and PDF")
        for view, tolerance in (
            ("plan", args.plan_tolerance_mm),
            ("elevation", args.elevation_tolerance_mm),
        ):
            diff = pathset_diff_mm(
                dwg["variants"][code]["views"][view]["paths_mm"],
                pdf["variants"][code]["views"][view]["paths_mm"],
            )
            record = {
                "model_code": code,
                "view": view,
                "dwg_source_kind": "native_dwg",
                "dwg_sha256": expected_hash,
                "pdf_sha256": EXPECTED["pdf"],
                "tolerance_mm": tolerance,
                "tolerance_basis": (
                    "published_vector_plan_curve_within_0.75_mm"
                    if view == "plan"
                    else "published_vector_elevation_curve_within_6_mm_or_0.67_percent_of_900_mm"
                ),
                **diff,
                "pass": diff["hausdorff_mm"] <= tolerance,
            }
            records.append(record)
    family_identity = {}
    for view in ("plan", "elevation"):
        family_identity[view] = pathset_diff_mm(
            dwg["variants"]["WFA"]["views"][view]["paths_mm"],
            dwg["variants"]["WFB"]["views"][view]["paths_mm"],
        )
    passed = all(record["pass"] for record in records)
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "mechanical_native_dwg_vs_official_vector_pdf",
        "scope": SCOPE,
        "dwg_linework_register": relative(args.dwg),
        "dwg_linework_register_sha256": sha256(args.dwg),
        "pdf_linework_register": relative(args.pdf),
        "pdf_linework_register_sha256": sha256(args.pdf),
        "records": records,
        "wfa_wfb_native_family_geometry_diff": family_identity,
        "identity_gates": {
            "wfa_2d_hash_match": True,
            "wfb_2d_hash_match": True,
            "pdf_hash_match": True,
            "wfa_3d_used_as_wfb": False,
            "project_review_model": "WFB",
            "blue_line_source_kind": "native_dwg",
        },
        "pass": passed,
    }
    write_json(args.output.resolve(), payload)
    if not passed:
        failures = [record for record in records if not record["pass"]]
        raise RuntimeError(f"PDF/DWG linework tolerance gate failed: {failures}")
    print(f"DWG/PDF linework verification passed: {relative(args.output)}")


if __name__ == "__main__":
    main()
