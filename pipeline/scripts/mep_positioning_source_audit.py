#!/usr/bin/env python3
"""Audit external kitchen MEP source drawings against current IFC inventories."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


EXPECTED_IFC_SHA256 = "7521c09991f3d0c7b7d91ca2324fd55ad961d8e32e9e3e9a9777a4cc19b06e81"
OLD_PDF_SHA256 = "4176a27a58f8c926ca3d66cbcb6015e679d7af0977fb8a26976980bf677a3488"
LATEST_PDF_SHA256 = "ec9f67ebe25fd6fe4be5ada495f919d6d92ff5e501f79d33a257eed59e4310a5"
EXPECTED_REQUIREMENT_IDS = {
    "MEP-SOURCE-001",
    "PLUM-SOURCE-B01",
    "PLUM-SOURCE-C01",
    "PLUM-SOURCE-D01",
    "PLUM-SOURCE-PLAN01",
    "ELEC-SOURCE-B01",
    "ELEC-SOURCE-C01",
    "ELEC-SOURCE-PLAN01",
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    source_dir = root.parent / "图纸"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument("--plum", type=Path, default=root / "build/plum/plum-report.json")
    parser.add_argument("--elec", type=Path, default=root / "build/elec/elec-existing-candidate.json")
    parser.add_argument("--register", type=Path, default=root / "pipeline/decisions/mep-positioning-source.csv")
    parser.add_argument("--old-pdf", type=Path, default=source_dir / "珠海-陈总(水电图)2025.5.29.pdf")
    parser.add_argument("--latest-pdf", type=Path, default=source_dir / "珠海-陈总(水电图)2025.12.23.pdf")
    parser.add_argument("--output", type=Path, default=root / "build/mep-positioning/source-audit.json")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def pdf_page_count(path: Path) -> int:
    output = subprocess.check_output(["pdfinfo", str(path)], text=True)
    for line in output.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"pdfinfo did not report a page count for {path}")


def read_register(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ids = [row["requirement_id"] for row in rows]
    if set(ids) != EXPECTED_REQUIREMENT_IDS or len(ids) != len(set(ids)):
        raise RuntimeError("MEP source register is incomplete or has duplicate IDs")
    required = ("source_sha256", "requirement", "basis", "confidence", "review_required", "status", "stop_condition")
    if any(not row[field] for row in rows for field in required):
        raise RuntimeError("MEP source register contains an incomplete evidence row")
    return rows


def main() -> int:
    args = parse_args()
    ifc_hash = sha256(args.ifc)
    old_hash = sha256(args.old_pdf)
    latest_hash = sha256(args.latest_pdf)
    if ifc_hash != EXPECTED_IFC_SHA256:
        raise RuntimeError(f"formal IFC hash changed: {ifc_hash}")
    if old_hash != OLD_PDF_SHA256 or latest_hash != LATEST_PDF_SHA256:
        raise RuntimeError("water/electric source PDF hash changed; re-review revisions")

    plum = read_json(args.plum)
    elec = read_json(args.elec)
    rows = read_register(args.register)
    if plum["source"]["ifc_sha256"] != ifc_hash or elec["source"]["sha256"] != ifc_hash:
        raise RuntimeError("PLUM or ELEC inventory is stale relative to the formal IFC")

    e301 = elec["sheets"]["E-301"]
    e303 = elec["sheets"]["E-303"]
    report = {
        "source_ifc_sha256": ifc_hash,
        "source_drawings": {
            "older": {"path": str(args.old_pdf), "sha256": old_hash, "pages": pdf_page_count(args.old_pdf)},
            "latest": {"path": str(args.latest_pdf), "sha256": latest_hash, "pages": pdf_page_count(args.latest_pdf)},
            "mechanical_page_difference": [6],
            "difference_review": "latest page 6 adds a red circle around the hood outlet; pages 1-5 and 7-10 render identically at 90 dpi",
        },
        "source_requirements": rows,
        "current_ifc_inventory": {
            "plum_registered_terminals": plum["inventory"]["p201_registered_terminal_count"],
            "plum_service_demand_candidates": plum["inventory"]["p201_service_demand_candidate_count"],
            "plum_non_service_components": plum["inventory"]["p201_non_service_component_count"],
            "plum_existing_location_objects": plum["inventory"]["p202_existing_object_count"],
            "lights": len(e301["lights"]),
            "sockets": len(e303["sockets"]),
            "typed_equipment": len(e303["typed_equipment"]),
            "unresolved_electrical_proxies": len(e303["proxy_handoffs"]),
            "switch_instances": len(elec["sheets"]["E-302"]["instances"]),
            "network_instances": len(elec["sheets"]["E-304"]["instances"]),
        },
        "gates": {
            "latest_kitchen_source_registered": True,
            "current_ifc_reports_match": True,
            "whole_home_water_positioning_complete": False,
            "whole_home_electrical_positioning_complete": False,
            "automatic_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"current_ifc_inventory": report["current_ifc_inventory"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
