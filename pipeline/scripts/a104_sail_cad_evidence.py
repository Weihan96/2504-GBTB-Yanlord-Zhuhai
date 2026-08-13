#!/usr/bin/env python3
"""Mechanically audit the official Rimadesio Sail monorail CAD against M05/M06.

The official CAD is a product-family drawing with several generic examples. It
is evidence for the Sail monorail system, but it is not a project shop drawing.
This audit therefore records the generic dimensions and proves why they cannot
be copied into IfcDoor.OverallWidth/OverallHeight for the project occurrences.
It never writes the IFC.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

import ifcopenshell

from a103_wall_plan_candidate import geometry_settings, world_bbox_mm
from e304_router_cad_evidence import parse_dxf, text_value


ROOT = Path(__file__).resolve().parents[2]
EXPECTED_DWG_SHA256 = "4eeaeb27b0668bb22b9b4fb191091e07eabbb2b259d90cedcdc45584a1f68e52"
EXPECTED_GENERIC_DIMENSIONS = {
    "panel_width_mm": [1000],
    "opening_width_mm": [976, 989, 1978],
    "rail_width_mm": [2011, 2037, 4022],
    "panel_height_mm": [2670],
    "opening_height_mm": [2666, 2678],
}
M05_GUID = "2D5BPoo2XFSvhTdfPenCh7"
M06_GUID = "0zjVS5FBbBewgUkk0fdfiv"
GROUP_NAME = "M05/M06 Rimadesio Sail 门组"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dwg",
        type=Path,
        default=ROOT / "drawings/evidence/RIMADESIO-official-Sail-monorotaia.dwg",
    )
    parser.add_argument(
        "--dxf",
        type=Path,
        default=ROOT / "drawings/evidence/RIMADESIO-Sail-monorotaia.dxf",
    )
    parser.add_argument(
        "--ifc",
        type=Path,
        default=ROOT / "2504 GBTB Yanlord Zhuhai.ifc",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "drawings/evidence/RIMADESIO-Sail-monorotaia-mechanical-audit.json",
    )
    parser.add_argument("--expected-ifc-sha256")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean_mtext(value: str) -> str:
    return re.sub(r"\\P", " ", value, flags=re.IGNORECASE).replace("  ", " ").strip()


def values_for(pattern: str, texts: list[str]) -> list[int]:
    values: set[int] = set()
    compiled = re.compile(pattern, re.IGNORECASE)
    for value in texts:
        match = compiled.search(value)
        if match:
            values.add(int(match.group(1)))
    return sorted(values)


def generic_dimensions(texts: list[str]) -> dict[str, list[int]]:
    return {
        "panel_width_mm": values_for(r"PANEL\s+W\s+(\d+)", texts),
        "opening_width_mm": values_for(r"OPENING\s+W\s+(\d+)", texts),
        "rail_width_mm": values_for(r"RAIL\s+W(?:\s+min)?\s+(\d+)", texts),
        "panel_height_mm": values_for(r"PANEL\s+H.*=\s*(\d+)\s*$", texts),
        "opening_height_mm": values_for(r"H\s+VANO\s+(\d+)", texts),
    }


def exact_group_members(model: ifcopenshell.file) -> list[str]:
    groups = [group for group in model.by_type("IfcGroup") if str(group.Name or "") == GROUP_NAME]
    if len(groups) != 1:
        raise RuntimeError(f"expected one {GROUP_NAME!r}, found {len(groups)}")
    return sorted(
        member.GlobalId
        for relation in (groups[0].IsGroupedBy or ())
        for member in relation.RelatedObjects
    )


def product_record(settings: Any, product: ifcopenshell.entity_instance) -> dict[str, Any]:
    bbox = world_bbox_mm(settings, product)
    return {
        "global_id": product.GlobalId,
        "ifc_class": product.is_a(),
        "name": str(product.Name or ""),
        "tag": str(product.Tag or ""),
        "overall_width_mm": product.OverallWidth,
        "overall_height_mm": product.OverallHeight,
        "world_bbox_mm": bbox,
    }


def main() -> None:
    args = parse_args()
    if not args.dxf.exists():
        raise RuntimeError(
            f"verified DXF conversion is missing: {args.dxf}; restore the committed AutoCAD conversion"
        )
    dwg_sha = sha256(args.dwg)
    if dwg_sha != EXPECTED_DWG_SHA256:
        raise RuntimeError(f"official Sail DWG hash drifted: {dwg_sha}")
    ifc_sha = sha256(args.ifc)
    if args.expected_ifc_sha256 and ifc_sha != args.expected_ifc_sha256:
        raise RuntimeError(
            f"formal IFC SHA-256 mismatch: expected {args.expected_ifc_sha256}, found {ifc_sha}"
        )

    entities = parse_dxf(args.dxf)
    texts = sorted(
        {
            clean_mtext(text_value(entity))
            for entity in entities
            if entity["type"] in {"TEXT", "MTEXT"} and clean_mtext(text_value(entity))
        }
    )
    dimensions = generic_dimensions(texts)
    if dimensions != EXPECTED_GENERIC_DIMENSIONS:
        raise RuntimeError(
            "Sail generic dimension extraction drifted: "
            + json.dumps(dimensions, ensure_ascii=False, sort_keys=True)
        )

    model = ifcopenshell.open(args.ifc)
    if model.schema != "IFC4":
        raise RuntimeError(f"expected IFC4, found {model.schema}")
    members = exact_group_members(model)
    if members != sorted([M05_GUID, M06_GUID]):
        raise RuntimeError(f"Sail group membership drifted: {members}")
    settings = geometry_settings()
    m05 = product_record(settings, model.by_guid(M05_GUID))
    m06 = product_record(settings, model.by_guid(M06_GUID))
    if any(
        value is not None
        for value in (
            m05["overall_width_mm"],
            m05["overall_height_mm"],
            m06["overall_width_mm"],
            m06["overall_height_mm"],
        )
    ):
        raise RuntimeError("M05/M06 nominal dimensions were filled without a project shop drawing")

    m05_dims = m05["world_bbox_mm"]["dimensions_mm"]
    m06_dims = m06["world_bbox_mm"]["dimensions_mm"]
    comparisons = {
        "m05_panel_width": {
            "ifc_observed_mm": m05_dims[0],
            "official_generic_examples_mm": dimensions["panel_width_mm"],
            "numeric_match": abs(m05_dims[0] - 1000.0) <= 0.1,
            "project_nominal_dimension_proven": False,
        },
        "m05_panel_height": {
            "ifc_observed_mm": m05_dims[2],
            "official_generic_examples_mm": dimensions["panel_height_mm"],
            "numeric_match": any(abs(m05_dims[2] - value) <= 0.1 for value in dimensions["panel_height_mm"]),
            "project_nominal_dimension_proven": False,
        },
        "m06_rail_long_axis": {
            "ifc_observed_mm": max(m06_dims[:2]),
            "official_generic_examples_mm": dimensions["rail_width_mm"],
            "numeric_match": any(abs(max(m06_dims[:2]) - value) <= 0.1 for value in dimensions["rail_width_mm"]),
            "project_nominal_dimension_proven": False,
        },
    }

    counts = Counter(entity["type"] for entity in entities)
    report = {
        "mode": "read-only-official-cad-to-ifc-audit",
        "source": {
            "official_dwg": str(args.dwg.resolve()),
            "official_dwg_sha256": dwg_sha,
            "verified_conversion_dxf": str(args.dxf.resolve()),
            "verified_conversion_dxf_sha256": sha256(args.dxf),
            "formal_ifc": str(args.ifc.resolve()),
            "formal_ifc_sha256": ifc_sha,
            "ifc_schema": model.schema,
        },
        "cad_inventory": {
            "entity_count": len(entities),
            "entity_type_counts": dict(sorted(counts.items())),
            "distinct_annotation_count": len(texts),
            "generic_dimensions": dimensions,
        },
        "ifc_occurrences": {"M05": m05, "M06": m06},
        "ifc_group": {"name": GROUP_NAME, "members": members},
        "comparisons": comparisons,
        "evidence_boundary": {
            "proves": [
                "Rimadesio Sail MONOROTAIA is a single-track sliding system family",
                "the official CAD contains generic one-panel and two-panel installation examples",
                "M05 and M06 are already the exact two-member project Sail group",
            ],
            "does_not_prove": [
                "project opening width or height",
                "project panel nominal width or height",
                "project rail nominal length",
                "sliding direction, host opening, rail substrate, closure detail, or installation tolerance",
            ],
            "formal_ifc_write_allowed": False,
            "dimension_closeout_rule": "retain blank OverallWidth/OverallHeight until a project order, door schedule, or manufacturer shop drawing identifies M05/M06",
        },
        "gates": {
            "official_dwg_hash_verified": True,
            "generic_dimensions_extracted": True,
            "ifc_group_exact": True,
            "ifc_nominal_dimensions_still_blank": True,
            "mechanical_pass": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "dwg_sha256": dwg_sha,
                "dxf_sha256": report["source"]["verified_conversion_dxf_sha256"],
                "ifc_sha256": ifc_sha,
                "entity_count": len(entities),
                "generic_dimensions": dimensions,
                "formal_ifc_write_allowed": False,
                "mechanical_pass": True,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
