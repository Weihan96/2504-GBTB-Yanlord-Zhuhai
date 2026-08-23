#!/usr/bin/env python3
"""Audit official antoniolupi Street DXF clusters against the custom project STREET Body."""

from __future__ import annotations

import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, relative, sha256, write_json
from street_h_linework import bounds, entity_path, first, records


PRODUCT_DIR = ROOT / "output/review/highpoly-types/street"
SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_ZIP = SOURCE_DIR / "ANTONIOLUPI-official-Street-2D-CAD.zip"
SOURCE_DXF = SOURCE_DIR / "AL_Street.dxf"
SOURCE_TECHNICAL_PDF = SOURCE_DIR / "ANTONIOLUPI-official-Street-technical.pdf"
SOURCE_CATALOGUE = SOURCE_DIR / "ANTONIOLUPI-official-Street-LevantoRed-extract.pdf"
PAGE_EVIDENCE = SOURCE_DIR / "official-product-page-evidence.json"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dxf-configuration-audit.json"
PRODUCT_PAGE = "https://www.antoniolupi.it/en/products/sinks/street"
DRAWING_UNITS_TO_MM = 10.0
EXPECTED = {
    "zip": "36161ac81bade86b7ec0419c0d9cf0c49d52ae0ecb9dced39006b9eddcb10e62",
    "dxf": "72b0a893545c4d9ff2a1d3fc79b7618eb10c00c651533eed36071e4130277c7e",
    "technical_pdf": "0b883bbf6df617331cc2408da4860b31160af67a6b6924b3caaee0ec389d2946",
    "catalogue": "6dd192c3da88c955c2357fc9a7e1e97442838f43393e4b343376d8b5dbacc2da",
    "page_evidence": "6db178b9bb6a405ebdc54401b73d1ff510e9a9c0dc27f1c18c85807579f499d0",
}
PROJECT_BODY_BOUNDS_MM = [1000.0, 470.0, 250.0]
SCOPE = (
    "official antoniolupi Street family CAD evidence for a domestic-custom project top; "
    "no exact 1000 x 470 x 250 mm official configuration is present and no project shop drawing is claimed"
)


def normalised_insert(entities, blocks, handle: str, block_name: str) -> dict:
    insertions = [item for item in entities if item["type"] == "INSERT" and first(item, "5") == handle]
    if len(insertions) != 1 or first(insertions[0], "2") != block_name:
        raise RuntimeError(f"Street insert {handle}/{block_name} identity drifted")
    source_entities = [item for item in blocks[block_name] if first(item, "8", "0") == "0"]
    raw_paths = [path for item in source_entities if (path := entity_path(item))]
    native = bounds(raw_paths)
    paths = [
        [
            [
                round((point[0] - native["minimum"][0]) * DRAWING_UNITS_TO_MM, 6),
                round((point[1] - native["minimum"][1]) * DRAWING_UNITS_TO_MM, 6),
            ]
            for point in path
        ]
        for path in raw_paths
    ]
    return {
        "native_insert_handle": handle,
        "native_block_name": block_name,
        "native_entity_counts": dict(sorted(Counter(item["type"] for item in source_entities if entity_path(item)).items())),
        "native_bounds_drawing_units": native,
        "bounds_mm": bounds(paths),
        "path_count": len(paths),
        "paths_mm": paths,
    }


def label(entities, handle: str, expected: str) -> dict:
    matches = [item for item in entities if item["type"] == "MTEXT" and first(item, "5") == handle]
    if len(matches) != 1 or first(matches[0], "1", "").lower() != expected.lower():
        raise RuntimeError(f"Street label {handle} identity drifted")
    return {"native_label_handle": handle, "label": first(matches[0], "1")}


def all_depth_47_widths(entities, blocks) -> list[float]:
    widths = set()
    for item in entities:
        if item["type"] != "INSERT":
            continue
        block_name = first(item, "2")
        raw_paths = [
            path
            for source in blocks.get(block_name, [])
            if first(source, "8", "0") == "0" and (path := entity_path(source))
        ]
        if not raw_paths:
            continue
        size = bounds(raw_paths)["size"]
        if abs(size[1] - 47.0) <= 1e-6:
            widths.add(round(size[0] * DRAWING_UNITS_TO_MM, 6))
    return sorted(widths)


def main() -> None:
    if (
        sha256(SOURCE_ZIP) != EXPECTED["zip"]
        or sha256(SOURCE_DXF) != EXPECTED["dxf"]
        or sha256(SOURCE_TECHNICAL_PDF) != EXPECTED["technical_pdf"]
        or sha256(SOURCE_CATALOGUE) != EXPECTED["catalogue"]
        or sha256(PAGE_EVIDENCE) != EXPECTED["page_evidence"]
    ):
        raise RuntimeError("Street official evidence hash mismatch")
    entities, blocks = records(SOURCE_DXF)
    described = {
        **label(entities, "1388", "street240 prof. 40 + street4054 prof. 40"),
        "plan": normalised_insert(entities, blocks, "1376", "*U4"),
        "front": normalised_insert(entities, blocks, "1369", "*U0"),
        "configuration_bounds_mm": [1080.0, 400.0, 250.0],
        "matches_project_body": False,
    }
    nearest_depth = {
        **label(entities, "13A2", "street147 prof. 47 + street4754 prof. 40"),
        "plan": normalised_insert(entities, blocks, "13A0", "*U15"),
        "configuration_bounds_mm": [1080.0, 470.0, 250.0],
        "matches_project_body": False,
    }
    widths = all_depth_47_widths(entities, blocks)
    if described["plan"]["bounds_mm"]["size"] != [1080.0, 400.0]:
        raise RuntimeError("Street described cluster plan dimensions drifted")
    if described["front"]["bounds_mm"]["size"] != [1080.0, 250.0]:
        raise RuntimeError("Street described cluster front dimensions drifted")
    if nearest_depth["plan"]["bounds_mm"]["size"] != [1080.0, 470.0]:
        raise RuntimeError("Street nearest-depth cluster dimensions drifted")
    if 1000.0 in widths:
        raise RuntimeError("official DXF unexpectedly contains an exact 1000 x 470 mm plan cluster")
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/street_linework.py",
        "manufacturer": "antoniolupi",
        "family": "Street",
        "project_ifc_type_name": "STREET",
        "project_ifc_type_description": "antoniolupi street240 prof. 40 + street4054 prof. 40",
        "product_page": PRODUCT_PAGE,
        "source_kind": "native_dxf_configuration_audit",
        "source_zip": relative(SOURCE_ZIP),
        "source_zip_sha256": sha256(SOURCE_ZIP),
        "source_dxf": relative(SOURCE_DXF),
        "source_dxf_sha256": sha256(SOURCE_DXF),
        "source_technical_pdf": relative(SOURCE_TECHNICAL_PDF),
        "source_technical_pdf_sha256": sha256(SOURCE_TECHNICAL_PDF),
        "native_drawing_units_to_mm": DRAWING_UNITS_TO_MM,
        "project_body_bounds_mm": PROJECT_BODY_BOUNDS_MM,
        "ifc_description_cluster": described,
        "nearest_official_depth_cluster": nearest_depth,
        "all_official_depth_470_plan_widths_mm": widths,
        "exact_1000_x_470_plan_cluster_present": False,
        "exact_project_configuration_match": False,
        "official_dxf_paths_used_as_project_representation": False,
        "drawing_geometry_source": "geometry_derived_simplified_proxy",
        "scope": SCOPE,
        "pass": True,
    }
    write_json(DEFAULT_OUTPUT, payload)
    print(json.dumps({"output": relative(DEFAULT_OUTPUT), "project": PROJECT_BODY_BOUNDS_MM, "described": described["configuration_bounds_mm"], "nearest": nearest_depth["configuration_bounds_mm"], "pass": True}, indent=2))


if __name__ == "__main__":
    main()
