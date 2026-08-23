#!/usr/bin/env python3
"""Extract and reject the antoniolupi Street parent-top DXF for STREET-H."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, relative, sha256, write_json


PRODUCT_DIR = ROOT / "output/review/highpoly-types/street-h"
SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_ZIP = SOURCE_DIR / "ANTONIOLUPI-official-Street-2D-CAD.zip"
SOURCE_DXF = SOURCE_DIR / "AL_Street.dxf"
SOURCE_TECHNICAL_PDF = SOURCE_DIR / "ANTONIOLUPI-official-Street-technical.pdf"
DEFAULT_OUTPUT = PRODUCT_DIR / "official-native-dxf-linework.json"
PRODUCT_PAGE = "https://www.antoniolupi.it/en/products/sinks/street"
EXPECTED_SHA256 = {
    "zip": "36161ac81bade86b7ec0419c0d9cf0c49d52ae0ecb9dced39006b9eddcb10e62",
    "dxf": "72b0a893545c4d9ff2a1d3fc79b7618eb10c00c651533eed36071e4130277c7e",
    "technical_pdf": "0b883bbf6df617331cc2408da4860b31160af67a6b6924b3caaee0ec389d2946",
}
EXACT_LABEL_HANDLE = "1388"
EXACT_LABEL = "street240 prof. 40 + street4054 prof. 40"
VIEW_INSERTS = {
    "plan": {"handle": "1376", "block": "*U4", "expected_paths": 7},
    "front": {"handle": "1369", "block": "*U0", "expected_paths": 1},
}
DRAWING_UNITS_TO_MM = 10.0
EXPECTED_COMPLETE_TOP_BOUNDS_MM = {"plan": [1080.0, 400.0], "front": [1080.0, 250.0]}
PROJECT_COMPONENT_BOUNDS_MM = [300.0, 150.0, 100.0]
SCOPE = (
    "official complete antoniolupi Street family-top CAD evidence only; the project STREET-H type is an "
    "isolated repeated sink-holder subcomponent, not the complete catalogue top and not a project shop drawing"
)


def records(path: Path) -> tuple[list[dict], dict[str, list[dict]]]:
    lines = path.read_text(encoding="latin1").splitlines()
    if len(lines) % 2:
        raise RuntimeError("Street DXF has an odd group-code line count")
    pairs = [(lines[index].strip(), lines[index + 1].strip()) for index in range(0, len(lines), 2)]
    section = None
    current = None
    entities: list[dict] = []
    blocks: dict[str, list[dict]] = {}
    current_block = None
    for code, value in pairs:
        if code == "0" and value == "SECTION":
            section = "PENDING"
            continue
        if section == "PENDING" and code == "2":
            section = value
            continue
        if code == "0" and value == "ENDSEC":
            if current is not None:
                if section == "BLOCKS" and current_block and current["type"] not in {"BLOCK", "ENDBLK"}:
                    blocks.setdefault(current_block, []).append(current)
                elif section == "ENTITIES":
                    entities.append(current)
            current = None
            current_block = None
            section = None
            continue
        if section not in {"BLOCKS", "ENTITIES"}:
            continue
        if code == "0":
            if current is not None:
                if section == "BLOCKS" and current_block and current["type"] not in {"BLOCK", "ENDBLK"}:
                    blocks.setdefault(current_block, []).append(current)
                elif section == "ENTITIES":
                    entities.append(current)
            current = {"type": value, "pairs": []}
            if section == "BLOCKS" and value == "ENDBLK":
                current_block = None
            continue
        if current is None:
            continue
        current["pairs"].append((code, value))
        if section == "BLOCKS" and current["type"] == "BLOCK" and code == "2":
            current_block = value
            blocks.setdefault(current_block, [])
    return entities, blocks


def first(entity: dict, code: str, default=None):
    return next((value for key, value in entity["pairs"] if key == code), default)


def vertices(entity: dict) -> list[list[float]]:
    result: list[list[float]] = []
    for code, value in entity["pairs"]:
        if code == "10":
            result.append([float(value), 0.0])
        elif code == "20" and result:
            result[-1][1] = float(value)
    return result


def circle_path(entity: dict, segments: int = 64) -> list[list[float]]:
    cx, cy, radius = float(first(entity, "10")), float(first(entity, "20")), float(first(entity, "40"))
    return [
        [cx + radius * math.cos(2.0 * math.pi * index / segments), cy + radius * math.sin(2.0 * math.pi * index / segments)]
        for index in range(segments + 1)
    ]


def entity_path(entity: dict) -> list[list[float]] | None:
    kind = entity["type"]
    if kind == "LINE":
        return [
            [float(first(entity, "10")), float(first(entity, "20"))],
            [float(first(entity, "11")), float(first(entity, "21"))],
        ]
    if kind == "LWPOLYLINE":
        points = vertices(entity)
        if any(code == "42" and abs(float(value)) > 1e-12 for code, value in entity["pairs"]):
            raise RuntimeError("selected Street DXF cluster contains an unsupported bulged polyline")
        if int(first(entity, "70", "0")) & 1 and points and points[0] != points[-1]:
            points.append(points[0])
        return points
    if kind == "CIRCLE":
        return circle_path(entity)
    return None


def bounds(paths: list[list[list[float]]]) -> dict:
    points = [point for path in paths for point in path]
    minimum = [min(point[axis] for point in points) for axis in range(2)]
    maximum = [max(point[axis] for point in points) for axis in range(2)]
    return {
        "minimum": [round(value, 6) for value in minimum],
        "maximum": [round(value, 6) for value in maximum],
        "size": [round(maximum[axis] - minimum[axis], 6) for axis in range(2)],
    }


def insert_record(entities: list[dict], handle: str, block: str) -> dict:
    matches = [item for item in entities if item["type"] == "INSERT" and first(item, "5") == handle]
    if len(matches) != 1 or first(matches[0], "2") != block:
        raise RuntimeError(f"Street DXF insert {handle}/{block} identity drifted")
    return matches[0]


def extract_view(entities: list[dict], blocks: dict[str, list[dict]], view: str, selector: dict) -> dict:
    insertion = insert_record(entities, selector["handle"], selector["block"])
    block_entities = [item for item in blocks[selector["block"]] if first(item, "8", "0") == "0"]
    raw_paths = [path for item in block_entities if (path := entity_path(item))]
    native = bounds(raw_paths)
    normalized = [
        [
            [
                round((point[0] - native["minimum"][0]) * DRAWING_UNITS_TO_MM, 6),
                round((point[1] - native["minimum"][1]) * DRAWING_UNITS_TO_MM, 6),
            ]
            for point in path
        ]
        for path in raw_paths
    ]
    metric = bounds(normalized)
    if len(normalized) != selector["expected_paths"] or metric["size"] != EXPECTED_COMPLETE_TOP_BOUNDS_MM[view]:
        raise RuntimeError(f"Street official {view} cluster geometry drifted: {len(normalized)} / {metric['size']}")
    return {
        "view": view,
        "source_dxf": relative(SOURCE_DXF),
        "source_dxf_sha256": sha256(SOURCE_DXF),
        "source_kind": "native_dxf",
        "native_insert_handle": selector["handle"],
        "native_block_name": selector["block"],
        "native_insert_point": [float(first(insertion, "10")), float(first(insertion, "20"))],
        "native_entity_counts": dict(sorted(Counter(item["type"] for item in block_entities if entity_path(item)).items())),
        "native_drawing_units_to_mm": DRAWING_UNITS_TO_MM,
        "native_bounds_drawing_units": native,
        "bounds_mm": metric,
        "path_count": len(normalized),
        "point_count": sum(len(path) for path in normalized),
        "paths_mm": normalized,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE_DXF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    source = args.source.resolve()
    if (
        sha256(SOURCE_ZIP) != EXPECTED_SHA256["zip"]
        or sha256(source) != EXPECTED_SHA256["dxf"]
        or sha256(SOURCE_TECHNICAL_PDF) != EXPECTED_SHA256["technical_pdf"]
    ):
        raise RuntimeError("antoniolupi Street official archive hash mismatch")
    entities, blocks = records(source)
    labels = [
        item for item in entities
        if item["type"] == "MTEXT" and first(item, "5") == EXACT_LABEL_HANDLE
    ]
    if len(labels) != 1 or first(labels[0], "1", "").lower() != EXACT_LABEL:
        raise RuntimeError("exact Street240 + Street4054 DXF label drifted")
    views = {view: extract_view(entities, blocks, view, selector) for view, selector in VIEW_INSERTS.items()}
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/street_h_linework.py",
        "manufacturer": "antoniolupi",
        "family": "Street",
        "project_ifc_type_name": "STREET-H",
        "project_ifc_type_description": "Sink holder",
        "product_page": PRODUCT_PAGE,
        "source_kind": "native_dxf",
        "source_zip": relative(SOURCE_ZIP),
        "source_zip_sha256": sha256(SOURCE_ZIP),
        "source_dxf": relative(source),
        "source_dxf_sha256": sha256(source),
        "source_technical_pdf": relative(SOURCE_TECHNICAL_PDF),
        "source_technical_pdf_sha256": sha256(SOURCE_TECHNICAL_PDF),
        "dxf_version": "AC1015",
        "dxf_header_insunits": 4,
        "unit_cross_check": {
            "native_cluster_plan_bounds_drawing_units": [108.0, 40.0],
            "native_cluster_front_bounds_drawing_units": [108.0, 25.0],
            "technical_pdf_page": 1,
            "technical_pdf_complete_top_dimensions_cm": {"minimum_length": 108.0, "depth": 40.0, "height": 25.0},
            "resolved_native_drawing_units_to_mm": DRAWING_UNITS_TO_MM,
            "resolved_complete_top_bounds_mm": EXPECTED_COMPLETE_TOP_BOUNDS_MM,
            "pass": True,
        },
        "selected_parent_family_cluster": {
            "label_handle": EXACT_LABEL_HANDLE,
            "label": EXACT_LABEL,
            "material_variant": "marble",
            "views_present": ["plan", "front"],
            "side_view_present": False,
        },
        "project_component_cross_check": {
            "project_STREET_H_body_bounds_mm": PROJECT_COMPONENT_BOUNDS_MM,
            "project_instance_count": 2,
            "official_cluster_scope": "complete configurable STREET240 + STREET4054 top",
            "exact_STREET_H_component_match": False,
            "official_parent_paths_used_as_STREET_H_representation": False,
            "reason": "The native DXF geometry is the complete parent sink top and is at least 1080 x 400 x 250 mm. The repeated project STREET-H Body is a 300 x 150 x 100 mm sink-holder subcomponent.",
            "pass": True,
        },
        "scope": SCOPE,
        "views": views,
        "pass": True,
    }
    write_json(args.output.resolve(), payload)
    print(json.dumps({"output": relative(args.output), "views": {view: item["bounds_mm"] for view, item in views.items()}, "pass": True}, indent=2))


if __name__ == "__main__":
    main()
