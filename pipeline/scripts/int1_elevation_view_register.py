#!/usr/bin/env python3
"""Compile the official D1 elevation-index views into a controlled CSV register."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

from e304_router_cad_evidence import entity_by_handle, first, parse_dxf, point, sha256, text_value


EXPECTED_PLAN_DWG_SHA256 = "ba355f6a90732ad07f843d59e8bab5e1da9daffe7aed74889d1d265b2ce22d7e"
EXPECTED_PLAN_DXF_SHA256 = "706483de83d90e526d7d7cf4f0b097902a4f97868c16aba5f6fac9d9e942f0f0"
EXPECTED_ELEVATION_DWG_SHA256 = "80c0e01a3e713941f729c3b01750181c18ea4d0ff4a5ae8d65800686f19646dc"

SOURCE_VIEWPORT_HANDLE = "243CB4"
SOURCE_SCALE = "1:50@A2"
INDEX_TITLE = "KP-01 立面索引图"
UNINDEXED_SCOPE = "I-503/玄关：01 DWG立面索引未发现官方视图；不得伪称来自DWG"

DIRECTION_BY_BLOCK = {
    "*U517": "+Y",
    "*U518": "+X",
    "*U516": "-Y",
    "*U519": "-X",
}

SHEET_TITLES = {
    "EL-01": "客餐厅立面图",
    "EL-02": "客餐厅立面图",
    "EL-03": "厨房立面图",
    "EL-04": "次卧一立面图",
    "EL-05": "主卧立面图",
    "EL-06": "主卫立面图",
    "EL-07": "次卧二立面图",
    "EL-08": "公卫立面图",
    "EL-09": "次卧三立面图",
}

# The source point is the resolved leader/arrow target in source model space, not
# the elevation marker circle insertion point. Connector inserts own the long
# leader geometry for A2-A5, A7-A8, A10 and A12; singleton markers use the arrow
# target inside their anonymous block for A1, A6, A9 and A11.
ANCHORS = {
    "A1":  {"handle": "243CB7", "block": "*U517", "source": (126967.662, -49873.101), "ifc": (-2285.579, 127.511), "space": "R14 次卧", "gid": "0WyQ2Z9pX5qgwfdTOZwAkw"},
    "A2":  {"handle": "243C70", "block": "*U520", "source": (131074.203, -49373.632), "ifc": (1820.962, 626.980), "space": "R20 客厅", "gid": "2wgBPVUpv2DvcZCfbe6fdv"},
    "A3":  {"handle": "243CC6", "block": "*U520", "source": (131074.203, -54241.016), "ifc": (1820.962, -4240.404), "space": "R04 中厨", "gid": "2fhEbDfK1EkhJwlPikNm$b"},
    "A4":  {"handle": "243D0A", "block": "*U521", "source": (127587.878, -53173.679), "ifc": (-1665.364, -3173.067), "space": "R07 餐厅", "gid": "2JnxeB$or6D8nffRg7GQbf"},
    "A5":  {"handle": "243D4E", "block": "*U522", "source": (124003.449, -51782.222), "ifc": (-5249.793, -1781.610), "space": "R09 主卧", "gid": "3gHz6U6BfFXgV6PnRzfOf$"},
    "A6":  {"handle": "243DD6", "block": "*U517", "source": (123564.439, -46965.612), "ifc": (-5688.803, 3035.000), "space": "R13 主卫湿区飘窗", "gid": "2QIFZrZIr8GQwOQOo0KVFP"},
    "A7":  {"handle": "243D92", "block": "*U523", "source": (123041.344, -48478.855), "ifc": (-6211.897, 1521.757), "space": "R12 主卫湿区", "gid": "3a4COIs5X7lgDirMBDT4Vs"},
    "A8":  {"handle": "243DE5", "block": "*U524", "source": (125885.604, -46450.394), "ifc": (-3367.638, 3550.218), "space": "R15 次卧飘窗", "gid": "004pVXwHv1W8qavLfulRuQ"},
    "A9":  {"handle": "243E6D", "block": "*U517", "source": (128011.070, -46492.886), "ifc": (-1242.172, 3507.726), "space": "R18 客卫飘窗", "gid": "3I284j2FzFG9RguqlzYOH2"},
    "A10": {"handle": "243E29", "block": "*U525", "source": (127923.918, -47310.877), "ifc": (-1329.323, 2689.735), "space": "R17 客卫", "gid": "3bF9ub5u5FFgmXQtSgfWnW"},
    "A11": {"handle": "243E7C", "block": "*U516", "source": (127795.830, -48990.124), "ifc": (-1457.412, 1010.488), "space": "R16 客卫干区", "gid": "25Lcecvwn9lg6g3QKgRlJU"},
    "A12": {"handle": "25F70A", "block": "*U643", "source": (134277.118, -48277.902), "ifc": (5023.876, 1722.710), "space": "R22 书房", "gid": "0XmeOOraz9tP2_CetH7YYh"},
}

VIEW_SPECS = [
    ("01", "EL-01", "A1",  "+Y", "243CB7"),
    ("02", "EL-01", "A2",  "+Y", "243C78"),
    ("03", "EL-01", "A2",  "+X", "243C87"),
    ("04", "EL-02", "A2",  "-Y", "243C96"),
    ("05", "EL-02", "A2",  "-X", "243CA5"),
    ("06", "EL-03", "A3",  "+Y", "243CCE"),
    ("07", "EL-03", "A3",  "+X", "243CDD"),
    ("08", "EL-03", "A3",  "-Y", "243CEC"),
    ("09", "EL-03", "A3",  "-X", "243CFB"),
    ("10", "EL-04", "A4",  "+Y", "243D12"),
    ("11", "EL-04", "A4",  "+X", "243D21"),
    ("12", "EL-04", "A4",  "-Y", "243D30"),
    ("13", "EL-04", "A4",  "-X", "243D3F"),
    ("14", "EL-05", "A5",  "+Y", "243D56"),
    ("15", "EL-05", "A5",  "+X", "243D65"),
    ("16", "EL-05", "A5",  "-Y", "243D74"),
    ("17", "EL-05", "A5",  "-X", "243D83"),
    ("18", "EL-06", "A6",  "+Y", "243DD6"),
    ("19", "EL-06", "A7",  "+Y", "243D9A"),
    ("20", "EL-06", "A7",  "-Y", "243DB8"),
    ("21", "EL-06", "A7",  "+X", "243DA9"),
    ("22", "EL-06", "A7",  "-X", "243DC7"),
    ("23", "EL-07", "A8",  "+Y", "243DED"),
    ("24", "EL-07", "A8",  "+X", "243DFC"),
    ("25", "EL-07", "A8",  "-Y", "243E0B"),
    ("26", "EL-07", "A8",  "-X", "243E1A"),
    ("27", "EL-08", "A9",  "+Y", "243E6D"),
    ("28", "EL-08", "A10", "+Y", "243E31"),
    ("29", "EL-08", "A10", "+X", "243E40"),
    ("30", "EL-08", "A11", "-Y", "243E7C"),
    ("31", "EL-08", "A10", "-Y", "243E4F"),
    ("32", "EL-08", "A10", "-X", "243E5E"),
    ("33", "EL-09", "A12", "+Y", "25F712"),
    ("34", "EL-09", "A12", "+X", "25F721"),
    ("35", "EL-09", "A12", "-Y", "25F730"),
    ("36", "EL-09", "A12", "-X", "25F73F"),
]

REQUIRED_COLUMNS = [
    "view_id", "sheet_id", "official_title", "anchor_id", "direction",
    "source_handle", "source_viewport_handle", "source_scale",
    "source_anchor_handle", "source_x_mm", "source_y_mm", "ifc_x_mm",
    "ifc_y_mm", "space_reference", "space_global_id", "review_status",
    "source_locator",
]
EXTRA_COLUMNS = ["direction_block", "duplicate_source_handles", "anchor_method", "scope_exclusion"]


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    source_dir = root.parent / "图纸" / "矩阵纵横"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dwg", type=Path, default=source_dir / "01 成品房(D1户型-115)平面系统图.dwg")
    parser.add_argument("--plan-dxf", type=Path, default=root / "tmp/dwg/d1-handover-plan.dxf")
    parser.add_argument("--elevation-dwg", type=Path, default=source_dir / "02 成品房(D1户型-115)立面图.dwg")
    parser.add_argument("--output", type=Path, default=root / "pipeline/decisions/int1-elevation-view-register.csv")
    return parser.parse_args()


def verify_hash(path: Path, expected: str, label: str) -> str:
    if not path.is_file():
        raise RuntimeError(f"{label} is missing: {path}")
    actual = sha256(path)
    if actual != expected:
        raise RuntimeError(f"{label} SHA-256 mismatch: expected {expected}, got {actual}")
    return actual


def owned_attributes(entities: list[dict[str, Any]]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for entity in entities:
        if entity["type"] != "ATTRIB":
            continue
        owner = first(entity, "330")
        if owner:
            result.setdefault(owner, []).append(text_value(entity))
    return result


def extract_index_candidates(entities: list[dict[str, Any]]) -> list[dict[str, Any]]:
    attributes = owned_attributes(entities)
    candidates = []
    for entity in entities:
        block = first(entity, "2")
        if entity["type"] != "INSERT" or block not in DIRECTION_BY_BLOCK:
            continue
        handle = first(entity, "5")
        values = attributes.get(handle, [])
        if len(values) != 2 or not values[0].isdigit() or values[1] not in SHEET_TITLES:
            continue
        candidates.append({
            "view_id": values[0].zfill(2),
            "sheet_id": values[1],
            "direction": DIRECTION_BY_BLOCK[block],
            "direction_block": block,
            "source_handle": handle,
            "paper_point": point(entity),
        })
    return candidates


def deduplicate(candidates: list[dict[str, Any]]) -> tuple[dict[str, dict[str, Any]], dict[str, list[str]]]:
    unique: dict[str, dict[str, Any]] = {}
    duplicates: dict[str, list[str]] = {}
    for candidate in candidates:
        view_id = candidate["view_id"]
        if view_id not in unique:
            unique[view_id] = candidate
            continue
        current = unique[view_id]
        same = all(candidate[key] == current[key] for key in ("sheet_id", "direction", "direction_block", "paper_point"))
        if not same:
            raise RuntimeError(f"view {view_id} has conflicting DXF index instances")
        duplicates.setdefault(view_id, []).append(candidate["source_handle"])
    return unique, duplicates


def review_status(anchor_id: str) -> str:
    if anchor_id in {"A1", "A4", "A12"}:
        return "review_required_current_space_differs_from_official_title"
    return "mapped_to_current_space"


def compile_rows(entities: list[dict[str, Any]]) -> tuple[list[dict[str, str]], dict[str, list[str]]]:
    viewport = entity_by_handle(entities, SOURCE_VIEWPORT_HANDLE)
    if viewport["type"] != "VIEWPORT" or first(viewport, "67", "0") != "1":
        raise RuntimeError(f"source viewport #{SOURCE_VIEWPORT_HANDLE} is not the paper-space viewport")
    scale = float(first(viewport, "45")) / float(first(viewport, "41"))
    if abs(scale - 50.0) > 1e-9:
        raise RuntimeError(f"source viewport scale changed: expected 1:50, got 1:{scale}")

    candidates, duplicates = deduplicate(extract_index_candidates(entities))
    expected_ids = {view_id for view_id, *_ in VIEW_SPECS}
    if set(candidates) != expected_ids:
        raise RuntimeError(f"DXF elevation index changed: expected {sorted(expected_ids)}, got {sorted(candidates)}")

    rows = []
    for view_id, sheet_id, anchor_id, direction, expected_handle in VIEW_SPECS:
        candidate = candidates[view_id]
        expected_block = next(block for block, value in DIRECTION_BY_BLOCK.items() if value == direction)
        if candidate["sheet_id"] != sheet_id or candidate["direction"] != direction:
            raise RuntimeError(f"view {view_id} sheet/direction changed in DXF")
        if candidate["source_handle"] != expected_handle:
            raise RuntimeError(f"view {view_id} canonical source handle changed")
        if candidate["direction_block"] != expected_block:
            raise RuntimeError(f"view {view_id} direction block changed")

        anchor = ANCHORS[anchor_id]
        anchor_entity = entity_by_handle(entities, anchor["handle"])
        if anchor_entity["type"] != "INSERT" or first(anchor_entity, "2") != anchor["block"]:
            raise RuntimeError(f"anchor {anchor_id} source entity changed")
        source_x, source_y = anchor["source"]
        ifc_x, ifc_y = anchor["ifc"]
        duplicate_handles = duplicates.get(view_id, [])
        locator = (
            f"01 DWG > {INDEX_TITLE} > VIEWPORT #{SOURCE_VIEWPORT_HANDLE} > "
            f"INSERT #{expected_handle} > resolved leader target {anchor_id}; "
            f"02 DWG > {sheet_id} {SHEET_TITLES[sheet_id]}"
        )
        rows.append({
            "view_id": view_id,
            "sheet_id": sheet_id,
            "official_title": SHEET_TITLES[sheet_id],
            "anchor_id": anchor_id,
            "direction": direction,
            "source_handle": expected_handle,
            "source_viewport_handle": SOURCE_VIEWPORT_HANDLE,
            "source_scale": SOURCE_SCALE,
            "source_anchor_handle": anchor["handle"],
            "source_x_mm": f"{source_x:.3f}",
            "source_y_mm": f"{source_y:.3f}",
            "ifc_x_mm": f"{ifc_x:.3f}",
            "ifc_y_mm": f"{ifc_y:.3f}",
            "space_reference": anchor["space"],
            "space_global_id": anchor["gid"],
            "review_status": review_status(anchor_id),
            "source_locator": locator,
            "direction_block": candidate["direction_block"],
            "duplicate_source_handles": ";".join(duplicate_handles),
            "anchor_method": "resolved_leader_or_arrow_target_not_marker_circle",
            "scope_exclusion": UNINDEXED_SCOPE,
        })
    return rows, duplicates


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS + EXTRA_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    try:
        hashes = {
            "plan_dwg": verify_hash(args.plan_dwg, EXPECTED_PLAN_DWG_SHA256, "01 plan DWG"),
            "plan_dxf": verify_hash(args.plan_dxf, EXPECTED_PLAN_DXF_SHA256, "01 verification DXF"),
            "elevation_dwg": verify_hash(args.elevation_dwg, EXPECTED_ELEVATION_DWG_SHA256, "02 elevation DWG"),
        }
        entities = parse_dxf(args.plan_dxf)
        rows, duplicates = compile_rows(entities)
        write_csv(args.output, rows)
        report = {
            "schema_version": "int1-elevation-view-register/v1",
            "output": str(args.output.resolve()),
            "hashes": hashes,
            "view_count": len(rows),
            "sheet_count": len({row["sheet_id"] for row in rows}),
            "directions": sorted({row["direction"] for row in rows}),
            "duplicates_removed": duplicates,
            "anchor_policy": "resolved leader/arrow target; marker circle insertion is not an anchor",
            "unindexed_scopes": [UNINDEXED_SCOPE],
            "formal_ifc_write": False,
            "blender_required": False,
        }
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 0
    except (OSError, RuntimeError, ValueError, KeyError) as error:
        print(f"int1 elevation view register failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
