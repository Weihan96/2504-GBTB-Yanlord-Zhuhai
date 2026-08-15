#!/usr/bin/env python3
"""Register the owner-confirmed APP-019 ice-cream-machine placeholder."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "pipeline/decisions"
SOURCE_ID = "OWNER-INPUT-APP019-20260815"
EVIDENCE_RELATIVE = "drawings/evidence/OWNER-INPUT-APP-019-20260815.md"
EVIDENCE_PATH = ROOT / EVIDENCE_RELATIVE
EQUIPMENT_ID = "APP-019"


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    source_fields, sources = read_csv(DECISIONS / "source-evidence-register.csv")
    source_row = {field: "" for field in source_fields}
    source_row.update({
        "source_id": SOURCE_ID,
        "discipline": "ELEC/PLUM/INT1",
        "sheet_id": "E-303/P-201/P-202/INT1/S-701",
        "decision_scope": "APP-019 冰淇淋机设备类别与数量",
        "source_kind": "owner_confirmation",
        "source_document": EVIDENCE_RELATIVE,
        "local_path": EVIDENCE_RELATIVE,
        "sha256": sha256(EVIDENCE_PATH),
        "locator": "业主新增一台冰淇淋机",
        "evidence": "业主确认新增冰淇淋机 1 台",
        "proves": "项目需要登记一台冰淇淋机",
        "does_not_prove": "已采购、品牌型号、位置、安装形式、功率、插头插座、回路、给排水、尺寸或接口中心",
        "status": "confirmed_owner_input_identity_only",
        "confidence": "1.00",
        "review_required": "yes",
        "formal_ifc_write_allowed": "no",
        "revision": "2026-08-15",
        "publication_date": "2026-08-15",
        "notes": "只登记设备类别与数量；所有安装参数保持 unknown。",
    })
    source_by_id = {row["source_id"]: row for row in sources}
    source_by_id[SOURCE_ID] = source_row
    write_csv(DECISIONS / "source-evidence-register.csv", source_fields, list(source_by_id.values()))

    master_fields, masters = read_csv(DECISIONS / "equipment-register.csv")
    master_row = {field: "" for field in master_fields}
    master_row.update({
        "equipment_id": EQUIPMENT_ID,
        "domain": "APPLIANCE",
        "category": "甜品设备",
        "item_name": "冰淇淋机",
        "quantity": "1",
        "procurement_status": "candidate",
        "decision_status": "partial",
        "storage_location_candidate": "待确认位置",
        "use_location_candidate": "待确认位置",
        "schedule_included": "yes",
        "selector_kind": "logical_input",
        "selector_value": EQUIPMENT_ID,
        "source_ids": SOURCE_ID,
        "identity_basis": "业主确认新增一台冰淇淋机；未提供品牌、型号、位置或安装资料",
        "confidence": "1.00",
        "human_review_required": "yes",
        "legacy_kind": "appliance",
        "legacy_id": EQUIPMENT_ID,
        "notes": "候选设备，非已选/已购；不计入 NS-01/NS-02 负荷，待位置和准确型号后再决定插座、回路与给排水。",
    })
    master_by_id = {row["equipment_id"]: row for row in masters}
    master_by_id[EQUIPMENT_ID] = master_row
    write_csv(DECISIONS / "equipment-register.csv", master_fields, list(master_by_id.values()))

    requirement_fields, requirements = read_csv(DECISIONS / "equipment-installation-requirements.csv")
    requirements = [row for row in requirements if row["equipment_id"] != EQUIPMENT_ID]
    definitions = [
        ("rated_power", "ELEC", "", "W", "准确型号名牌功率"),
        ("simultaneous_group", "ELEC", "", "", "使用位置和同时工况"),
        ("water_required", "PLUM", "", "", "准确型号安装图"),
        ("drain_required", "PLUM", "", "", "准确型号安装图"),
        ("gas_required", "GAS", "", "", "准确型号安装图"),
        ("ventilation_required", "HVAC/INT1", "", "", "准确型号散热和净距要求"),
        ("installation_mode", "INT1", "unknown", "", "台面、落地或嵌入形式"),
        ("appliance_plug_rating_a", "ELEC", "unknown", "A", "实物插头或厂家电气页"),
        ("wall_socket_rating_a", "ELEC", "unknown", "A", "按实购插头和支路计算匹配"),
        ("branch_breaker_rating_a", "ELEC", "unknown", "A", "名牌负荷、同时工况和回路计算"),
        ("dedicated_branch_circuit", "ELEC", "unknown", "", "名牌负荷、位置和共路工况"),
        ("product_width", "INT1", "unknown", "mm", "准确型号尺寸图"),
        ("product_height", "INT1", "unknown", "mm", "准确型号尺寸图"),
        ("product_depth", "INT1", "unknown", "mm", "准确型号尺寸图"),
        ("service_clearance", "INT1", "unknown", "mm", "厂家散热、操作和检修要求"),
        ("interface_center_coordinates", "MULTI", "unknown", "", "厂家安装图与项目定位"),
    ]
    for index, (key, discipline, value, unit, required_evidence) in enumerate(definitions, start=1):
        row = {field: "" for field in requirement_fields}
        row.update({
            "requirement_id": f"REQ-APP019-{index:03d}",
            "equipment_id": EQUIPMENT_ID,
            "discipline": discipline,
            "parameter_key": key,
            "value_text": value,
            "unit": unit,
            "value_origin": "pending",
            "status": "pending",
            "source_id": SOURCE_ID,
            "source_locator": "业主仅确认设备类别与数量",
            "blocks_release": "yes",
            "notes": f"保持 unknown；关闭证据：{required_evidence}。",
        })
        requirements.append(row)
    write_csv(DECISIONS / "equipment-installation-requirements.csv", requirement_fields, requirements)


if __name__ == "__main__":
    main()
