#!/usr/bin/env python3
"""Validate and synchronize the owner input workbook and its SSOT projections.

The default mode is a read-only dry run. ``--apply`` updates the owner input
registers and PM status block, then rebuilds the workbook's three read-only
views from canonical CSVs. Read-only workbook views never write back to those
canonical CSVs. The script never writes IFC.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import posixpath
import re
import sys
import tempfile
import zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Optional
from xml.etree import ElementTree as ET

from equipment_ssot import apply_owner_appliances, projections


DECISION_HEADERS = [
    "input_id", "workstream", "priority", "blocks_release", "question",
    "candidate_value", "user_value", "unit", "status", "evidence_reference",
    "source_basis", "sync_target", "notes",
]
APPLIANCE_HEADERS = [
    "appliance_id", "appliance_name", "category", "storage_location_candidate",
    "use_location_candidate", "storage_location_confirmed", "use_location_confirmed",
    "quantity", "rated_power_w", "simultaneous_group", "water_required",
    "drain_required", "gas_required", "ventilation_required", "model",
    "evidence_reference", "status", "notes",
]
CLOSEOUT_HEADERS = [
    "input_id", "closeout_kind", "responsible_party", "required_evidence",
    "automatic_close_allowed", "notes",
]
DECISION_DISPLAY_HEADERS = [
    "ID", "专业", "优先级", "阻塞无保留发布", "需要你确认", "常见候选（不等于确认）",
    "你的确认值", "单位", "状态", "证据/链接", "现有依据", "同步目标", "备注",
]
APPLIANCE_DISPLAY_HEADERS = [
    "设备ID", "设备名称", "类别", "候选存放位置", "候选使用位置", "确认存放位置",
    "确认使用位置", "数量", "额定功率(W)", "同时使用组", "需给水", "需排水",
    "需燃气", "需通风", "型号", "证据/链接", "状态", "备注",
]
DECISION_STATUSES = {"待填写", "采用候选", "自定义确认", "需证据", "暂缓", "不适用"}
APPLIANCE_STATUSES = {"待填写", "部分确认", "已确认", "暂缓", "不适用"}
PM_START = "<!-- OWNER_INPUTS:START -->"
PM_END = "<!-- OWNER_INPUTS:END -->"
DECISION_USER_FIELDS = {"user_value", "status", "evidence_reference", "notes"}
APPLIANCE_USER_FIELDS = {
    "storage_location_confirmed", "use_location_confirmed", "quantity",
    "rated_power_w", "simultaneous_group", "water_required", "drain_required",
    "gas_required", "ventilation_required", "model", "evidence_reference",
    "status", "notes",
}
EXPECTED_WORKBOOK_SHEETS = [
    "使用说明", "设计决策", "家电清单", "设备主表", "安装条件", "证据索引",
]
READONLY_VIEW_SPECS = {
    "设备主表": {
        "argument": "equipment_register",
        "primary_key": "equipment_id",
    },
    "安装条件": {
        "argument": "installation_requirements",
        "primary_key": "requirement_id",
    },
    "证据索引": {
        "argument": "evidence_register",
        "primary_key": "source_id",
    },
}
MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
PACKAGE_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "output/forms/滨海湾施工输入清单.xlsx")
    parser.add_argument("--decisions", type=Path, default=root / "pipeline/decisions/owner-input-register.csv")
    parser.add_argument("--appliances", type=Path, default=root / "pipeline/decisions/appliance-input-register.csv")
    parser.add_argument("--closeout-rules", type=Path, default=root / "pipeline/decisions/owner-input-closeout-rules.csv")
    parser.add_argument("--equipment-register", type=Path, default=root / "pipeline/decisions/equipment-register.csv")
    parser.add_argument("--installation-requirements", type=Path, default=root / "pipeline/decisions/equipment-installation-requirements.csv")
    parser.add_argument("--evidence-register", type=Path, default=root / "pipeline/decisions/source-evidence-register.csv")
    parser.add_argument("--pm", type=Path, default=root / "drawings/滨海湾装修施工图深化工作管理.md")
    parser.add_argument("--report", type=Path, help="Optional user-owned JSON preview path")
    parser.add_argument("--open-items", type=Path, help="Optional user-owned Markdown open-items path")
    parser.add_argument("--apply", action="store_true")
    return parser.parse_args()


def column_index(reference: str) -> int:
    letters = "".join(character for character in reference if character.isalpha())
    value = 0
    for character in letters.upper():
        value = value * 26 + ord(character) - ord("A") + 1
    return value - 1


def read_xlsx(
    path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, list[list[str]]], list[str]]:
    with zipfile.ZipFile(path) as archive:
        strings = _shared_strings(archive)
        paths = _sheet_paths(archive)
        missing = set(EXPECTED_WORKBOOK_SHEETS) - set(paths)
        if missing:
            raise ValueError(f"workbook missing sheets: {', '.join(sorted(missing))}")
        decisions = _table(
            _sheet_rows(archive, paths["设计决策"], strings),
            DECISION_HEADERS, DECISION_DISPLAY_HEADERS, "设计决策",
        )
        appliances = _table(
            _sheet_rows(archive, paths["家电清单"], strings),
            APPLIANCE_HEADERS, APPLIANCE_DISPLAY_HEADERS, "家电清单",
        )
        readonly_rows = {
            sheet_name: _sheet_rows(archive, paths[sheet_name], strings)
            for sheet_name in READONLY_VIEW_SPECS
        }
    return decisions, appliances, readonly_rows, list(paths)


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    path = "xl/sharedStrings.xml"
    if path not in archive.namelist():
        return []
    root = ET.fromstring(archive.read(path))
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    return ["".join(node.text or "" for node in item.findall(".//x:t", ns)) for item in root]


def _sheet_paths(archive: zipfile.ZipFile) -> dict[str, str]:
    main_ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    rel_ns = {"r": "http://schemas.openxmlformats.org/package/2006/relationships"}
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {node.attrib["Id"]: node.attrib["Target"] for node in relationships.findall("r:Relationship", rel_ns)}
    result: dict[str, str] = {}
    rel_key = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
    for sheet in workbook.findall(".//x:sheet", main_ns):
        target = targets[sheet.attrib[rel_key]].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        result[sheet.attrib["name"]] = target
    return result


def _sheet_rows(archive: zipfile.ZipFile, path: str, strings: list[str]) -> list[list[str]]:
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    root = ET.fromstring(archive.read(path))
    rows: list[list[str]] = []
    for row in root.findall(".//x:sheetData/x:row", ns):
        values: dict[int, str] = {}
        for cell in row.findall("x:c", ns):
            index = column_index(cell.attrib["r"])
            cell_type = cell.attrib.get("t", "")
            if cell_type == "inlineStr":
                value = "".join(node.text or "" for node in cell.findall(".//x:t", ns))
            else:
                node = cell.find("x:v", ns)
                value = node.text if node is not None and node.text is not None else ""
                if cell_type == "s" and value:
                    value = strings[int(value)]
            values[index] = value
        if values:
            rows.append([values.get(index, "") for index in range(max(values) + 1)])
    return rows


def _table(rows: list[list[str]], headers: list[str], display_headers: list[str], sheet_name: str) -> list[dict[str, str]]:
    if not rows or [value.strip() for value in rows[0]][: len(display_headers)] != display_headers:
        raise ValueError(f"{sheet_name} headers changed; expected {display_headers}")
    records: list[dict[str, str]] = []
    for values in rows[1:]:
        padded = values + [""] * (len(headers) - len(values))
        record = {header: str(padded[index]).strip() for index, header in enumerate(headers)}
        if any(record.values()):
            records.append(record)
    return records


def read_csv(path: Path, headers: list[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != headers:
            raise ValueError(f"{path} headers changed; expected {headers}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def read_canonical_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        headers = list(reader.fieldnames or [])
        if not headers:
            raise ValueError(f"{path} has no header")
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]
    return headers, rows


def write_csv(path: Path, headers: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def readonly_view_parity(
    sheet_name: str,
    workbook_rows: list[list[str]],
    canonical_headers: list[str],
    canonical_rows: list[dict[str, str]],
    primary_key: str,
) -> dict[str, Any]:
    actual_headers = [value.strip() for value in workbook_rows[1]] if len(workbook_rows) > 1 else []
    columns_match = actual_headers == canonical_headers
    actual_records: list[dict[str, str]] = []
    if actual_headers:
        for values in workbook_rows[2:]:
            padded = values + [""] * (len(actual_headers) - len(values))
            record = {
                header: str(padded[index]).strip()
                for index, header in enumerate(actual_headers)
            }
            if any(record.values()):
                actual_records.append(record)

    expected_keys = [row.get(primary_key, "") for row in canonical_rows]
    actual_keys = [row.get(primary_key, "") for row in actual_records]
    duplicate_keys = sorted(
        key for key, count in Counter(actual_keys).items() if key and count > 1
    )
    expected_key_set = set(expected_keys)
    actual_key_set = set(actual_keys)
    missing_keys = sorted(expected_key_set - actual_key_set)
    extra_keys = sorted(actual_key_set - expected_key_set)
    keys_match = (
        not duplicate_keys
        and "" not in actual_key_set
        and actual_key_set == expected_key_set
        and len(actual_keys) == len(expected_keys)
    )

    canonical_by_key = {row[primary_key]: row for row in canonical_rows}
    actual_by_key = {row.get(primary_key, ""): row for row in actual_records}
    value_mismatches: list[dict[str, Any]] = []
    cell_mismatch_count = 0
    row_value_mismatch_count = 0
    if columns_match:
        for key in sorted(expected_key_set & actual_key_set):
            changed_fields = [
                header
                for header in canonical_headers
                if actual_by_key[key].get(header, "") != canonical_by_key[key].get(header, "")
            ]
            if changed_fields:
                row_value_mismatch_count += 1
                cell_mismatch_count += len(changed_fields)
                if len(value_mismatches) < 50:
                    value_mismatches.append({"id": key, "changed_fields": changed_fields})
    values_match = columns_match and not missing_keys and not extra_keys and cell_mismatch_count == 0
    row_count_match = len(actual_records) == len(canonical_rows)
    all_match = columns_match and row_count_match and keys_match and values_match
    return {
        "sheet": sheet_name,
        "primary_key": primary_key,
        "canonical_row_count": len(canonical_rows),
        "workbook_row_count": len(actual_records),
        "row_count_match": row_count_match,
        "canonical_columns": canonical_headers,
        "workbook_columns": actual_headers,
        "columns_match": columns_match,
        "keys_match": keys_match,
        "missing_keys": missing_keys,
        "extra_keys": extra_keys,
        "duplicate_keys": duplicate_keys,
        "blank_key_count": actual_keys.count(""),
        "values_match": values_match,
        "row_value_mismatch_count": row_value_mismatch_count,
        "cell_mismatch_count": cell_mismatch_count,
        "value_mismatches": value_mismatches,
        "all_match": all_match,
    }


def readonly_views_parity(
    readonly_rows: dict[str, list[list[str]]],
    workbook_sheet_names: list[str],
    canonical_tables: dict[str, tuple[list[str], list[dict[str, str]]]],
) -> dict[str, Any]:
    views = {
        sheet_name: readonly_view_parity(
            sheet_name,
            readonly_rows[sheet_name],
            canonical_tables[sheet_name][0],
            canonical_tables[sheet_name][1],
            spec["primary_key"],
        )
        for sheet_name, spec in READONLY_VIEW_SPECS.items()
    }
    workbook_structure_match = workbook_sheet_names == EXPECTED_WORKBOOK_SHEETS
    return {
        "workbook_sheet_names": workbook_sheet_names,
        "expected_workbook_sheet_names": EXPECTED_WORKBOOK_SHEETS,
        "workbook_structure_match": workbook_structure_match,
        "views": views,
        "all_match": workbook_structure_match and all(
            view["all_match"] for view in views.values()
        ),
    }


def _column_reference(index: int) -> str:
    result = ""
    value = index + 1
    while value:
        value, remainder = divmod(value - 1, 26)
        result = chr(ord("A") + remainder) + result
    return result


def _cell(
    row_number: int, column_index_value: int, value: str, style: Optional[str],
) -> ET.Element:
    attributes = {
        "r": f"{_column_reference(column_index_value)}{row_number}",
        "t": "inlineStr",
    }
    if style is not None:
        attributes["s"] = style
    cell = ET.Element(f"{{{MAIN_NS}}}c", attributes)
    inline = ET.SubElement(cell, f"{{{MAIN_NS}}}is")
    text = ET.SubElement(inline, f"{{{MAIN_NS}}}t")
    if value != value.strip() or "\n" in value:
        text.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    text.text = value
    return cell


def _table_path_for_sheet(
    archive: zipfile.ZipFile, sheet_path: str,
) -> str:
    relationship_path = posixpath.join(
        posixpath.dirname(sheet_path), "_rels", posixpath.basename(sheet_path) + ".rels",
    )
    relationship_root = ET.fromstring(archive.read(relationship_path))
    table_relationships = [
        relationship
        for relationship in relationship_root.findall(f"{{{PACKAGE_REL_NS}}}Relationship")
        if relationship.attrib.get("Type", "").endswith("/table")
    ]
    if len(table_relationships) != 1:
        raise ValueError(f"{sheet_path}: expected exactly one table relationship")
    target = table_relationships[0].attrib["Target"]
    if target.startswith("/"):
        return target.lstrip("/")
    return posixpath.normpath(posixpath.join(posixpath.dirname(sheet_path), target))


def _rebuilt_sheet_xml(
    source: bytes, headers: list[str], rows: list[dict[str, str]], sheet_name: str,
) -> bytes:
    root = ET.fromstring(source)
    sheet_data = root.find(f"{{{MAIN_NS}}}sheetData")
    if sheet_data is None:
        raise ValueError(f"{sheet_name}: sheetData missing")
    existing_rows = list(sheet_data)
    if len(existing_rows) < 3:
        raise ValueError(f"{sheet_name}: title/header/data style templates missing")
    title_row, header_template, data_template = existing_rows[:3]
    header_cells = list(header_template)
    data_cells = list(data_template)
    header_styles = [cell.attrib.get("s") for cell in header_cells]
    data_styles = [cell.attrib.get("s") for cell in data_cells]
    header_fallback = header_styles[-1] if header_styles else None
    data_fallback = data_styles[-1] if data_styles else None

    for child in existing_rows:
        sheet_data.remove(child)
    sheet_data.append(title_row)

    header_attributes = dict(header_template.attrib)
    header_attributes["r"] = "2"
    header_row = ET.Element(f"{{{MAIN_NS}}}row", header_attributes)
    for index, header in enumerate(headers):
        style = header_styles[index] if index < len(header_styles) else header_fallback
        header_row.append(_cell(2, index, header, style))
    sheet_data.append(header_row)

    data_attributes = dict(data_template.attrib)
    for offset, record in enumerate(rows, start=3):
        row_attributes = dict(data_attributes)
        row_attributes["r"] = str(offset)
        row = ET.Element(f"{{{MAIN_NS}}}row", row_attributes)
        for index, header in enumerate(headers):
            style = data_styles[index] if index < len(data_styles) else data_fallback
            row.append(_cell(offset, index, record.get(header, ""), style))
        sheet_data.append(row)

    merge_cells = root.find(f"{{{MAIN_NS}}}mergeCells")
    if merge_cells is not None and len(merge_cells):
        list(merge_cells)[0].set("ref", f"A1:{_column_reference(len(headers)-1)}1")
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def _rebuilt_table_xml(source: bytes, headers: list[str], row_count: int) -> bytes:
    root = ET.fromstring(source)
    table_ref = f"A2:{_column_reference(len(headers)-1)}{row_count+2}"
    root.set("ref", table_ref)
    auto_filter = root.find(f"{{{MAIN_NS}}}autoFilter")
    if auto_filter is not None:
        auto_filter.set("ref", table_ref)
    table_columns = root.find(f"{{{MAIN_NS}}}tableColumns")
    if table_columns is None:
        raise ValueError("tableColumns missing")
    for child in list(table_columns):
        table_columns.remove(child)
    table_columns.set("count", str(len(headers)))
    for index, header in enumerate(headers, start=1):
        ET.SubElement(
            table_columns,
            f"{{{MAIN_NS}}}tableColumn",
            {"id": str(index), "name": header},
        )
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def _updated_instruction_summary_xml(
    source: bytes, summary_values: dict[str, int],
) -> bytes:
    root = ET.fromstring(source)
    for reference, value in summary_values.items():
        cell = root.find(f".//{{{MAIN_NS}}}c[@r='{reference}']")
        if cell is None:
            raise ValueError(f"使用说明: summary cell {reference} missing")
        for child in list(cell):
            cell.remove(child)
        cell.set("t", "n")
        value_node = ET.SubElement(cell, f"{{{MAIN_NS}}}v")
        value_node.text = str(value)
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def rebuild_readonly_views(
    workbook_path: Path,
    canonical_tables: dict[str, tuple[list[str], list[dict[str, str]]]],
    summary_values: dict[str, int],
) -> None:
    with zipfile.ZipFile(workbook_path, "r") as source_archive:
        sheet_paths = _sheet_paths(source_archive)
        replacements: dict[str, bytes] = {
            sheet_paths["使用说明"]: _updated_instruction_summary_xml(
                source_archive.read(sheet_paths["使用说明"]), summary_values,
            ),
        }
        for sheet_name in READONLY_VIEW_SPECS:
            sheet_path = sheet_paths[sheet_name]
            headers, rows = canonical_tables[sheet_name]
            replacements[sheet_path] = _rebuilt_sheet_xml(
                source_archive.read(sheet_path), headers, rows, sheet_name,
            )
            table_path = _table_path_for_sheet(source_archive, sheet_path)
            replacements[table_path] = _rebuilt_table_xml(
                source_archive.read(table_path), headers, len(rows),
            )

        temporary = tempfile.NamedTemporaryFile(
            prefix=workbook_path.stem + "-", suffix=".xlsx",
            dir=workbook_path.parent, delete=False,
        )
        temporary_path = Path(temporary.name)
        temporary.close()
        try:
            with zipfile.ZipFile(temporary_path, "w") as target_archive:
                for info in source_archive.infolist():
                    target_archive.writestr(info, replacements.get(info.filename, source_archive.read(info.filename)))
            os.replace(temporary_path, workbook_path)
        finally:
            if temporary_path.exists():
                temporary_path.unlink()


def validate_unique(rows: list[dict[str, str]], key: str, label: str) -> list[str]:
    ids = [row[key] for row in rows]
    errors = [f"{label}: blank {key}"] if any(not value for value in ids) else []
    duplicates = sorted(value for value, count in Counter(ids).items() if count > 1)
    if duplicates:
        errors.append(f"{label}: duplicate IDs: {', '.join(duplicates)}")
    return errors


def validate_protected_fields(
    rows: list[dict[str, str]],
    baseline: list[dict[str, str]],
    key: str,
    user_fields: set[str],
    label: str,
) -> list[str]:
    baseline_by_id = {row[key]: row for row in baseline}
    errors: list[str] = []
    for row in rows:
        previous = baseline_by_id.get(row[key])
        if previous is None:
            continue
        changed = sorted(
            field for field in row
            if field not in user_fields and row[field] != previous.get(field, "")
        )
        if changed:
            errors.append(
                f"{label} {row[key]}: protected fields changed: {', '.join(changed)}"
            )
    return errors


def validate_decisions(rows: list[dict[str, str]], baseline: list[dict[str, str]]) -> list[str]:
    errors = validate_unique(rows, "input_id", "decisions")
    expected_ids = {row["input_id"] for row in baseline}
    actual_ids = {row["input_id"] for row in rows}
    if actual_ids != expected_ids:
        errors.append(f"decisions: IDs changed; missing={sorted(expected_ids-actual_ids)} extra={sorted(actual_ids-expected_ids)}")
    errors.extend(validate_protected_fields(
        rows, baseline, "input_id", DECISION_USER_FIELDS, "decision",
    ))
    for row in rows:
        status = row["status"]
        if status not in DECISION_STATUSES:
            errors.append(f"{row['input_id']}: invalid status {status!r}")
        if status == "采用候选" and not row["candidate_value"]:
            errors.append(f"{row['input_id']}: candidate value is blank")
        if status == "自定义确认" and not row["user_value"]:
            errors.append(f"{row['input_id']}: custom confirmation requires user_value")
    return errors


def validate_appliances(rows: list[dict[str, str]], baseline: list[dict[str, str]]) -> list[str]:
    errors = validate_unique(rows, "appliance_id", "appliances")
    expected_ids = {row["appliance_id"] for row in baseline}
    actual_ids = {row["appliance_id"] for row in rows}
    if actual_ids != expected_ids:
        errors.append(f"appliances: IDs changed; missing={sorted(expected_ids-actual_ids)} extra={sorted(actual_ids-expected_ids)}")
    errors.extend(validate_protected_fields(
        rows, baseline, "appliance_id", APPLIANCE_USER_FIELDS, "appliance",
    ))
    for row in rows:
        if row["status"] not in APPLIANCE_STATUSES:
            errors.append(f"{row['appliance_id']}: invalid status {row['status']!r}")
        try:
            quantity = int(float(row["quantity"]))
            if quantity < 0 or (quantity == 0 and row["status"] != "不适用"):
                raise ValueError
            row["quantity"] = str(quantity)
        except ValueError:
            errors.append(
                f"{row['appliance_id']}: quantity must be a positive integer; "
                "zero is allowed only for a not-applicable alias"
            )
        if row["rated_power_w"]:
            try:
                rated_power = float(row["rated_power_w"])
                if rated_power < 0 or (
                    rated_power == 0
                    and row["gas_required"] != "yes"
                    and row["status"] != "不适用"
                ):
                    raise ValueError
            except ValueError:
                errors.append(
                    f"{row['appliance_id']}: rated_power_w must be positive; "
                    "zero is allowed only for a gas appliance with no mains load "
                    "or a not-applicable alias"
                )
        if row["status"] == "已确认":
            for field in ("storage_location_confirmed", "use_location_confirmed", "rated_power_w", "model"):
                if not row[field]:
                    errors.append(f"{row['appliance_id']}: confirmed row requires {field}")
    return errors


def validate_closeout_rules(
    rows: list[dict[str, str]], decisions: list[dict[str, str]],
) -> list[str]:
    errors = validate_unique(rows, "input_id", "closeout rules")
    decision_ids = {row["input_id"] for row in decisions}
    rule_ids = {row["input_id"] for row in rows}
    if rule_ids != decision_ids:
        errors.append(
            "closeout rules: IDs changed; "
            f"missing={sorted(decision_ids-rule_ids)} extra={sorted(rule_ids-decision_ids)}"
        )
    for row in rows:
        if row["automatic_close_allowed"] not in {"yes", "no"}:
            errors.append(
                f"{row['input_id']}: automatic_close_allowed must be yes or no"
            )
        for field in ("closeout_kind", "responsible_party", "required_evidence"):
            if not row[field]:
                errors.append(f"{row['input_id']}: closeout rule requires {field}")
    return errors


def decision_closeout_status(
    decision: dict[str, str], closeout_rule: dict[str, str],
) -> str:
    if decision["status"] == "不适用":
        return "not_applicable"
    if (
        decision["status"] == "自定义确认"
        and closeout_rule["closeout_kind"] == "human_design_selection"
        and decision["user_value"]
        and decision["evidence_reference"]
    ):
        return "closed_by_human_confirmation"
    if decision["status"] in {"采用候选", "自定义确认"}:
        return "decision_confirmed_evidence_pending"
    return "open"


def normalized_inputs(
    decisions: list[dict[str, str]], appliances: list[dict[str, str]],
    closeout_rules: Optional[list[dict[str, str]]] = None,
) -> dict[str, Any]:
    closeout_by_id = {
        row["input_id"]: row for row in (closeout_rules or [])
    }
    normalized_decisions = []
    for row in decisions:
        effective_value = ""
        if row["status"] == "采用候选":
            effective_value = row["candidate_value"]
        elif row["status"] == "自定义确认":
            effective_value = row["user_value"]
        normalized_decisions.append({
            "input_id": row["input_id"],
            "status": row["status"],
            "decision_status": row["status"],
            "closeout_status": (
                decision_closeout_status(row, closeout_by_id[row["input_id"]])
                if row["input_id"] in closeout_by_id
                else "not_evaluated"
            ),
            "candidate_value": row["candidate_value"],
            "effective_value": effective_value,
            "evidence_reference": row["evidence_reference"] if effective_value else "",
            "sync_target": row["sync_target"],
        })

    normalized_appliances = []
    for row in appliances:
        accepts_partial_fields = row["status"] in {"部分确认", "已确认"}
        effective = {
            field: row[field] if accepts_partial_fields else ""
            for field in (
                "storage_location_confirmed", "use_location_confirmed", "quantity",
                "rated_power_w", "simultaneous_group", "water_required",
                "drain_required", "gas_required", "ventilation_required", "model",
                "evidence_reference",
            )
        }
        normalized_appliances.append({
            "appliance_id": row["appliance_id"],
            "status": row["status"],
            "candidate_storage_location": row["storage_location_candidate"],
            "candidate_use_location": row["use_location_candidate"],
            "effective": effective,
        })
    return {"decisions": normalized_decisions, "appliances": normalized_appliances}


def diff_rows(old: list[dict[str, str]], new: list[dict[str, str]], key: str) -> list[dict[str, Any]]:
    old_by_id = {row[key]: row for row in old}
    changes = []
    for row in new:
        previous = old_by_id.get(row[key])
        if previous is None:
            changes.append({"id": row[key], "kind": "added", "changed_fields": sorted(row)})
        else:
            fields = [field for field in row if row[field] != previous.get(field, "")]
            if fields:
                changes.append({"id": row[key], "kind": "updated", "changed_fields": fields})
    return changes


def summary(
    decisions: list[dict[str, str]], appliances: list[dict[str, str]],
    closeout_rules: list[dict[str, str]],
) -> dict[str, Any]:
    closeout_by_id = {row["input_id"]: row for row in closeout_rules}
    closeout_status_by_id = {
        row["input_id"]: decision_closeout_status(
            row, closeout_by_id[row["input_id"]],
        )
        for row in decisions
    }
    release_open = [
        row["input_id"]
        for row in decisions
        if row["blocks_release"] == "yes"
        and closeout_status_by_id[row["input_id"]]
        not in {"closed_by_human_confirmation", "not_applicable"}
    ]
    open_closeout = [closeout_by_id[input_id] for input_id in release_open]
    appliance_open = [row["appliance_id"] for row in appliances if row["status"] not in {"已确认", "不适用"}]
    load_by_group: dict[str, float] = {}
    unknown_power: dict[str, list[str]] = {}
    for row in appliances:
        group = row["simultaneous_group"] or "UNASSIGNED"
        if row["rated_power_w"]:
            load_by_group[group] = load_by_group.get(group, 0.0) + float(row["rated_power_w"]) * int(row["quantity"])
        elif row["status"] != "不适用":
            unknown_power.setdefault(group, []).append(row["appliance_id"])
    return {
        "decision_status_counts": dict(sorted(Counter(row["status"] for row in decisions).items())),
        "decision_closeout_status_counts": dict(sorted(Counter(closeout_status_by_id.values()).items())),
        "appliance_status_counts": dict(sorted(Counter(row["status"] for row in appliances).items())),
        "release_blocking_open_count": len(release_open),
        "release_blocking_open_ids": release_open,
        "open_closeout_kind_counts": dict(sorted(Counter(row["closeout_kind"] for row in open_closeout).items())),
        "open_responsible_party_counts": dict(sorted(Counter(row["responsible_party"] for row in open_closeout).items())),
        "automatic_close_open_count": sum(row["automatic_close_allowed"] == "yes" for row in open_closeout),
        "human_or_external_closeout_open_count": sum(row["automatic_close_allowed"] == "no" for row in open_closeout),
        "appliance_open_count": len(appliance_open),
        "appliance_open_ids": appliance_open,
        "known_connected_load_w_by_simultaneous_group": dict(sorted(load_by_group.items())),
        "unknown_power_ids_by_simultaneous_group": dict(sorted(unknown_power.items())),
    }


def open_items_markdown(
    decisions: list[dict[str, str]], appliances: list[dict[str, str]],
    closeout_rules: list[dict[str, str]], data: dict[str, Any],
) -> str:
    closeout_by_id = {row["input_id"]: row for row in closeout_rules}
    lines = [
        "# 业主输入开放项", "",
        f"生成时间：{datetime.now().astimezone().isoformat(timespec='seconds')}", "",
        f"- 无保留发布阻塞输入：{data['release_blocking_open_count']} 项",
        f"- 可由本地脚本自动关闭：{data['automatic_close_open_count']} 项",
        f"- 必须由人审、现场、厂家或主管方关闭：{data['human_or_external_closeout_open_count']} 项",
        f"- 未完全确认家电：{data['appliance_open_count']} 项",
        "- 候选值只有在状态改为“采用候选”后才视为确认。", "",
        "- “采用候选”只确认设计选择，不会自动关闭施工证据、厂家接口或回路计算。", "",
        "## 无保留发布阻塞输入", "",
        "| ID | 专业 | 关闭类型 | 责任方 | 所需证据 | 状态 | 已填写值 |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in decisions:
        if row["input_id"] in data["release_blocking_open_ids"]:
            rule = closeout_by_id[row["input_id"]]
            values = [
                row["input_id"], row["workstream"], rule["closeout_kind"],
                rule["responsible_party"], rule["required_evidence"],
                row["status"], row["user_value"],
            ]
            lines.append("| " + " | ".join(value.replace("|", "／").replace("\n", " ") for value in values) + " |")
    lines.extend(["", "## 家电信息缺口", "", "| ID | 设备 | 使用点 | 缺少 |", "| --- | --- | --- | --- |"])
    for row in appliances:
        if row["appliance_id"] not in data["appliance_open_ids"]:
            continue
        missing = [field for field in ("storage_location_confirmed", "use_location_confirmed", "rated_power_w", "model") if not row[field]]
        lines.append(f"| {row['appliance_id']} | {row['appliance_name']} | {row['use_location_confirmed'] or row['use_location_candidate']} | {', '.join(missing)} |")
    return "\n".join(lines) + "\n"


def pm_block(data: dict[str, Any]) -> str:
    timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
    return "\n".join([
        PM_START, "### 1.4 用户输入同步状态", "",
        f"- 最近同步：`{timestamp}`",
        f"- 无保留发布阻塞输入仍开放：**{data['release_blocking_open_count']}** 项；不阻止明确披露未决项的待复核候选版。",
        f"- 其中本地脚本可自动关闭：**{data['automatic_close_open_count']}** 项；须由人审、现场、厂家或主管方关闭：**{data['human_or_external_closeout_open_count']}** 项。",
        f"- 家电条目尚未完全确认：**{data['appliance_open_count']}** 项。",
        "- “采用候选”只关闭设计选择，不自动关闭施工证据；设计状态与施工关闭状态分别计算。",
        "- 同步脚本只更新决策登记与本状态块，不直接写正式 IFC。", PM_END,
    ])


def update_pm(path: Path, block: str) -> None:
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(re.escape(PM_START) + r".*?" + re.escape(PM_END), re.DOTALL)
    if pattern.search(text):
        text = pattern.sub(block, text)
    else:
        marker = "## 2. 数据边界与执行规则"
        if marker not in text:
            raise ValueError("PM insertion marker not found")
        text = text.replace(marker, block + "\n\n" + marker, 1)
    path.write_text(text, encoding="utf-8")


def load_readonly_canonical_tables(
    args: argparse.Namespace,
) -> dict[str, tuple[list[str], list[dict[str, str]]]]:
    return {
        sheet_name: read_canonical_csv(getattr(args, spec["argument"]))
        for sheet_name, spec in READONLY_VIEW_SPECS.items()
    }


def main() -> int:
    args = parse_args()
    baseline_decisions = read_csv(args.decisions, DECISION_HEADERS)
    baseline_appliances = read_csv(args.appliances, APPLIANCE_HEADERS)
    closeout_rules = read_csv(args.closeout_rules, CLOSEOUT_HEADERS)
    canonical_tables = load_readonly_canonical_tables(args)
    readonly_parity: Optional[dict[str, Any]] = None
    if args.input.suffix.lower() == ".xlsx":
        decisions, appliances, readonly_rows, workbook_sheet_names = read_xlsx(args.input)
        readonly_parity = readonly_views_parity(
            readonly_rows, workbook_sheet_names, canonical_tables,
        )
    elif args.input.suffix.lower() == ".csv":
        decisions = read_csv(args.input, DECISION_HEADERS)
        appliances = baseline_appliances
    else:
        raise ValueError("input must be .xlsx or .csv")

    errors = (
        validate_decisions(decisions, baseline_decisions)
        + validate_appliances(appliances, baseline_appliances)
        + validate_closeout_rules(closeout_rules, baseline_decisions)
    )
    if errors:
        print("Input validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 2

    data = summary(decisions, appliances, closeout_rules)
    report = {
        "mode": "apply" if args.apply else "dry-run",
        "input": str(args.input.resolve()),
        "formal_ifc_write": False,
        "decision_changes": diff_rows(baseline_decisions, decisions, "input_id"),
        "appliance_changes": diff_rows(baseline_appliances, appliances, "appliance_id"),
        "readonly_view_parity": readonly_parity,
        "summary": data,
        "normalized_inputs": normalized_inputs(decisions, appliances, closeout_rules),
    }
    if args.apply:
        default_appliances = Path(__file__).resolve().parents[2] / "pipeline/decisions/appliance-input-register.csv"
        if args.appliances.resolve() == default_appliances.resolve():
            apply_owner_appliances(Path(__file__).resolve().parents[2], appliances)
            projections(Path(__file__).resolve().parents[2])
        write_csv(args.decisions, DECISION_HEADERS, decisions)
        if args.appliances.resolve() != default_appliances.resolve():
            write_csv(args.appliances, APPLIANCE_HEADERS, appliances)
        update_pm(args.pm, pm_block(data))
        if args.input.suffix.lower() == ".xlsx":
            canonical_tables = load_readonly_canonical_tables(args)
            rebuild_readonly_views(
                args.input,
                canonical_tables,
                {
                    "B13": len(decisions),
                    "B14": data["release_blocking_open_count"],
                    "B15": len(appliances),
                    "B16": sum(row["status"] == "已确认" for row in appliances),
                },
            )
            _, _, rebuilt_rows, rebuilt_sheet_names = read_xlsx(args.input)
            report["readonly_view_parity_after_apply"] = readonly_views_parity(
                rebuilt_rows, rebuilt_sheet_names, canonical_tables,
            )
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if args.open_items:
        args.open_items.parent.mkdir(parents=True, exist_ok=True)
        args.open_items.write_text(
            open_items_markdown(decisions, appliances, closeout_rules, data),
            encoding="utf-8",
        )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if args.apply and args.input.suffix.lower() == ".xlsx":
        if not report["readonly_view_parity_after_apply"]["all_match"]:
            return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
