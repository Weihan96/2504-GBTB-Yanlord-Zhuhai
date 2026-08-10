#!/usr/bin/env python3
"""Validate and synchronize the owner input workbook into project registers.

The default mode is a read-only dry run. ``--apply`` updates only the owner
input CSV registers and the generated PM status block. It never writes IFC.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET


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


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "output/forms/滨海湾施工输入清单.xlsx")
    parser.add_argument("--decisions", type=Path, default=root / "pipeline/decisions/owner-input-register.csv")
    parser.add_argument("--appliances", type=Path, default=root / "pipeline/decisions/appliance-input-register.csv")
    parser.add_argument("--pm", type=Path, default=root / "drawings/滨海湾装修施工图深化工作管理.md")
    parser.add_argument("--report", type=Path, default=root / "build/owner-inputs/sync-preview.json")
    parser.add_argument("--open-items", type=Path, default=root / "build/owner-inputs/open-inputs.md")
    parser.add_argument("--apply", action="store_true")
    return parser.parse_args()


def column_index(reference: str) -> int:
    letters = "".join(character for character in reference if character.isalpha())
    value = 0
    for character in letters.upper():
        value = value * 26 + ord(character) - ord("A") + 1
    return value - 1


def read_xlsx(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    with zipfile.ZipFile(path) as archive:
        strings = _shared_strings(archive)
        paths = _sheet_paths(archive)
        missing = {"设计决策", "家电清单"} - set(paths)
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
    return decisions, appliances


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


def write_csv(path: Path, headers: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate_unique(rows: list[dict[str, str]], key: str, label: str) -> list[str]:
    ids = [row[key] for row in rows]
    errors = [f"{label}: blank {key}"] if any(not value for value in ids) else []
    duplicates = sorted(value for value, count in Counter(ids).items() if count > 1)
    if duplicates:
        errors.append(f"{label}: duplicate IDs: {', '.join(duplicates)}")
    return errors


def validate_decisions(rows: list[dict[str, str]], baseline: list[dict[str, str]]) -> list[str]:
    errors = validate_unique(rows, "input_id", "decisions")
    expected_ids = {row["input_id"] for row in baseline}
    actual_ids = {row["input_id"] for row in rows}
    if actual_ids != expected_ids:
        errors.append(f"decisions: IDs changed; missing={sorted(expected_ids-actual_ids)} extra={sorted(actual_ids-expected_ids)}")
    for row in rows:
        status = row["status"]
        if status not in DECISION_STATUSES:
            errors.append(f"{row['input_id']}: invalid status {status!r}")
        if status == "采用候选" and not row["candidate_value"]:
            errors.append(f"{row['input_id']}: candidate value is blank")
        if status == "自定义确认" and not row["user_value"]:
            errors.append(f"{row['input_id']}: custom confirmation requires user_value")
    return errors


def validate_appliances(rows: list[dict[str, str]]) -> list[str]:
    errors = validate_unique(rows, "appliance_id", "appliances")
    for row in rows:
        if row["status"] not in APPLIANCE_STATUSES:
            errors.append(f"{row['appliance_id']}: invalid status {row['status']!r}")
        try:
            quantity = int(float(row["quantity"]))
            if quantity < 1:
                raise ValueError
            row["quantity"] = str(quantity)
        except ValueError:
            errors.append(f"{row['appliance_id']}: quantity must be a positive integer")
        if row["rated_power_w"]:
            try:
                if float(row["rated_power_w"]) <= 0:
                    raise ValueError
            except ValueError:
                errors.append(f"{row['appliance_id']}: rated_power_w must be positive")
        if row["status"] == "已确认":
            for field in ("storage_location_confirmed", "use_location_confirmed", "rated_power_w", "model"):
                if not row[field]:
                    errors.append(f"{row['appliance_id']}: confirmed row requires {field}")
    return errors


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


def summary(decisions: list[dict[str, str]], appliances: list[dict[str, str]]) -> dict[str, Any]:
    release_open = [row["input_id"] for row in decisions if row["blocks_release"] == "yes" and row["status"] not in {"采用候选", "自定义确认", "不适用"}]
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
        "appliance_status_counts": dict(sorted(Counter(row["status"] for row in appliances).items())),
        "release_blocking_open_count": len(release_open),
        "release_blocking_open_ids": release_open,
        "appliance_open_count": len(appliance_open),
        "appliance_open_ids": appliance_open,
        "known_connected_load_w_by_simultaneous_group": dict(sorted(load_by_group.items())),
        "unknown_power_ids_by_simultaneous_group": dict(sorted(unknown_power.items())),
    }


def open_items_markdown(decisions: list[dict[str, str]], appliances: list[dict[str, str]], data: dict[str, Any]) -> str:
    lines = [
        "# 业主输入开放项", "",
        f"生成时间：{datetime.now().astimezone().isoformat(timespec='seconds')}", "",
        f"- 无保留发布阻塞输入：{data['release_blocking_open_count']} 项",
        f"- 未完全确认家电：{data['appliance_open_count']} 项",
        "- 候选值只有在状态改为“采用候选”后才视为确认。", "",
        "## 无保留发布阻塞输入", "",
        "| ID | 专业 | 问题 | 候选 | 状态 | 已填写值 |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in decisions:
        if row["input_id"] in data["release_blocking_open_ids"]:
            values = [row["input_id"], row["workstream"], row["question"], row["candidate_value"], row["status"], row["user_value"]]
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
        f"- 家电条目尚未完全确认：**{data['appliance_open_count']}** 项。",
        "- 候选值不等于确认；只有“采用候选”或“自定义确认”的设计决策才关闭输入项。",
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


def main() -> int:
    args = parse_args()
    baseline_decisions = read_csv(args.decisions, DECISION_HEADERS)
    baseline_appliances = read_csv(args.appliances, APPLIANCE_HEADERS)
    if args.input.suffix.lower() == ".xlsx":
        decisions, appliances = read_xlsx(args.input)
    elif args.input.suffix.lower() == ".csv":
        decisions = read_csv(args.input, DECISION_HEADERS)
        appliances = baseline_appliances
    else:
        raise ValueError("input must be .xlsx or .csv")

    errors = validate_decisions(decisions, baseline_decisions) + validate_appliances(appliances)
    if errors:
        print("Input validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 2

    data = summary(decisions, appliances)
    report = {
        "mode": "apply" if args.apply else "dry-run",
        "input": str(args.input.resolve()),
        "formal_ifc_write": False,
        "decision_changes": diff_rows(baseline_decisions, decisions, "input_id"),
        "appliance_changes": diff_rows(baseline_appliances, appliances, "appliance_id"),
        "summary": data,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    args.open_items.parent.mkdir(parents=True, exist_ok=True)
    args.open_items.write_text(open_items_markdown(decisions, appliances, data), encoding="utf-8")
    if args.apply:
        write_csv(args.decisions, DECISION_HEADERS, decisions)
        write_csv(args.appliances, APPLIANCE_HEADERS, appliances)
        update_pm(args.pm, pm_block(data))
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
