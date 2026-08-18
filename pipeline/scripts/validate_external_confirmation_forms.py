#!/usr/bin/env python3
"""Mechanical guardrails for the nine current external confirmation forms."""

from __future__ import annotations

import csv
import hashlib
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FORM_DIR = ROOT / "output" / "forms" / "对外确认表"
CURRENT = [
    "01-门组做法确认表-发全屋定制.md",
    "02-日立空调接口确认表-发空调厂家.md",
    "03-给排水设备确认表-发设备与施工方.md",
    "04-燃气消防确认表-发主管单位.md",
    "05-弱电现场记录表-发现场负责人.md",
    "06-门组厂家复核表-发Rimadesio与Poliform.md",
    "07-智能面板电气接口确认表-发JINK与电气方.md",
    "08-定制家具与墙脚节点确认表-发全屋定制.md",
    "09-厨房设备与柜体深化确认表-发橱柜设备方.md",
]
ALLOWED_STATUSES = {
    "confirmed owner decision",
    "official evidence",
    "official/research conclusion",
    "research conclusion",
    "external signoff pending",
    "site measurement pending",
    "authority signoff pending",
    "project internal pending",
    "owner preference pending",
    "unknown",
}
SOURCE_IDS = {
    "01-门组做法确认表-发全屋定制.md": "OUTBOUND-FORM-DOOR-20260817",
    "02-日立空调接口确认表-发空调厂家.md": "OUTBOUND-FORM-HVAC-20260817",
    "03-给排水设备确认表-发设备与施工方.md": "OUTBOUND-FORM-PLUM-20260817",
    "04-燃气消防确认表-发主管单位.md": "OUTBOUND-FORM-GASFIRE-20260817",
    "05-弱电现场记录表-发现场负责人.md": "OUTBOUND-FORM-SITE-20260817",
    "06-门组厂家复核表-发Rimadesio与Poliform.md": "OUTBOUND-FORM-DOOR-VENDOR-20260817",
    "07-智能面板电气接口确认表-发JINK与电气方.md": "OUTBOUND-FORM-SMART-PANEL-20260817",
    "08-定制家具与墙脚节点确认表-发全屋定制.md": "OUTBOUND-FORM-JOINERY-DETAIL-20260817",
    "09-厨房设备与柜体深化确认表-发橱柜设备方.md": "OUTBOUND-FORM-KITCHEN-20260818",
}
HISTORICAL_EXTERNAL_SOURCE_IDS = {
    "A104-SHOP-REQUEST-20260815",
    "RCP1-PLUM-REQUEST-20260815",
    "GAS-CONSULTATION-TEMPLATE-001",
    "A106-CONSULTATION-PACK-20260815",
    "E304-SITE-CHECKLIST-20260815",
}


def read_csv(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle)}


def cells(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def validate_tables(name: str, text: str, errors: list[str]) -> None:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if not line.startswith("|"):
            continue
        header = cells(line)
        if "当前状态" not in header:
            continue
        status_index = header.index("当前状态")
        known_index = next((header.index(item) for item in header if item.startswith("已知")), None)
        action_index = next((header.index(item) for item in header if item.startswith("只需")), None)
        responsibility_index = header.index("责任方") if "责任方" in header else None
        row_index = index + 2
        while row_index < len(lines) and lines[row_index].startswith("|"):
            row = cells(lines[row_index])
            if len(row) != len(header):
                errors.append(f"{name}:{row_index + 1}: table column count differs from header")
                row_index += 1
                continue
            status = row[status_index]
            if status not in ALLOWED_STATUSES:
                errors.append(f"{name}:{row_index + 1}: unsupported status {status!r}")
            if known_index is not None and not row[known_index]:
                errors.append(f"{name}:{row_index + 1}: missing known-information context")
            if action_index is not None and not row[action_index]:
                errors.append(f"{name}:{row_index + 1}: missing remaining action")
            if responsibility_index is not None and not row[responsibility_index]:
                errors.append(f"{name}:{row_index + 1}: missing responsible party")
            row_index += 1


def require(text: str, name: str, tokens: list[str], errors: list[str]) -> None:
    for token in tokens:
        if token not in text:
            errors.append(f"{name}: missing regression token {token!r}")


def main() -> int:
    errors: list[str] = []
    actual = sorted(path.name for path in FORM_DIR.glob("*.md"))
    if actual != sorted(CURRENT):
        errors.append(f"current Markdown set differs: expected {CURRENT}, found {actual}")
    if (ROOT / "drawings" / "evidence" / "EXT-外部信息最短清单-20260817.md").exists():
        errors.append("deleted external-information index has reappeared")

    texts: dict[str, str] = {}
    for name in CURRENT:
        path = FORM_DIR / name
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        texts[name] = text
        if "待填写" in text:
            errors.append(f"{name}: contains contextless 待填写")
        if "外部信息最短清单" in text:
            errors.append(f"{name}: references deleted duplicate dispatch index")
        validate_tables(name, text, errors)

    require(texts.get(CURRENT[0], ""), CURRENT[0], ["已完成", "不再作为厂家或业主填写入口", "向东／图纸右侧滑开", "不要求“两道水平横档”"], errors)
    require(texts.get(CURRENT[1], ""), CURRENT[1], ["现场提供的保温管", "未给材料、导热系数、防火等级或厚度", "旧 IFC 紫色管线不是最终带保温外径", "一只线控器最多控制 6 台室内机"], errors)
    require(texts.get(CURRENT[2], ""), CURRENT[2], ["两个 16A 插座或两个回路", "已被当前方案替代", "网络案例中的“四分铝塑套管＋二分 PE 管”不是项目规格"], errors)
    require(texts.get(CURRENT[3], ""), CURRENT[3], ["ER9EPA33MP", "未获批准", "厨房探测器点位已经确认"], errors)
    require(texts.get(CURRENT[4], ""), CURRENT[4], ["现有约测 110 mm，待带尺复核", "525 mm", "H+350 mm", "(4600.016, -735.369) mm", "玄关高柜右下柜格", "五孔插座", "入户临时置物位"], errors)
    require(texts.get(CURRENT[5], ""), CURRENT[5], ["只答“是／否”", "项目加工图", "地面无通长下轨"], errors)
    require(texts.get(CURRENT[6], ""), CURRENT[6], ["至少 4 键", "客厅双控、书房双控、餐厅双控、照明总控", "真实物理接线双控", "氛围 LED、重点射灯"], errors)
    require(texts.get(CURRENT[7], ""), CURRENT[7], ["350–450 mm", "owner_reference_candidate_not_final", "业主偏好"], errors)
    require(texts.get(CURRENT[8], ""), CURRENT[8], ["LS33R6VB9W/01", "895×345×873 mm", "不得把 APP-015 重复计算", "灶台背板为台面同材大板，远端墙为浅置物架"], errors)

    sources = read_csv(ROOT / "pipeline" / "decisions" / "source-evidence-register.csv", "source_id")
    owner_inputs = read_csv(ROOT / "pipeline" / "decisions" / "owner-input-register.csv", "input_id")
    equipment = read_csv(ROOT / "pipeline" / "decisions" / "equipment-register.csv", "equipment_id")
    requirements = read_csv(ROOT / "pipeline" / "decisions" / "equipment-installation-requirements.csv", "requirement_id")
    for name, source_id in SOURCE_IDS.items():
        if source_id not in sources:
            errors.append(f"source register missing {source_id}")
            continue
        expected = hashlib.sha256((FORM_DIR / name).read_bytes()).hexdigest()
        if sources[source_id]["sha256"] != expected:
            errors.append(f"{source_id}: stored hash differs from current Markdown")
        if "现行唯一可编辑源" not in sources[source_id]["notes"]:
            errors.append(f"{source_id}: not marked as current unique editable source")
    for source_id in HISTORICAL_EXTERNAL_SOURCE_IDS:
        row = sources.get(source_id)
        if row is None:
            errors.append(f"historical source register missing {source_id}")
            continue
        if row["status"] != "superseded_historical_evidence_current_form_linked":
            errors.append(f"{source_id}: historical request is not retired from current filling entry")
        if "现行填写入口" not in row["notes"]:
            errors.append(f"{source_id}: historical request does not point to its current form")

    ssot_checks = {
        "E304-CABINET-DIMENSIONS": ["110", "525", "H+350", "4600.016", "只缺准确净宽"],
        "RCP1-HVAC-PORTS": ["未给保温材料", "珠海露点", "旧 IFC 紫色管线"],
        "INT1-ENTRY-PARCEL-FUNCTION": ["临时置物", "不是专用快递柜"],
        "INT1-ENTRY-PARCEL-LAYOUT": ["I-503 前", "不要求门开／门关"],
        "APP014-REPLACEABLE-SLEEVE": ["连续可抽换", "网络案例规格不作为项目规格"],
    }
    for item_id, tokens in ssot_checks.items():
        row = owner_inputs.get(item_id)
        if row is None:
            errors.append(f"owner input missing {item_id}")
            continue
        joined = " ".join(row.values())
        for token in tokens:
            if token not in joined:
                errors.append(f"{item_id}: missing SSOT token {token!r}")

    if "at least 4" not in equipment["CTRL-ENTRY-A"]["variant"] and "minimum 4" not in equipment["CTRL-ENTRY-A"]["variant"]:
        errors.append("CTRL-ENTRY-A equipment row does not preserve minimum four functions")
    if "3-Key/3-Relay（研究优选" in equipment["CTRL-ENTRY-A"]["model"]:
        errors.append("CTRL-ENTRY-A still promotes incompatible 3-key candidate")
    for req_id in [f"REQ-APP014-SLEEVE-20260818-00{index}" for index in range(1, 6)]:
        if req_id not in requirements:
            errors.append(f"installation requirements missing {req_id}")
    if "EXT-EXTERNAL-INFO-MINIMUM-20260817" in sources:
        errors.append("deleted duplicate dispatch index remains in source register")

    history = (ROOT / "drawings" / "evidence" / "E304-现场最短取证清单-20260815.md").read_text(encoding="utf-8")
    if "不是当前填写入口" not in history or "当前唯一可编辑入口" not in history:
        errors.append("historical E304 checklist is not clearly retired")

    if errors:
        print("external confirmation form validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 2
    print(f"external confirmation form validation passed: {len(CURRENT)} current forms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
