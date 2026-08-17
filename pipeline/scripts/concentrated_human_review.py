#!/usr/bin/env python3
"""Compile canonical decision and requirement blockers into two derived review views."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


CONFIRMED_DECISION_STATUSES = {"采用候选", "自定义确认", "已确认", "confirmed"}
NOT_APPLICABLE_DECISION_STATUSES = {"不适用"}
KNOWN_DECISION_STATUSES = CONFIRMED_DECISION_STATUSES | NOT_APPLICABLE_DECISION_STATUSES | {"待填写", "需证据", "暂缓"}
CLOSED_CLOSEOUT_STATUSES = {"verified", "not_applicable"}
EVIDENCE_GATED_CLOSEOUT_KINDS = {
    "authority_product_evidence",
    "human_material_selection",
    "human_product_selection",
    "product_interface_evidence",
    "site_measurement",
    "site_or_drawing_evidence",
    "site_test",
}


INT1_PACKAGES = (
    {
        "package_id": "INT1-GEBERIT-FINISH-INTERFACE",
        "title": "Geberit 完成面与安装标高",
        "equipment_ids": {"SAN-001", "SAN-003", "SAN-014"},
        "root_cause": "标准产品资料不能替代本项目包管厚度、面板标高和地砖排版",
        "affected_sheets": ["I-502", "P-202", "D-602"],
        "responsible_party": "建筑设计/给排水设计/现场",
        "required_evidence": "卫生间完成面立面、地砖排版和可定位到实际构造的安装节点",
    },
    {
        "package_id": "INT1-FOSTER-KITCHEN-SINK",
        "title": "Foster 厨房水槽加工与柜内净空",
        "equipment_ids": {"SAN-KIT-001"},
        "root_cause": "公开资料和 IFC 不能关闭项目安装身份、存水弯净空与模板加工公差",
        "affected_sheets": ["I-501", "P-201", "P-202"],
        "responsible_party": "橱柜深化/台面加工/给排水设计",
        "required_evidence": "随货开孔模板、实物复核和含柜内存水弯检修净空的厨房加工图",
    },
    {
        "package_id": "INT1-ISLAND-CUSTOM-SINK",
        "title": "岛台定制石材水槽加工尺寸",
        "equipment_ids": {"SAN-022"},
        "root_cause": "定制石材水槽尚无可施工外尺寸、盆腔尺寸和深度",
        "affected_sheets": ["I-501", "P-201", "P-202"],
        "responsible_party": "石材加工/家具深化/给排水设计",
        "required_evidence": "盖章加工图，明确成品外尺寸、盆腔平面、深度及与排水接口的关系",
    },
    {
        "package_id": "INT1-CUSTOM-FLOOR-DRAIN",
        "title": "厕所定制水母地漏／中央集水器与线性排水渠",
        "equipment_ids": {"DRAIN-CUSTOM-001"},
        "root_cause": "业主只确认定制渠道和组件方向，尚无逐房间数量、准确组件、加工尺寸、排水与防水接口或与吉博力构件的分配关系",
        "affected_sheets": ["P-202", "I-502", "D-602", "S-701"],
        "responsible_party": "定制地漏商家/给排水设计/防水/全屋定制",
        "required_evidence": "逐卫生间盖章 shop drawing 与组件清单，明确数量和房间映射、长宽深、材质表面、过滤、水封、设计流量、出水口、完成标高、坡向、防水法兰、清扫检修及与吉博力构件的沿用／替换／连接关系",
    },
    {
        "package_id": "INT1-ISLAND-WORKTOP",
        "title": "岛台台面厚度、边型与统一开孔表",
        "equipment_ids": {"FURIFC-003"},
        "root_cause": "定制台面厚度、边型、设备开孔和加工图尚未统一",
        "affected_sheets": ["I-501", "E-303"],
        "responsible_party": "家具深化/石材加工/电气设计",
        "required_evidence": "岛台盖章加工图及水槽、龙头、插座等全部实物模板汇总开孔表",
    },
    {
        "package_id": "INT1-MASTER-BATH-STREET",
        "title": "主卫 antoniolupi Street 定制台盆",
        "equipment_ids": {"SAN-021"},
        "root_cause": "国内定制 Street 的成品尺寸、盆数和盆宽不能由 IFC 包围盒推定",
        "affected_sheets": ["I-502", "P-201", "P-202"],
        "responsible_party": "卫浴加工/建筑设计/给排水设计",
        "required_evidence": "主卫立面与盖章加工图，明确长深高、盆数、盆宽和接口开孔",
    },
    {
        "package_id": "INT1-GUEST-BATH-SORGENTE",
        "title": "客卫 Falper Sorgente 定制台盆",
        "equipment_ids": {"SAN-018"},
        "root_cause": "国内定制 Sorgente 尚无项目成品尺寸和盆体几何",
        "affected_sheets": ["I-502", "P-201", "P-202"],
        "responsible_party": "卫浴加工/建筑设计/给排水设计",
        "required_evidence": "客卫干区立面与盖章加工图，明确最终尺寸、盆体、固定和排水定位",
    },
    {
        "package_id": "INT1-GUEST-BATH-MIRROR-CABINET",
        "title": "客卫干区全高镜柜",
        "equipment_ids": {"FUR-BATHG-MIRROR-001"},
        "root_cause": "视觉参考不能证明项目尺寸、柜深、门扇开启与内部构造",
        "affected_sheets": ["I-502", "I-504"],
        "responsible_party": "家具深化/建筑设计/业主",
        "required_evidence": "客卫干区立面和镜柜加工图，明确宽高深、分格开启、层板通风与满载条件",
    },
    {
        "package_id": "INT1-JINK-CONTROL-PANELS",
        "title": "JINK EGG 控制面板嵌装适配",
        "equipment_ids": {"CTRL-ENTRY-A", "CTRL-MASTER-A", "CTRL-MASTER-B"},
        "root_cause": "候选面板族外形不能替代准确 SKU、端子图、底盒净深和本项目平嵌收口验证",
        "affected_sheets": ["E-302", "I-504"],
        "responsible_party": "电气设计/设备方/建筑设计",
        "required_evidence": "准确 SKU、厂家端子图、底盒与饰面节点及样板安装复核",
    },
    {
        "package_id": "INT1-SIEMENS-LAUNDRY-STACK",
        "title": "西门子洗衣机与热泵干衣机叠放连接",
        "equipment_ids": {"APP-017", "APP-017-WASHER", "APP-017-DRYER"},
        "root_cause": "说明书订货号 17008829 已确认，但抽板配置、WTZ27510 商业型号映射和两台 E-Nr./FD 逐型号兼容尚未书面关闭",
        "affected_sheets": ["I-503", "E-303", "P-201", "S-701"],
        "responsible_party": "西门子客服/设备方/橱柜深化",
        "required_evidence": "西门子书面确认本两台设备的原厂连接件商业型号、抽板配置和逐型号兼容性",
    },
    {
        "package_id": "INT1-ICE-CREAM-MACHINE-SELECTION",
        "title": "冰淇淋机选型与安装接口",
        "equipment_ids": {"APP-019"},
        "root_cause": "业主只确认新增一台冰淇淋机，未确认位置、安装形式或准确型号",
        "affected_sheets": ["I-501", "E-303", "E-304", "P-201", "P-202"],
        "responsible_party": "业主/商家/设备厂家",
        "required_evidence": "确认使用位置与安装形式，并提供准确型号铭牌或官方安装图，含功率、插头、给排水、通风、机身尺寸、检修净空和接口定位",
    },
    {
        "package_id": "INT1-KITCHEN-GAS-ALARM",
        "title": "中厨燃气报警器准入与平嵌条件",
        "equipment_ids": {"SENSOR-GAS-R04"},
        "root_cause": "研究候选不能替代珠海燃气公司的准入、认证、联动阀及平嵌书面许可",
        "affected_sheets": ["A-106", "E-303", "I-501"],
        "responsible_party": "燃气公司/设备方/电气设计",
        "required_evidence": "燃气公司书面回复及最终产品安装图，含供电、联动、标高、禁距与收口许可",
    },
    {
        "package_id": "INT1-KITCHEN-FIRE-SENSOR",
        "title": "中厨火灾探测器最终类型与装饰收口",
        "equipment_ids": {"SENSOR-001"},
        "root_cause": "已确认点位不等于已确认感温、感烟或复合产品及其平嵌装饰条件",
        "affected_sheets": ["A-106", "E-301", "I-501"],
        "responsible_party": "消防/设备方/建筑设计",
        "required_evidence": "消防或设备方确认最终探测类型、准确型号、供电通信、安装禁距和装饰收口许可",
    },
)


DISCIPLINE_PACKAGES = {
    "ARCH": {
        "package_id": "REQ-ARCH-PROJECT-DIMENSIONS",
        "title": "门窗项目尺寸与开启关系",
        "root_cause": "产品系列或 IFC 现状不能替代项目门洞尺寸与开启证据",
        "affected_sheets": ["A-104", "I-504"],
        "responsible_party": "建筑设计/全屋定制/现场",
        "required_evidence": "项目订单、门窗表、厂家 shop drawing 或可定位现场证据",
    },
    "ELEC": {
        "package_id": "REQ-ELEC-FINAL-POWER",
        "title": "设备最终功率与电源接口",
        "root_cause": "候选功率和未选型号不能关闭回路容量或精确电源接口",
        "affected_sheets": ["E-303", "E-304"],
        "responsible_party": "设备方/电气设计/业主",
        "required_evidence": "最终型号铭牌、官方说明书和项目回路计算或接线接口图",
    },
    "HVAC/PLUM": {
        "package_id": "REQ-HVAC-FINAL-INTERFACES",
        "title": "空调最终型号与冷媒冷凝水接口",
        "root_cause": "候选机位不能替代最终设备尺寸、管径、坡度与接口坐标",
        "affected_sheets": ["M-401"],
        "responsible_party": "空调厂家/机电设计",
        "required_evidence": "最终型号厂家安装图或 BIM/DWG 接口图，含管径、坡度和接口方向",
    },
    "MULTI": {
        "package_id": "REQ-MULTI-INTERFACE-COORDINATION",
        "title": "跨专业隐藏接口与检修协调",
        "root_cause": "单专业产品资料不能关闭隐藏接口、检修和相邻构造关系",
        "affected_sheets": ["I-502", "P-202", "E-303"],
        "responsible_party": "建筑设计/机电设计/设备方",
        "required_evidence": "跨专业综合节点，标注隐藏接口坐标、检修路径及完成面关系",
    },
    "PLUM": {
        "package_id": "REQ-PLUM-FINAL-ROUGHINS",
        "title": "给排水最终粗装接口",
        "root_cause": "产品族资料和现状端点不能关闭项目给排水坐标、管径与坡度",
        "affected_sheets": ["P-201", "P-202", "I-501", "I-502"],
        "responsible_party": "设备方/给排水设计/现场",
        "required_evidence": "最终产品安装图、项目定位尺寸和现场管路复核记录",
    },
    "PROCUREMENT": {
        "package_id": "REQ-PROCUREMENT-IDENTITY-QUANTITY",
        "title": "采购身份、配置与数量差额",
        "root_cause": "系列意向或部分收货证据不能证明最终配置和完整数量",
        "affected_sheets": ["A-001", "I-501", "I-502"],
        "responsible_party": "业主/采购/设备方",
        "required_evidence": "最终订单、精确 SKU、数量清单和收货核对照片",
    },
    "SAFETY": {
        "package_id": "REQ-SAFETY-PRODUCT-SPEC",
        "title": "安全构造与产品规格",
        "root_cause": "视觉参考或普通材料描述不能证明安全构造符合项目要求",
        "affected_sheets": ["I-502", "I-504"],
        "responsible_party": "建筑设计/家具深化/产品方",
        "required_evidence": "产品安全规格、适用标准和项目构造节点",
    },
    "STRUCTURE": {
        "package_id": "REQ-STRUCTURE-SUPPORT-FIXING",
        "title": "设备、石材与家具支撑固定",
        "root_cause": "产品或造型意向不能替代实际墙体、跨度和满载条件下的支撑设计",
        "affected_sheets": ["I-501", "I-502", "I-504"],
        "responsible_party": "结构/建筑设计/加工方",
        "required_evidence": "按实际基层、跨度、开孔和满载条件完成的支撑固定节点或盖章加工图",
    },
    "WFIN": {
        "package_id": "REQ-WFIN-MATERIAL-WATERPROOFING",
        "title": "饰面、防水与密封体系",
        "root_cause": "材料意向或设备自带构造不能替代项目饰面、防水和密封系统",
        "affected_sheets": ["D-601", "D-602", "I-501", "I-502"],
        "responsible_party": "材料设计/建筑设计/加工方",
        "required_evidence": "最终材料体系、样板确认和与设备接口一致的防水密封节点",
    },
}


WORKSTREAM_SHEETS = {
    "A-103": ["A-103"],
    "A-104": ["A-104"],
    "E-302": ["E-302"],
    "E-303": ["E-303"],
    "E-304": ["E-304"],
    "A-106": ["A-106"],
    "WFIN": ["D-601"],
    "DET1": ["D-601", "D-602"],
    "RCP1": ["RCP-101", "M-401"],
    "PLUM": ["P-201", "P-202"],
    "INT1": ["I-501", "I-502", "I-503", "I-504"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument("--owner-input-register", type=Path, default=Path("pipeline/decisions/owner-input-register.csv"))
    parser.add_argument("--owner-input-closeout-rules", type=Path, default=Path("pipeline/decisions/owner-input-closeout-rules.csv"))
    parser.add_argument("--equipment-register", type=Path, default=Path("pipeline/decisions/equipment-register.csv"))
    parser.add_argument("--requirements", type=Path, default=Path("pipeline/decisions/equipment-installation-requirements.csv"))
    parser.add_argument("--drawing-register", type=Path, default=Path("pipeline/decisions/drawing-register.csv"))
    parser.add_argument("--reviewed-release", type=Path, default=Path("build/release/reviewed-candidate-current.json"))
    parser.add_argument("--construction-release", type=Path, default=Path("build/release/construction-candidate-current.json"))
    parser.add_argument("--json-output", type=Path, default=Path("build/release/concentrated-human-review-current.json"))
    parser.add_argument("--markdown-output", type=Path, default=Path("build/release/concentrated-human-review-current.md"))
    return parser.parse_args()


def resolve(root: Path, path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path, required: set[str], id_field: str) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"missing canonical CSV: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise ValueError(f"{path.name} missing fields: {', '.join(missing)}")
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    ids = [row[id_field] for row in rows]
    if any(not item for item in ids):
        raise ValueError(f"{path.name} contains blank {id_field}")
    duplicates = sorted(item for item, count in Counter(ids).items() if count > 1)
    if duplicates:
        raise ValueError(f"{path.name} duplicate {id_field}: {', '.join(duplicates)}")
    return rows


def optional_release_context(path: Path, formal_ifc_hash: str) -> dict[str, Any]:
    if not path.is_file():
        return {
            "path": str(path),
            "available": False,
            "current": False,
            "ifc_sha256": "",
            "result": {},
            "error": "current release report is not available",
        }
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return {
            "path": str(path),
            "available": True,
            "current": False,
            "ifc_sha256": "",
            "result": {},
            "error": f"cannot read release report: {exc}",
        }
    if not isinstance(value, dict):
        return {
            "path": str(path),
            "available": True,
            "current": False,
            "ifc_sha256": "",
            "result": {},
            "error": "release report is not a JSON object",
        }
    value_hash = value.get("source", {}).get("formal_ifc_sha256")
    valid_hash = isinstance(value_hash, str) and len(value_hash) == 64
    current = valid_hash and value_hash == formal_ifc_hash
    return {
        "path": str(path),
        "available": True,
        "current": current,
        "ifc_sha256": value_hash if valid_hash else "",
        "result": value.get("result", {}),
        "error": "" if current else "release report does not match the formal IFC",
    }


def affected_owner_sheets(workstream: str, drawing_ids: set[str]) -> list[str]:
    result: list[str] = []
    for token in (part.strip() for part in workstream.split("/") if part.strip()):
        for sheet in WORKSTREAM_SHEETS.get(token, [token] if token in drawing_ids else []):
            if sheet in drawing_ids and sheet not in result:
                result.append(sheet)
    if not result:
        raise ValueError(f"unmapped owner workstream: {workstream}")
    return result


def closeout_status(decision: dict[str, str], rule: dict[str, str]) -> str:
    status = decision["status"]
    if status not in KNOWN_DECISION_STATUSES:
        raise ValueError(f"unknown decision status {status!r} for {decision['input_id']}")
    automatic = rule["automatic_close_allowed"].lower()
    if automatic not in {"yes", "no"}:
        raise ValueError(f"invalid automatic_close_allowed for {decision['input_id']}: {automatic!r}")
    if status in NOT_APPLICABLE_DECISION_STATUSES:
        return "not_applicable"
    if status == "暂缓":
        return "deferred_blocker"
    if status == "待填写":
        return "decision_pending"
    if status == "需证据":
        return "decision_pending_evidence_required"
    if status == "采用候选" and automatic == "no":
        return "decision_confirmed_evidence_pending"
    if rule["closeout_kind"] in EVIDENCE_GATED_CLOSEOUT_KINDS and automatic == "no":
        return "decision_confirmed_evidence_pending"
    if (
        status in CONFIRMED_DECISION_STATUSES
        and rule["closeout_kind"] == "human_design_selection"
        and decision["user_value"]
        and decision["evidence_reference"]
    ):
        return "verified"
    if status in CONFIRMED_DECISION_STATUSES and automatic == "yes" and decision["evidence_reference"]:
        return "verified"
    return "decision_confirmed_evidence_pending"


def owner_review_items(
    decisions: list[dict[str, str]],
    rules: list[dict[str, str]],
    drawing_ids: set[str],
) -> list[dict[str, Any]]:
    decision_ids = {row["input_id"] for row in decisions}
    rule_ids = {row["input_id"] for row in rules}
    if decision_ids != rule_ids:
        raise ValueError(
            "owner input/rule ID mismatch; missing rules="
            f"{sorted(decision_ids - rule_ids)}, orphan rules={sorted(rule_ids - decision_ids)}"
        )
    rules_by_id = {row["input_id"]: row for row in rules}
    items: list[dict[str, Any]] = []
    for decision in decisions:
        rule = rules_by_id[decision["input_id"]]
        items.append(
            {
                "input_id": decision["input_id"],
                "workstream": decision["workstream"],
                "priority": decision["priority"],
                "blocks_release": decision["blocks_release"].lower() == "yes",
                "question": decision["question"],
                "candidate_value": decision["candidate_value"],
                "user_value": decision["user_value"],
                "decision_status": decision["status"],
                "closeout_status": closeout_status(decision, rule),
                "closeout_kind": rule["closeout_kind"],
                "responsible_party": rule["responsible_party"],
                "required_evidence": rule["required_evidence"],
                "evidence_reference": decision["evidence_reference"],
                "affected_sheets": affected_owner_sheets(decision["workstream"], drawing_ids),
            }
        )
    return items


def make_requirement_package(
    definition: dict[str, Any],
    rows: list[dict[str, str]],
    equipment: dict[str, dict[str, str]],
    drawing_ids: set[str],
    root_cause_id: str,
) -> dict[str, Any]:
    missing_sheets = sorted(set(definition["affected_sheets"]) - drawing_ids)
    if missing_sheets:
        raise ValueError(f"package {definition['package_id']} references unknown sheets: {missing_sheets}")
    equipment_ids = sorted({row["equipment_id"] for row in rows})
    missing_equipment = sorted(set(equipment_ids) - set(equipment))
    if missing_equipment:
        raise ValueError(f"package {definition['package_id']} references unknown equipment: {missing_equipment}")
    return {
        "package_id": definition["package_id"],
        "title": definition["title"],
        "root_cause_id": root_cause_id,
        "root_cause": definition["root_cause"],
        "requirement_ids": sorted(row["requirement_id"] for row in rows),
        "equipment_ids": equipment_ids,
        "equipment": [
            {
                "equipment_id": item,
                "item_name": equipment[item]["item_name"],
                "manufacturer": equipment[item]["manufacturer"],
                "model": equipment[item]["model"],
            }
            for item in equipment_ids
        ],
        "affected_sheets": list(definition["affected_sheets"]),
        "responsible_party": definition["responsible_party"],
        "required_evidence": definition["required_evidence"],
    }


def requirement_review_packages(
    requirements: list[dict[str, str]],
    equipment_rows: list[dict[str, str]],
    drawing_ids: set[str],
) -> list[dict[str, Any]]:
    equipment = {row["equipment_id"]: row for row in equipment_rows}
    blockers = [row for row in requirements if row["blocks_release"].lower() == "yes"]
    packages: list[dict[str, Any]] = []
    mapped: set[str] = set()

    int1_rows = [
        row for row in blockers
        if "INT1" in {part.strip() for part in row["discipline"].upper().split("/")}
    ]
    for definition in INT1_PACKAGES:
        # Stable root-cause packages are equipment-owned, so capture every
        # blocking discipline for that equipment instead of losing compound
        # values such as ELEC/INT1, ELEC/GAS, or ELEC/PLUM/INT1.
        rows = [row for row in blockers if row["equipment_id"] in definition["equipment_ids"]]
        if not rows:
            raise ValueError(f"INT1 package has no current requirements: {definition['package_id']}")
        packages.append(
            make_requirement_package(definition, rows, equipment, drawing_ids, definition["package_id"])
        )
        mapped.update(row["requirement_id"] for row in rows)

    unmapped_int1 = sorted(row["requirement_id"] for row in int1_rows if row["requirement_id"] not in mapped)
    if unmapped_int1:
        raise ValueError(f"unmapped INT1 root causes: {', '.join(unmapped_int1)}")

    for discipline, definition in DISCIPLINE_PACKAGES.items():
        rows = [
            row
            for row in blockers
            if row["discipline"].upper() == discipline and row["requirement_id"] not in mapped
        ]
        if not rows:
            continue
        packages.append(
            make_requirement_package(definition, rows, equipment, drawing_ids, f"DISCIPLINE-{discipline}")
        )
        mapped.update(row["requirement_id"] for row in rows)

    unmapped = sorted(row["requirement_id"] for row in blockers if row["requirement_id"] not in mapped)
    if unmapped:
        details = [
            f"{row['requirement_id']}:{row['discipline']}:{row['parameter_key']}"
            for row in blockers
            if row["requirement_id"] in unmapped
        ]
        raise ValueError("unmapped release-blocking root causes: " + ", ".join(details))
    if len(mapped) != len(blockers):
        raise ValueError("release-blocking requirement mapping is not one-to-one")
    return packages


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# 集中人审清单（派生视图）",
        "",
        "> 本报告只读现有决策与设备单一真源，不创建新决策真源；禁止写入正式 IFC。",
        "",
        f"- 生成时间（UTC）：{report['generated_at_utc']}",
        f"- 正式 IFC SHA-256：`{report['formal_ifc']['sha256']}`",
        f"- 业主输入：{report['summary']['owner_input_count']} 条；待关闭：{report['summary']['owner_closeout_open_count']} 条",
        f"- release-blocking requirement：{report['summary']['blocking_requirement_count']} 条；聚合为 {report['summary']['requirement_review_package_count']} 包",
        f"- INT1：{report['summary']['int1_requirement_count']} 条；稳定聚合为 {report['summary']['int1_review_package_count']} 包",
        "",
        "## 业主输入关闭状态",
        "",
        "| ID | 决策状态 | 关闭状态 | 责任方 | 所需证据 | 影响图纸 |",
        "|---|---|---|---|---|---|",
    ]
    for item in report["decision_review_items"]:
        lines.append(
            "| {input_id} | {decision_status} | {closeout_status} | {responsible_party} | {required_evidence} | {sheets} |".format(
                **item,
                sheets=", ".join(item["affected_sheets"]),
            )
        )
    lines.extend(["", "## Requirement 根因审核包", ""])
    for package in report["requirement_review_packages"]:
        lines.extend(
            [
                f"### {package['package_id']} · {package['title']}",
                "",
                f"- 根因：{package['root_cause']}",
                f"- Requirement：{', '.join(package['requirement_ids'])}",
                f"- 设备：{', '.join(package['equipment_ids'])}",
                f"- 影响图纸：{', '.join(package['affected_sheets'])}",
                f"- 责任方：{package['responsible_party']}",
                f"- 所需证据：{package['required_evidence']}",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    paths = {
        "ifc": resolve(root, args.ifc),
        "owner_inputs": resolve(root, args.owner_input_register),
        "closeout_rules": resolve(root, args.owner_input_closeout_rules),
        "equipment": resolve(root, args.equipment_register),
        "requirements": resolve(root, args.requirements),
        "drawings": resolve(root, args.drawing_register),
        "reviewed_release": resolve(root, args.reviewed_release),
        "construction_release": resolve(root, args.construction_release),
        "json_output": resolve(root, args.json_output),
        "markdown_output": resolve(root, args.markdown_output),
    }
    try:
        if not paths["ifc"].is_file():
            raise ValueError(f"formal IFC does not exist: {paths['ifc']}")
        formal_ifc_hash = sha256(paths["ifc"])
        decisions = read_csv(
            paths["owner_inputs"],
            {"input_id", "workstream", "priority", "blocks_release", "question", "candidate_value", "user_value", "status", "evidence_reference"},
            "input_id",
        )
        rules = read_csv(
            paths["closeout_rules"],
            {"input_id", "closeout_kind", "responsible_party", "required_evidence", "automatic_close_allowed"},
            "input_id",
        )
        equipment = read_csv(
            paths["equipment"],
            {"equipment_id", "item_name", "manufacturer", "model"},
            "equipment_id",
        )
        requirements = read_csv(
            paths["requirements"],
            {"requirement_id", "equipment_id", "discipline", "parameter_key", "status", "blocks_release"},
            "requirement_id",
        )
        drawings = read_csv(paths["drawings"], {"sheet_number", "title", "status"}, "sheet_number")
        drawing_ids = {row["sheet_number"] for row in drawings}
        reviewed_release = optional_release_context(paths["reviewed_release"], formal_ifc_hash)
        construction_release = optional_release_context(paths["construction_release"], formal_ifc_hash)

        decision_items = owner_review_items(decisions, rules, drawing_ids)
        requirement_packages = requirement_review_packages(requirements, equipment, drawing_ids)
        blockers = [row for row in requirements if row["blocks_release"].lower() == "yes"]
        int1 = [row for row in blockers if row["discipline"].upper() == "INT1"]
        closeout_counts = Counter(item["closeout_status"] for item in decision_items)
        review_items = [
            {"review_item_kind": "owner_input", **item}
            for item in decision_items
            if item["blocks_release"] and item["closeout_status"] not in CLOSED_CLOSEOUT_STATUSES
        ] + [
            {"review_item_kind": "requirement_package", **package}
            for package in requirement_packages
        ]
        construction_release_ready = bool(
            construction_release["current"]
            and construction_release["result"].get("construction_release_candidate_ready")
        ) and not review_items
        report = {
            "schema_version": 1,
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "generated_view": True,
            "read_only_sources": True,
            "source_ifc_sha256": formal_ifc_hash,
            "canonical_source_of_truth": [str(paths[key]) for key in ("owner_inputs", "closeout_rules", "equipment", "requirements", "drawings")],
            "formal_ifc": {
                "path": str(paths["ifc"]),
                "sha256": formal_ifc_hash,
                "write_prohibited": True,
            },
            "release_context": {
                "reviewed_candidate": reviewed_release,
                "construction_release_candidate": construction_release,
            },
            "summary": {
                "owner_input_count": len(decision_items),
                "owner_closeout_status_counts": dict(sorted(closeout_counts.items())),
                "owner_closeout_open_count": sum(item["closeout_status"] not in CLOSED_CLOSEOUT_STATUSES for item in decision_items),
                "blocking_requirement_count": len(blockers),
                "int1_requirement_count": len(int1),
                "requirement_review_package_count": len(requirement_packages),
                "int1_review_package_count": sum(package["root_cause_id"].startswith("INT1-") for package in requirement_packages),
                "unmapped_requirement_count": 0,
                "open_root_review_item_count": len(review_items),
                "unmapped_blocker_count": 0,
                "review_item_count": len(review_items),
            },
            "decision_review_items": decision_items,
            "requirement_review_packages": requirement_packages,
            "review_items": review_items,
            "gates": {
                "current_release_reports_match_formal_ifc": (
                    reviewed_release["current"] and construction_release["current"]
                ),
                "owner_rules_complete": True,
                "all_release_blocking_requirements_mapped_once": True,
                "all_blocking_requirements_mapped": True,
                "int1_requirements_mapped_to_stable_packages": True,
                "formal_ifc_write_prohibited": True,
                "construction_release_ready": construction_release_ready,
            },
        }
        atomic_write(paths["json_output"], json.dumps(report, ensure_ascii=False, indent=2) + "\n")
        atomic_write(paths["markdown_output"], markdown(report))
        print(json.dumps({
            "status": "ok",
            "json_output": str(paths["json_output"]),
            "markdown_output": str(paths["markdown_output"]),
            "summary": report["summary"],
            "gates": report["gates"],
        }, ensure_ascii=False, indent=2))
        return 0
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        print(json.dumps({"status": "fail_closed", "error": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
