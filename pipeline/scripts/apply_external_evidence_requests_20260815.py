#!/usr/bin/env python3
"""Register the 2026-08-15 external evidence request packages in project SSOT."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "pipeline/decisions"


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
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def joined(value: str, addition: str) -> str:
    return ";".join(dict.fromkeys([part for part in (value or "").split(";") if part] + [addition]))


PACKAGES = {
    "A104-SHOP-REQUEST-20260815": {
        "path": "drawings/evidence/A104-M05-M07-shop-drawing-request-20260815.md",
        "discipline": "A-104/INT1/ELEC",
        "sheet_id": "A-104/E-302",
        "decision_scope": "M05/M06 Sail 与 M07 Pivot 项目 shop drawing 索取",
        "locator": "M05/M06；M07；证据边界",
        "evidence": "官方 CAD、正式 IFC 与开发商 DXF 的机械证据边界及项目索取字段",
        "proves": "通用产品构造、当前项目身份/空间关系和不得用通用尺寸替代订单尺寸的关闭规则",
        "does_not_prove": "项目名义宽高、Sail 滑向、Pivot 左右手、宿主、安装中心或厂家签认",
    },
    "RCP1-PLUM-REQUEST-20260815": {
        "path": "drawings/evidence/RCP1-PLUM-厂家接口索取表-20260815.md",
        "discipline": "RCP1/PLUM/INT1/ELEC",
        "sheet_id": "RCP1/P-201/P-202/I-501/I-502/I-503/E-303",
        "decision_scope": "HVAC、给排水和 APP-016 厂家接口最短索取表",
        "locator": "RCP1；PLUM；APP-016",
        "evidence": "现有官方资料已关闭项与逐设备剩余厂家/现场字段的精确分流",
        "proves": "内部可确认接口已吸收，剩余阻断已转换为责任人、证据和受影响图纸",
        "does_not_prove": "缺失的厂家端口坐标、A05/A06 管径、APP-016 功率/法兰/排水或定制加工尺寸",
    },
    "E304-SITE-CHECKLIST-20260815": {
        "path": "drawings/evidence/E304-现场最短取证清单-20260815.md",
        "discipline": "ELEC",
        "sheet_id": "E-304",
        "decision_scope": "弱电箱净尺寸、温升和网线通断现场取证",
        "locator": "弱电箱与柜格净尺寸；温升试验；网线通断与 PoE",
        "evidence": "带单位、方法、照片和 blocker 对应关系的可填写现场表",
        "proves": "现场团队已有无需 Blender 的最短可执行关闭清单",
        "does_not_prove": "任何尚未填写的尺寸、温升、通断、PoE 功率或最终端口表",
    },
    "A106-CONSULTATION-PACK-20260815": {
        "path": "drawings/evidence/A106-燃气消防最终咨询包-20260815.md",
        "discipline": "ELEC/GAS/FIRE",
        "sheet_id": "A-106/E-303/I-501",
        "decision_scope": "珠海燃气公司与厨房消防设备方最终咨询包",
        "locator": "珠海燃气公司；厨房消防探测器",
        "evidence": "燃气模板附件清单、书面回填字段和消防类型/产品一次回复字段",
        "proves": "候选未升级为批准产品，厨房探测器点位与最终类型/产品责任边界分开",
        "does_not_prove": "燃气公司批准、消防设备选型、供电通信、安装净距或验收完成",
    },
}


def main() -> None:
    source_fields, sources = read_csv(DECISIONS / "source-evidence-register.csv")
    source_by_id = {row["source_id"]: row for row in sources}
    for source_id, definition in PACKAGES.items():
        relative = definition["path"]
        source_by_id[source_id] = {
            "source_id": source_id,
            "discipline": definition["discipline"],
            "sheet_id": definition["sheet_id"],
            "decision_scope": definition["decision_scope"],
            "source_kind": "project_external_evidence_request_package",
            "source_document": relative,
            "source_url": "",
            "local_path": relative,
            "sha256": sha256(ROOT / relative),
            "locator": definition["locator"],
            "evidence": definition["evidence"],
            "proves": definition["proves"],
            "does_not_prove": definition["does_not_prove"],
            "status": "verified_request_package_pending_external_response",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "",
            "model_scope": "",
            "revision": "2026-08-15",
            "publication_date": "2026-08-15",
            "legacy_targets": "",
            "legacy_projection_json": "",
            "notes": "请求包只定义最短关闭证据，不把待回复字段标为已确认。",
        }
    write_csv(DECISIONS / "source-evidence-register.csv", source_fields, list(source_by_id.values()))

    owner_fields, owners = read_csv(DECISIONS / "owner-input-register.csv")
    owner_updates = {
        "A104-M05-M06-DIMENSIONS": (
            "A104-SHOP-REQUEST-20260815",
            "已读取官方 Sail CAD 并提取通用构造；项目洞口、门板、轨道、滑向、基层和公差只接受订单或项目 shop drawing。",
        ),
        "A104-M07-EVIDENCE": (
            "A104-SHOP-REQUEST-20260815",
            "官方 Pivot 证明衣柜与墙间构造；开发商 DXF 只关闭现有门套/墙线关系。Master A 必须避开活动门扇和门套，左右手、名义宽高与净距仍待项目 shop drawing。",
        ),
        "RCP1-HVAC-PORTS": (
            "RCP1-PLUM-REQUEST-20260815",
            "A01–A04 管径/坡度已由官方型号族关闭；A05 几何已关闭但管径/端口未闭；A06 产品与全部接口未闭。索取表已按对象列出精确字段。",
        ),
        "PLUM-ROUGHINS": (
            "RCP1-PLUM-REQUEST-20260815",
            "洗碗机与净饮机通用官方接口已内部吸收；吉博力、Foster、洗烘和定制盆剩余标高、方向、中心、防水、检修及 shop drawing 已逐台分流。",
        ),
        "APP016-DISPOSER-DATA": (
            "RCP1-PLUM-REQUEST-20260815",
            "未找到可核验的勒科斯 F50 厂家官方规格或安装图；只向商家索取同机铭牌与安装图，不估算功率、法兰或排水。",
        ),
        "E304-CABINET-DIMENSIONS": (
            "E304-SITE-CHECKLIST-20260815",
            "现场清单已明确箱体/柜格 W×H×D、插头/弯线占用、单位和照片要求；待填写实测值。",
        ),
        "E304-CABINET-VENTILATION": (
            "E304-SITE-CHECKLIST-20260815",
            "现场清单已定义柜门关闭、全负载、0/60/120 min 三点温度和照片；手感温度不关闭该门。",
        ),
        "E304-CABLE-CONTINUITY": (
            "E304-SITE-CHECKLIST-20260815",
            "现场清单已定义逐根 1–8 线序、两端标签、端口与 PoE 记录；照片标签不替代测线。",
        ),
        "A106-FIRE-TYPE": (
            "A106-CONSULTATION-PACK-20260815",
            "厨房点位已确认；感温/感烟/复合、准确产品、认证、供电通信、净距和装饰件适配由消防/设备方一次回复。",
        ),
        "A106-GAS-ALARM": (
            "A106-CONSULTATION-PACK-20260815",
            "正式咨询包已就绪并引用保留模板；所有型号仍是咨询候选，不代表珠海燃气公司批准。",
        ),
    }
    for row in owners:
        update = owner_updates.get(row["input_id"])
        if update:
            evidence_id, notes = update
            row["evidence_reference"] = joined(row["evidence_reference"], evidence_id)
            row["notes"] = notes
    write_csv(DECISIONS / "owner-input-register.csv", owner_fields, owners)

    close_fields, closeouts = read_csv(DECISIONS / "owner-input-closeout-rules.csv")
    closeout_updates = {
        "A104-M05-M06-DIMENSIONS": "按 A104-SHOP-REQUEST-20260815 提交项目订单/shop drawing：洞口、门板、轨道、滑向、基层、收口和公差",
        "A104-M07-EVIDENCE": "按 A104-SHOP-REQUEST-20260815 提交 Pivot 项目 shop drawing：名义宽高、左右手、开启包络、门套/衣柜净距和 Master A 避让",
        "RCP1-HVAC-PORTS": "按 RCP1-PLUM-REQUEST-20260815 的 A01–A06 行提交准确接口图和设计参数",
        "PLUM-ROUGHINS": "按 RCP1-PLUM-REQUEST-20260815 逐台提交厂家图、现场标高/方向、防水检修和定制 shop drawing",
        "APP016-DISPOSER-DATA": "勒科斯 F50 同机铭牌照片和官方安装图，包含功率、控制/电源、法兰/开孔、机身净距和排水接口",
        "E304-CABINET-DIMENSIONS": "填写 E304-SITE-CHECKLIST-20260815 的弱电箱/柜格净尺寸并附带尺照片",
        "E304-CABINET-VENTILATION": "完成 E304-SITE-CHECKLIST-20260815 的 120 min 柜门关闭全负载温升试验",
        "E304-CABLE-CONTINUITY": "填写 E304-SITE-CHECKLIST-20260815 的逐根 1–8 线序、端点、标签和照片",
        "A106-FIRE-TYPE": "按 A106-CONSULTATION-PACK-20260815 由消防/设备方书面确认类型、型号、认证、供电通信和安装净距",
        "A106-GAS-ALARM": "发送 A106-燃气公司咨询模板并按 A106-CONSULTATION-PACK-20260815 收齐书面准入、联动、供电、净距与验收回复",
    }
    for row in closeouts:
        if row["input_id"] in closeout_updates:
            row["required_evidence"] = closeout_updates[row["input_id"]]
            row["automatic_close_allowed"] = "no"
    write_csv(DECISIONS / "owner-input-closeout-rules.csv", close_fields, closeouts)

    req_fields, requirements = read_csv(DECISIONS / "equipment-installation-requirements.csv")
    for row in requirements:
        if row["equipment_id"] == "APP-016" and row["parameter_key"] in {"rated_power", "switch_method", "flange_size", "drain_connection"}:
            row["source_id"] = "RCP1-PLUM-REQUEST-20260815"
            row["source_locator"] = "APP-016｜勒科斯 F50 垃圾处理器"
            row["notes"] = "准确值保持 unknown；只接受同机铭牌与厂家安装图，不按同类产品估算。"
    write_csv(DECISIONS / "equipment-installation-requirements.csv", req_fields, requirements)

    drawing_fields, drawings = read_csv(DECISIONS / "drawing-register.csv")
    notes = {
        "A-104": "M05/M06 官方 Sail CAD 只关闭单轨产品族；通用尺寸不写入项目。M07 官方 Pivot 与开发商 DXF 关闭衣柜—墙—既有门套空间关系，左右手、项目宽高、宿主、净距和三樘 shop drawing 仍待厂家。索取表：drawings/evidence/A104-M05-M07-shop-drawing-request-20260815.md。",
        "E-302": "Entry/Master 功能关系已关闭；Master A 只能落固定阻燃可检修见光板并避开 M07 活动门扇/门套。准确 SKU、负载/浪涌、底盒和 M07 shop drawing 仍阻断发布。",
        "E-304": "两台 AP 的 PoE 星型拓扑已关闭；弱电箱/柜格净尺寸、120 min 关门温升、逐根网线通断及 PoE 端口记录已有可直接填写的现场清单，待现场返回。",
        "A-106": "厨房火灾探测点位已确认；感温/感烟/复合与产品由消防/设备方确认。燃气公司咨询模板及最终附件/回填清单已就绪，所有报警器型号仍为咨询候选，未获批准。",
        "P-201": "已确认产品的官方介质需求继续采用；吉博力、Foster、洗烘、定制盆和 APP-016 缺项已转为逐台厂家/现场索取表，未提供的压力、中心、管径和路线保持 unknown。",
        "P-202": "现状排水保持只读；最终标高、方向、坡度、防水、检修和定制 shop drawing 按 RCP1-PLUM-厂家接口索取表关闭，不从 IFC proxy 或效果图推测。",
        "M-401": "A01–A04 官方管径/坡度已关闭；A05 仍缺准确冷媒/排水接口，A06 缺产品与全部接口；六台精确端口坐标均未批准，按厂家索取表关闭后再冻结路线。",
    }
    for row in drawings:
        if row["sheet_number"] in notes:
            row["notes"] = notes[row["sheet_number"]]
    write_csv(DECISIONS / "drawing-register.csv", drawing_fields, drawings)


if __name__ == "__main__":
    main()
