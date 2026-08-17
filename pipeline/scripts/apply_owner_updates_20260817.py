#!/usr/bin/env python3
"""Apply the 2026-08-17 kitchen, product-sheet, and custom-drain owner updates."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "pipeline/decisions"


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def upsert(rows: list[dict[str, str]], key: str, record: dict[str, str]) -> None:
    for index, row in enumerate(rows):
        if row[key] == record[key]:
            rows[index] = {field: record.get(field, "") for field in row}
            return
    rows.append(record)


def sha256(relative: str) -> str:
    digest = hashlib.sha256()
    with (ROOT / relative).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def update_sources() -> None:
    path = DECISIONS / "source-evidence-register.csv"
    fields, rows = read_csv(path)
    by_id = {row["source_id"]: row for row in rows}
    old = by_id["OWNER-WFIN-KITCHEN-PUTTY-20260816"]
    old.update({
        "status": "superseded_owner_design_scope",
        "does_not_prove": "当前中厨墙面方案；该全墙全生堂决定已由 2026-08-17 台面同材大板＋浅置物架方案取代",
        "notes": "保留为决策演变历史；不得继续投影为当前 WFIN、I-501、DET1 或 S-701 要求。",
    })
    request = by_id["RCP1-PLUM-REQUEST-20260815"]
    request.update({
        "sha256": sha256("drawings/evidence/RCP1-PLUM-厂家接口索取表-20260815.md"),
        "evidence": "逐台列出日立 HVAC、吉博力、Foster、APP-017、定制盆、厕所定制地漏和 APP-016 仍缺的接口、加工与现场证据",
        "proves": "当前可直接发送的最短厂家/现场索取范围，并明确定制地漏不自动替代已登记吉博力构件",
        "does_not_prove": "任何尚未返回的接口中心、管径、标高、路线、定制地漏房间分配或加工尺寸",
        "notes": "2026-08-17 增补 DRAIN-CUSTOM-001 逐房间 shop drawing 与吉博力分配/连接门。",
    })
    a104_request = by_id["A104-SHOP-REQUEST-20260815"]
    a104_request.update({
        "sha256": sha256("drawings/evidence/A104-M05-M07-shop-drawing-request-20260815.md"),
        "locator": "M05 格栅滑门门扇＋M06 配套顶部单轨轨道；M07；证据边界",
        "evidence": "以全屋定制可直接理解的名称说明 M05/M06 属于同一套 Rimadesio Sail MONOROTAIA 门组，并列出官方 CAD、正式 IFC 与开发商 DXF 的机械证据边界及项目索取字段",
        "notes": "2026-08-17 将抽象的‘M05/M06 对应关系’改为门扇、轨道和同一门组的白话说明；请求包仍只定义最短关闭证据。",
    })

    additions = [
        {
            "source_id": "VENDOR-WFIN-QUANSHENGTANG-KITCHEN-20260817",
            "discipline": "WFIN/DET1/INT1",
            "sheet_id": "WFIN/D-601/D-602/I-501/S-701",
            "decision_scope": "全生堂关于厨房墙面适用性的商家回复",
            "source_kind": "user_provided_vendor_chat",
            "source_document": "drawings/evidence/VENDOR-WFIN-全生堂厨房建议贴砖-20260817.md",
            "source_url": "",
            "local_path": "drawings/evidence/VENDOR-WFIN-全生堂厨房建议贴砖-20260817.md",
            "sha256": sha256("drawings/evidence/VENDOR-WFIN-全生堂厨房建议贴砖-20260817.md"),
            "locator": "业主询问‘厨房可以用吗’；商家回复‘厨房建议贴砖便于打理’",
            "evidence": "商家未推荐全生堂用于本项目中厨，并从日常打理角度明确建议厨房贴砖",
            "proves": "不得把中厨全墙全生堂写成厂家认可方案",
            "does_not_prove": "准确产品性能、证书、项目基层、质保或本项目必须采用小规格瓷砖",
            "status": "verified_user_provided_vendor_chat_scope_only",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "全生堂",
            "model_scope": "古法糯米石灰／准确产品未核验",
            "revision": "2026-08-17",
            "publication_date": "2026-08-17",
            "legacy_targets": "",
            "legacy_projection_json": "",
            "notes": "商家一般宣传不替代检测报告、准确包装或项目施工系统。",
        },
        {
            "source_id": "OWNER-WFIN-KITCHEN-SLAB-SHELF-20260817",
            "discipline": "WFIN/DET1/INT1",
            "sheet_id": "WFIN/D-601/D-602/I-501/S-701",
            "decision_scope": "中厨台面同材大板与浅置物架墙方向",
            "source_kind": "owner_confirmation",
            "source_document": "drawings/evidence/OWNER-WFIN-中厨同材大板浅置物架-20260817.md",
            "source_url": "",
            "local_path": "drawings/evidence/OWNER-WFIN-中厨同材大板浅置物架-20260817.md",
            "sha256": sha256("drawings/evidence/OWNER-WFIN-中厨同材大板浅置物架-20260817.md"),
            "locator": "业主确认灶台后背板区用厨房台面同款材质大板，远端墙面用浅置物架",
            "evidence": "灶台操作墙台面同材大板＋远端墙浅置物架的设计方向",
            "proves": "当前中厨墙面设计方向及旧全生堂方案已被取代",
            "does_not_prove": "准确墙段、大板材料与加工尺寸、置物架宽深高/承载/固定、架后背衬、防水和收口",
            "status": "confirmed_owner_design_direction_detail_pending",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "",
            "model_scope": "",
            "revision": "2026-08-17",
            "publication_date": "2026-08-17",
            "legacy_targets": "",
            "legacy_projection_json": "",
            "notes": "方向进入 I-501/WFIN/DET1；未提供的数据保持 unknown，不新增正式 IFC 材料或几何。",
        },
        {
            "source_id": "OWNER-DRAWING-MATERIALIZED-ORTHO-20260817",
            "discipline": "INT1/QA01",
            "sheet_id": "A-001/I-501/I-502/I-503/I-504",
            "decision_scope": "材质化正投影产品页交付要求",
            "source_kind": "owner_confirmation",
            "source_document": "drawings/evidence/OWNER-DRAWING-材质化正投影产品页-20260817.md",
            "source_url": "",
            "local_path": "drawings/evidence/OWNER-DRAWING-材质化正投影产品页-20260817.md",
            "sha256": sha256("drawings/evidence/OWNER-DRAWING-材质化正投影产品页-20260817.md"),
            "locator": "业主确认该效果为项目最终需完成的出图内容",
            "evidence": "最终图纸应包含平面、正立面/展开立面、材质纹理、尺寸链和统一版式的材质化正投影产品页",
            "proves": "交付内容和产品页/施工协调图双轨边界",
            "does_not_prove": "最终材料、柜体 shop drawing、设备尺寸、现场完成面尺寸或厂家加工尺寸",
            "status": "confirmed_owner_deliverable_scope_inputs_pending",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "",
            "model_scope": "",
            "revision": "2026-08-17",
            "publication_date": "2026-08-17",
            "legacy_targets": "",
            "legacy_projection_json": "",
            "notes": "可先生成带状态标识的候选页；不得用示意纹理冒充已选材料或替代加工图。",
        },
        {
            "source_id": "OWNER-PLUM-CUSTOM-DRAIN-20260817",
            "discipline": "PLUM/DET1/INT1",
            "sheet_id": "P-202/I-502/D-602/S-701",
            "decision_scope": "厕所定制地漏渠道与组合方向",
            "source_kind": "owner_confirmation_with_external_profile_reference",
            "source_document": "drawings/evidence/OWNER-PLUM-厕所定制地漏渠道-20260817.md",
            "source_url": "https://www.xiaohongshu.com/user/profile/61137f3c000000000100af65",
            "local_path": "drawings/evidence/OWNER-PLUM-厕所定制地漏渠道-20260817.md",
            "sha256": sha256("drawings/evidence/OWNER-PLUM-厕所定制地漏渠道-20260817.md"),
            "locator": "小红书主页 ID 61137f3c000000000100af65；业主指定定制组合",
            "evidence": "定制渠道已选；方向为水母地漏/中央集水器＋托克乐思网＋线性地漏排水渠",
            "proves": "业主选择的定制渠道和组件方向",
            "does_not_prove": "商家主体、准确型号、数量/房间映射、尺寸、流量、水封、排水/防水中心、完成标高或 shop drawing",
            "status": "confirmed_owner_vendor_direction_shop_drawing_pending",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "",
            "model_scope": "custom floor drain assembly",
            "revision": "2026-08-17",
            "publication_date": "2026-08-17",
            "legacy_targets": "",
            "legacy_projection_json": "",
            "notes": "只保存稳定主页，不保存临时 xsec_token；不自动替代已购 Geberit CleanLine50。",
        },
        {
            "source_id": "EXT-EXTERNAL-INFO-MINIMUM-20260817",
            "discipline": "PM/QA01/ALL",
            "sheet_id": "A-001/PM/QA01",
            "decision_scope": "当前外部信息最短清单与分发入口",
            "source_kind": "project_external_evidence_request_package",
            "source_document": "drawings/evidence/EXT-外部信息最短清单-20260817.md",
            "source_url": "",
            "local_path": "drawings/evidence/EXT-外部信息最短清单-20260817.md",
            "sha256": sha256("drawings/evidence/EXT-外部信息最短清单-20260817.md"),
            "locator": "先发出的 4 份现成文件；厂家／现场／主管部门／业主偏好；当前结论",
            "evidence": "把 33 项无保留发布阻断和 19 个稳定证据包压缩为可直接分发的四类最短清单，并链接 4 份现成关闭包",
            "proves": "当前缺失证据、责任方、受影响图纸和可执行分发入口已完整整理",
            "does_not_prove": "清单所索取的任何尚未返回的厂家、现场、主管部门或业主偏好证据",
            "status": "verified_request_package_pending_external_response",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "",
            "model_scope": "",
            "revision": "2026-08-17",
            "publication_date": "2026-08-17",
            "legacy_targets": "",
            "legacy_projection_json": "",
            "notes": "这是外部取证协调入口，不是返回证据本身；不得据此关闭任何项目尺寸或接口。",
        },
    ]
    for record in additions:
        upsert(rows, "source_id", record)
    write_csv(path, fields, rows)


def update_owner_inputs() -> None:
    path = DECISIONS / "owner-input-register.csv"
    fields, rows = read_csv(path)
    replacements = [
        {
            "input_id": "WFIN-KITCHEN-PUTTY-COVERAGE", "workstream": "WFIN/DET1/I-501", "priority": "P1", "blocks_release": "yes",
            "question": "中厨 R04 墙面采用什么材料以及是否保留瓷砖", "candidate_value": "中厨 R04 全部室内垂直墙面采用全生堂腻子；不用瓷砖",
            "user_value": "旧方案已撤回；当前改为灶台操作墙台面同材大板＋远端墙浅置物架", "unit": "", "status": "不适用",
            "evidence_reference": "VENDOR-WFIN-QUANSHENGTANG-KITCHEN-20260817;OWNER-WFIN-KITCHEN-SLAB-SHELF-20260817",
            "source_basis": "业主于 2026-08-16 纠正此前墙砖决定",
            "sync_target": "WFIN-R05;D-601;D-602;I-501;S-701", "notes": "本条为历史决定的受控替代；不得继续投影为当前中厨全墙材料。",
        },
        {
            "input_id": "WFIN-KITCHEN-PUTTY-SYSTEM", "workstream": "WFIN/DET1/I-501", "priority": "P1", "blocks_release": "yes",
            "question": "中厨全生堂腻子的准确产品、颜色表面、涂层、厚度、基层底涂、防潮耐污和收口系统", "candidate_value": "全生堂品牌与腻子类别已确认；准确产品和完整施工系统待厂家资料",
            "user_value": "旧中厨全生堂体系已撤回，不再适用当前厨房方案", "unit": "", "status": "不适用",
            "evidence_reference": "VENDOR-WFIN-QUANSHENGTANG-KITCHEN-20260817;OWNER-WFIN-KITCHEN-SLAB-SHELF-20260817",
            "source_basis": "品牌和材料类别由业主确认，但准确产品及施工体系没有证据",
            "sync_target": "WFIN-R05;D-601;D-602;I-501;S-701", "notes": "关闭原全生堂厨房系统索取项；不影响今后在其他合适空间重新选择该材料。",
        },
        {
            "input_id": "WFIN-KITCHEN-SLAB-SHELF-SCOPE", "workstream": "WFIN/DET1/I-501", "priority": "P1", "blocks_release": "no",
            "question": "中厨墙面采用什么总体设计方向", "candidate_value": "灶台操作墙台面同材大板；远端墙浅置物架",
            "user_value": "灶台后背板区用厨房台面同款材质大板；远端墙面用浅置物架", "unit": "", "status": "自定义确认",
            "evidence_reference": "OWNER-WFIN-KITCHEN-SLAB-SHELF-20260817", "source_basis": "业主于 2026-08-17 直接确认",
            "sync_target": "WFIN-R05;D-601;D-602;I-501;S-701", "notes": "只关闭设计方向；置物架不替代架后墙面完成层。",
        },
        {
            "input_id": "WFIN-KITCHEN-SLAB-SHELF-DETAIL", "workstream": "WFIN/DET1/I-501", "priority": "P1", "blocks_release": "yes",
            "question": "中厨大板和浅置物架如何形成可加工、可清洁且与防水连续的完整节点", "candidate_value": "由 I-501 分墙段并由台面/大板供应商与全屋定制联合出 shop drawing",
            "user_value": "", "unit": "mm", "status": "需证据", "evidence_reference": "OWNER-WFIN-KITCHEN-SLAB-SHELF-20260817",
            "source_basis": "方向已定，但准确墙段、材料、尺寸、承载、固定和背衬均无证据",
            "sync_target": "WFIN-R05;D-601;D-602;I-501;S-701", "notes": "不得用房间级 WFIN 分段反推加工范围或接口中心。",
        },
        {
            "input_id": "DRAWING-MATERIALIZED-ORTHO-SCOPE", "workstream": "A-001/I-501/I-502/I-503/I-504/QA01", "priority": "P2", "blocks_release": "no",
            "question": "最终是否需要材质化正投影产品页", "candidate_value": "平面＋正/展开立面＋材质纹理＋尺寸链＋统一版式",
            "user_value": "这是项目最终要做的内容", "unit": "", "status": "自定义确认",
            "evidence_reference": "OWNER-DRAWING-MATERIALIZED-ORTHO-20260817", "source_basis": "业主于 2026-08-17 确认交付要求",
            "sync_target": "A-001;I-501;I-502;I-503;I-504;QA01", "notes": "产品页与施工协调图并行，不能替代厂家 shop drawing。",
        },
        {
            "input_id": "DRAWING-MATERIALIZED-ORTHO-INPUTS", "workstream": "I-501/I-502/I-503/I-504/QA01", "priority": "P2", "blocks_release": "yes",
            "question": "最终材质化正投影产品页还缺哪些可发布输入", "candidate_value": "最终柜体 shop drawing、现场完成面、材料编号/纹理、设备尺寸及大板加工节点",
            "user_value": "", "unit": "", "status": "需证据", "evidence_reference": "OWNER-DRAWING-MATERIALIZED-ORTHO-20260817",
            "source_basis": "现有 Bonsai 平立面可作为几何底座，但最终材质和加工信息未齐",
            "sync_target": "I-501;I-502;I-503;I-504;QA01", "notes": "允许先生成带 confirmed/tentative/unknown 标识的候选页。",
        },
        {
            "input_id": "PLUM-CUSTOM-DRAIN-VENDOR", "workstream": "P-202/I-502/D-602/S-701", "priority": "P1", "blocks_release": "no",
            "question": "厕所定制地漏由谁深化以及采用什么组合方向", "candidate_value": "指定主页商家；水母地漏/中央集水器＋托克乐思网＋线性地漏排水渠",
            "user_value": "找该小红书主页商家定制；采用水母地漏/中央集水器＋托克乐思网＋线性地漏排水渠", "unit": "", "status": "自定义确认",
            "evidence_reference": "OWNER-PLUM-CUSTOM-DRAIN-20260817", "source_basis": "业主于 2026-08-17 直接确认",
            "sync_target": "P-202;I-502;D-602;S-701", "notes": "渠道和组合方向已定；商家主体及准确产品仍待订单/shop drawing。",
        },
        {
            "input_id": "PLUM-CUSTOM-DRAIN-SHOP-DRAWING", "workstream": "P-202/I-502/D-602/S-701", "priority": "P0", "blocks_release": "yes",
            "question": "厕所定制地漏如何与项目卫生间、坡面、防水和排水接口闭合", "candidate_value": "由定制商家按每个卫生间提交盖章/签认 shop drawing 与组件清单",
            "user_value": "", "unit": "mm", "status": "需证据", "evidence_reference": "OWNER-PLUM-CUSTOM-DRAIN-20260817",
            "source_basis": "现阶段只有渠道和组合方向，没有准确型号、数量、房间映射或安装接口",
            "sync_target": "P-202;I-502;D-602;S-701", "notes": "不自动替代已购 Geberit CleanLine50，不修改正式 IFC 地漏、坡面或防水几何。",
        },
    ]
    for record in replacements:
        upsert(rows, "input_id", record)
    write_csv(path, fields, rows)


def update_closeout_rules() -> None:
    path = DECISIONS / "owner-input-closeout-rules.csv"
    fields, rows = read_csv(path)
    records = [
        {"input_id": "WFIN-KITCHEN-PUTTY-COVERAGE", "closeout_kind": "human_design_selection", "responsible_party": "业主/室内设计", "required_evidence": "业主确认旧全生堂中厨方案已由台面同材大板＋浅置物架方案取代", "automatic_close_allowed": "no", "notes": "已由 OWNER-WFIN-KITCHEN-SLAB-SHELF-20260817 关闭并保留历史。"},
        {"input_id": "WFIN-KITCHEN-PUTTY-SYSTEM", "closeout_kind": "human_design_selection", "responsible_party": "业主/室内设计", "required_evidence": "确认不再把全生堂产品系统作为中厨当前材料索取项", "automatic_close_allowed": "no", "notes": "已关闭为不适用当前中厨方案；不再等待全生堂厨房系统资料。"},
        {"input_id": "WFIN-KITCHEN-SLAB-SHELF-SCOPE", "closeout_kind": "human_design_selection", "responsible_party": "业主/室内设计", "required_evidence": "业主确认灶台操作墙台面同材大板＋远端墙浅置物架", "automatic_close_allowed": "no", "notes": "设计方向已关闭，准确构造由独立深化项关闭。"},
        {"input_id": "WFIN-KITCHEN-SLAB-SHELF-DETAIL", "closeout_kind": "product_interface_evidence", "responsible_party": "室内设计/全屋定制/台面与大板供应商/现场", "required_evidence": "I-501 分墙段平立面及联合 shop drawing：大板材料厚度板幅拼缝耐热开孔收口；浅置物架位置宽深高层数承载固定基层背衬和清洁方式；防潮防水连续节点", "automatic_close_allowed": "no", "notes": "不得用常见置物架尺寸或房间级 WFIN 分段代替项目加工图。"},
        {"input_id": "DRAWING-MATERIALIZED-ORTHO-SCOPE", "closeout_kind": "human_design_selection", "responsible_party": "业主/室内设计", "required_evidence": "业主确认材质化正投影产品页为最终交付内容", "automatic_close_allowed": "no", "notes": "交付范围已关闭。"},
        {"input_id": "DRAWING-MATERIALIZED-ORTHO-INPUTS", "closeout_kind": "product_interface_evidence", "responsible_party": "全屋定制/材料供应商/设备方/现场/室内设计", "required_evidence": "最终柜体 shop drawing、现场完成面尺寸、材料编号与纹理/样板、设备外形与开启包络、大板范围厚度拼缝边缘及开孔", "automatic_close_allowed": "no", "notes": "未齐前允许输出明确标记未知项的候选产品页，不得冒充最终材料或加工图。"},
        {"input_id": "PLUM-CUSTOM-DRAIN-VENDOR", "closeout_kind": "human_design_selection", "responsible_party": "业主/给排水设计", "required_evidence": "业主确认定制主页渠道与水母地漏/中央集水器＋托克乐思网＋线性排水渠组合方向", "automatic_close_allowed": "no", "notes": "渠道和方向已关闭。"},
        {"input_id": "PLUM-CUSTOM-DRAIN-SHOP-DRAWING", "closeout_kind": "product_interface_evidence", "responsible_party": "定制地漏商家/给排水设计/防水/现场", "required_evidence": "商家主体与订单；逐卫生间数量/房间映射/组件清单；准确型号材质表面；总长宽深、排水口径方向、水封、流量、找平层/完成面范围、防水法兰、坡向、检修清洁、排水中心和安装公差的 shop drawing", "automatic_close_allowed": "no", "notes": "shop drawing 未签认前不替代 Geberit、不改正式 IFC 地漏与坡面。"},
    ]
    for record in records:
        upsert(rows, "input_id", record)
    write_csv(path, fields, rows)


def update_equipment() -> None:
    path = DECISIONS / "equipment-register.csv"
    fields, rows = read_csv(path)
    record = {
        "equipment_id": "DRAIN-CUSTOM-001", "domain": "DRAINAGE", "category": "custom_floor_drain_assembly",
        "item_name": "厕所定制水母地漏/中央集水器＋托克乐思网＋线性排水渠组合", "manufacturer": "", "model": "",
        "variant": "owner_selected_custom_vendor_direction", "quantity": "", "procurement_status": "selected", "decision_status": "partial",
        "storage_location_candidate": "", "use_location_candidate": "卫生间；具体房间和数量待 shop drawing", "storage_location_confirmed": "", "use_location_confirmed": "",
        "schedule_included": "yes", "selector_kind": "logical_input", "selector_value": "DRAIN-CUSTOM-001", "ifc_class": "", "ifc_type_name": "",
        "ifc_type_global_id": "", "ifc_global_ids": "", "source_ids": "OWNER-PLUM-CUSTOM-DRAIN-20260817",
        "identity_basis": "业主确认小红书主页定制渠道与组件组合方向；公开检索不能替代订单或产品证据", "confidence": "1.00",
        "human_review_required": "yes", "legacy_kind": "owner_custom_drain", "legacy_id": "DRAIN-CUSTOM-001",
        "notes": "渠道和方向已定，非已采购/已加工；房间、数量、品牌型号、尺寸、接口和防水均待 shop drawing；不自动替代已购 Geberit CleanLine50。",
    }
    upsert(rows, "equipment_id", record)
    write_csv(path, fields, rows)


def update_requirements() -> None:
    path = DECISIONS / "equipment-installation-requirements.csv"
    fields, rows = read_csv(path)
    confirmed = [
        ("REQ-DRAINCUST-001", "vendor_profile", "https://www.xiaohongshu.com/user/profile/61137f3c000000000100af65", "稳定主页 URL；账号主体待订单核验"),
        ("REQ-DRAINCUST-002", "assembly_direction", "水母地漏/中央集水器＋托克乐思网＋线性地漏排水渠", "业主确认组合方向；不是准确产品型号"),
    ]
    for req_id, key, value, note in confirmed:
        upsert(rows, "requirement_id", {
            "requirement_id": req_id, "equipment_id": "DRAIN-CUSTOM-001", "discipline": "PLUM/DET1/INT1", "parameter_key": key,
            "value_text": value, "value_number": "", "unit": "", "datum": "", "value_origin": "user_input", "status": "confirmed",
            "source_id": "OWNER-PLUM-CUSTOM-DRAIN-20260817", "source_locator": "业主 2026-08-17 确认", "blocks_release": "no", "notes": note,
        })
    pending = [
        ("REQ-DRAINCUST-003", "component_models", "关闭证据：逐组件准确品牌/型号/材质/表面和订单组件表。"),
        ("REQ-DRAINCUST-004", "quantity", "关闭证据：逐卫生间数量表。"),
        ("REQ-DRAINCUST-005", "room_mapping", "关闭证据：主卫/客卫等逐空间映射。"),
        ("REQ-DRAINCUST-006", "shop_drawing", "关闭证据：商家签认项目 shop drawing。"),
        ("REQ-DRAINCUST-007", "overall_dimensions", "关闭证据：总长、宽、深及组件包络。"),
        ("REQ-DRAINCUST-008", "outlet_nominal_diameter", "关闭证据：准确排水口径和连接制式。"),
        ("REQ-DRAINCUST-009", "outlet_orientation", "关闭证据：出水方向与现场连接剖面。"),
        ("REQ-DRAINCUST-010", "water_seal_depth", "关闭证据：同型号水封构造与深度。"),
        ("REQ-DRAINCUST-011", "design_flow_rate", "关闭证据：厂家设计/测试流量及适用边界。"),
        ("REQ-DRAINCUST-012", "waterproof_flange_and_build_up", "关闭证据：防水法兰/翼环、找平层与完成面适用范围。"),
        ("REQ-DRAINCUST-013", "finish_level_and_slope_interface", "关闭证据：完成标高、坡向、拼接与公差。"),
        ("REQ-DRAINCUST-014", "service_and_cleaning_access", "关闭证据：毛发网取出、清洁、检修和更换方式。"),
        ("REQ-DRAINCUST-015", "interface_center_coordinates", "关闭证据：厂家安装图与项目定位；不得推测中心坐标。"),
    ]
    for req_id, key, note in pending:
        upsert(rows, "requirement_id", {
            "requirement_id": req_id, "equipment_id": "DRAIN-CUSTOM-001", "discipline": "PLUM/DET1/INT1", "parameter_key": key,
            "value_text": "unknown", "value_number": "", "unit": "", "datum": "", "value_origin": "pending", "status": "pending",
            "source_id": "OWNER-PLUM-CUSTOM-DRAIN-20260817", "source_locator": "当前仅确认渠道与组合方向", "blocks_release": "yes", "notes": note,
        })
    write_csv(path, fields, rows)


def update_drawing_register() -> None:
    path = DECISIONS / "drawing-register.csv"
    fields, rows = read_csv(path)
    by_id = {row["sheet_number"]: row for row in rows}
    by_id["A-001"]["notes"] = "M072 自动汇总全部图号、阶段、版本、说明和图例；最终交付增加材质化正投影产品页要求：平面＋正/展开立面＋材质纹理＋尺寸链＋统一版式，与施工协调图并行且不替代厂家加工图"
    by_id["I-501"]["notes"] = "70 个既有对象协调包络入图；APP-009/010 两台未采购 SJ85ZX26MC 候选的官方柜孔、G3/4 冷水与 Ø38 排水约束已作为 2 条无定位接口表项入图，不生成粗装 XYZ、阀门、软管路径或开孔位置；EL-03 视图 06～09 及 EL-P01/P02 公共空间折线展开和贯穿长立面已写入正式 IFC 的 Bonsai 原生 Drawing；新增 1:30 无纹理 SVG 与手机 PNG；其余厂家安装图、燃气实测和五金运动包络未关闭；中厨墙面更新为灶台操作墙台面同材大板＋远端墙浅置物架，准确墙段、材料、加工和架后背衬待深化；最终需输出材质化正投影厨房产品页"
    by_id["I-502"]["notes"] = "33 个既有对象协调包络入图；EL-06 视图 18～22 与 EL-08 视图 27～32 已写入正式 IFC；原 9 个深层 BRep、公共空间新增 2 个对象及全量重编发现的 BED02 现均有纯折线 ELEVATION_VIEW，详细 Body 未替换；洁具粗装图、五金运动包络和节点未关闭；厕所定制地漏渠道与水母地漏/中央集水器＋托克乐思网＋线性排水渠方向已确认，准确房间、数量和 shop drawing 待商家；最终需输出材质化正投影卫生间产品页"
    by_id["S-701"]["notes"] = "53 条家具产品身份、家电使用/存放、门窗五金、既有暖通类型、定制地漏方向和墙面材料系统记录按已确认/候选/未决分栏；中厨采用灶台操作墙台面同材大板＋远端墙浅置物架方向；厕所定制地漏渠道已定但 shop drawing 未齐；不是下单表"
    write_csv(path, fields, rows)


def main() -> None:
    update_sources()
    update_owner_inputs()
    update_closeout_rules()
    update_equipment()
    update_requirements()
    update_drawing_register()
    print("Applied 2026-08-17 kitchen, product-sheet, and custom-drain owner updates.")


if __name__ == "__main__":
    main()
