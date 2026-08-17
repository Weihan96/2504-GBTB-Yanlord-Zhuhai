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
        "notes": "2026-08-17 将抽象的‘M05/M06 对应关系’改为门扇、轨道和同一门组的白话说明；MS 已在 01 表确认项目方向并纠正 M05 为向东／图纸右侧滑、M07 为图纸右侧停靠；厂家签认改由 06 表及项目 shop drawing 关闭。",
    })
    outbound_forms = {
        "OUTBOUND-FORM-DOOR-20260817": {
            "filename": "01-门组做法确认表-发全屋定制.md",
            "discipline": "A-104/INT1/E-302",
            "sheet_id": "A-104/E-302",
            "decision_scope": "MS 门组设计决定回复",
            "source_kind": "owner_completed_confirmation_markdown",
            "locator": "门组确认 D01-D13；回复人 MS；回复日期 2026-08-17",
            "evidence": "MS 确认 M05/M06 系统、第二种官方 DWG 配置、白橡木／浅色木饰面、2000 mm 顶轨、向东／图纸右侧滑开；确认 M07 为 950×2400 mm Poliform Pivot、与衣柜背面新增护墙板齐平、向主卧内开并停靠图纸右侧；确认 Master A 首选安装面及后备位置",
            "proves": "业主设计方向与优先安装关系已经确认，可进入项目 SSOT 和候选图纸",
            "does_not_prove": "Rimadesio／Poliform 或全屋定制已经签认、准确下单尺寸、原厂饰面编号、下部防摆／定位构件、安装基层、现场净距或正式 IFC 已修改",
            "status": "verified_owner_response_vendor_followup_pending",
            "notes": "保留 MS 原始回答；原表中的向左／西文字示意被 D06 回答明确纠正，当前值以结构化 SSOT 和 06 厂家复核表为准。",
        },
        "OUTBOUND-FORM-HVAC-20260817": {
            "filename": "02-日立空调接口确认表-发空调厂家.md",
            "locator": "Markdown 日立空调 H01-H09、PC-P1HEQ P01-P07 及附件依据",
            "evidence": "可直接在 Codex 中编辑和回填的 Markdown；把已知型号与仍需厂家确认的接口、兼容性和逐点控制映射分开",
            "proves": "HVAC 对外问询已形成封闭式确认项，且未虚构接口中心、最终路线或控制映射",
            "does_not_prove": "厂家已经确认、现场安装条件已经复核或正式 IFC 已修改",
            "notes": "2026-08-17 改为唯一 Markdown 可编辑源；原 XLSX 退出正式输出目录。",
        },
        "OUTBOUND-FORM-PLUM-20260817": {
            "filename": "03-给排水设备确认表-发设备与施工方.md",
            "locator": "Markdown 西门子洗烘、吉博力与卫浴、厨房设备、VVD、定制盆与地漏及附件依据",
            "evidence": "可直接在 Codex 中编辑和回填的 Markdown；逐台区分已确认产品条件与仍缺 shop drawing、粗装接口和现场核对",
            "proves": "给排水和厨房设备问询已收敛为封闭式确认项，且定制产品继续保留 shop drawing 门",
            "does_not_prove": "厂家已经确认、定制加工尺寸已冻结、现场粗装条件已复核或正式 IFC 已修改",
            "notes": "2026-08-17 改为唯一 Markdown 可编辑源；原 XLSX 退出正式输出目录。",
        },
        "OUTBOUND-FORM-GASFIRE-20260817": {
            "filename": "04-燃气消防确认表-发主管单位.md",
            "locator": "Markdown 燃气 R01-R08、消防 F01-F07 及附件依据",
            "evidence": "可直接在 Codex 中编辑和回填的 Markdown；候选设备、主管部门准入、探测器类型和联动边界分别列示",
            "proves": "燃气和消防咨询已形成可直接签认的封闭式问题，候选报警器未升级为已批准",
            "does_not_prove": "燃气公司、消防或设备方已经批准任何候选产品、点位或联动方式",
            "notes": "2026-08-17 改为唯一 Markdown 可编辑源；原 XLSX 退出正式输出目录。",
        },
        "OUTBOUND-FORM-SITE-20260817": {
            "filename": "05-弱电现场记录表-发现场负责人.md",
            "locator": "Markdown 弱电箱、温升、网线、门口设备、PC-P1HEQ 和快递现场记录",
            "evidence": "可直接在 Codex 中编辑和回填的 Markdown；狄耐克品牌已预填，只要求现场补准确型号、端子线缆、物业兼容与保留／迁移／接入信息",
            "proves": "弱电与门口设备问询已删除可由现场照片回答的对讲品牌问题",
            "does_not_prove": "狄耐克室内机准确型号、端子接法、物业系统兼容、门铃身份或外部已经回复",
            "notes": "2026-08-17 改为唯一 Markdown 可编辑源；按现场照片预填 DNAKE／狄耐克，IP/MAC 不进入对外表。",
        },
        "OUTBOUND-FORM-DOOR-VENDOR-20260817": {
            "filename": "06-门组厂家复核表-发Rimadesio与Poliform.md",
            "discipline": "A-104/INT1",
            "sheet_id": "A-104/I-504/S-701",
            "decision_scope": "Rimadesio Sail 与 Poliform Pivot 厂家复核",
            "locator": "Rimadesio R01-R06；Poliform P01-P05",
            "evidence": "把 MS 已确认的设计方向转换为厂家只需回答是／否的复核项，并要求回传项目加工图",
            "proves": "门组厂家问询已与业主设计决定分开，且当前东向滑开、右侧停靠和护墙板齐平关系已进入问询",
            "does_not_prove": "厂家已经确认、项目加工图已经返回、准确下单尺寸和安装接口已经冻结",
            "notes": "发送时随附 Rimadesio Sail 官方资料、MONOROTAIA DWG、Poliform Architectural PDF 和 A-104 候选图。",
        },
        "OUTBOUND-FORM-SMART-PANEL-20260817": {
            "filename": "07-智能面板电气接口确认表-发JINK与电气方.md",
            "discipline": "ELEC/NETWORK/INT1",
            "sheet_id": "E-302/E-304/I-503",
            "decision_scope": "JINK 智能面板产品与电气接口复核",
            "locator": "J01-J10",
            "evidence": "逐项询问 Entry A、Master A、Master B 的准确 SKU、直接负载／场景语义、负载与浪涌、接线、Matter 和见光板安装条件",
            "proves": "智能面板剩余产品与电气接口问题已有可直接发送的封闭式表格",
            "does_not_prove": "JINK 已经回复、准确 SKU 已选择、各回路负载已验算或现场底盒净距已复核",
            "notes": "不重复询问三块面板的角色；只关闭产品、接线、负载和安装兼容。",
        },
        "OUTBOUND-FORM-JOINERY-DETAIL-20260817": {
            "filename": "08-定制家具与墙脚节点确认表-发全屋定制.md",
            "discipline": "INT1/DET1/WFIN",
            "sheet_id": "I-502/I-504/D-601/S-701",
            "decision_scope": "客卫服务塔与干区墙脚节点复核",
            "locator": "客卫 F01-F04；干区墙脚 B01-B05",
            "evidence": "把客卫纸巾／垃圾服务塔及干区齐平宽踢脚＋阴影缝转成全屋定制可签认的节点问题",
            "proves": "此前没有对外关闭路径的定制家具和墙脚节点已有直接发送表格",
            "does_not_prove": "全屋定制已经回复、节点尺寸已冻结或实物样板已经批准",
            "notes": "湿区不直接套用干区木作节点；准确尺寸和材料仍由节点图及样板关闭。",
        },
    }
    for source_id, metadata in outbound_forms.items():
        relative_path = f"output/forms/对外确认表/{metadata['filename']}"
        row = by_id.get(source_id, {field: "" for field in fields})
        row.update({
            "source_id": source_id,
            "discipline": metadata.get("discipline", row.get("discipline", "")),
            "sheet_id": metadata.get("sheet_id", row.get("sheet_id", "")),
            "decision_scope": metadata.get("decision_scope", row.get("decision_scope", "")),
            "source_kind": metadata.get("source_kind", "project_external_confirmation_markdown"),
            "source_document": metadata["filename"],
            "local_path": relative_path,
            "sha256": sha256(relative_path),
            "locator": metadata["locator"],
            "evidence": metadata["evidence"],
            "proves": metadata["proves"],
            "does_not_prove": metadata["does_not_prove"],
            "status": metadata.get("status", "verified_confirmation_markdown_pending_external_response"),
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "revision": "2026-08-17",
            "publication_date": "2026-08-17",
            "notes": metadata["notes"],
        })
        by_id[source_id] = row
        upsert(rows, "source_id", row)
    site_checklist = by_id["E304-SITE-CHECKLIST-20260815"]
    site_checklist.update({
        "sha256": sha256("drawings/evidence/E304-现场最短取证清单-20260815.md"),
        "evidence": "六张现场表的测量、单位、照片和回传规则；狄耐克对讲品牌已关闭，只补准确型号、端子线缆和物业系统证据",
        "does_not_prove": "狄耐克准确型号、物业兼容、门铃身份、外部已经回复或项目值已经写入正式 IFC",
        "notes": "2026-08-17 删除重复询问对讲品牌；IP/MAC 不进入对外表。",
    })
    by_id["GAS-CONSULTATION-TEMPLATE-001"].update({
        "sha256": sha256("drawings/evidence/A106-燃气公司咨询模板.md"),
        "notes": "2026-08-17 将可填写表引用改为唯一 Markdown 源；其余准入边界不变。",
    })
    by_id["A106-CONSULTATION-PACK-20260815"].update({
        "sha256": sha256("drawings/evidence/A106-燃气消防最终咨询包-20260815.md"),
        "notes": "2026-08-17 将可填写表引用改为唯一 Markdown 源；候选产品仍须主管单位签认。",
    })

    additions = [
        {
            "source_id": "OWNER-WFIN-BASEBOARD-SHADOW-GAP-20260817",
            "discipline": "WFIN/DET1/INT1",
            "sheet_id": "WFIN/D-601/I-502/I-503/I-504/S-701",
            "decision_scope": "干区齐平宽踢脚与上下阴影缝候选",
            "source_kind": "owner_shared_design_reference",
            "source_document": "drawings/evidence/OWNER-WFIN-干区齐平踢脚阴影缝候选-20260817.md",
            "source_url": "https://www.xiaohongshu.com/explore/68773092000000001202076c",
            "local_path": "drawings/evidence/OWNER-WFIN-干区齐平踢脚阴影缝候选-20260817.md",
            "sha256": sha256("drawings/evidence/OWNER-WFIN-干区齐平踢脚阴影缝候选-20260817.md"),
            "locator": "齐平墙面同色宽踢脚、上下约 10 mm 阴影缝、墙下口型材和通缝关系",
            "evidence": "业主要求记录该做法；参考案例说明齐平宽踢脚、上下阴影缝和墙门柜通缝的设计方向及施工前置条件",
            "proves": "可作为卧室、走廊等干区的候选方向，并须与墙面找平、门套、隐形门和柜体统一深化",
            "does_not_prove": "本项目最终采用、准确高度厚度材料颜色、阴影缝公差、湿区适用性或已完成实物样板",
            "status": "owner_reference_candidate_not_final",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "",
            "model_scope": "",
            "revision": "2026-08-17",
            "publication_date": "2025-08-05",
            "legacy_targets": "",
            "legacy_projection_json": "",
            "notes": "中厨 VVD Pewter 60 mm 踢脚保持独立；湿区暂不采用上下双缝。",
        },
        {
            "source_id": "OWNER-E304-DNAKE-PHOTO-20260817",
            "discipline": "ELEC/INT1/NETWORK",
            "sheet_id": "E-304/I-503",
            "decision_scope": "既有可视对讲室内机品牌与版本信息",
            "source_kind": "owner_provided_site_photo",
            "source_document": "E304-DNAKE-indoor-monitor-version-photo-20260817.jpg",
            "source_url": "",
            "local_path": "drawings/evidence/E304-DNAKE-indoor-monitor-version-photo-20260817.jpg",
            "sha256": sha256("drawings/evidence/E304-DNAKE-indoor-monitor-version-photo-20260817.jpg"),
            "locator": "设备正面 DNAKE／狄耐克标识与版本信息页",
            "evidence": "既有可视对讲室内机品牌为 DNAKE／狄耐克；系统版本 1.6.0 20210615；应用版本 1.1.0 20210615 16M",
            "proves": "既有室内机的品牌和照片拍摄时可见的软件版本",
            "does_not_prove": "准确型号、背部端子、既有线缆接法、物业系统兼容性、保留迁移接入方案或最终安装坐标",
            "status": "verified_brand_and_software_version_only",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "DNAKE／狄耐克",
            "model_scope": "unknown",
            "revision": "2026-08-17",
            "publication_date": "2026-08-17",
            "legacy_targets": "",
            "legacy_projection_json": "",
            "notes": "原始照片中的网络地址不投影到 SSOT、对外表或图纸；准确型号继续保持 unknown。",
        },
        {
            "source_id": "OWNER-E304-DNAKE-NOTE-20260817",
            "discipline": "ELEC/INT1/NETWORK",
            "sheet_id": "E-304/I-503",
            "decision_scope": "狄耐克既有可视对讲证据边界说明",
            "source_kind": "project_evidence_boundary_note",
            "source_document": "OWNER-E304-狄耐克可视对讲室内机-20260817.md",
            "source_url": "",
            "local_path": "drawings/evidence/OWNER-E304-狄耐克可视对讲室内机-20260817.md",
            "sha256": sha256("drawings/evidence/OWNER-E304-狄耐克可视对讲室内机-20260817.md"),
            "locator": "品牌、可见软件版本、未证明事项和最短后续取证",
            "evidence": "把现场照片可确认与不可确认的信息拆开，并明确网络地址不进入 SSOT、对外表或图纸",
            "proves": "本次狄耐克现场照片的受控投影边界和剩余取证范围",
            "does_not_prove": "准确型号、端子接法、物业系统兼容、保留迁移接入结论或最终坐标",
            "status": "verified_project_evidence_boundary_note",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "DNAKE／狄耐克",
            "model_scope": "unknown",
            "revision": "2026-08-17",
            "publication_date": "2026-08-17",
            "legacy_targets": "",
            "legacy_projection_json": "",
            "notes": "原始现场照片由 OWNER-E304-DNAKE-PHOTO-20260817 单独登记。",
        },
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
            "input_id": "A104-M05-M06-DIMENSIONS", "workstream": "A-104/INT1", "priority": "P0", "blocks_release": "yes",
            "question": "请 Rimadesio／全屋定制确认项目方已定的 Sail 门组做法并回传项目加工图",
            "candidate_value": "M05/M06 为同一套 Sail MONOROTAIA 单轨单扇；官方 DWG 第二种配置；1000×2400 mm 白橡木／浅色木门扇；2000 mm 暗藏顶轨；无通长地轨；向东／图纸右侧滑开",
            "user_value": "MS 于 2026-08-17 确认：采用官方 DWG 第二种配置，顶部轨道嵌入天花，门垛墙与门平行，白橡木／浅色木饰面；M05 向东／图纸右侧滑开；轨道末端保留检修条件。", "unit": "", "status": "自定义确认",
            "evidence_reference": "OUTBOUND-FORM-DOOR-20260817;RIMADESIO-SAIL-PRODUCT-PDF-001;RIMADESIO-SAIL-MONOROTAIA-001;OUTBOUND-FORM-DOOR-VENDOR-20260817",
            "source_basis": "01 表记录业主设计决定；官方资料证明产品族和通用构造；准确下单尺寸、饰面编号和安装接口仍由厂家加工图关闭。",
            "sync_target": "A104 M05-M06 dimensions;I-504;S-701;P0 IDS",
            "notes": "D03 所述‘2 道水平横档’未获业主确认，不进入当前设计要求。向东／图纸右侧是对原向左／西候选的明确纠正；未修改正式 IFC。",
        },
        {
            "input_id": "A104-M07-EVIDENCE", "workstream": "A-104/E-302/INT1", "priority": "P0", "blocks_release": "yes",
            "question": "请 Poliform／全屋定制确认项目方已定的 Pivot 门组做法并回传项目加工图",
            "candidate_value": "Poliform Pivot 墙装单扇，950×2400 mm；关闭时与 Senzafine 背面新增护墙板齐平；向主卧内开并停靠图纸右侧固定墙／次卧侧",
            "user_value": "MS 于 2026-08-17 确认：采用 Poliform Pivot 墙装单扇及官方门框、顶／地轴、金属门扇框系统；950×2400 mm；与衣柜背面新增护墙板齐平；向主卧内开，停靠图纸右侧固定墙／次卧侧。", "unit": "", "status": "自定义确认",
            "evidence_reference": "OUTBOUND-FORM-DOOR-20260817;POLIFORM-ARCHITECTURAL-PDF-001;OUTBOUND-FORM-DOOR-VENDOR-20260817",
            "source_basis": "01 表记录业主设计决定；Poliform 官方资料证明 Pivot 系统范围；准确下单尺寸、轴位、基层、净距和收口仍由项目加工图关闭。",
            "sync_target": "A104-R03;E302-MASTER-SIDE;I-504;S-701;P0 IDS",
            "notes": "与衣柜正面齐平和停靠图纸左侧的旧候选已被明确纠正；产品候选不得写成已采购。",
        },
        {
            "input_id": "E302-MASTER-SIDE", "workstream": "E-302/A-104/INT1", "priority": "P0", "blocks_release": "yes",
            "question": "请现场复核 Master A 首选墙面及后备安装面的实际净宽、门套／见光板构造与底盒净深",
            "candidate_value": "首选主卧门与主卫门同时打开后两门之间的剩余固定墙面；如无可用固定墙面，后备为左门套／见光板侧",
            "user_value": "MS 于 2026-08-17 确认：Master A 优先安装在主卧门和主卫门都打开后两者之间的剩余墙面；若没有剩余墙面，则安装在左门套／见光板侧。", "unit": "mm", "status": "自定义确认",
            "evidence_reference": "OUTBOUND-FORM-DOOR-20260817;OUTBOUND-FORM-SITE-20260817",
            "source_basis": "01 表关闭业主位置偏好；候选图和现场记录表负责验证门扇全开时的实际可安装面。",
            "sync_target": "E302 master control;A104-R03;I-504",
            "notes": "不安装在 Pivot 活动门扇、门套五金或开启包络内；接口中心坐标保持 unknown，直到现场完成面测量。",
        },
        {
            "input_id": "WFIN-DRY-BASEBOARD-SHADOW-GAP-SCOPE", "workstream": "WFIN/DET1/INT1", "priority": "P1", "blocks_release": "no",
            "question": "卧室、走廊等干区踢脚是否采用齐平宽踢脚＋上下阴影缝方向", "candidate_value": "墙面同色齐平宽踢脚；参考案例上下各约 10 mm 阴影缝；门套、隐形门和柜体边界连续通缝",
            "user_value": "记录为干区候选，暂不最终冻结", "unit": "", "status": "采用候选",
            "evidence_reference": "OWNER-WFIN-BASEBOARD-SHADOW-GAP-20260817", "source_basis": "业主于 2026-08-17 分享参考并要求记录；参考正文给出设计方向与施工前置条件",
            "sync_target": "WFIN-R01;D-601;I-502;I-503;I-504;S-701", "notes": "仅适用于干区候选；中厨 VVD Pewter 60 mm 踢脚保持独立，湿区暂不采用上下双缝。",
        },
        {
            "input_id": "WFIN-DRY-BASEBOARD-SHADOW-GAP-DETAIL", "workstream": "WFIN/DET1/INT1", "priority": "P1", "blocks_release": "yes",
            "question": "干区齐平宽踢脚＋上下阴影缝如何形成可施工、可清洁并与门柜连续的节点", "candidate_value": "由 D-601 输出 1:5 墙脚节点与门套／隐形门／柜体通缝展开，并经实物样板确认",
            "user_value": "", "unit": "mm", "status": "需证据", "evidence_reference": "OWNER-WFIN-BASEBOARD-SHADOW-GAP-20260817",
            "source_basis": "参考案例没有提供本项目高度、厚度、材料、公差或湿区适用证据",
            "sync_target": "WFIN-R01;D-601;I-502;I-503;I-504;S-701", "notes": "需关闭分房间适用范围、踢脚高度厚度材料、上下缝宽与公差、基层型材、转角、门柜通缝、耐撞拖地水和积灰清洁。",
        },
        {
            "input_id": "E304-VIDEO-INTERCOM-DOORBELL", "workstream": "E-304", "priority": "P0", "blocks_release": "yes",
            "question": "请确认现有狄耐克可视对讲室内机的准确型号、背部端子与线缆、物业系统兼容性，以及保留／迁移／接入方案；门铃另行核对",
            "candidate_value": "可视对讲室内机品牌已确认为 DNAKE／狄耐克；开发商参考点按底边 1400 mm 协调；准确型号和系统接口待确认；门铃参考点按底边 1300 mm 协调",
            "user_value": "门禁面板是狄耐克的，具体型号没看到", "unit": "", "status": "需证据",
            "evidence_reference": "OWNER-E304-DNAKE-PHOTO-20260817;build/elec/elec-developer-control-reference.json",
            "source_basis": "现场照片确认 DNAKE／狄耐克品牌及软件版本；开发商图证明既有对讲和门铃参考点",
            "sync_target": "E-304;I-503", "notes": "不再询问可视对讲品牌；准确型号、端子、物业接口、保留迁移接入和最终坐标仍待现场／物业证据。照片中的 IP/MAC 不进入对外资料。",
        },
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
        {"input_id": "A104-M05-M06-DIMENSIONS", "closeout_kind": "vendor_shop_drawing_and_compatibility", "responsible_party": "Rimadesio供货安装方/全屋定制", "required_evidence": "在 06-门组厂家复核表 R01-R06 逐项回答是／否，并回传标有实际下单尺寸、原厂饰面编号、2000 mm 顶轨、向东／图纸右侧滑向、下部防摆／定位构件、吊顶固定、收口和检修路径的项目加工图。", "automatic_close_allowed": "no", "notes": "01 表只证明 MS 的设计决定，不证明厂家兼容或项目加工图；不再询问 2 道水平横档。"},
        {"input_id": "A104-M07-EVIDENCE", "closeout_kind": "vendor_shop_drawing_and_compatibility", "responsible_party": "Poliform供货安装方/全屋定制", "required_evidence": "在 06-门组厂家复核表 P01-P05 逐项回答是／否，并回传标有 950×2400 mm 名义尺寸、顶／地轴、与衣柜背面新增护墙板齐平、向主卧内开、图纸右侧停靠、基层、净距和收口的项目加工图。", "automatic_close_allowed": "no", "notes": "01 表只证明 MS 的设计决定；产品是否支持及准确下单尺寸仍由厂家签认。"},
        {"input_id": "E302-MASTER-SIDE", "closeout_kind": "site_or_drawing_evidence", "responsible_party": "建筑设计/全屋定制/现场", "required_evidence": "门扇全开状态下，首选两门之间固定墙面的可安装净宽、底边标高和底盒净深；如首选无可用墙面，再测左门套／见光板侧；同时证明不落入 M07 活动门扇和五金范围。", "automatic_close_allowed": "no", "notes": "业主位置优先级已确认；只剩现场几何和基层构造，不再问业主选择墙侧。"},
        {"input_id": "WFIN-DRY-BASEBOARD-SHADOW-GAP-SCOPE", "closeout_kind": "human_design_selection", "responsible_party": "业主/室内设计", "required_evidence": "业主确认该参考只作为卧室、走廊等干区候选，并明确中厨和湿区不套用", "automatic_close_allowed": "no", "notes": "候选方向已记录，最终采用仍由节点和样板关闭。"},
        {"input_id": "WFIN-DRY-BASEBOARD-SHADOW-GAP-DETAIL", "closeout_kind": "material_and_mockup_evidence", "responsible_party": "室内设计/墙面施工/全屋定制/现场", "required_evidence": "D-601 1:5 节点；分房间适用表；踢脚高度厚度材料和颜色；上下阴影缝宽与公差；墙下口型材；门套/隐形门/柜体/转角通缝展开；耐撞、拖地水和积灰清洁实物样板", "automatic_close_allowed": "no", "notes": "不得按小红书图片直接下单或写正式 IFC。"},
        {"input_id": "E304-VIDEO-INTERCOM-DOORBELL", "closeout_kind": "site_and_system_evidence", "responsible_party": "物业/门禁维护方/弱电设计/现场", "required_evidence": "狄耐克室内机背面或工程信息页的准确型号；背部端子、线缆和接法照片；物业系统兼容与移动限制书面确认；保留／迁移／接入结论；门铃实物身份及门套完成面后的最终定位", "automatic_close_allowed": "no", "notes": "DNAKE／狄耐克品牌和软件版本已由现场照片关闭，不再重复询问；开发商图仍只证明既有参考点和标注高度。"},
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
    by_id = {row["equipment_id"]: row for row in rows}
    by_id["DW-M05"].update({
        "variant": "Sail MONOROTAIA single-track single-leaf; official DWG second configuration; white oak/light wood owner direction",
        "decision_status": "partial",
        "source_ids": "IFC-FORMAL-001;RIMADESIO-SAIL-PRODUCT-PDF-001;RIMADESIO-SAIL-MONOROTAIA-001;OUTBOUND-FORM-DOOR-20260817;OUTBOUND-FORM-DOOR-VENDOR-20260817",
        "identity_basis": "正式 IFC 证明既有对象与 1000×38×2400 mm 几何包络；官方资料证明 Sail 产品族；MS 确认第二种 DWG 配置、白橡木／浅色木和向东／图纸右侧滑向",
        "notes": "业主设计方向已确认；准确下单尺寸、原厂饰面编号、吊顶固定、下部防摆／定位构件和收口仍待厂家项目加工图。D03 的 2 道水平横档不再作为项目要求。",
    })
    by_id["DW-M06"].update({
        "variant": "Sail MONOROTAIA concealed ceiling top track for DW-M05",
        "decision_status": "partial",
        "source_ids": "IFC-FORMAL-001;RIMADESIO-SAIL-MONOROTAIA-001;OUTBOUND-FORM-DOOR-20260817;OUTBOUND-FORM-DOOR-VENDOR-20260817",
        "identity_basis": "正式 IFC 证明既有轨道对象及早期 1200 mm 包络；MS 确认采用 2000 mm 暗藏顶轨、无通长地轨并保留检修条件",
        "notes": "M06 是 M05 同一门组的顶轨，不是第二樘门；2000 mm 为项目名义值，最终长度、吊顶固定和下部防摆／定位构件仍待厂家项目加工图，不写正式 IFC。",
    })
    by_id["DW-M07"].update({
        "manufacturer": "Poliform", "model": "Pivot",
        "variant": "wall-mounted single Pivot integrated with Senzafine rear added wall panel",
        "decision_status": "partial",
        "source_ids": "IFC-FORMAL-001;POLIFORM-ARCHITECTURAL-PDF-001;OUTBOUND-FORM-DOOR-20260817;OUTBOUND-FORM-DOOR-VENDOR-20260817",
        "identity_basis": "正式 IFC 证明既有门对象与空间位置；Poliform 官方资料证明 Pivot 系统；MS 确认 950×2400 mm、衣柜背面新增护墙板齐平、向主卧内开并停靠图纸右侧",
        "notes": "业主设计方向已确认，产品仍未记为已采购；准确下单尺寸、顶／地轴、安装基层、门后净距和收口待厂家项目加工图。",
    })
    by_id["CTRL-MASTER-A"].update({
        "use_location_candidate": "首选主卧门与主卫门全开后两门之间的剩余固定墙；后备左门套／见光板侧",
        "source_ids": "OWNER-INBOX-20260814-001;OUTBOUND-FORM-DOOR-20260817;TB-JINK-EGG-001;JINK-EGG-OFFICIAL-001;TB-PANASURFACE-SWITCH-INSET-001;ZOYLIGHT-D1-20250815-001;ELEC-RESEARCH-20260815-001;OUTBOUND-FORM-SMART-PANEL-20260817",
        "identity_basis": "业主确认 Master A 独立角色及安装面优先级；直接负载与纯场景方案尚未选择",
        "notes": "首选与后备安装面已定，但须在两门全开状态实测固定墙／见光板净宽、底盒净深和五金避让；接口中心坐标保持 unknown。",
    })
    records = [{
        "equipment_id": "DRAIN-CUSTOM-001", "domain": "DRAINAGE", "category": "custom_floor_drain_assembly",
        "item_name": "厕所定制水母地漏/中央集水器＋托克乐思网＋线性排水渠组合", "manufacturer": "", "model": "",
        "variant": "owner_selected_custom_vendor_direction", "quantity": "", "procurement_status": "selected", "decision_status": "partial",
        "storage_location_candidate": "", "use_location_candidate": "卫生间；具体房间和数量待 shop drawing", "storage_location_confirmed": "", "use_location_confirmed": "",
        "schedule_included": "yes", "selector_kind": "logical_input", "selector_value": "DRAIN-CUSTOM-001", "ifc_class": "", "ifc_type_name": "",
        "ifc_type_global_id": "", "ifc_global_ids": "", "source_ids": "OWNER-PLUM-CUSTOM-DRAIN-20260817",
        "identity_basis": "业主确认小红书主页定制渠道与组件组合方向；公开检索不能替代订单或产品证据", "confidence": "1.00",
        "human_review_required": "yes", "legacy_kind": "owner_custom_drain", "legacy_id": "DRAIN-CUSTOM-001",
        "notes": "渠道和方向已定，非已采购/已加工；房间、数量、品牌型号、尺寸、接口和防水均待 shop drawing；不自动替代已购 Geberit CleanLine50。",
    }, {
        "equipment_id": "ACCESS-INTERCOM-001", "domain": "NETWORK", "category": "video_intercom_indoor_monitor",
        "item_name": "入户既有可视对讲室内机", "manufacturer": "DNAKE／狄耐克", "model": "",
        "variant": "existing_indoor_monitor_exact_model_unknown", "quantity": "1", "procurement_status": "existing", "decision_status": "partial",
        "storage_location_candidate": "", "use_location_candidate": "入户既有对讲参考点；开发商标注底边 1400 mm，现场待复核", "storage_location_confirmed": "", "use_location_confirmed": "",
        "schedule_included": "yes", "selector_kind": "logical_input", "selector_value": "ACCESS-INTERCOM-001", "ifc_class": "", "ifc_type_name": "",
        "ifc_type_global_id": "", "ifc_global_ids": "", "source_ids": "OWNER-E304-DNAKE-PHOTO-20260817",
        "identity_basis": "业主现场照片正面可见 DNAKE／狄耐克标识和版本信息；照片没有型号铭牌",
        "confidence": "1.00", "human_review_required": "yes", "legacy_kind": "existing_video_intercom", "legacy_id": "DEV-INTERCOM-001",
        "notes": "品牌已确认；准确型号保持 unknown。端子、线缆、物业系统兼容、保留／迁移／接入方案及最终安装坐标均未关闭；照片中的网络地址不投影。",
    }]
    for record in records:
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
    intercom_confirmed = [
        ("REQ-INTERCOM-001", "manufacturer", "DNAKE／狄耐克", "设备正面品牌标识可见。"),
        ("REQ-INTERCOM-002", "system_version", "1.6.0 20210615", "仅记录照片拍摄时可见软件版本，不用于反推硬件型号。"),
        ("REQ-INTERCOM-003", "application_version", "1.1.0 20210615 16M", "仅记录照片拍摄时可见应用版本，不用于反推硬件型号。"),
    ]
    for req_id, key, value, note in intercom_confirmed:
        upsert(rows, "requirement_id", {
            "requirement_id": req_id, "equipment_id": "ACCESS-INTERCOM-001", "discipline": "ELEC/INT1/NETWORK", "parameter_key": key,
            "value_text": value, "value_number": "", "unit": "", "datum": "", "value_origin": "site_observed", "status": "confirmed",
            "source_id": "OWNER-E304-DNAKE-PHOTO-20260817", "source_locator": "现场照片正面／版本信息页", "blocks_release": "no", "notes": note,
        })
    intercom_pending = [
        ("REQ-INTERCOM-004", "exact_model", "关闭证据：背面／侧面型号标签或工程信息页；不得凭外观或软件版本猜测。"),
        ("REQ-INTERCOM-005", "terminal_diagram_and_existing_wiring", "关闭证据：端子与线缆近照、端子定义和现有接法。"),
        ("REQ-INTERCOM-006", "property_system_compatibility", "关闭证据：物业或门禁维护方书面确认系统角色、兼容与移动限制。"),
        ("REQ-INTERCOM-007", "keep_move_integrate_strategy", "关闭证据：保留、迁移或接入装修后系统的签认结论。"),
        ("REQ-INTERCOM-008", "final_panel_coordinates_and_mounting", "关闭证据：门套完成面后的安装面、净距和底边标高；接口中心坐标不得推测。"),
    ]
    for req_id, key, note in intercom_pending:
        upsert(rows, "requirement_id", {
            "requirement_id": req_id, "equipment_id": "ACCESS-INTERCOM-001", "discipline": "ELEC/INT1/NETWORK", "parameter_key": key,
            "value_text": "unknown", "value_number": "", "unit": "", "datum": "", "value_origin": "pending", "status": "pending",
            "source_id": "OWNER-E304-DNAKE-PHOTO-20260817", "source_locator": "现有照片未显示该信息", "blocks_release": "yes", "notes": note,
        })
    write_csv(path, fields, rows)


def update_drawing_register() -> None:
    path = DECISIONS / "drawing-register.csv"
    fields, rows = read_csv(path)
    by_id = {row["sheet_number"]: row for row in rows}
    by_id["A-001"]["notes"] = "M072 自动汇总全部图号、阶段、版本、说明和图例；最终交付增加材质化正投影产品页要求：平面＋正/展开立面＋材质纹理＋尺寸链＋统一版式，与施工协调图并行且不替代厂家加工图"
    by_id["A-104"]["notes"] = "MS 已确认：M05/M06 为 Sail MONOROTAIA 单轨单扇、官方 DWG 第二种配置、白橡木/浅色木、1000×2400 mm 门扇、2000 mm 暗藏顶轨、无通长地轨、向东/图纸右侧滑开；M07 为 Poliform Pivot 950×2400 mm，与 Senzafine 背面新增护墙板齐平、向主卧内开并停靠图纸右侧。候选图投影业主决定；准确下单尺寸、饰面编号、轴位、基层、净距与收口仍待 06 厂家复核表和项目加工图；正式 IFC 未修改。"
    by_id["E-302"]["notes"] = "Entry A/Master A/Master B 的 4/2/3 键角色保持；MS 已确认 Master A 首选主卧门与主卫门全开后两门之间的剩余固定墙，后备为左门套/见光板侧。候选图显示该优先级但不生成接口中心坐标；实际净宽、底盒净深、M07 开启包络、准确 SKU、端子、负载/浪涌和平嵌适配仍阻断最终发布。"
    by_id["I-501"]["notes"] = "70 个既有对象协调包络入图；APP-009/010 两台未采购 SJ85ZX26MC 候选的官方柜孔、G3/4 冷水与 Ø38 排水约束已作为 2 条无定位接口表项入图，不生成粗装 XYZ、阀门、软管路径或开孔位置；EL-03 视图 06～09 及 EL-P01/P02 公共空间折线展开和贯穿长立面已写入正式 IFC 的 Bonsai 原生 Drawing；新增 1:30 无纹理 SVG 与手机 PNG；其余厂家安装图、燃气实测和五金运动包络未关闭；中厨墙面更新为灶台操作墙台面同材大板＋远端墙浅置物架，准确墙段、材料、加工和架后背衬待深化；最终需输出材质化正投影厨房产品页"
    by_id["I-502"]["notes"] = "33 个既有对象协调包络入图；EL-06 视图 18～22 与 EL-08 视图 27～32 已写入正式 IFC；原 9 个深层 BRep、公共空间新增 2 个对象及全量重编发现的 BED02 现均有纯折线 ELEVATION_VIEW，详细 Body 未替换；洁具粗装图、五金运动包络和节点未关闭；厕所定制地漏渠道与水母地漏/中央集水器＋托克乐思网＋线性排水渠方向已确认，准确房间、数量和 shop drawing 待商家；最终需输出材质化正投影卫生间产品页"
    by_id["I-504"]["notes"] = "25 个既有对象协调包络入图；EL-01/02/04/05/07/09 及公共空间 EL-P01/P02 已写入正式 IFC 的 Bonsai 原生 Drawing；SIS04、HIMA01、BED01 等复杂家具使用受控折线立面。MS 已确认 M07 与 Senzafine 背面新增护墙板齐平、向主卧内开并停靠图纸右侧；Master A 首选两门全开后的中间固定墙、后备左门套/见光板侧。候选图只表达协调关系，准确护墙板厚度、门后净距、底盒和加工节点仍待厂家/现场。"
    by_id["E-304"]["notes"] = "两个吸顶 AP 已采用 CAT6 星型回弱电箱 PoE 交换机的拓扑；既有可视对讲室内机品牌已由现场照片确认为 DNAKE／狄耐克，准确型号、端子、物业兼容和保留／迁移／接入方案仍待现场／物业确认；既有门铃仍为开发商参考点。实购 AP 功耗、PoE 总预算、弱电箱净尺寸/柜门温升、网线通断、门禁接口和端接顺序仍待厂家/物业/现场证据。"
    by_id["S-701"]["notes"] = "家具产品身份、家电使用/存放、门窗五金、既有暖通类型、定制地漏方向和墙面材料系统按已确认/候选/未决分栏；M05/M06 的第二种官方 DWG 配置、白橡木/浅色木、东向/图纸右侧滑动与 M07 的护墙板齐平、图纸右侧停靠已按 MS 回复登记，厂家项目加工图仍未齐；中厨采用台面同材大板＋远端墙浅置物架；厕所定制地漏渠道已定但 shop drawing 未齐；不是下单表。"
    by_id["D-601"]["notes"] = "材料与收口变量节点待复核；新增干区墙面同色齐平宽踢脚＋上下阴影缝候选，须以 1:5 墙脚节点、分房间适用表、门套/隐形门/柜体通缝展开和实物样板关闭；中厨 VVD Pewter 60 mm 踢脚保持独立，湿区暂不套用"
    write_csv(path, fields, rows)


def update_finish_details() -> None:
    path = DECISIONS / "wfin-open-issues.csv"
    fields, rows = read_csv(path)
    upsert(rows, "issue_id", {
        "issue_id": "WFIN-R06",
        "scope": "dry-area-flush-wide-baseboard-shadow-gap",
        "object_guid": "",
        "current_evidence": "业主要求记录干区墙面同色齐平宽踢脚候选；参考案例为上下各约 10 mm 阴影缝、墙下口型材和门墙柜连续通缝。中厨 VVD Pewter 60 mm 踢脚保持独立，湿区暂不采用上下双缝。",
        "required_action_or_decision": "由 D-601 输出 1:5 墙脚节点、分房间适用表、门套/隐形门/柜体/转角通缝展开，并以实物样板关闭材料、尺寸、公差、耐撞、拖地水和积灰清洁。",
        "basis": "OWNER-WFIN-BASEBOARD-SHADOW-GAP-20260817；参考案例只支持候选方向，不支持项目加工尺寸或湿区适用性。",
        "confidence": "1.00",
        "review_required": "yes",
        "status": "candidate_pending_detail_and_mockup",
        "stop_condition": "节点和样板未批准前不得把参考案例尺寸写成项目尺寸、不得扩展到中厨或湿区、不得写正式 IFC 材料或几何。",
    })
    write_csv(path, fields, rows)

    path = DECISIONS / "det1-detail-review.csv"
    fields, rows = read_csv(path)
    by_id = {row["node_id"]: row for row in rows}
    row = by_id["D601-N04"]
    row.update({
        "confirmed_evidence": "中厨已确认 VVD 材料组合及 Pewter 60 mm 踢脚；业主另要求记录卧室、走廊等干区墙面同色齐平宽踢脚＋上下阴影缝候选，门套、隐形门和柜体连续通缝；准确节点未定",
        "variable_parameters": "房间/墙段｜踢脚高度厚度材料｜上下阴影缝宽与公差｜墙下口型材｜门柜转角通缝｜墙顶收口｜材料分界｜密封/清洁策略",
        "unresolved_material_or_product": "Tadelakt 与大白墙系统仍待选；干区踢脚材料、颜色和截面未定；中厨只待项目样板批次和加工节点",
        "unresolved_construction": "干区 1:5 墙脚节点、分房间适用、基层型材、门柜转角通缝、耐撞、拖地水和积灰清洁样板未定；湿区排除边界待确认",
        "source_references": "wfin-open-issues.csv WFIN-R01/R02/R03/R05/R06；OWNER-WFIN-BASEBOARD-SHADOW-GAP-20260817；OWNER-WFIN-VVD-SINKS-20260817；施工图检查清单 D-601",
    })
    write_csv(path, fields, rows)


def update_door_coordination() -> None:
    path = DECISIONS / "equipment-installation-requirements.csv"
    fields, rows = read_csv(path)
    by_id = {row["requirement_id"]: row for row in rows}
    by_id["REQ-0215"].update({
        "value_origin": "user_input", "status": "confirmed", "source_id": "OUTBOUND-FORM-DOOR-20260817",
        "source_locator": "D01-D04；MS／2026-08-17", "blocks_release": "no",
        "notes": "MS 确认项目名义门扇宽度；准确下单尺寸仍由 06 表和厂家项目加工图关闭，未写入正式 IFC。",
    })
    by_id["REQ-0216"].update({
        "value_origin": "user_input", "status": "confirmed", "source_id": "OUTBOUND-FORM-DOOR-20260817",
        "source_locator": "D01-D04；MS／2026-08-17", "blocks_release": "no",
        "notes": "MS 确认项目名义门扇高度；D03 的 2 道水平横档未确认并从项目要求删除；准确下单尺寸仍待厂家项目加工图。",
    })
    by_id["REQ-0217"].update({
        "value_text": "Sail MONOROTAIA single-track single-leaf; slide east/drawing-right",
        "value_origin": "user_input", "status": "confirmed", "source_id": "OUTBOUND-FORM-DOOR-20260817",
        "source_locator": "D01-D02、D06；MS／2026-08-17", "blocks_release": "no",
        "notes": "MS 明确纠正旧向左／西候选；当前项目方向为向东／图纸右侧滑开。厂家兼容与完整开启包络仍由 06 表及加工图关闭。",
    })
    by_id["REQ-0220"].update({
        "value_origin": "user_input", "status": "confirmed", "source_id": "OUTBOUND-FORM-DOOR-20260817",
        "source_locator": "D04；MS／2026-08-17", "blocks_release": "no",
        "notes": "MS 确认 2000 mm 项目名义顶轨；不沿用正式 IFC 的 1200 mm 早期示意。准确下单长度和固定方式待厂家加工图。",
    })
    by_id["REQ-0222"].update({
        "value_text": "Sail MONOROTAIA concealed ceiling track; east/drawing-right slide; no continuous floor track; lower anti-sway/positioning component per Rimadesio final shop drawing",
        "value_origin": "research_conclusion", "status": "candidate",
        "source_id": "OUTBOUND-FORM-DOOR-VENDOR-20260817",
        "source_locator": "D02、D04-D07；厂家复核 R01-R06", "blocks_release": "yes",
        "notes": "MS 已确认同一门组、暗藏顶轨、无通长地轨、向东／图纸右侧滑开和检修要求；下部防摆／定位构件、吊顶固定与收口仍由 Rimadesio 项目加工图确认。",
    })
    by_id["REQ-0225"].update({
        "value_origin": "user_input", "status": "confirmed", "source_id": "OUTBOUND-FORM-DOOR-20260817",
        "source_locator": "D08-D09；MS／2026-08-17", "blocks_release": "no",
        "notes": "MS 确认项目名义宽度；准确下单尺寸仍由 06 表和厂家项目加工图关闭。",
    })
    by_id["REQ-0226"].update({
        "value_origin": "user_input", "status": "confirmed", "source_id": "OUTBOUND-FORM-DOOR-20260817",
        "source_locator": "D08-D09；MS／2026-08-17", "blocks_release": "no",
        "notes": "MS 确认项目名义高度；准确下单尺寸仍由 06 表和厂家项目加工图关闭。",
    })
    by_id["REQ-0227"].update({
        "value_text": "single Pivot; opens inward to master; parks at drawing-right fixed wall/guest-bedroom side",
        "value_origin": "user_input", "status": "confirmed", "source_id": "OUTBOUND-FORM-DOOR-20260817",
        "source_locator": "D08-D12；MS／2026-08-17", "blocks_release": "no",
        "notes": "MS 明确纠正旧图纸左侧停靠候选；当前为向主卧内开并停靠图纸右侧固定墙／次卧侧。厂家兼容、轴位和净距仍由 06 表及加工图关闭。",
    })
    confirmed = [
        ("REQ-DOOR-OWNER-001", "DW-M05", "selected_official_configuration", "official DWG second configuration; concealed ceiling top track; jamb wall parallel to door travel", "D03；MS／2026-08-17", "不推测原厂配置代码；由厂家在 06 表填写准确配置／加工图编号。"),
        ("REQ-DOOR-OWNER-002", "DW-M05", "finish_direction", "white oak / light wood", "D03；MS／2026-08-17", "业主饰面方向已定；准确原厂饰面编号和样板仍待厂家。"),
        ("REQ-DOOR-OWNER-003", "DW-M06", "service_access", "removable end access allowing trolley adjustment and replacement without dismantling main ceiling", "D07；MS／2026-08-17", "具体检修口尺寸和构造由项目加工图关闭。"),
        ("REQ-DOOR-OWNER-004", "DW-M07", "closed_alignment", "align with added wall panel on rear of Senzafine wardrobe", "D10；MS／2026-08-17", "已纠正旧衣柜正面齐平候选；护墙板厚度和收口待加工图。"),
        ("REQ-DOOR-OWNER-005", "CTRL-MASTER-A", "preferred_installation_surface", "primary: remaining fixed wall between master-bedroom and master-bath doors when both open; fallback: left jamb/reveal-panel side", "D13；MS／2026-08-17", "最终安装坐标、净宽、标高和底盒净深保持 unknown，待现场实测。"),
    ]
    for req_id, equipment_id, key, value, locator, notes in confirmed:
        upsert(rows, "requirement_id", {
            "requirement_id": req_id, "equipment_id": equipment_id, "discipline": "ARCH/INT1/ELEC", "parameter_key": key,
            "value_text": value, "value_number": "", "unit": "", "datum": "", "value_origin": "user_input", "status": "confirmed",
            "source_id": "OUTBOUND-FORM-DOOR-20260817", "source_locator": locator, "blocks_release": "no", "notes": notes,
        })
    write_csv(path, fields, rows)

    path = DECISIONS / "s701-schedule-review.csv"
    fields, rows = read_csv(path)
    by_id = {row["schedule_id"]: row for row in rows}
    by_id["S701-M05"].update({
        "confirmed_scope": "身份/定位；项目名义尺寸 1000×2400 mm；官方 DWG 第二种配置；白橡木／浅色木；向东／图纸右侧滑开",
        "candidate_or_observed_scope": "Sail MONOROTAIA 单轨单扇；正式 IFC 仍为早期几何观察",
        "unresolved_for_release": "厂家项目加工图、准确下单尺寸、原厂饰面编号、完整开启包络、安装/收口及现场复核",
        "evidence_reference": "IFC-FORMAL-001;RIMADESIO-SAIL-MONOROTAIA-001;OUTBOUND-FORM-DOOR-20260817;OUTBOUND-FORM-DOOR-VENDOR-20260817",
    })
    by_id["S701-M06"].update({
        "confirmed_scope": "身份/定位；项目名义顶轨 2000 mm；暗藏顶轨、无通长地轨；向东／图纸右侧滑开；保留末端检修",
        "candidate_or_observed_scope": "Sail MONOROTAIA 同门组顶轨；下部防摆／定位构件按最终 Rimadesio 项目加工图",
        "unresolved_for_release": "厂家项目加工图、准确下单长度、下部防摆／定位构件、吊顶固定、安装/收口及现场复核",
        "evidence_reference": "IFC-FORMAL-001;RIMADESIO-SAIL-MONOROTAIA-001;OUTBOUND-FORM-DOOR-20260817;OUTBOUND-FORM-DOOR-VENDOR-20260817",
    })
    by_id["S701-M07"].update({
        "confirmed_scope": "Poliform Pivot 墙装单扇；项目名义尺寸 950×2400 mm；与衣柜背面新增护墙板齐平；向主卧内开并停靠图纸右侧",
        "candidate_or_observed_scope": "官方 Pivot 门框、顶／地轴和金属门扇框；正式 IFC OperationType 仍为 NOTDEFINED",
        "unresolved_for_release": "厂家项目加工图、准确下单尺寸、顶／地轴、门框护墙板收口、门后净距、基层与现场复核",
        "evidence_reference": "IFC-FORMAL-001;POLIFORM-ARCHITECTURAL-PDF-001;OUTBOUND-FORM-DOOR-20260817;OUTBOUND-FORM-DOOR-VENDOR-20260817",
    })
    write_csv(path, fields, rows)


def main() -> None:
    update_sources()
    update_owner_inputs()
    update_closeout_rules()
    update_equipment()
    update_requirements()
    update_finish_details()
    update_door_coordination()
    update_drawing_register()
    print("Applied 2026-08-17 kitchen, product-sheet, custom-drain, and DNAKE intercom owner updates.")


if __name__ == "__main__":
    main()
