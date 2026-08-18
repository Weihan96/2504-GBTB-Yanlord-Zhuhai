#!/usr/bin/env python3
"""Absorb the 2026-08-18 HVAC, plumbing, and kitchen owner responses."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "pipeline/decisions"
FORM_DIR = "output/forms/对外确认表"


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def index(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in rows}


def upsert(rows: list[dict[str, str]], key: str, record: dict[str, str], fields: list[str]) -> None:
    clean = {field: record.get(field, "") for field in fields}
    for position, row in enumerate(rows):
        if row[key] == record[key]:
            rows[position] = clean
            return
    rows.append(clean)


def merge_ids(*values: str) -> str:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        for item in value.split(";"):
            item = item.strip()
            if item and item not in seen:
                seen.add(item)
                result.append(item)
    return ";".join(result)


def append_once(value: str, addition: str) -> str:
    """Append one normalized sentence block and remove prior repeated copies."""
    base = value
    while addition in base:
        base = base.replace(addition, "")
    return " ".join(part for part in (base.strip(), addition.strip()) if part)


def append_semicolon_once(value: str, addition: str) -> str:
    """Append one semicolon-delimited clause while collapsing prior copies."""
    clause = addition.strip().strip("；")
    parts = [
        part.strip()
        for part in value.split("；")
        if part.strip() and part.strip() != clause
    ]
    return "；".join([*parts, clause])


def sha256(relative: str) -> str:
    digest = hashlib.sha256()
    with (ROOT / relative).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_record(fields: list[str], **values: str) -> dict[str, str]:
    row = {field: "" for field in fields}
    row.update({
        "confidence": "1.00",
        "review_required": "yes",
        "formal_ifc_write_allowed": "no",
        "revision": "2026-08-18",
        "publication_date": "2026-08-18",
    })
    row.update(values)
    return row


def update_sources() -> None:
    path = DECISIONS / "source-evidence-register.csv"
    fields, rows = read_csv(path)
    by_id = index(rows, "source_id")

    form_updates = {
        "OUTBOUND-FORM-HVAC-20260817": {
            "file": "02-日立空调接口确认表-发空调厂家.md",
            "status": "verified_partial_owner_response_external_followup_pending",
            "locator": "H01-H03 由 MS 回复；H04-H09、P01-P07 待厂家；回复日期 2026-08-18",
            "evidence": "MS 确认 A01-A04 型号族、已登记管径和接管侧，并要求核实保温厚度及模型匹配；其余型号、接口和线控问题仍空白",
            "proves": "业主认可 A01-A04 当前型号族和接管侧方向，并提出保温厚度核验任务",
            "does_not_prove": "厂家已回复、保温厚度已确定、最终管线与保温层已建模，或 PC-P1HEQ 已完成兼容和控制映射",
        },
        "OUTBOUND-FORM-PLUM-20260817": {
            "file": "03-给排水设备确认表-发设备与施工方.md",
            "status": "verified_owner_response_external_followup_pending",
            "locator": "S01-S07、G01-G07、K01-K06、V01-V11、C01-C06；MS／2026-08-18",
            "evidence": "MS 确认洗烘、VVD、定制地漏和石材盆方向，否决指定 Geberit 安装项、F50、trolley 和悬挂式油烟机，并提出 Sistema 7、Eco kits 与厨房布局深化任务",
            "proves": "业主产品和设计方向可进入 SSOT，并可删除已由官方资料或业主回答关闭的重复问题",
            "does_not_prove": "厂家兼容、项目加工图、接口中心、定制地漏防水节点、油烟机风管或 APP-016 最终选型",
        },
    }
    for source_id, metadata in form_updates.items():
        row = by_id[source_id]
        relative = f"{FORM_DIR}/{metadata['file']}"
        row.update({
            "sha256": sha256(relative),
            "locator": metadata["locator"],
            "evidence": metadata["evidence"],
            "proves": metadata["proves"],
            "does_not_prove": metadata["does_not_prove"],
            "status": metadata["status"],
            "revision": "2026-08-18",
            "publication_date": "2026-08-18",
            "notes": "原始 MS 回答保留；未填写项继续由外部方关闭，不将业主回答冒充厂家签认。",
        })
    by_id["EXT-EXTERNAL-INFO-MINIMUM-20260817"].update({
        "sha256": sha256("drawings/evidence/EXT-外部信息最短清单-20260817.md"),
        "status": "superseded_project_external_closeout_index",
        "does_not_prove": "当前派件内容、外部方已回复或任何项目接口已冻结",
        "notes": "保留为 2026-08-17 历史状态；不得继续作为现行发件清单。",
    })
    by_id["OWNER-WFIN-VVD-SINKS-20260817"].update({
        "sha256": sha256("drawings/evidence/OWNER-WFIN-VVD厨房材质与水槽-20260817.md"),
        "evidence": "VVD 材料与两只项目水槽决定；2026-08-18 补充无 trolley、无悬挂式 hood、Sistema 7 四门玻璃墙柜内整合油烟机、灯光供货协作及 Eco kits 盆下收纳",
        "proves": "业主采用的厨房产品语言、材料、水槽和本轮修正方向",
        "does_not_prove": "已下单、项目加工尺寸、油烟机／风管／灯光／收纳接口或厂家签认",
        "status": "verified_owner_design_direction_shop_drawing_pending",
        "revision": "2026-08-18",
        "notes": "2026-08-18 决定取代同文件中 2021 产品页的 trolley 和 suspended hood 示例；保留原产品页只作演变追溯。",
    })

    local_sources = [
        source_record(
            fields,
            source_id="OWNER-RESPONSE-HVAC-PLUM-20260818",
            discipline="HVAC/PLUM/ELEC/INT1/WFIN",
            sheet_id="RCP1/M-401/P-201/P-202/E-303/I-501/I-503/S-701",
            decision_scope="MS 第二轮 HVAC 与给排水回复结构化吸收",
            source_kind="project_owner_response_normalization",
            source_document="OWNER-RESPONSE-HVAC-PLUM-20260818.md",
            local_path="drawings/evidence/OWNER-RESPONSE-HVAC-PLUM-20260818.md",
            sha256=sha256("drawings/evidence/OWNER-RESPONSE-HVAC-PLUM-20260818.md"),
            locator="可直接写入的决定；不能关闭的事项；证据边界",
            evidence="把 02／03 表的 MS 回复拆为 confirmed、research conclusion 与 pending，并明确正式 IFC 未修改",
            proves="本轮业主回复的受控投影范围和后续外部证据门",
            does_not_prove="任何厂家、施工方或主管部门已签认",
            status="verified_owner_response_projection_boundary",
            notes="原始逐项回答仍以 OUTBOUND-FORM-HVAC-20260817 和 OUTBOUND-FORM-PLUM-20260817 为准。",
        ),
        source_record(
            fields,
            source_id="RCP1-HITACHI-INSULATION-BOUNDARY-20260818",
            discipline="HVAC/PLUM",
            sheet_id="RCP1/M-401",
            decision_scope="日立风管机保温厚度与模型匹配边界",
            source_kind="project_evidence_boundary_note",
            source_document="RCP1-HITACHI-保温厚度证据边界-20260818.md",
            source_url="https://www.hisensehitachi.com/upload/accessory/20236/p1h2esc855pnidai1v7k10ha17rp4.pdf",
            local_path="drawings/evidence/RCP1-HITACHI-保温厚度证据边界-20260818.md",
            sha256=sha256("drawings/evidence/RCP1-HITACHI-保温厚度证据边界-20260818.md"),
            locator="P02010Q pp.5-6, 11；与当前模型的关系",
            evidence="官方要求冷媒管和冷凝水管保温，但保温管由现场提供且资料未给厚度；正式 IFC 仅有旧路线骨架且没有最终保温层",
            proves="保温是必需条件，厚度和最终模型匹配仍须厂家／机电设计关闭",
            does_not_prove="任何保温厚度、材料、防火等级、最终外径、路线或接口中心",
            status="verified_official_requirement_thickness_unknown",
            manufacturer="日立／海信日立",
            model_scope="RPIZ-22FSLN5QD/P; RPIZ-22FSLN5QDF/P",
        ),
        source_record(
            fields,
            source_id="APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818",
            discipline="ELEC/PLUM/INT1",
            sheet_id="E-303/P-201/I-501/S-701",
            decision_scope="APP-016 Franke Slim 50 CN 替代研究",
            source_kind="project_research_boundary_note",
            source_document="APP-016-Franke-Slim50-CN研究-20260818.md",
            source_url="https://www.franke.com/cn/zh/home-solutions/产品中心/waste-management/product-detail-page.html/134.0721.210.html",
            local_path="drawings/evidence/APP-016-Franke-Slim50-CN研究-20260818.md",
            sha256=sha256("drawings/evidence/APP-016-Franke-Slim50-CN研究-20260818.md"),
            locator="当前决定；可核验研究候选；证据边界",
            evidence="F50 已否决；Franke 中国准确 SKU 134.0721.210 的官方公开参数可核验，但官方可访问资料未给完整安装包络",
            proves="Franke Slim 50 CN 可作为更强官方证据的研究候选",
            does_not_prove="与 F50 同高、最终选型、额定输入功率、法兰、排水、洗碗机支管或拆换净空",
            status="official_research_candidate_project_selection_pending",
            manufacturer="Franke",
            model_scope="Slim 50 CN / LD370-E01B / 134.0721.210",
        ),
        source_record(
            fields,
            source_id="MOLTENI-SISTEMA7-ECOKITS-RESEARCH-20260818",
            discipline="INT1/WFIN/ELEC/PLUM",
            sheet_id="I-501/P-201/E-303/S-701",
            decision_scope="Sistema 7 柜内油烟机与 Eco kits 盆下收纳",
            source_kind="project_research_boundary_note",
            source_document="MOLTENI-Sistema7-EcoKits研究-20260818.md",
            source_url="https://molteni.it/en/download/document/3c1c2c00b23e4dbb8ab7c89fa68936d8008a7a56",
            local_path="drawings/evidence/MOLTENI-Sistema7-EcoKits研究-20260818.md",
            sha256=sha256("drawings/evidence/MOLTENI-Sistema7-EcoKits研究-20260818.md"),
            locator="Sistema 7 油烟机整合；Eco kits 盆下收纳；IFC 语义纠正",
            evidence="Molteni 官方产品族支持 Sistema 7 整合油烟机及 Eco kits 盆下抽屉／垃圾桶／托盘；正式 IFC 对象名支持玻璃门柜语义",
            proves="项目方向具有官方产品族依据，可进入候选深化和封闭式厂家复核",
            does_not_prove="项目四门柜、油烟机、风管、600 mm 盆下收纳或加工接口已完成厂家签认",
            status="verified_official_family_research_project_shop_drawing_pending",
            manufacturer="Molteni&C / Dada",
            model_scope="Sistema 7; Eco kits",
        ),
        source_record(
            fields,
            source_id="OUTBOUND-FORM-KITCHEN-20260818",
            discipline="INT1/PLUM/ELEC/WFIN",
            sheet_id="I-501/P-201/E-303/S-701",
            decision_scope="厨房设备与柜体项目深化复核",
            source_kind="project_external_confirmation_markdown",
            source_document="09-厨房设备与柜体深化确认表-发橱柜设备方.md",
            local_path=f"{FORM_DIR}/09-厨房设备与柜体深化确认表-发橱柜设备方.md",
            sha256=sha256(f"{FORM_DIR}/09-厨房设备与柜体深化确认表-发橱柜设备方.md"),
            locator="K01-K07",
            evidence="只询问现有岛台服务路径、上方干区插座、Sistema 7 柜内油烟机、玻璃门柜、Eco kits、柜内灯和 APP-016 安装证据",
            proves="已回答和已有官方证据的问题已删除，剩余问题均可由项目加工图或同型号资料关闭",
            does_not_prove="外部方已经回复或任何项目加工尺寸已冻结",
            status="verified_confirmation_markdown_pending_external_response",
        ),
    ]
    for record in local_sources:
        upsert(rows, "source_id", record, fields)
    write_csv(path, fields, rows)


def update_owner_inputs() -> None:
    path = DECISIONS / "owner-input-register.csv"
    fields, rows = read_csv(path)
    by_id = index(rows, "input_id")
    by_id["RCP1-HVAC-PORTS"].update({
        "user_value": "MS 确认 A01-A04 当前型号族、裸管规格和接管侧；要求核实保温厚度与模型匹配",
        "status": "需证据",
        "evidence_reference": "OUTBOUND-FORM-HVAC-20260817;RCP1-HITACHI-INSULATION-BOUNDARY-20260818",
        "source_basis": "机位和路线骨架已确认但精确端口坐标尚无批准证据",
        "notes": "A01-A04 业主方向已关闭；A05/A06、全部精确端口、保温厚度、最终路线和模型匹配继续由日立／机电设计关闭。",
    })
    by_id["PLUM-ROUGHINS"].update({
        "user_value": "MS 已确认洗烘、VVD、定制地漏和石材盆方向；Foster 官方尺寸已内部核验；剩余厨房项目加工接口见 09 表",
        "evidence_reference": "OUTBOUND-FORM-PLUM-20260817;OWNER-RESPONSE-HVAC-PLUM-20260818;OUTBOUND-FORM-KITCHEN-20260818",
        "notes": "不再重复询问 Foster、洗烘电气语义或已否决 Geberit；只收准确墙排、柜体服务路径、定制地漏／石材盆／油烟机加工图和洁具接口。",
    })
    by_id["APP016-DISPOSER-DATA"].update({
        "question": "补勒科斯 F50 订单/铭牌/厂家安装图",
        "candidate_value": "保持候选，不按同类产品估算",
        "user_value": "勒科斯 F50 否决；寻找高度接近但品质／品牌更好的产品",
        "status": "需证据",
        "evidence_reference": "OUTBOUND-FORM-PLUM-20260817;APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818;OUTBOUND-FORM-KITCHEN-20260818",
        "source_basis": "仅有业主候选名称",
        "notes": "Franke 只是研究候选，不是已选／已购；只接受最终同型号铭牌和完整尺寸／接口图。",
    })
    by_id["APP017-SIEMENS-STACK"].update({
        "user_value": "WG54M7D20W 下＋WQ55M7U20W 上；WTZ27510 采用方向；两个侧置 10A 插座共用 C16；中厨蝴蝶门；洗衣和干衣冷凝水接专用墙排",
        "status": "自定义确认",
        "evidence_reference": merge_ids(by_id["APP017-SIEMENS-STACK"]["evidence_reference"], "OUTBOUND-FORM-PLUM-20260817;OWNER-RESPONSE-HVAC-PLUM-20260818"),
        "notes": "业主方向已关闭；WTZ27510 与两台 E-Nr./FD 兼容、准确双口墙排产品／接口、蝴蝶门通风和整机抽出节点仍待外部书面确认。",
    })

    additions = [
        {
            "input_id": "PLUM-GEBERIT-CURRENT-USE-20260818", "workstream": "PLUM/INT1/S-701", "priority": "P1", "blocks_release": "no",
            "question": "已登记 Geberit 排水组件是否继续用于当前方案", "candidate_value": "按采购历史保留并按项目适用性决定",
            "user_value": "152.464.00.1、154.150.00.1、388.013.00.2、154.446.KS.1 当前方案不采用；224.212.00.2、146.361.00.1 保留",
            "unit": "", "status": "自定义确认", "evidence_reference": "OUTBOUND-FORM-PLUM-20260817;OWNER-RESPONSE-HVAC-PLUM-20260818",
            "source_basis": "MS 逐项回复 G01-G07 和 C02", "sync_target": "equipment SSOT;P-201;P-202;I-502;S-701",
            "notes": "采购／收货记录继续保留历史；不采用不等于未购买。115.770.11.5 为 Sigma01 冲水面板，仍需补足数量和安装节点。",
        },
        {
            "input_id": "APP017-WALL-DRAIN-20260818", "workstream": "PLUM/INT1/I-503", "priority": "P0", "blocks_release": "yes",
            "question": "洗衣机与干衣机冷凝水专用墙排的准确产品和项目接口", "candidate_value": "可检修双接口墙排；洗衣和干衣冷凝水分别接入",
            "user_value": "两台设备均接洗衣机专用墙排，干衣机接其烘干机排水孔位", "unit": "mm", "status": "需证据",
            "evidence_reference": "OUTBOUND-FORM-PLUM-20260817;OWNER-RESPONSE-HVAC-PLUM-20260818",
            "source_basis": "MS 确认使用意图，但 388.013.00.2 被否决且它是地漏，不是已确认的双接口墙排",
            "sync_target": "P-201;I-503;S-701", "notes": "关闭证据：准确品牌／型号、双接口图、排水标高、管径、防回流、防臭、可清洁和柜内检修剖面。",
        },
        {
            "input_id": "INT1-VVD-HOOD-TROLLEY-20260818", "workstream": "INT1/WFIN/ELEC", "priority": "P1", "blocks_release": "no",
            "question": "VVD 岛台 trolley 与油烟机采用什么方向", "candidate_value": "按 2021 产品页示意", "user_value": "无 trolley；无悬挂式 hood；油烟机整合在 Sistema 7 四门玻璃墙柜内",
            "unit": "", "status": "自定义确认", "evidence_reference": "OUTBOUND-FORM-PLUM-20260817;MOLTENI-SISTEMA7-ECOKITS-RESEARCH-20260818",
            "source_basis": "MS 明确否决产品页示意；Molteni 官方产品族支持柜内油烟机可行性",
            "sync_target": "I-501;E-303;S-701", "notes": "方向已定；准确四门柜、油烟机、风管、净距和检修由 09 表项目加工图关闭。",
        },
        {
            "input_id": "INT1-KITCHEN-ISLAND-LAYOUT-20260818", "workstream": "INT1/PLUM/ELEC", "priority": "P1", "blocks_release": "yes",
            "question": "两台洗碗机相邻还是分置水槽两侧", "candidate_value": "保持正式 IFC 当前相邻布局：600 石材盆＋洗碗机＋洗碗机",
            "user_value": "希望比较动线；现状相邻可获得最大连续台面", "unit": "mm", "status": "需证据",
            "evidence_reference": "IFC-FORMAL-001;OUTBOUND-FORM-PLUM-20260817;OUTBOUND-FORM-KITCHEN-20260818",
            "source_basis": "现有 IFC 四个 600 mm 模块和两个相邻 BD600DW 可机械确认；相邻布置保留连续台面并减少改管",
            "sync_target": "I-501;P-201;E-303;S-701", "notes": "由 09 表关闭第二台洗碗机的软管长度、柜孔和检修路径；未签认前不改正式 IFC。",
        },
        {
            "input_id": "SAN022-ECO-KIT-20260818", "workstream": "INT1/PLUM/WFIN", "priority": "P1", "blocks_release": "yes",
            "question": "SAN-022 600 mm 石材盆下方如何收纳并避让管线", "candidate_value": "参考 Molteni Eco kits 的盆下抽屉、垃圾桶和清洁用品托盘",
            "user_value": "抽拉式柜门；内设垃圾桶和海绵收纳", "unit": "mm", "status": "需证据",
            "evidence_reference": "OUTBOUND-FORM-PLUM-20260817;MOLTENI-SISTEMA7-ECOKITS-RESEARCH-20260818;OUTBOUND-FORM-KITCHEN-20260818",
            "source_basis": "业主确认功能；Molteni 官方页面支持产品族做法，但没有本项目 600 mm 加工尺寸",
            "sync_target": "I-501;P-201;S-701", "notes": "关闭证据：盆体、抽屉、桶体、存水弯／排水、溢水、支撑和检修联合剖面。",
        },
    ]
    for record in additions:
        upsert(rows, "input_id", record, fields)
    write_csv(path, fields, rows)


def update_closeout_rules() -> None:
    path = DECISIONS / "owner-input-closeout-rules.csv"
    fields, rows = read_csv(path)
    by_id = index(rows, "input_id")
    by_id["RCP1-HVAC-PORTS"].update({
        "required_evidence": "日立／机电设计提交 A05/A06 准确型号和 A01-A06 项目接口图；书面给出冷媒／冷凝水保温材料与厚度、完成外径、送回风、端口、冷凝水排放和检修条件；证明最终路线与吊顶净高相容。",
        "notes": "A01-A04 型号族、裸管规格和接管侧已由 MS 接受；P02010Q 未给保温厚度，旧 IFC 管线不能关闭最终模型匹配。",
    })
    by_id["PLUM-ROUGHINS"].update({
        "required_evidence": "提交准确墙排、柜体服务路径、洁具接口、定制地漏／石材盆／油烟机联合 shop drawing；说明书未给的中心坐标由项目图关闭。",
        "notes": "Foster 官方参数、洗烘电气语义和被否决 Geberit 不再重复询问；厨房项目深化问题直接见 09 表。",
    })
    by_id["APP016-DISPOSER-DATA"].update({
        "required_evidence": "最终所选垃圾处理器同型号铭牌与官方／供货商签认尺寸图，包含额定输入、插头／控制、总高直径、法兰／开孔、排水、洗碗机支管、空气开关和拆换包络。",
        "notes": "F50 已否决；Franke 134.0721.210 只作研究候选，不能用其他地区 SKU 或同类机尺寸关闭。",
    })
    by_id["APP017-SIEMENS-STACK"].update({
        "required_evidence": "西门子书面确认 WTZ27510 与两台实购 E-Nr./FD 兼容；给排水／柜体方提交准确双接口墙排产品图及洗衣／干衣冷凝水独立接入口、阀门、两个侧置 10A 插座、蝴蝶门、通风和整机抽出联合剖面。",
        "notes": "上下顺序、两个 10A 插座、共用 C16 和侧置检修已由 MS 确认，不再重复询问。",
    })
    additions = [
        {"input_id": "PLUM-GEBERIT-CURRENT-USE-20260818", "closeout_kind": "human_design_selection", "responsible_party": "业主/给排水设计", "required_evidence": "MS 对 G01-G07、C02 的逐项确认", "automatic_close_allowed": "no", "notes": "当前用途决定已关闭；采购历史保留。"},
        {"input_id": "APP017-WALL-DRAIN-20260818", "closeout_kind": "product_interface_evidence", "responsible_party": "西门子/给排水/家政柜", "required_evidence": "准确墙排品牌型号、双接口图、排水标高管径、防回流防臭、清洁和柜内检修剖面", "automatic_close_allowed": "no", "notes": "388.013.00.2 已否决且不是可直接替代的墙排。"},
        {"input_id": "INT1-VVD-HOOD-TROLLEY-20260818", "closeout_kind": "human_design_selection", "responsible_party": "业主/室内设计", "required_evidence": "MS 确认无 trolley、无悬挂式 hood、油烟机整合 Sistema 7 四门玻璃墙柜", "automatic_close_allowed": "no", "notes": "方向已关闭，项目接口由 09 表关闭。"},
        {"input_id": "INT1-KITCHEN-ISLAND-LAYOUT-20260818", "closeout_kind": "vendor_shop_drawing_and_compatibility", "responsible_party": "全屋定制/西门子/给排水/电气", "required_evidence": "在 09 表 K01 答是／否并提交两台进水排水插座柜孔软管检修路径的岛台平面剖面", "automatic_close_allowed": "no", "notes": "保持相邻是研究结论；厂家图返回前不改正式 IFC。"},
        {"input_id": "SAN022-ECO-KIT-20260818", "closeout_kind": "vendor_shop_drawing_and_compatibility", "responsible_party": "全屋定制/石材/给排水", "required_evidence": "09 表 K05 联合剖面：盆体、抽屉、垃圾桶、海绵托盘、存水弯排水、溢水、支撑和检修", "automatic_close_allowed": "no", "notes": "Eco kits 只证明产品族做法。"},
    ]
    for record in additions:
        upsert(rows, "input_id", record, fields)
    write_csv(path, fields, rows)


def update_appliance_inputs() -> None:
    path = DECISIONS / "appliance-input-register.csv"
    fields, rows = read_csv(path)
    by_id = index(rows, "appliance_id")
    by_id["APP-016"].update({
        "model": "Franke Slim 50 CN / LD370-E01B / 134.0721.210（研究候选）",
        "evidence_reference": "F50 已由 MS 否决；APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818；最终产品未选",
        "status": "部分确认",
        "notes": "Franke 仅为研究候选，非已选／已购；额定输入、总高直径、法兰／开孔、排水／洗碗机支管、控制和拆换净空保持 unknown。",
    })
    by_id["APP-017"].update({
        "evidence_reference": "MS 确认上下顺序、侧置 10A 插座、共用 C16、蝴蝶门和专用墙排方向；准确兼容和墙排接口待外部图",
        "notes": "WG54M7D20W 下＋WQ55M7U20W 上；两个侧置 10A 插座共用一路 C16 RCBO；中厨蝴蝶门。WTZ27510 逐 E-Nr./FD 兼容、准确双接口墙排和项目接口中心仍待确认。",
    })
    write_csv(path, fields, rows)


def update_equipment() -> None:
    path = DECISIONS / "equipment-register.csv"
    fields, rows = read_csv(path)
    by_id = index(rows, "equipment_id")

    app16 = by_id["APP-016"]
    app16.update({
        "manufacturer": "Franke（研究候选）",
        "model": "Slim 50 CN / LD370-E01B / 134.0721.210（研究候选）",
        "variant": "F50 rejected; final product not selected",
        "procurement_status": "candidate",
        "decision_status": "partial",
        "source_ids": merge_ids(app16["source_ids"], "OUTBOUND-FORM-PLUM-20260817;APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818;OUTBOUND-FORM-KITCHEN-20260818"),
        "identity_basis": "MS 否决勒科斯 F50；Franke 中国官方资料支持准确 SKU 和部分性能，但未给完整安装包络",
        "confidence": "1.00",
        "notes": "未选定、未采购；不得把 Franke 研究候选写成最终产品，额定输入、机身高度、法兰、排水和接口中心继续 unknown。",
    })
    app17 = by_id["APP-017"]
    app17.update({
        "source_ids": merge_ids(app17["source_ids"], "OUTBOUND-FORM-PLUM-20260817;OWNER-RESPONSE-HVAC-PLUM-20260818"),
        "identity_basis": "MS 确认精确双机上下顺序、WTZ27510 采用方向、侧置 10A 插座、共用 C16、蝴蝶门和专用墙排；官方资料核实设备参数",
        "notes": "业主方向已确认；WTZ27510 逐 E-Nr./FD 兼容、准确双接口墙排、接口中心、蝴蝶门通风和整机抽出节点仍待外部关闭。",
    })
    for equipment_id in ["APP-017-WASHER", "APP-017-DRYER"]:
        row = by_id[equipment_id]
        row["source_ids"] = merge_ids(row["source_ids"], "OUTBOUND-FORM-PLUM-20260817;OWNER-RESPONSE-HVAC-PLUM-20260818")
        row["notes"] = append_once(
            row["notes"],
            "MS 已确认上下叠放与专用墙排方向；准确墙排和接口中心仍待项目图。",
        )

    san22 = by_id["SAN-022"]
    san22.update({
        "item_name": "600 mm Travertino titan marble island sink with pull-out waste storage",
        "model": "Custom Travertino titan marble sink",
        "variant": "600 mm module; Molteni Eco kits product-language reference",
        "source_ids": merge_ids(san22["source_ids"], "OUTBOUND-FORM-PLUM-20260817;MOLTENI-SISTEMA7-ECOKITS-RESEARCH-20260818;OUTBOUND-FORM-KITCHEN-20260818"),
        "identity_basis": "正式 IFC 证明岛台 SIN01 与 600 mm 协调包络；MS 确认 Travertino titan marble 600 模数石材盆和盆下抽拉垃圾／海绵收纳",
        "notes": "材料、宽度和收纳功能已确认；IFC 其余包围盒不是加工尺寸。Eco kits 只作产品语言参考，盆体、抽屉、桶体、排水、溢水、支撑和检修待联合 shop drawing。",
    })
    vvd = by_id["KIT-VVD-FINISH-001"]
    vvd.update({
        "variant": "project palette; no trolley; Sistema 7 four-door glass wall unit integrated hood direction",
        "source_ids": merge_ids(vvd["source_ids"], "OUTBOUND-FORM-PLUM-20260817;MOLTENI-SISTEMA7-ECOKITS-RESEARCH-20260818;OUTBOUND-FORM-KITCHEN-20260818"),
        "identity_basis": "MS 确认 VVD 材料、无 trolley、无悬挂 hood、Sistema 7 四门玻璃墙柜内整合油烟机；Molteni 官方资料支持产品族可行性",
        "notes": "非已下单；四门柜、油烟机、风管、材料样板、灯光、石材和收纳联合加工图仍待外部关闭。",
    })
    for equipment_id in ["DRAIN-GEB-001", "DRAIN-GEB-002", "DRAIN-GEB-003", "DRAIN-GEB-004", "SAN-001"]:
        row = by_id[equipment_id]
        row.update({
            "decision_status": "superseded",
            "schedule_included": "no",
            "source_ids": merge_ids(row["source_ids"], "OUTBOUND-FORM-PLUM-20260817;OWNER-RESPONSE-HVAC-PLUM-20260818"),
            "human_review_required": "no",
            "notes": append_once(
                row["notes"],
                "MS 于 2026-08-18 确认当前方案不采用；采购／收货和官方参数保留为历史证据，不再生成当前安装要求。",
            ),
        })
    san3 = by_id["SAN-003"]
    san3.update({
        "source_ids": merge_ids(san3["source_ids"], "OWNER-RESPONSE-HVAC-PLUM-20260818"),
        "identity_basis": append_semicolon_once(
            san3["identity_basis"],
            "本轮向 MS 解释其为 Sigma01 双冲水面板而非隐藏阀件",
        ),
        "notes": "115.770.11.5 为正面操作的 Sigma01 白色亮面双冲水面板；面板开口兼作水箱维护入口。当前只证明 1 件收货，第二件和项目完成面节点仍待关闭。",
    })
    write_csv(path, fields, rows)


def update_requirements() -> None:
    path = DECISIONS / "equipment-installation-requirements.csv"
    fields, rows = read_csv(path)
    by_id = index(rows, "requirement_id")

    for req_id in ["REQ-0169", "REQ-0600", "REQ-0601", "REQ-0602"]:
        by_id[req_id].update({
            "source_id": "APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818",
            "source_locator": "最终产品未选；研究候选官方资料不完整",
            "notes": "F50 已否决；Franke 134.0721.210 仅研究候选。最终同型号铭牌／安装图返回前保持 unknown，不从其他 SKU 外推。",
        })
    by_id["REQ-0666"].update({
        "value_origin": "user_input",
        "status": "candidate",
        "source_id": "OUTBOUND-FORM-PLUM-20260817",
        "source_locator": "S02；MS／2026-08-18",
        "notes": "MS 确认采用方向；仍须西门子书面证明与两台实购 E-Nr./FD 兼容，不把业主选择冒充厂家兼容。",
    })
    by_id["REQ-0679"].update({
        "value_text": "direct drain to dedicated laundry wall-drain dryer port",
        "value_origin": "user_input",
        "status": "confirmed",
        "source_id": "OUTBOUND-FORM-PLUM-20260817",
        "source_locator": "S06；MS／2026-08-18",
        "blocks_release": "no",
        "notes": "使用方向已确认；准确墙排产品、标高、接口和防回流／防臭由 APP017-WALL-DRAIN-20260818 关闭。",
    })
    for row in rows:
        if row["equipment_id"] in {"DRAIN-GEB-001", "DRAIN-GEB-002", "DRAIN-GEB-003", "DRAIN-GEB-004", "SAN-001"} and row["status"] in {"pending", "not_applicable"}:
            row.update({
                "value_text": "not_applicable_current_scheme",
                "value_number": "",
                "status": "not_applicable",
                "blocks_release": "no",
                "source_id": "OUTBOUND-FORM-PLUM-20260817",
                "source_locator": "G01-G04、C02；MS／2026-08-18",
                "notes": append_once(
                    row["notes"],
                    "当前方案不采用该组件，本项目发布阻断取消；若未来恢复使用须重新启用全部项目接口门。",
                ),
            })
    by_id["REQ-0511"].update({
        "value_text": "Travertino titan marble",
        "value_origin": "user_input",
        "status": "confirmed",
        "source_id": "OUTBOUND-FORM-PLUM-20260817",
        "source_locator": "V10、C06；MS／2026-08-18",
        "notes": "材料方向已确认；准确荒料、表面和批次仍由 VVD 材料样板门关闭。",
    })
    by_id["REQ-VVD-005"].update({
        "value_text": "pull-out drawers; snack; column doors; glass-door wall unit; side panels; extension doors; no trolley",
        "source_id": "OUTBOUND-FORM-PLUM-20260817",
        "source_locator": "V03、V05、V07；MS／2026-08-18",
        "notes": "产品页 trolley 被本项目明确删除；GlobalId 1C2NYt_qT35P4jTGTCHJyh 是玻璃门柜，不是无门开放柜。",
    })
    by_id["REQ-VVD-006"].update({
        "value_text": "metal lacquer side panels; horizontal finger recess; skirting 60 mm; no suspended hood/support",
        "source_id": "OUTBOUND-FORM-PLUM-20260817",
        "source_locator": "V04；MS／2026-08-18",
        "notes": "删除 2021 产品页悬挂 hood/support；油烟机改由 Sistema 7 四门墙柜项目加工图关闭。",
    })
    by_id["REQ-VVD-008"].update({
        "value_text": "all base units internal LED; Sistema 7 all sides; project lighting designer supplies luminaires to cabinet factory",
        "source_id": "OUTBOUND-FORM-PLUM-20260817",
        "source_locator": "V09；MS／2026-08-18",
        "notes": "只确认范围和供货协作；准确灯具、功率、色温、驱动器、回路、接线和检修继续由 REQ-VVD-011 关闭。",
    })

    additions = [
        ("REQ-HVAC-20260818-001", "HVAC-001", "HVAC", "pipe_insulation_required", "yes", "", "", "official_model_family", "confirmed", "RCP1-HVAC-IFACE-001", "P02010Q p.11", "no", "官方要求冷媒管和冷凝水管保温。"),
        ("REQ-HVAC-20260818-002", "HVAC-001", "HVAC", "pipe_insulation_thickness", "unknown", "", "mm", "pending", "pending", "RCP1-HITACHI-INSULATION-BOUNDARY-20260818", "P02010Q 未给厚度", "yes", "由日立／机电设计按防结露、防火和吊顶净高书面确定。"),
        ("REQ-HVAC-20260818-003", "HVAC-002", "HVAC", "pipe_insulation_required", "yes", "", "", "official_model_family", "confirmed", "RCP1-HVAC-IFACE-002", "P02010Q p.11", "no", "官方要求冷媒管和冷凝水管保温。"),
        ("REQ-HVAC-20260818-004", "HVAC-002", "HVAC", "pipe_insulation_thickness", "unknown", "", "mm", "pending", "pending", "RCP1-HITACHI-INSULATION-BOUNDARY-20260818", "P02010Q 未给厚度", "yes", "由日立／机电设计按防结露、防火和吊顶净高书面确定。"),
        ("REQ-APP017-OWNER-20260818-001", "APP-017", "INT1", "stack_order", "WG54M7D20W lower; WQ55M7U20W upper", "", "", "user_input", "confirmed", "OUTBOUND-FORM-PLUM-20260817", "S01；MS／2026-08-18", "no", "上下顺序已关闭。"),
        ("REQ-APP017-OWNER-20260818-002", "APP-017", "INT1", "cabinet_front", "butterfly doors; full appliance extraction retained", "", "", "user_input", "confirmed", "OUTBOUND-FORM-PLUM-20260817", "S07；MS／2026-08-18", "no", "准确五金、开启包络、通风和收口由柜体加工图关闭。"),
        ("REQ-APP017-OWNER-20260818-003", "APP-017", "PLUM", "dedicated_wall_drain_product", "unknown", "", "", "pending", "pending", "OWNER-RESPONSE-HVAC-PLUM-20260818", "S05-S07；准确产品未给", "yes", "须为洗衣排水和干衣冷凝水给出明确接口、标高和检修。"),
        ("REQ-APP014-OWNER-20260818-001", "APP-014", "ELEC/PLUM/INT1", "socket_service_location", "above appliance in accessible dry service compartment; not directly above wet interfaces", "", "", "user_input", "confirmed", "OUTBOUND-FORM-PLUM-20260817", "K02；MS／2026-08-18", "no", "准确标高、分隔和柜体剖面由 09 表 K02 关闭。"),
        ("REQ-APP016-RESEARCH-20260818-001", "APP-016", "ELEC/PLUM", "research_candidate_sku", "Franke Slim 50 CN / LD370-E01B / 134.0721.210", "", "", "official_exact_model", "candidate", "APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818", "Franke 中国官方页面", "no", "研究候选，不是项目已选。"),
        ("REQ-APP016-RESEARCH-20260818-002", "APP-016", "ELEC", "motor_rating", "", "0.5", "HP", "official_exact_model", "candidate", "APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818", "Franke 中国官方页面", "no", "只属于研究候选；不得换算成项目额定输入功率。"),
        ("REQ-VVD-20260818-001", "KIT-VVD-FINISH-001", "INT1/WFIN", "trolley_included", "no", "", "", "user_input", "confirmed", "OUTBOUND-FORM-PLUM-20260817", "V03；MS／2026-08-18", "no", "删除产品页 trolley。"),
        ("REQ-VVD-20260818-002", "KIT-VVD-FINISH-001", "INT1/HVAC/ELEC", "hood_configuration", "integrated in project Sistema 7 four-door glass wall unit; no suspended VVD hood", "", "", "user_input", "confirmed", "OUTBOUND-FORM-PLUM-20260817", "V04；MS／2026-08-18", "no", "只关闭设计方向。"),
        ("REQ-VVD-20260818-003", "KIT-VVD-FINISH-001", "INT1/HVAC/ELEC", "integrated_hood_shop_drawing", "unknown", "", "", "pending", "pending", "MOLTENI-SISTEMA7-ECOKITS-RESEARCH-20260818", "官方产品族可行；项目配置未知", "yes", "须给四门分格、油烟机、风管、散热隔油、防火净距和拆换路径。"),
        ("REQ-VVD-20260818-004", "KIT-VVD-FINISH-001", "INT1", "glass_door_cabinet_global_id", "1C2NYt_qT35P4jTGTCHJyh; glass doors; not open column", "", "", "site_observed", "confirmed", "MOLTENI-SISTEMA7-ECOKITS-RESEARCH-20260818", "正式 IFC 对象名 Glass Cabinet Frame LP", "no", "只修正协调语义，不修改正式 IFC。"),
        ("REQ-SAN022-20260818-001", "SAN-022", "INT1", "under_sink_storage_function", "pull-out waste bin and sponge/cleaning storage", "", "", "user_input", "confirmed", "OUTBOUND-FORM-PLUM-20260817", "C06；MS／2026-08-18", "no", "功能已确认。"),
        ("REQ-SAN022-20260818-002", "SAN-022", "INT1", "under_sink_storage_product_language", "Molteni Eco kits family reference", "", "", "research_conclusion", "candidate", "MOLTENI-SISTEMA7-ECOKITS-RESEARCH-20260818", "Molteni 官方 Eco kits", "no", "不证明本项目加工尺寸。"),
        ("REQ-SAN022-20260818-003", "SAN-022", "INT1/PLUM", "sink_storage_plumbing_shop_drawing", "unknown", "", "", "pending", "pending", "OUTBOUND-FORM-KITCHEN-20260818", "K05", "yes", "须联合盆体、抽屉、桶体、存水弯、排水、溢水、支撑和检修。"),
        ("REQ-APP009-20260818-001", "APP-009", "INT1/PLUM/ELEC", "island_layout_relationship", "adjacent to APP-010 on one side of SAN-022; existing IFC candidate", "", "", "research_conclusion", "candidate", "OUTBOUND-FORM-KITCHEN-20260818", "K01", "yes", "须由柜体图证明软管、柜孔和检修路径。"),
        ("REQ-APP010-20260818-001", "APP-010", "INT1/PLUM/ELEC", "island_layout_relationship", "adjacent to APP-009; one 600 mm module beyond SAN-022", "", "", "research_conclusion", "candidate", "OUTBOUND-FORM-KITCHEN-20260818", "K01", "yes", "须由柜体图证明软管、柜孔和检修路径。"),
    ]
    for req_id, equipment_id, discipline, parameter_key, value_text, value_number, unit, origin, status, source_id, locator, blocks, notes in additions:
        record = {
            "requirement_id": req_id, "equipment_id": equipment_id, "discipline": discipline, "parameter_key": parameter_key,
            "value_text": value_text, "value_number": value_number, "unit": unit, "datum": "", "value_origin": origin,
            "status": status, "source_id": source_id, "source_locator": locator, "blocks_release": blocks, "notes": notes,
        }
        upsert(rows, "requirement_id", record, fields)
    write_csv(path, fields, rows)


def update_reviews_and_registers() -> None:
    path = DECISIONS / "rcp1-hvac-remodel-review.csv"
    fields, rows = read_csv(path)
    by_id = index(rows, "review_id")
    by_id["RCP1-HVAC-PIPE-001"].update({
        "current_fact": "五个管道对象只表达拆改前／旧方案底稿；P02010Q 要求冷媒管与冷凝水管保温但未给现场保温管厚度；正式 IFC 没有最终改造管线或保温层",
        "candidate_action": "按厂家／机电设计书面厚度生成最终完成外径后，再机械检查吊顶净高、洞口、坡度、碰撞、吊架和检修",
        "basis": "正式 IFC 与旧 blend 分支审计；RCP1-HITACHI-INSULATION-BOUNDARY-20260818",
        "required_input": "A01-A06 项目接口图、冷媒／冷凝水保温材料与厚度、完成外径、最终排放点和安装净距",
        "stop_condition": "不得把旧紫色管线或裸管外径标作最终保温完成路线；厚度与端口未知时不得发布最终路由",
    })
    by_id["RCP1-HVAC-ROUTE-001"].update({
        "current_fact": "装修后的冷媒液管、冷媒气管、冷凝水和全部保温完成外径尚未设计；A01-A04 仅关闭型号族、裸管规格和接管侧方向",
        "basis": "现有五条紫色管线仅为旧方案底稿；OUTBOUND-FORM-HVAC-20260817；RCP1-HITACHI-INSULATION-BOUNDARY-20260818",
        "required_input": "A01-A06 端口、室外机／立管接口、冷凝水排放点、保温厚度和完成外径",
    })
    write_csv(path, fields, rows)

    path = DECISIONS / "c003-fixed-surface-service-origin-reset-targets.csv"
    fields, rows = read_csv(path)
    for row in rows:
        if row["global_id"] == "1C2NYt_qT35P4jTGTCHJyh":
            row["basis"] = "Glass Cabinet Frame LP；正式 IFC 对象名与 MS 回复共同确认这是玻璃门柜，不是无门开放柱；目标仍只用于原点整数化，不能解释为柜体加工尺寸或安装连接点"
    write_csv(path, fields, rows)

    path = DECISIONS / "s701-schedule-review.csv"
    fields, rows = read_csv(path)
    by_id = index(rows, "schedule_id")
    by_id["S701-APP-016"].update({
        "confirmed_scope": "使用位置=西厨水槽柜；勒科斯 F50 已否决",
        "candidate_or_observed_scope": "Franke Slim 50 CN／LD370-E01B／134.0721.210 研究候选",
        "unresolved_for_release": "最终选型、同型号铭牌、额定输入、总高直径、法兰／开孔、排水／洗碗机支管、控制和拆换净空",
        "evidence_reference": "OUTBOUND-FORM-PLUM-20260817;APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818;OUTBOUND-FORM-KITCHEN-20260818",
        "review_status": "research_candidate_project_selection_pending",
    })
    by_id["S701-APP-017"].update({
        "confirmed_scope": "WG54M7D20W 下＋WQ55M7U20W 上；侧置两个 10A 插座共用 C16；中厨蝴蝶门；专用墙排方向",
        "candidate_or_observed_scope": "WTZ27510 带抽板连接件采用方向；650×800×1900 mm 项目协调净空",
        "unresolved_for_release": "WTZ27510 逐 E-Nr./FD 兼容；准确双接口墙排、阀门、接口中心、蝴蝶门通风与整机抽出联合剖面",
        "evidence_reference": merge_ids(by_id["S701-APP-017"]["evidence_reference"], "OUTBOUND-FORM-PLUM-20260817;OWNER-RESPONSE-HVAC-PLUM-20260818"),
        "review_status": "owner_direction_confirmed_vendor_interfaces_pending",
    })
    by_id["S701-DRAIN-CUSTOM-001"].update({
        "confirmed_scope": "主卫和客卫两个淋浴区均采用定制水母地漏／中央集水器＋托克乐思网＋线性排水渠；盖板打开后毛发网可取出",
        "candidate_or_observed_scope": "指定小红书主页商家定制；已登记 Geberit 线性排水组件当前不采用",
        "unresolved_for_release": "逐房间数量和组件型号、项目 shop drawing、排水／防水／完成面／找坡／清洁检修接口；装修宝典准确页面",
        "evidence_reference": "OWNER-PLUM-CUSTOM-DRAIN-20260817;OUTBOUND-FORM-PLUM-20260817;OWNER-RESPONSE-HVAC-PLUM-20260818",
        "review_status": "owner_scope_confirmed_shop_drawing_pending",
    })
    write_csv(path, fields, rows)

    path = DECISIONS / "drawing-register.csv"
    fields, rows = read_csv(path)
    by_id = index(rows, "sheet_number")
    by_id["M-401"]["notes"] = "A01-A04 型号族、裸管规格和接管侧已由 MS 接受；P02010Q 要求冷媒管／冷凝水管保温但未给现场保温管厚度。现有紫色管线只是旧方案底稿且无最终保温层，不能证明模型相符；A05/A06、全部精确端口、保温厚度、完成外径和最终路线仍待日立／机电设计。"
    by_id["P-201"]["notes"] = "Foster 1014850 官方参数已机械关闭；洗烘改为专用双接口墙排方向但准确产品、标高和接口待图；VVD 600 mm 石材盆采用盆下抽拉垃圾／海绵收纳；两台洗碗机保持相邻布局候选；APP-016 F50 已否决、Franke 134.0721.210 仅研究候选。未提供的中心、管径和路线继续 unknown。"
    by_id["P-202"]["notes"] = "主卫和客卫淋浴区均采用定制水母地漏／中央集水器＋托克乐思网＋线性排水渠；已登记 Geberit 154.446.KS.1／154.150.00.1 当前不采用。准确组件、排水口、水封、防水翼环、完成标高、找坡和清洁检修由逐房间 shop drawing 关闭；不修改正式 IFC。"
    by_id["E-303"]["notes"] = append_once(
        by_id["E-303"]["notes"],
        "WS7FSB0C1C 插座采用上方可触及干区方向；准确标高与湿区分隔待柜体剖面。APP-016 F50 已否决，Franke 134.0721.210 仅研究候选，额定输入和电源接口保持 unknown。",
    )
    by_id["I-501"]["notes"] = append_once(
        by_id["I-501"]["notes"],
        "本轮确认无 trolley、无悬挂式 VVD hood；油烟机整合 Sistema 7 四门玻璃墙柜，GlobalId 1C2NYt_qT35P4jTGTCHJyh 按玻璃门柜协调；SAN-022 为 600 mm Travertino titan marble 石材盆并采用盆下抽拉垃圾／海绵收纳。两台 SJ85ZX26MC 保持相邻布局研究候选；所有项目接口由 09 表加工图关闭。",
    )
    by_id["I-503"]["notes"] = append_once(
        by_id["I-503"]["notes"],
        "MS 已确认 WTZ27510 采用方向、中厨蝴蝶门和专用墙排；准确 E-Nr./FD 兼容、双接口墙排、蝴蝶门通风及整机抽出剖面仍待外部关闭。",
    )
    by_id["S-701"]["notes"] = append_once(
        by_id["S-701"]["notes"],
        "F50 已否决；Franke 134.0721.210 只作 APP-016 研究候选。VVD 无 trolley／无悬挂 hood，Sistema 7 柜内油烟机和 Eco kits 盆下收纳待项目加工图；被否决 Geberit 排水组件只保留采购历史。",
    )
    write_csv(path, fields, rows)


def main() -> None:
    update_sources()
    update_owner_inputs()
    update_closeout_rules()
    update_appliance_inputs()
    update_equipment()
    update_requirements()
    update_reviews_and_registers()
    print("Applied 2026-08-18 owner HVAC, plumbing, and kitchen updates.")


if __name__ == "__main__":
    main()
