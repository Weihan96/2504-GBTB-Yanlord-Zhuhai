#!/usr/bin/env python3
"""Reconcile the nine current external Markdown forms back into project SSOT."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "pipeline" / "decisions"
FORMS = ROOT / "output" / "forms" / "对外确认表"

FORM_SOURCES = {
    "OUTBOUND-FORM-DOOR-20260817": (
        "01-门组做法确认表-发全屋定制.md",
        "D01-D13 已完成业主决定记录；不再作为外部填写入口",
        "MS 已确认的 M05/M06/M07 系统、方向、对齐、检修和 Master A 安装面优先级",
        "厂家签认、项目下单尺寸、饰面编号、基层、公差或正式 IFC 已修改",
        "verified_owner_decision_record_vendor_followup_separated",
    ),
    "OUTBOUND-FORM-HVAC-20260817": (
        "02-日立空调接口确认表-发空调厂家.md",
        "A01-A04 型号族、裸管径、官方保温边界及 PC-P1HEQ 通用接线已预填；只保留外部接口",
        "现行表已消费业主回复、P02010Q 和 P02037Q，不再询问官方保温棉厚度",
        "A05/A06 最终型号、准确接口中心、保温计算、项目接线和外部签认",
        "verified_prefilled_confirmation_pending_external_interfaces",
    ),
    "OUTBOUND-FORM-PLUM-20260817": (
        "03-给排水设备确认表-发设备与施工方.md",
        "洗烘、现行吉博力、定制地漏、石材盆及 APP-014 套管边界已预填",
        "旧 MS 回复已拆为 confirmed owner decision、official/research conclusion、external pending 和 unknown",
        "准确双口墙排、项目节点、定制产品加工图或任何接口中心已冻结",
        "verified_prefilled_confirmation_pending_external_interfaces",
    ),
    "OUTBOUND-FORM-GASFIRE-20260817": (
        "04-燃气消防确认表-发主管单位.md",
        "ER9EPA33MP 燃气事实、报警联动决定和消防点位已预填；候选产品保持候选",
        "只向燃气／消防主管方询问准入、兼容、联动、净距和验收",
        "任何候选报警器或探测器已获主管单位批准",
        "verified_prefilled_authority_confirmation_pending",
    ),
    "OUTBOUND-FORM-SITE-20260817": (
        "05-弱电现场记录表-发现场负责人.md",
        "弱电箱约 110 mm、距墙 525 mm、H+350、CAD 坐标、柜格位置、五孔插座及五张照片已预填",
        "只剩准确净宽高、有效净深、柜格 W×H×D、温升检修、线序及设备现场接口",
        "现场准确尺寸、120 分钟温升、线缆通断、DNAKE 型号或 PC-P1HEQ 映射已完成",
        "verified_prefilled_site_record_pending_measurement",
    ),
    "OUTBOUND-FORM-DOOR-VENDOR-20260817": (
        "06-门组厂家复核表-发Rimadesio与Poliform.md",
        "业主既定门组方案转换为厂家是／否复核并要求回项目加工图",
        "厂家问询与 01 业主决定记录分开；不再把加工尺寸变成业主选择题",
        "厂家已签认或项目加工尺寸已冻结",
        "verified_closed_question_vendor_confirmation_pending",
    ),
    "OUTBOUND-FORM-SMART-PANEL-20260817": (
        "07-智能面板电气接口确认表-发JINK与电气方.md",
        "Entry A、Master A、Master B 数量、键序、负载和双控关系已预填",
        "只询问准确 SKU、直接负载／场景语义、接线、负载能力、协议和安装节点",
        "准确产品已选购、接线已签认或现场底盒已复核",
        "verified_prefilled_product_interface_confirmation_pending",
    ),
    "OUTBOUND-FORM-JOINERY-DETAIL-20260817": (
        "08-定制家具与墙脚节点确认表-发全屋定制.md",
        "服务塔候选尺寸和干区踢脚参考候选状态已预填",
        "服务塔与踢脚线均未冒充已下单或已采用；只保留设计、加工图、样板及真正业主偏好",
        "节点尺寸已冻结、样板已批准或业主已决定采用踢脚候选",
        "verified_prefilled_candidate_detail_confirmation_pending",
    ),
    "OUTBOUND-FORM-KITCHEN-20260818": (
        "09-厨房设备与柜体深化确认表-发橱柜设备方.md",
        "VVD 材料、Foster 参数、准确油烟机、厨房墙面、设备电气语义及候选边界已预填",
        "只保留项目加工图、服务路径、厂家准确 APP-016 证据和联合节点",
        "外部方已签认、项目加工接口已冻结或候选产品已最终选择",
        "verified_prefilled_kitchen_shopdrawing_confirmation_pending",
    ),
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def update_by_id(rows: list[dict[str, str]], key: str, item_id: str, values: dict[str, str]) -> None:
    for row in rows:
        if row[key] == item_id:
            row.update(values)
            return
    raise KeyError(f"missing {item_id} in {key}")


def append_if_missing(rows: list[dict[str, str]], key: str, row: dict[str, str]) -> None:
    if not any(existing[key] == row[key] for existing in rows):
        rows.append(row)


def upsert_full(rows: list[dict[str, str]], key: str, row: dict[str, str]) -> None:
    for index, existing in enumerate(rows):
        if existing[key] == row[key]:
            rows[index] = {field: row.get(field, "") for field in existing}
            return
    rows.append(row)


def merge_ids(value: str, *extra: str) -> str:
    result = [item for item in value.split(";") if item]
    for item in extra:
        if item and item not in result:
            result.append(item)
    return ";".join(result)


def append_sentence_once(value: str, sentence: str) -> str:
    """Normalize a generated note so rerunning reconciliation is idempotent."""
    base = value
    while sentence in base:
        base = base.replace(sentence, "")
    return " ".join(part for part in (base.strip(), sentence) if part)


def form_hash(name: str) -> str:
    return hashlib.sha256((FORMS / name).read_bytes()).hexdigest()


def reconcile_source_register() -> None:
    path = DECISIONS / "source-evidence-register.csv"
    fields, rows = read_csv(path)
    rows = [row for row in rows if row["source_id"] != "EXT-EXTERNAL-INFO-MINIMUM-20260817"]
    for source_id, (name, evidence, proves, does_not_prove, status) in FORM_SOURCES.items():
        update_by_id(rows, "source_id", source_id, {
            "source_document": name,
            "local_path": f"output/forms/对外确认表/{name}",
            "sha256": form_hash(name),
            "locator": "现行 Markdown；逐行已知信息／当前状态／只缺什么／责任方／回复附件",
            "evidence": evidence,
            "proves": proves,
            "does_not_prove": does_not_prove,
            "status": status,
            "notes": "现行唯一可编辑源；历史问询仅保留为证据，不得作为当前填写入口。",
        })

    historical_requests = {
        "A104-SHOP-REQUEST-20260815": "output/forms/对外确认表/06-门组厂家复核表-发Rimadesio与Poliform.md",
        "RCP1-PLUM-REQUEST-20260815": "output/forms/对外确认表/02-日立空调接口确认表-发空调厂家.md 与 03-给排水设备确认表-发设备与施工方.md",
        "GAS-CONSULTATION-TEMPLATE-001": "output/forms/对外确认表/04-燃气消防确认表-发主管单位.md",
        "A106-CONSULTATION-PACK-20260815": "output/forms/对外确认表/04-燃气消防确认表-发主管单位.md",
    }
    for source_id, current_entry in historical_requests.items():
        update_by_id(rows, "source_id", source_id, {
            "status": "superseded_historical_evidence_current_form_linked",
            "notes": f"保留为历史索取／咨询证据；现行填写入口仅为 {current_entry}。",
        })

    history_name = "E304-现场最短取证清单-20260815.md"
    update_by_id(rows, "source_id", "E304-SITE-CHECKLIST-20260815", {
        "source_document": history_name,
        "local_path": f"drawings/evidence/{history_name}",
        "sha256": hashlib.sha256((ROOT / "drawings" / "evidence" / history_name).read_bytes()).hexdigest(),
        "locator": "历史任务与现行 05 表边界说明",
        "evidence": "保留 2026-08-15 现场任务历史；明确弱电箱已有证据并指向现行 05 表",
        "proves": "历史文件不是当前填写入口，已知约测、CAD、照片和剩余实测边界均被保留",
        "does_not_prove": "任何剩余现场尺寸、温升、通断、设备接口已经完成",
        "status": "superseded_historical_evidence_current_form_linked",
        "notes": "现行填写入口仅为 output/forms/对外确认表/05-弱电现场记录表-发现场负责人.md。",
    })

    sleeve_path = ROOT / "drawings" / "evidence" / "OWNER-APP-014-净水管可抽换套管候选-20260817.md"
    sleeve = {field: "" for field in fields}
    sleeve.update({
        "source_id": "OWNER-APP014-REPLACEABLE-SLEEVE-20260817",
        "discipline": "PLUM/INT1",
        "sheet_id": "P-201/I-501/S-701",
        "decision_scope": "APP-014 连续可抽换净水管套管候选",
        "source_kind": "owner_reference_candidate_boundary",
        "source_document": sleeve_path.name,
        "local_path": str(sleeve_path.relative_to(ROOT)),
        "sha256": hashlib.sha256(sleeve_path.read_bytes()).hexdigest(),
        "locator": "业主候选方向、官方资料边界与外部关闭条件",
        "evidence": "业主希望预埋连续可抽换套管并在装修后穿 PE 管；网络案例规格不是项目规格",
        "proves": "套管作为项目候选方向；西门子说明书已给 2–3 m 服务距离、约 3 m 随机管及不小于 80 mm 柜侧孔",
        "does_not_prove": "随机管准确材料／外径、允许套管、项目套管规格、弯曲半径、路线、防水节点或抽换可行性",
        "status": "verified_owner_candidate_external_compatibility_pending",
        "confidence": "1.00",
        "review_required": "yes",
        "formal_ifc_write_allowed": "no",
        "manufacturer": "Siemens",
        "model_scope": "WS7FSB0C1C",
        "revision": "2026-08-17",
        "publication_date": "2026-08-17",
        "notes": "不得把网络案例中的四分铝塑套管和二分 PE 管写成项目规格。",
    })
    append_if_missing(rows, "source_id", sleeve)

    rcp_dir = ROOT / "drawings" / "evidence" / "OWNER-RCP-天花深化参考-20260817"
    rcp_sources = [
        {
            "source_id": "OWNER-RCP-WOOD-SLAT-PDF-20260817",
            "source_kind": "owner_design_reference_candidate",
            "source_document": "kitchen-wood-slat-ceiling-reference.pdf",
            "locator": "2 页 A3 厨房窄条木吊顶概念及构造参考",
            "evidence": "窄条木纹模块、黑色设备槽、可拆检修路径及厨房材料样板的设计方向",
            "proves": "RCP1／A-106 可深化研究的业主参考意向",
            "does_not_prove": "项目已采用、准确模数、槽宽、基层、防火、风量、检修或任何施工做法已签认",
            "status": "verified_reference_candidate_not_construction_detail",
            "notes": "只作设计意向与深化问题来源；不得据此修改正式 IFC 或冻结风口／吊顶节点。",
        },
        {
            "source_id": "OWNER-RCP-REFERENCE-README-20260817",
            "source_kind": "evidence_bundle_index_and_boundary",
            "source_document": "README.md",
            "locator": "天花参考目录说明、三条稳定小红书链接及证据边界",
            "evidence": "区分公开案例、业主设计意向、厂家资料、shop drawing 和现场证据",
            "proves": "该目录中资料的现行证据边界及进入施工深化前必须关闭的事项",
            "does_not_prove": "任何案例适配性、产品性能、准确尺寸、项目加工图或已批准施工节点",
            "status": "verified_evidence_boundary_index",
            "notes": "后续选定厂家时另行登记正式资料和项目节点。",
        },
    ]
    for data in rcp_sources:
        source_path = rcp_dir / data["source_document"]
        record = {field: "" for field in fields}
        record.update({
            "discipline": "RCP1/INT1",
            "sheet_id": "RCP1/A-106",
            "decision_scope": "厨房窄条木吊顶、包梁、风口与设备槽深化参考",
            "local_path": str(source_path.relative_to(ROOT)),
            "sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "revision": "2026-08-18",
            "publication_date": "2026-08-18",
        })
        record.update(data)
        upsert_full(rows, "source_id", record)
    write_csv(path, fields, rows)


def reconcile_electrical_evidence() -> None:
    path = DECISIONS / "elec-source-evidence.csv"
    fields, rows = read_csv(path)
    source = next(row for row in read_csv(DECISIONS / "source-evidence-register.csv")[1]
                  if row["source_id"] == "OWNER-APP014-REPLACEABLE-SLEEVE-20260817")
    item = {field: "" for field in fields}
    item.update({
        "evidence_id": source["source_id"],
        "discipline": source["discipline"],
        "sheet_id": source["sheet_id"],
        "decision_scope": source["decision_scope"],
        "source_kind": source["source_kind"],
        "source_document": source["local_path"],
        "source_sha256": source["sha256"],
        "source_locator": source["locator"],
        "evidence": source["evidence"],
        "proves": source["proves"],
        "does_not_prove": source["does_not_prove"],
        "status": source["status"],
        "confidence": source["confidence"],
        "review_required": source["review_required"],
        "formal_ifc_write_allowed": "no",
        "notes": source["notes"],
    })
    append_if_missing(rows, "evidence_id", item)
    write_csv(path, fields, rows)


def reconcile_owner_inputs() -> None:
    path = DECISIONS / "owner-input-register.csv"
    fields, rows = read_csv(path)
    weak_sources = "OWNER-INBOX-20260814-001;E304-CAD-001;E304-CAD-002;E304-PHOTO-001;E304-PHOTO-002;E304-PHOTO-003;E304-PHOTO-004;E304-PHOTO-005;OUTBOUND-FORM-SITE-20260817"
    update_by_id(rows, "input_id", "E302-ENTRY-SIDE", {
        "status": "自定义确认",
        "user_value": "Entry A；功能和键序已定，准确安装面、标高与底盒仍待现场",
        "notes": "不再询问 Entry A／B 选择；只由 05／07 表关闭现场安装与产品接口。",
    })
    update_by_id(rows, "input_id", "E304-CABINET-DIMENSIONS", {
        "question": "弱电箱和所在柜格尚缺哪些准确净尺寸",
        "candidate_value": "现有约测净深 110 mm；距入户墙 525 mm；右下柜格；底边 H+350 mm；CAD XY 约 (4600.016,-735.369) mm；箱内五孔插座",
        "user_value": "已纠正 525 mm；现有约测净深约 110 mm，待带尺复核",
        "status": "需证据",
        "evidence_reference": weak_sources,
        "source_basis": "业主约测、开发商 E-2／CAD 机械换算和五张现场照片",
        "notes": "只缺准确净宽、净高、插头和网线弯曲后的有效净深、所在柜格 W×H×D 及开门检修净空。",
    })
    update_by_id(rows, "input_id", "E304-CABINET-VENTILATION", {
        "status": "需证据",
        "evidence_reference": "E304-PHOTO-001;E304-PHOTO-002;E304-PHOTO-003;E304-USER-002;OUTBOUND-FORM-SITE-20260817",
        "notes": "通风孔和现状箱内运行已确认；只缺柜门关闭 0/60/120 min 温升、告警掉线、散热与拆换空间。",
    })
    update_by_id(rows, "input_id", "E304-CABLE-CONTINUITY", {
        "status": "需证据",
        "evidence_reference": "E304-PHOTO-005;OUTBOUND-FORM-SITE-20260817",
        "notes": "照片标签候选为客厅、次卧、主卧；只缺逐根 1–8 线序、两端标签、端口和 PoE 记录。",
    })
    update_by_id(rows, "input_id", "RCP1-HVAC-PORTS", {
        "question": "日立 A01-A06 仍缺哪些项目接口、保温计算和外部签认",
        "candidate_value": "A01-A04 型号族／裸管径／接管侧及保温必需性已知；A05 候选、A06 占位；准确接口中心保持 unknown",
        "user_value": "MS 已接受 A01-A04 型号族、裸管规格和接管侧；官方资料未给保温材料、导热系数、防火或厚度",
        "status": "需证据",
        "notes": "机电设计按珠海露点、管温、材料导热系数、防火和净空计算保温厚度／完成外径，再由日立／安装方签认兼容；旧 IFC 紫色管线不得冒充最终带保温管径。",
    })
    update_by_id(rows, "input_id", "INT1-ENTRY-PARCEL-FUNCTION", {
        "question": "入户是否需要临时置物以便取件时腾出双手",
        "candidate_value": "玄关柜台面、抽拉板、托盘或换鞋凳的一部分；不是专用快递柜",
        "user_value": "入户设计要考虑拿快递；需求已确认",
        "notes": "功能已关闭，不再问业主；只需在 I-503 形成时验证不挡门、通道、Entry 面板和强弱电箱检修。",
    })
    update_by_id(rows, "input_id", "INT1-ENTRY-PARCEL-LAYOUT", {
        "question": "入户临时置物位如何与玄关平立面和检修范围协调",
        "candidate_value": "在玄关台面、抽拉板、托盘或换鞋凳中整合；I-503 前只验证避让关系",
        "source_basis": "业主功能确认；准确形式由 I-503 和必要现场证据关闭",
        "notes": "当前不要求门开／门关／柜门开三种状态填写整套宽深高；不是专用快递柜。",
    })
    update_by_id(rows, "input_id", "E303-NS02-FORM", {"status": "采用候选"})
    update_by_id(rows, "input_id", "APP004-COFFEE-FINAL", {
        "blocks_release": "no",
        "question": "咖啡机最终采购时确认 GS3 或 E1 Prima EXP 的实购版本",
        "candidate_value": "两款均可安装；装修按 2600 W 不利值、独立 C16 回路和可选给排水预留，不需现在二选一",
        "user_value": "咖啡机为选装；现在不二选一",
        "status": "采用候选",
        "notes": "选购前再匹配插头面板和准确接口；不得作为当前装修 release blocker。",
    })

    sleeve = {field: "" for field in fields}
    sleeve.update({
        "input_id": "APP014-REPLACEABLE-SLEEVE",
        "workstream": "P-201/I-501/S-701",
        "priority": "P1",
        "blocks_release": "yes",
        "question": "APP-014 连续可抽换净水管套管如何形成可施工节点",
        "candidate_value": "连续可抽换套管，装修完成后再穿 PE 净水管；网络案例规格不作为项目规格",
        "user_value": "采用连续可抽换套管作为候选方向",
        "unit": "mm",
        "status": "需证据",
        "evidence_reference": "OWNER-APP014-REPLACEABLE-SLEEVE-20260817;APP-014-SIEMENS-MANUAL-001;OUTBOUND-FORM-PLUM-20260817;OUTBOUND-FORM-KITCHEN-20260818",
        "source_basis": "业主候选方向与西门子官方服务距离／柜孔边界",
        "sync_target": "APP-014;P-201;I-501;S-701",
        "notes": "关闭证据：随机管材料／外径和允许套管的厂家确认；实际路线、弯曲半径、端部防水固定、渗漏可见性及全程抽换样板。",
    })
    append_if_missing(rows, "input_id", sleeve)
    write_csv(path, fields, rows)


def reconcile_closeout_rules() -> None:
    path = DECISIONS / "owner-input-closeout-rules.csv"
    fields, rows = read_csv(path)
    update_by_id(rows, "input_id", "E304-CABINET-DIMENSIONS", {
        "responsible_party": "现场负责人",
        "required_evidence": "05 表 N01-N05：准确净宽、净高、带尺复核净深、线缆弯曲后有效净深、柜格 W×H×D 和开门检修净空；每项附带尺照片",
        "notes": "约测 110 mm、距墙 525 mm、H+350、CAD XY、右下柜格和五孔插座已知，不得重新按完全未知询问。",
    })
    update_by_id(rows, "input_id", "E304-CABINET-VENTILATION", {
        "responsible_party": "现场负责人/弱电施工方",
        "required_evidence": "05 表柜门关闭 0/60/120 min 全负载温升、设备／PoE 负载、告警掉线、通风孔和检修照片",
    })
    update_by_id(rows, "input_id", "RCP1-HVAC-PORTS", {
        "responsible_party": "机电设计/日立技术方/安装方",
        "required_evidence": "机电设计按珠海露点、管温、材料导热系数、防火和净空计算保温厚度及完成外径；日立／安装方签认兼容；另提交 A05/A06 型号和 A01-A06 项目接口图",
        "notes": "P02010Q 只证明现场提供保温管，未给材料或厚度；旧 IFC 紫色管线不是最终带保温外径。",
    })
    update_by_id(rows, "input_id", "INT1-ENTRY-PARCEL-FUNCTION", {
        "responsible_party": "业主",
        "required_evidence": "MS 已确认入户需要临时置物以便取件时腾出双手",
        "automatic_close_allowed": "yes",
        "notes": "功能已关闭，不再询问。",
    })
    update_by_id(rows, "input_id", "INT1-ENTRY-PARCEL-LAYOUT", {
        "responsible_party": "室内设计/现场",
        "required_evidence": "I-503 平立面形成时证明候选台面、抽拉板、托盘或换鞋凳不挡入户门、通道、Entry 面板及强弱电箱检修；必要照片或简图",
        "notes": "I-503 前不强制填写门开／门关三套宽深高；不是专用快递柜。",
    })
    sleeve = {field: "" for field in fields}
    sleeve.update({
        "input_id": "APP014-REPLACEABLE-SLEEVE",
        "closeout_kind": "manufacturer_and_mockup_evidence",
        "responsible_party": "西门子技术方/给排水施工/橱柜方",
        "required_evidence": "随机管材料／外径及套管兼容签认；项目路线、弯曲半径、端部防水固定、渗漏可见节点和全程抽换样板",
        "automatic_close_allowed": "no",
        "notes": "网络案例尺寸不是项目规格；无准确资料时接口和套管尺寸保持 unknown。",
    })
    append_if_missing(rows, "input_id", sleeve)
    write_csv(path, fields, rows)


def reconcile_equipment() -> None:
    path = DECISIONS / "equipment-register.csv"
    fields, rows = read_csv(path)
    updates = {
        "CTRL-ENTRY-A": {
            "model": "JINK EGG family（研究候选；准确 SKU 未定）",
            "variant": "minimum 4 functions; direct-load/scene implementation pending",
            "identity_basis": "业主确认 Entry A 为一个至少 4 键逻辑面板及完整键序；JINK EGG 仅为候选产品族",
            "notes": "已确认由近至远：客厅双控、书房双控、餐厅双控、照明总控；总控不切连续负载。准确 SKU、继电器、接线、负载和底盒待外部签认。",
        },
        "CTRL-MASTER-A": {
            "model": "JINK EGG family（研究候选；准确 SKU 未定）",
            "variant": "minimum 2 functions; ambient LED then accent spot; direct-load/scene pending",
            "notes": "一个至少 2 键面板，氛围 LED 在前、重点射灯在后；直接负载或场景实现待产品／电气复核，安装坐标保持 unknown。",
        },
        "CTRL-MASTER-B": {
            "model": "JINK EGG family（研究候选；准确 SKU 未定）",
            "variant": "minimum 3 functions; true wired two-way required",
            "notes": "一个至少 3 键面板：客厅、书房、餐厅；必须与 Entry A 形成断网仍可用的真实物理接线双控，不以场景替代。",
        },
    }
    for equipment_id, values in updates.items():
        update_by_id(rows, "equipment_id", equipment_id, values)
    for row in rows:
        if row["equipment_id"] == "APP-014":
            row["source_ids"] = merge_ids(row["source_ids"], "OWNER-APP014-REPLACEABLE-SLEEVE-20260817")
            row["notes"] = append_sentence_once(
                row["notes"],
                "连续可抽换净水管套管为候选方向，厂家兼容、规格、路线和抽换样板未关闭。",
            )
    write_csv(path, fields, rows)


def reconcile_requirements() -> None:
    path = DECISIONS / "equipment-installation-requirements.csv"
    fields, rows = read_csv(path)
    additions = [
        ("REQ-APP014-SLEEVE-20260818-001", "replaceable_service_sleeve_direction", "continuous_pull-through sleeve; install PE tube after fit-out", "user_input", "candidate", "no", "业主候选方向；不是已批准施工规格。"),
        ("REQ-APP014-SLEEVE-20260818-002", "supplied_tube_material_and_od", "unknown", "pending", "pending", "yes", "由西门子准确型号资料或书面回复关闭。"),
        ("REQ-APP014-SLEEVE-20260818-003", "sleeve_compatibility_and_bend_radius", "unknown", "pending", "pending", "yes", "不得用网络案例尺寸代填。"),
        ("REQ-APP014-SLEEVE-20260818-004", "sleeve_route_end_waterproof_and_leak_visibility", "unknown", "pending", "pending", "yes", "由 P-201／I-501 联合节点关闭。"),
        ("REQ-APP014-SLEEVE-20260818-005", "full_route_pull_through_mockup", "required", "project_candidate", "pending", "yes", "施工前完成全程抽换样板并留照片。"),
    ]
    for req_id, requirement, text_value, basis, status, blocks, notes in additions:
        item = {field: "" for field in fields}
        item.update({
            "requirement_id": req_id,
            "equipment_id": "APP-014",
            "discipline": "PLUM/INT1",
            "parameter_key": requirement,
            "value_text": text_value,
            "value_origin": basis,
            "status": status,
            "source_id": "OWNER-APP014-REPLACEABLE-SLEEVE-20260817",
            "source_locator": "owner candidate and evidence boundary",
            "blocks_release": blocks,
            "notes": notes,
        })
        if any(row["requirement_id"] == req_id for row in rows):
            update_by_id(rows, "requirement_id", req_id, item)
        else:
            rows.append(item)
    write_csv(path, fields, rows)


def main() -> None:
    reconcile_source_register()
    reconcile_electrical_evidence()
    reconcile_owner_inputs()
    reconcile_closeout_rules()
    reconcile_equipment()
    reconcile_requirements()
    print(json.dumps({"status": "ok", "forms": len(FORM_SOURCES)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
