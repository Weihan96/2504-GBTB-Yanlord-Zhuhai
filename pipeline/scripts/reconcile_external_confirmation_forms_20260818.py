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
        "开发商既有 A01-A05 型号、QD/QDF 镜像关系、裸管径、华美 Class 1 现场保温规格、回风口检修和 PC-P1HEQ 已预填",
        "现行表已消费业主回复、P02010Q、P02037Q、现场铭牌及保温管照片；不再询问采购、A05 型号、统一保温厚度或面板兼容性",
        "A06 与现场实物的机位绑定、项目路线／支吊架／控制图、未提供的接口中心及按图施工和竣工记录",
        "verified_prefilled_confirmation_pending_external_interfaces",
    ),
    "OUTBOUND-FORM-PLUM-20260817": (
        "03-给排水设备确认表-发设备与施工方.md",
        "洗烘双存水弯候选、现行吉博力、定制地漏、石材盆、APP-014 套管及三款已购龙头证据边界已预填",
        "旧 MS 回复已拆成已由业主决定、已有官方资料、仍需外部回答和确实未知四类，不再把机器状态词显示给收件人",
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
        "弱电箱位置、约 110 mm 深度、三江品牌、五孔插座、有效安装包络原则、风道／风扇预留、TP-Link AP 和 DNAKE 280M-S3 已预填",
        "只剩准确净宽高、插线后有效深度、柜格 W×H×D、固定／进线／通风节点、线序、实际供电及 PC-P1HEQ 追线",
        "任何现场净尺寸、网线通断、散热节点、DNAKE 实际端子／供电或 PC-P1HEQ 五点映射已完成",
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
        "客卫台盆左侧全高服务塔和干区齐平同色宽踢脚＋上下双阴影缝的项目结论已预填",
        "位置、功能和双阴影缝方向已关闭；只保留项目出图、加工图、五金／型材信息与 1:1 样板复核",
        "准确加工尺寸、型材、基层、公差、清洁背衬或样板已经批准",
        "verified_project_direction_detail_confirmation_pending",
    ),
    "OUTBOUND-FORM-KITCHEN-20260818": (
        "09-厨房设备与柜体深化确认表-发橱柜设备方.md",
        "VVD 材料、Foster 参数、准确油烟机、厨房墙面、设备电气、小家电／已购工具收纳、冰淇淋机包络及 APP-016 比选边界已预填",
        "只保留项目加工图、服务路径、实物收纳尺寸、F50／替代品准确资料报价和联合节点",
        "外部方已签认、项目加工接口已冻结或候选产品已最终选择",
        "verified_prefilled_kitchen_shopdrawing_confirmation_pending",
    ),
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    for row in rows:
        overflow = row.pop(None, None)
        if not overflow:
            continue
        if path.name == "source-evidence-register.csv" and row.get("source_id") == "E304-USER-003":
            # A pre-existing row had one extra semantic cell before does_not_prove,
            # which shifted all following values. Normalize it once before rewrite.
            row["proves"] = f'{row["proves"]}；{row["does_not_prove"]}'
            row["does_not_prove"] = row["status"]
            row["status"] = row["confidence"]
            row["confidence"] = row["review_required"]
            row["review_required"] = row["formal_ifc_write_allowed"]
            row["formal_ifc_write_allowed"] = row["manufacturer"]
            row["manufacturer"] = row["model_scope"]
            row["model_scope"] = row["revision"]
            row["revision"] = row["publication_date"]
            row["publication_date"] = row["legacy_targets"]
            row["legacy_targets"] = row["legacy_projection_json"]
            row["legacy_projection_json"] = row["notes"]
            row["notes"] = overflow[0]
            continue
        raise ValueError(f"{path}: unexpected extra CSV columns in {row.get(fields[0], '<unknown>')}")
    return fields, rows


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


def source_record(fields: list[str], **values: str) -> dict[str, str]:
    record = {field: "" for field in fields}
    record.update(values)
    return record


def reconcile_source_register() -> None:
    path = DECISIONS / "source-evidence-register.csv"
    fields, rows = read_csv(path)
    rows = [row for row in rows if row["source_id"] != "EXT-EXTERNAL-INFO-MINIMUM-20260817"]
    for source_id, (name, evidence, proves, does_not_prove, status) in FORM_SOURCES.items():
        update_by_id(rows, "source_id", source_id, {
            "source_document": name,
            "local_path": f"output/forms/对外确认表/{name}",
            "sha256": form_hash(name),
            "locator": "现行 Markdown；逐行说明已经知道什么／现在要做什么／由谁回答／在哪里回复",
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

    local_evidence = [
        (
            "OWNER-HVAC-INSULATION-PHOTOS-20260819",
            "RCP1",
            "M-401/RCP1",
            "现场华美 Class 1 保温管规格",
            "owner_provided_site_photo_set",
            "HVAC-INS-20260819-P01-20x15.jpg",
            "P01/P02/P04/P06/P07/P13/P14/P17/P19 九张照片；部分印字须旋转 180° 阅读",
            "现场存在 6×15、10×15、13×9、13×15、16×9、16×15、16×20、20×15、25×15、32×9、32×15 mm 多种 ID×TK 保温套",
            "现场并非统一一种保温厚度；ID 为内径，TK 为单边壁厚",
            "各规格与 A01～A06 每一段管线的逐段对应、接缝连续性或最终施工路线",
            "verified_site_product_and_size_set_mapping_pending",
        ),
        (
            "OWNER-PLUM-GEBERIT-DOUBLE-TRAP-20260819",
            "PLUM/INT1",
            "P-201/I-503/S-701",
            "洗衣机与热泵干衣机双独立存水弯候选",
            "seller_product_image_owner_candidate",
            "PLUM-GEBERIT-TAOBAO-703917654224-SKU-5125125938318-double-trap.webp",
            "淘宝 item 703917654224／SKU 5125125938318 商家尺寸图",
            "商家提供双独立存水弯组合及部分尺寸；图中明确部分连接件不是吉博力原厂件",
            "业主选定的优先比选 SKU 和双独立水封方向",
            "吉博力原厂系统身份、与两台西门子设备的兼容、最终标高、管径、防串水或项目安装批准",
            "verified_seller_candidate_external_signoff_pending",
        ),
        (
            "OWNER-WIN-CRANK-RETROFIT-20260820",
            "A/DET1",
            "A-104/D-601/S-701",
            "既有窗手摇开窗器改造现场证据",
            "owner_site_evidence_note",
            "OWNER-WIN-手摇开窗器旧窗改造案例-20260820.md",
            "旧窗改造案例、22 mm 闭合边和厨房下轨固定照片",
            "现有窗框具备约 22 mm 闭合边及可研究的下轨固定位置",
            "手摇开窗器作为既有窗改造方向并进入项目节点深化",
            "最终产品、推拉力、开启角、固定基层、防水或厂家项目适配已经签认",
            "verified_owner_direction_detail_pending",
        ),
        (
            "OWNER-WIN-CRANK-PHOTO-01-20260820",
            "A/DET1",
            "A-104/D-601/S-701",
            "既有窗 22 mm 闭合边现场测量",
            "owner_site_photo",
            "OWNER-WIN-手摇开窗器现场-20260820/01-下轨与22mm闭合边测量.jpg",
            "现场照片；卷尺读数与窗框闭合边",
            "照片显示既有窗框约 22 mm 闭合边及下轨局部条件",
            "手摇开窗器安装包络研究的现场几何依据",
            "最终产品、固定基层、防水、受力、准确连续可用长度或厂家适配",
            "verified_site_dimension_detail_pending",
        ),
        (
            "OWNER-WIN-CRANK-PHOTO-02-20260820",
            "A/DET1",
            "A-104/D-601/S-701",
            "厨房既有手摇器下轨固定现场照片",
            "owner_site_photo",
            "OWNER-WIN-手摇开窗器现场-20260820/02-厨房手摇器下轨固定.jpg",
            "现场照片；厨房既有手摇器及下轨固定位置",
            "既有厨房窗可观察到手摇器下轨固定实例",
            "手摇开窗器改造方向的现场案例依据",
            "其他窗位可直接复制、最终产品、受力、防水或厂家项目适配",
            "verified_site_precedent_detail_pending",
        ),
    ]
    for source_id, discipline, sheet_id, scope, kind, name, locator, evidence, proves, does_not_prove, status in local_evidence:
        local_path = ROOT / "drawings" / "evidence" / name
        record = source_record(
            fields,
            source_id=source_id,
            discipline=discipline,
            sheet_id=sheet_id,
            decision_scope=scope,
            source_kind=kind,
            source_document=name,
            local_path=str(local_path.relative_to(ROOT)),
            sha256=hashlib.sha256(local_path.read_bytes()).hexdigest(),
            locator=locator,
            evidence=evidence,
            proves=proves,
            does_not_prove=does_not_prove,
            status=status,
            confidence="1.00",
            review_required="yes",
            formal_ifc_write_allowed="no",
            revision="2026-08-20",
            publication_date="2026-08-20",
            notes="由现行对外结论表消费；不得越过证据边界写入施工接口中心。",
        )
        upsert_full(rows, "source_id", record)

    hvac_photo_names = [
        "HVAC-INS-20260819-P01-20x15.jpg",
        "HVAC-INS-20260819-P02-16x9-13x9.jpg",
        "HVAC-INS-20260819-P04-10x15.jpg",
        "HVAC-INS-20260819-P06-6x15-16x20.jpg",
        "HVAC-INS-20260819-P07-32x9.jpg",
        "HVAC-INS-20260819-P13-16x15.jpg",
        "HVAC-INS-20260819-P14-13x15.jpg",
        "HVAC-INS-20260819-P17-25x15.jpg",
        "HVAC-INS-20260819-P19-32x15.jpg",
    ]
    for index, name in enumerate(hvac_photo_names, start=1):
        local_path = ROOT / "drawings" / "evidence" / name
        record = source_record(
            fields,
            source_id=f"OWNER-HVAC-INS-PHOTO-{index:02d}-20260819",
            discipline="RCP1",
            sheet_id="M-401/RCP1",
            decision_scope="现场华美 Class 1 保温管规格",
            source_kind="owner_provided_site_photo",
            source_document=name,
            local_path=str(local_path.relative_to(ROOT)),
            sha256=hashlib.sha256(local_path.read_bytes()).hexdigest(),
            locator=f"现行 02 表保温规格逐行配图；{name}",
            evidence="现场管套品牌和 ID×TK 印字；倒置文字按现行 02 表注明方向读取",
            proves="该照片中可读的现场保温管规格存在",
            does_not_prove="该规格与某一机位或整条管路逐段对应，或最终施工验收完成",
            status="verified_site_photo_mapping_pending",
            confidence="1.00",
            review_required="yes",
            formal_ifc_write_allowed="no",
            manufacturer="Huamei／华美",
            model_scope="Class 1 Rubber Foam",
            revision="2026-08-19",
            publication_date="2026-08-19",
            notes="汇总语义见 OWNER-HVAC-INSULATION-PHOTOS-20260819。",
        )
        upsert_full(rows, "source_id", record)

    live_sources = [
        source_record(
            fields,
            source_id="HUAMEI-CLASS1-OFFICIAL-20260819",
            discipline="RCP1",
            sheet_id="M-401/RCP1",
            decision_scope="华美 Class 1 橡塑保温产品性能和执行标准",
            source_kind="manufacturer_official_webpage",
            source_document="华美 Class 1 Rubber Foam 官方页",
            source_url="https://www.huameiworld.com/rubber-foam/class-1-rubber-foam.html",
            sha256="not_applicable_live_reference",
            locator="官方产品页；2026-08-19 读取",
            evidence="0℃ 平均导热系数≤0.034 W/(m·K)；燃烧性能 B1；列明 GB/T 6343、2406、8627、10294、17146、17794、8811、10808、6669、7762、16259 及 GB 8624-2012",
            proves="现场照片所示产品族的官方性能与标准边界",
            does_not_prove="本项目逐段采用规格、最终保温厚度或施工验收已经完成",
            status="registered_official_live_reference",
            confidence="0.95",
            review_required="yes",
            formal_ifc_write_allowed="no",
            manufacturer="Huamei／华美",
            model_scope="Class 1 Rubber Foam",
            revision="2026-08-19",
            publication_date="2026-08-19",
            notes="现场实物规格另由 OWNER-HVAC-INSULATION-PHOTOS-20260819 固化。",
        ),
        source_record(
            fields,
            source_id="TPLINK-XAP1500GE-OFFICIAL-20260819",
            discipline="NETWORK/ELEC",
            sheet_id="E-304/RCP1",
            decision_scope="嵌入式吸顶 AP 候选产品",
            source_kind="manufacturer_official_webpage",
            source_document="TP-Link TL-XAP1500GE-PoE/DC 易展版官方页",
            source_url="https://www.tp-link.com.cn/product_3166.html",
            sha256="not_applicable_live_reference",
            locator="官方产品页；2026-08-19 读取",
            evidence="Wi-Fi 6 AX1500；Ø184×40 mm；开孔 Ø155 mm；单千兆口；802.3at PoE；PoE 最大功耗 11.1 W",
            proves="候选 AP 的准确型号、外形、开孔、网络和 PoE 参数",
            does_not_prove="已经采购、现场网线通断、覆盖、控制器／网关或吊顶机械适配已经关闭",
            status="verified_official_candidate",
            confidence="1.00",
            review_required="yes",
            formal_ifc_write_allowed="no",
            manufacturer="TP-Link",
            model_scope="TL-XAP1500GE-PoE/DC 易展版",
            revision="2026-08-19",
            publication_date="2026-08-19",
            notes="替代此前 Huawei／Ruijie 优先研究方向；仍保持 candidate。",
        ),
        source_record(
            fields,
            source_id="DNAKE-280M-S3-OFFICIAL-20260819",
            discipline="NETWORK/ELEC/INT1",
            sheet_id="E-304/I-503",
            decision_scope="既有狄耐克可视对讲室内机型号识别",
            source_kind="manufacturer_official_webpage_with_site_photo_match",
            source_document="DNAKE 280M-S3 官方产品页",
            source_url="https://www.dnake-global.com/10-1-inch-linux-based-indoor-monitor-280m-s3-product/",
            sha256="not_applicable_live_reference",
            locator="官方产品页与业主现场正面照片外观匹配",
            evidence="10.1 英寸 Linux 室内机；约 270×168×15 mm；支持 802.3af PoE 或 DC 12 V／2 A",
            proves="现有室内机按 280M-S3 识别并可用于安装包络；型号不再重复询问",
            does_not_prove="本项目实际供电、背部端子、物业系统兼容、保留／迁移或最终坐标",
            status="research_identification_confirmed_site_interfaces_pending",
            confidence="0.95",
            review_required="yes",
            formal_ifc_write_allowed="no",
            manufacturer="DNAKE／狄耐克",
            model_scope="280M-S3",
            revision="2026-08-19",
            publication_date="2026-08-19",
            notes="实际供电和端子继续由现场复核；照片中的网络地址不得投影。",
        ),
    ]
    for record in live_sources:
        upsert_full(rows, "source_id", record)

    legacy_path = ROOT / "drawings" / "evidence" / "OWNER-LEGACY-精装房改造深化说明-20251216-对账.md"
    legacy_record = source_record(
        fields,
        source_id="OWNER-LEGACY-DESIGN-BRIEF-20251216",
        discipline="ARCH/INT1/WFIN/RCP1/ELEC/PLUM",
        sheet_id="A-104/A-105/A-106/I-501/I-502/I-503/I-504/D-601/E-303/M-401/S-701",
        decision_scope="2025-12 旧版精装房改造设计方向及方案演变",
        source_kind="legacy_owner_design_brief_audit",
        source_document=legacy_path.name,
        local_path=str(legacy_path.relative_to(ROOT)),
        sha256=hashlib.sha256(legacy_path.read_bytes()).hexdigest(),
        locator="46 页旧 Pages 临时导出逐页对账；原文件 SHA-256 见文内",
        evidence="旧稿功能分区、墙体、卫浴、材料、天花、空调、电气、设备和收口方向；逐项标明已吸收、同向深化、已替代和仍需项目深化",
        proves="项目早期设计意图和后续方案演变，可补充当前图纸深化输入",
        does_not_prove="旧稿示例产品已经选定／采购、图片尺寸可施工、厂家已批准、正式 IFC 已修改或旧决定优先于现行 SSOT",
        status="verified_legacy_brief_reconciled_current_ssot_controls",
        confidence="1.00",
        review_required="yes",
        formal_ifc_write_allowed="no",
        revision="2025-12-16",
        publication_date="2025-12-16",
        notes="原 .pages 文件保留在项目根目录但不作为现行填写入口；原文件 SHA-256 为 70ad2df264cdad9b75b9715ef4c2982ff0c4d4105fafb7d817826730718e47c9。",
    )
    upsert_full(rows, "source_id", legacy_record)

    update_by_id(rows, "source_id", "APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818", {
        "sha256": hashlib.sha256((ROOT / "drawings" / "evidence" / "APP-016-Franke-Slim50-CN研究-20260818.md").read_bytes()).hexdigest(),
        "evidence": "F50 保留为基准候选；Franke 中国准确 SKU 134.0721.210 的官方公开参数可核验，尺寸图确认安装面以下总高 365 mm",
        "proves": "Franke Slim 50 CN 可作为与 F50 比较的准确型号候选；365 mm 总高可进入项目联合剖面",
        "does_not_prove": "F50 被否决、Franke 已选／已购、最终额定输入、法兰、洗碗机支管、拆换净空、报价或项目适配已经关闭",
        "status": "official_comparison_candidate_project_section_and_selection_pending",
        "notes": "业主已纠正：F50 没有被否决；替代品只有在品质更好且价格可接受时才考虑。",
    })
    update_by_id(rows, "source_id", "OWNER-RESPONSE-HVAC-PLUM-20260818", {
        "sha256": hashlib.sha256((ROOT / "drawings" / "evidence" / "OWNER-RESPONSE-HVAC-PLUM-20260818.md").read_bytes()).hexdigest(),
        "notes": "原始逐项回答仍以 OUTBOUND-FORM-HVAC-20260817 和 OUTBOUND-FORM-PLUM-20260817 为准；后续业主纠正以现行 SSOT 和确认表为准。",
    })
    update_by_id(rows, "source_id", "OWNER-INT1-BATHG-TISSUE-WASTE-20260817", {
        "proves": "客卫台盆左侧全高服务塔的项目方向；400 mm 优选、350–450 mm 可调、250–280 mm 深候选控制值",
        "does_not_prove": "现场完成面尺寸、准确加工尺寸、五金、内胆、给排水避让或样板已经批准",
        "status": "owner_project_direction_confirmed_shop_drawing_pending",
        "notes": "不再作为外露架／服务塔二选一；项目先出 I-502，外部按图复核。",
    })
    update_by_id(rows, "source_id", "OWNER-WFIN-BASEBOARD-SHADOW-GAP-20260817", {
        "proves": "干区采用墙面同色齐平宽踢脚、上下约 10 mm 双阴影缝及墙门柜连续通缝的项目方向",
        "does_not_prove": "准确高度、厚度、材料、型材、公差、清洁背衬、湿区适用或样板已经批准",
        "status": "owner_project_direction_confirmed_detail_and_mockup_pending",
        "notes": "业主确认保留原参考做法；中厨 VVD Pewter 60 mm 踢脚和湿区另行处理。",
    })
    update_by_id(rows, "source_id", "OWNER-E304-DNAKE-PHOTO-20260817", {
        "does_not_prove": "单凭照片证明本项目实际供电、背部端子、物业系统兼容、保留迁移接入方案或最终坐标",
        "status": "verified_brand_software_and_280m_s3_visual_match",
        "model_scope": "280M-S3",
        "notes": "结合 DNAKE-280M-S3-OFFICIAL-20260819 识别准确型号；照片中的网络地址不投影。",
    })
    remote_owner_sources = [
        {
            "source_id": "OWNER-KITCHEN-H70FT-20260820",
            "discipline": "INT1/ELEC",
            "sheet_id": "I-501/E-303/S-701",
            "decision_scope": "惠人 H70FT 已购状态与岛台下存放决定",
            "source_kind": "owner_notion_purchase_and_storage_confirmation",
            "source_document": "业主 Notion 小家电存放安排＋Codex task 2026-08-20",
            "source_url": "https://app.notion.com/p/weihanshen/278fda08ab6c800884eecac88b3d871b?source=copy_link",
            "sha256": "not_applicable_owner_notion_live_page",
            "locator": "岛台下分组；业主补充准确型号 H70FT",
            "evidence": "H70FT 已购买；存放在岛台下；设备在台面使用而不是在关闭柜格内运行",
            "proves": "准确型号、采购事实和存放位置决定",
            "does_not_prove": "已经到货、实物 W×D×H、拆件高度、重量、电源线、附件包络或柜体加工尺寸",
            "status": "confirmed_purchase_storage_dimensions_pending",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "manufacturer": "HUROM／惠人",
            "model_scope": "H70FT",
            "revision": "2026-08-20",
            "publication_date": "2026-08-20",
            "notes": "购买不等于到货；存放位置不等于使用位置。",
        },
        {
            "source_id": "OWNER-KITCHEN-TOOLS-PURCHASED-20260820",
            "discipline": "INT1",
            "sheet_id": "I-501/S-701",
            "decision_scope": "八项已购厨房工具与 K-D01～K-D04 抽屉分组",
            "source_kind": "owner_notion_purchase_list_and_project_storage_decision",
            "source_document": "业主 Notion 已购厨房工具清单＋09 厨房结论表",
            "source_url": "https://app.notion.com/p/weihanshen/3c2fda08ab6c80ccad27e0786345dc42?source=copy_link",
            "sha256": "not_applicable_owner_notion_live_page",
            "locator": "8 项已购工具；K-D01～K-D04 收纳结论",
            "evidence": "8 项工具均按已购买处理；项目已决定长工具、酒具、锅具和易碎咖啡器具四类收纳",
            "proves": "采购事实、数量口径和抽屉功能分组",
            "does_not_prove": "已经到货、每件准确订单 SKU、完整实物包络、抽屉加工尺寸、导轨或摆样已经通过",
            "status": "confirmed_purchase_storage_groups_measurement_pending",
            "confidence": "1.00",
            "review_required": "yes",
            "formal_ifc_write_allowed": "no",
            "revision": "2026-08-20",
            "publication_date": "2026-08-20",
            "notes": "未取得可靠尺寸的项目只缺实物测量，不重新开放收纳分组决定。",
        },
    ]
    for values in remote_owner_sources:
        upsert_full(rows, "source_id", source_record(fields, **values))
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
    update_by_id(rows, "input_id", "E302-HVAC-CONTROL-PANELS", {
        "question": "5 个既有 PC-P1HEQ 控制点分别对应 A01–A06 的哪台或哪组室内机，以及装修后如何保留、迁移或合并",
        "candidate_value": "PC-P1HEQ 是开发商原配并在既有日立系统中运行，型号与兼容性已关闭；只待五点物理映射和项目控制策略",
        "user_value": "精装房原配面板，不再询问采购或兼容性；项目先画五点控制关系，现场追线并记录实际机组",
        "status": "需证据",
        "evidence_reference": "OWNER-HVAC-PC-P1HEQ-20260817;HITACHI-PC-P1HEQ-P02037Q-001;OUTBOUND-FORM-HVAC-20260817;build/elec/elec-developer-control-reference.json",
        "source_basis": "业主确认开发商原配；官方说明书确认通用安装；开发商图证明 5 个既有参考点",
        "notes": "兼容性和准确型号均已关闭。只由项目／现场关闭五点控制对象、保留／迁移／合并、端子图线长和最终位置；不得按相邻机位猜测系统归属。",
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
        "question": "弱电箱及玄关柜怎样预留被动风道和可加装风扇节点",
        "candidate_value": "低位进风＋高位排风；高位预留 120 mm 级低速直流排风扇的加固可拆安装面、低压电源、温控器与防护网",
        "user_value": "先把散热解决方案和固定／检修条件设计周全；温升只用于验收和决定是否立即装风扇",
        "status": "需证据",
        "evidence_reference": "E304-PHOTO-001;E304-PHOTO-002;E304-PHOTO-003;E304-USER-002;OUTBOUND-FORM-SITE-20260817",
        "source_basis": "业主要求先设计失效安全的散热节点；现场照片证明现状箱盖有通风孔",
        "notes": "不再把关门温升作为前置设计问题；只缺有效进排风口、加固安装板、电源／温控、拆换路线和实机安装后的验收记录。",
    })
    update_by_id(rows, "input_id", "E304-CABLE-CONTINUITY", {
        "status": "需证据",
        "evidence_reference": "E304-PHOTO-005;OUTBOUND-FORM-SITE-20260817",
        "notes": "照片标签候选为客厅、次卧、主卧；只缺逐根 1–8 线序、两端标签、端口和 PoE 记录。",
    })
    update_by_id(rows, "input_id", "RCP1-HVAC-PORTS", {
        "question": "日立 A01-A06 仍缺哪些项目接口、保温计算和外部签认",
        "candidate_value": "A01/A04=RPIZ-22FSLN5QD/P；A02/A03=RPIZ-22FSLN5QDF/P 镜像；A05=RPIZ-71FSLN5QD/P；A06 机位已定但现场型号映射待施工开放遮挡后记录；接口中心未由说明书提供",
        "user_value": "全部为开发商随精装交付的既有设备；正式 IFC 的 QD/QDF 镜像正确；回风口兼作检修口；项目画路线，安装后按支吊点标高验收",
        "status": "需证据",
        "evidence_reference": "HITACHI-P02010Q-001;OWNER-HVAC-INSULATION-PHOTOS-20260819;HUAMEI-CLASS1-OFFICIAL-20260819;OUTBOUND-FORM-HVAC-20260817",
        "source_basis": "日立官方安装资料、现场铭牌／保温管照片、华美官方资料和业主设计决定",
        "notes": "不再问采购、A05 型号、面板兼容或统一保温厚度。只缺 A06 机位映射、项目 M-401/RCP1/E-302、说明书未给的接口中心和按图施工／竣工记录；旧 IFC 紫色管线不得冒充最终带保温外径。",
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
    update_by_id(rows, "input_id", "APP016-DISPOSER-DATA", {
        "question": "F50 与可接受替代品如何完成高度、性能、价格和接口比选",
        "candidate_value": "F50 保留为基准候选；Franke Slim 50 CN 为高度 365 mm 的比较候选；两者均未最终选定",
        "user_value": "没有否决 F50；继续寻找品质更好且价格可接受的方案",
        "status": "需证据",
        "evidence_reference": "OUTBOUND-FORM-PLUM-20260817;APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818;OUTBOUND-FORM-KITCHEN-20260818",
        "source_basis": "业主纠正＋Franke 官方研究资料",
        "notes": "由项目先做含 Foster 200 mm 盆深的联合剖面；外部只提供候选准确铭牌、完整尺寸／接口图、报价和保修。不得写成 F50 已否决。",
    })
    update_by_id(rows, "input_id", "E304-AP-POWER", {
        "question": "主卧和次卧嵌入式吸顶 AP 如何供电并接回弱电箱",
        "candidate_value": "TP-Link TL-XAP1500GE-PoE/DC 易展版；Cat6 星型回弱电箱；802.3at PoE；AP 点不设 220 V",
        "user_value": "采用 TP-Link 嵌入式款式作为优先方向",
        "status": "采用候选",
        "evidence_reference": "TPLINK-XAP1500GE-OFFICIAL-20260819;OUTBOUND-FORM-SITE-20260817",
        "source_basis": "业主选型方向＋TP-Link 官方规格",
        "notes": "准确版本、控制器／网关、吊顶开孔与拆换空间、Cat6 通断和 PoE 预算仍由 E-304／RCP1 深化关闭；不是已采购。",
    })
    update_by_id(rows, "input_id", "E304-VIDEO-INTERCOM-DOORBELL", {
        "question": "现有 DNAKE 280M-S3 的实际供电、端子和物业系统保留／复装条件是什么",
        "candidate_value": "现有室内机按 DNAKE 280M-S3 识别；约 270×168×15 mm；官方支持 802.3af PoE 或 DC 12 V／2 A",
        "user_value": "品牌和型号已由照片／官方产品页匹配关闭；不再拆机读取背面铭牌",
        "status": "需证据",
        "evidence_reference": "OWNER-E304-DNAKE-PHOTO-20260817;DNAKE-280M-S3-OFFICIAL-20260819;OUTBOUND-FORM-SITE-20260817",
        "source_basis": "业主现场照片与狄耐克官方产品页的外观／规格匹配",
        "notes": "只缺背部实际端子与供电、物业系统兼容、装修期间拆下复装条件和最终安装坐标；照片中的 IP/MAC 不投影。",
    })
    update_by_id(rows, "input_id", "INT1-BATHG-TISSUE-WASTE-SCHEME", {
        "question": "客卫干区如何设置抽纸和隐藏垃圾收纳",
        "candidate_value": "台盆左侧全高服务塔；优选面宽 400 mm，可调 350–450 mm；深 250–280 mm、上限 300 mm；柜体不超过台盆前缘",
        "user_value": "采用服务塔方向，由项目先设计再交全屋定制复核",
        "status": "自定义确认",
        "notes": "功能、位置和优先尺寸已关闭；不得再作为业主二选一。项目先出 I-502 平立剖，外部只回加工图与五金／检修复核。",
    })
    update_by_id(rows, "input_id", "WFIN-DRY-BASEBOARD-SHADOW-GAP-SCOPE", {
        "question": "卧室、走廊等干区踢脚采用什么设计方向",
        "candidate_value": "墙面同色齐平宽踢脚；上下各约 10 mm 阴影缝；门套、隐形门和柜体边界连续通缝",
        "user_value": "保留原参考做法",
        "status": "自定义确认",
        "notes": "方向已关闭，不再提供 A/B 方案。中厨 VVD Pewter 60 mm 踢脚保持独立；湿区不套用；准确截面、型材、公差和清洁背衬由 D-601 与 1:1 样板关闭。",
    })
    update_by_id(rows, "input_id", "APP017-SIEMENS-STACK", {
        "candidate_value": "WG54M7D20W + WQ55M7U20W；有效净空 650×800×1900 mm；双独立存水弯候选 SKU 5125125938318",
        "user_value": "洗衣机下、干衣机上；WTZ27510 采用方向；两个侧置 10A 插座共用 C16；蝴蝶门；洗衣和干衣冷凝水分别进独立存水弯后下游合流",
        "evidence_reference": merge_ids(next(row["evidence_reference"] for row in rows if row["input_id"] == "APP017-SIEMENS-STACK"), "OWNER-PLUM-GEBERIT-DOUBLE-TRAP-20260819"),
        "notes": "双存水弯拓扑和优先 SKU 已关闭为项目候选；仍缺商家确认两只存水弯准确货号、非吉博力配件、接口、包络、防串水、标高、检修和双机排水试验。",
    })
    newly_absorbed_inputs = [
        {
            "input_id": "KITCHEN-H70FT-STORAGE",
            "workstream": "INT1/ELEC",
            "priority": "P1",
            "blocks_release": "yes",
            "question": "已购惠人 H70FT 如何在岛台下存放并在台面安全使用",
            "candidate_value": "岛台下直立格或全拉出托盘；柜内只存放；台面使用位置设可触及 10A 接地插座",
            "user_value": "准确型号 H70FT，已经购买，存放位置为岛台下",
            "unit": "mm",
            "status": "需证据",
            "evidence_reference": "OWNER-KITCHEN-H70FT-20260820;OUTBOUND-FORM-KITCHEN-20260818",
            "source_basis": "业主 Notion 存放安排、业主型号／采购确认和 09 表项目结论",
            "sync_target": "APP-020;I-501;E-303;S-701",
            "notes": "采购、型号和存放位置已关闭；只缺到货实测、取放包络、托盘承重、防污垫、湿附件通风和使用位置插座。",
        },
        {
            "input_id": "KITCHEN-TOOLS-DRAWERS",
            "workstream": "INT1",
            "priority": "P1",
            "blocks_release": "yes",
            "question": "八项已购厨房工具如何进入 K-D01～K-D04 抽屉加工图",
            "candidate_value": "K-D01 长工具；K-D02 Coravin 酒具；K-D03 奶锅／锅盖；K-D04 HARIO 易碎器具",
            "user_value": "八项均已购买；收纳分组按 09 表执行，不交给橱柜方重新分配",
            "unit": "mm",
            "status": "需证据",
            "evidence_reference": "OWNER-KITCHEN-TOOLS-PURCHASED-20260820;OUTBOUND-FORM-KITCHEN-20260818",
            "source_basis": "业主 Notion 已购清单和项目收纳决定",
            "sync_target": "KTOOL-001..008;I-501;S-701",
            "notes": "只缺实物 SKU／包络、摆样、抽屉内净尺寸、导轨和承重；不重新开放四类收纳决定。",
        },
        {
            "input_id": "INT1-WINDOW-TREATMENTS",
            "workstream": "INT1/ELEC/A-104",
            "priority": "P1",
            "blocks_release": "yes",
            "question": "全屋已确认窗饰方向如何形成逐窗报价和安装图",
            "candidate_value": "公共区 Luminette K5-501；卧室 Silhouette N31-204；厨卫 25 mm 灰黑铝百叶；夜帘 A／B 两案比价",
            "user_value": "日帘和百叶产品家族／颜色方向确认；夜帘保持 LightLock 优先与 Duolite 对照，尚未终选或下单",
            "unit": "mm",
            "status": "需证据",
            "evidence_reference": "OWNER-INT1-WINDOW-TREATMENTS-20260818",
            "source_basis": "业主窗饰设计方向和亨特道格拉斯官方产品研究",
            "sync_target": "WINTR-001..004;A-104;A-106;E-303;I-504;S-701",
            "notes": "只缺逐窗完成面复尺、分幅、收拢／操作侧、轨道、实体色卡、电源控制、报价和厂家安装图。",
        },
        {
            "input_id": "WIN-CRANK-RETROFIT",
            "workstream": "A-104/DET1",
            "priority": "P1",
            "blocks_release": "yes",
            "question": "既有窗批量加装手摇开窗器如何完成产品与样板验证",
            "candidate_value": "机械手摇款；参考厨房既有下轨固定；其余窗约 22 mm 闭合边仅作现场几何依据",
            "user_value": "采用手摇机械款方向，不增加电源、控制线或智能联动",
            "unit": "mm",
            "status": "需证据",
            "evidence_reference": "OWNER-WIN-CRANK-RETROFIT-20260820;OWNER-WIN-CRANK-PHOTO-01-20260820;OWNER-WIN-CRANK-PHOTO-02-20260820",
            "source_basis": "业主现场测量、厨房既有实例和业主机械款决定",
            "sync_target": "WIN-OPENER-001;A-104;D-601;S-701",
            "notes": "只缺准确产品、推拉力、开启角、孔距、基层／排水腔避让、防水保修和一樘非厨房窗样板。",
        },
        {
            "input_id": "E303-STRONGBOX-IDENTITY",
            "workstream": "E-303/I-503",
            "priority": "P2",
            "blocks_release": "no",
            "question": "既有强电配电箱铭牌能关闭哪些现状信息",
            "candidate_value": "三江电气；手写型号 P230、箱号 PX1；铭牌 380/220 V、63 A",
            "user_value": "按现场照片记录，不把 63 A 当成断路器整定或可增容结论",
            "unit": "A",
            "status": "自定义确认",
            "evidence_reference": "E303-PHOTO-001;E303-PHOTO-002;ELEC-BOX-EVIDENCE-20260820",
            "source_basis": "业主现场铭牌与正面照片",
            "sync_target": "ELEC-DB-EXISTING-001;E-303;I-503",
            "notes": "P230／PX1 为手写识读；如用于订货或更换须近距离复核。",
        },
        {
            "input_id": "E303-STRONGBOX-CAPACITY",
            "workstream": "E-303",
            "priority": "P0",
            "blocks_release": "yes",
            "question": "装修新增回路前如何核实现有配电箱真实容量和回路余量",
            "candidate_value": "铭牌 63 A 只描述箱体；不代表总开关整定、导线容量或剩余模数",
            "user_value": "由电气方打开箱门记录，不凭外壳铭牌判断增容",
            "unit": "A",
            "status": "需证据",
            "evidence_reference": "E303-PHOTO-001;E303-PHOTO-002;ELEC-BOX-EVIDENCE-20260820",
            "source_basis": "现场铭牌证据边界",
            "sync_target": "ELEC-DB-EXISTING-001;E-303",
            "notes": "只缺总开关／分路断路器型号整定、回路标签、导线截面、剩余模数和新增回路负荷计算。",
        },
    ]
    for item in newly_absorbed_inputs:
        upsert_full(rows, "input_id", {field: item.get(field, "") for field in fields})

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
    legacy_inputs = [
        {
            "input_id": "LEGACY-RCP-CEILING-ACCESS",
            "workstream": "RCP1/A-106/DET1",
            "priority": "P1",
            "blocks_release": "yes",
            "question": "木饰面／木条天花及拼缝隐藏检修如何形成可施工节点",
            "candidate_value": "厨房木饰面拼板；卫浴／阳台木条拼缝隐藏整板检修；其余双层石膏板；保留现有天花高低关系",
            "user_value": "旧稿设计方向继续保留；不把旧稿龙骨间距和图片做法直接当施工规格",
            "unit": "mm",
            "status": "需证据",
            "evidence_reference": "OWNER-LEGACY-DESIGN-BRIEF-20251216;OWNER-RCP-WOOD-SLAT-PDF-20260817;RCP1-CEILING-COVERING-001",
            "source_basis": "旧稿 p.14–15、21–23 与现行 RCP1 证据",
            "sync_target": "RCP1;A-106;D-601;S-701",
            "notes": "关闭材料防火防潮、基层／龙骨、固定、模数、液压五金、设备拆出路线和 1:1 样板；正式 IFC 不在本批修改。",
        },
        {
            "input_id": "LEGACY-RCP-CURTAIN-COVE-WINERACK",
            "workstream": "RCP1/INT1/ELEC",
            "priority": "P1",
            "blocks_release": "yes",
            "question": "灯槽、200 mm 客厅窗帘盒和约 1300 mm 酒架上方加固如何联合深化",
            "candidate_value": "反灯槽直边铝型材收口；客厅窗帘盒 200 mm；顶天立地酒架上方按最终满载和固定做结构加固",
            "user_value": "旧稿方向保留",
            "unit": "mm",
            "status": "需证据",
            "evidence_reference": "OWNER-LEGACY-DESIGN-BRIEF-20251216",
            "source_basis": "旧稿 p.15–16",
            "sync_target": "RCP1;A-106;I-501;E-301;S-701",
            "notes": "窗帘盒须与最终窗帘系统、电源、灯具、空调风口和检修协调；酒架不得只固定在饰面板。",
        },
        {
            "input_id": "LEGACY-INT1-BAYWINDOW-WOOD",
            "workstream": "INT1/WFIN/DET1",
            "priority": "P2",
            "blocks_release": "yes",
            "question": "客厅和次卧飘窗木作如何与大白墙、窗边防水和检修收口",
            "candidate_value": "次卧飘窗落地浅色木；客厅飘窗木盒与大白墙脱缝并齐平墙面",
            "user_value": "旧稿方向保留，准确节点待项目深化",
            "unit": "mm",
            "status": "需证据",
            "evidence_reference": "OWNER-LEGACY-DESIGN-BRIEF-20251216",
            "source_basis": "旧稿 p.12–13",
            "sync_target": "I-504;D-601;S-701",
            "notes": "关闭防潮、日晒、基层、伸缩、检修、窗边防水和准确完成面。",
        },
        {
            "input_id": "LEGACY-WFIN-BALCONY-DECK",
            "workstream": "WFIN/A-105/DET1",
            "priority": "P2",
            "blocks_release": "yes",
            "question": "阳台塑木地板方向如何满足完成面、排水和耐候要求",
            "candidate_value": "阳台塑木地板；准确产品和节点未定",
            "user_value": "旧稿材料方向保留",
            "unit": "mm",
            "status": "需证据",
            "evidence_reference": "OWNER-LEGACY-DESIGN-BRIEF-20251216",
            "source_basis": "旧稿 p.12",
            "sync_target": "A-105;D-601;S-701",
            "notes": "关闭产品、完成面高度、排水、检修、收边、防火和耐候；旧稿不是采购证据。",
        },
        {
            "input_id": "LEGACY-M401-L-SUPPLY-AIR",
            "workstream": "M-401/RCP1",
            "priority": "P1",
            "blocks_release": "yes",
            "question": "客厅 L 形转角送风口能否在满足性能和检修前提下实现",
            "candidate_value": "用 L 形转角风口弱化窗帘盒阳角；侧出风＋下回风，回风口兼检修口",
            "user_value": "旧稿造型方向保留；项目先计算并画图",
            "unit": "mm",
            "status": "需证据",
            "evidence_reference": "OWNER-LEGACY-DESIGN-BRIEF-20251216;OUTBOUND-FORM-HVAC-20260817",
            "source_basis": "旧稿 p.20 与现行 HVAC 结论",
            "sync_target": "M-401;RCP1;A-106",
            "notes": "关闭风量、静压、噪声、有效开口、型材、转角压损、风管和检修；不得只按效果图定尺寸。",
        },
    ]
    for item in legacy_inputs:
        upsert_full(rows, "input_id", {field: item.get(field, "") for field in fields})
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
        "responsible_party": "项目设计/全屋定制/弱电施工方",
        "required_evidence": "先提交低位进风、高位排风、120 mm 级风扇加固可拆安装面、低压电源／温控和拆换路线节点；实机安装后再用 0/60/120 min 温度、负载、告警掉线和照片验收",
        "notes": "温升不是设计前置问题；即使首次测试正常，也必须预留失效安全的被动风道和可加装风扇条件。",
    })
    update_by_id(rows, "input_id", "RCP1-HVAC-PORTS", {
        "responsible_party": "机电设计/日立技术方/安装方",
        "required_evidence": "项目先完成 M-401/RCP1/E-302，逐段标裸管、华美 Class 1 套管规格／完成外径、冷凝水路线坡度、回风口检修与五点控制映射；施工开放遮挡后补 A06 环境＋铭牌对应；日立／安装方只按图复核并提交竣工记录",
        "notes": "A01-A05、QD/QDF 镜像、既有采购状态、PC-P1HEQ 兼容和现场存在的多种保温规格已关闭；接口中心未由说明书提供时继续不猜。",
    })
    update_by_id(rows, "input_id", "E302-HVAC-CONTROL-PANELS", {
        "responsible_party": "机电设计/日立安装方/现场",
        "required_evidence": "5 个既有点位逐点控制对象和保留／迁移／合并表；项目端子图与线长校核；最终安装位置、端子和机组对应照片",
        "notes": "PC-P1HEQ 是开发商原配并在既有日立系统中运行，准确型号与兼容性已关闭；不再向厂家询问型号或兼容性。",
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
    update_by_id(rows, "input_id", "APP016-DISPOSER-DATA", {
        "responsible_party": "室内设计/橱柜方/设备供货方",
        "required_evidence": "项目先提交 Foster 838×500×200 mm 水槽柜联合剖面和允许的设备包络；供货方再为 F50 与替代候选提交准确铭牌、总高／机身高、法兰、排水／电气接口、报价和保修；橱柜方按图答可实施／冲突",
        "notes": "F50 未否决；Franke Slim 50 CN 高 365 mm 只作比较候选。不得让外部单位在没有项目剖面的情况下代做布局。",
    })
    update_by_id(rows, "input_id", "E304-AP-POWER", {
        "responsible_party": "网络设计/弱电施工/现场",
        "required_evidence": "按 TP-Link TL-XAP1500GE-PoE/DC 的 Ø184×40 mm、Ø155 mm 开孔和 802.3at PoE 参数完成吊顶开孔／拆换剖面、Cat6 通断、交换机端口和总 PoE 预算",
        "notes": "TP-Link 为项目优先候选，非已购；不再重复比较 Huawei／Ruijie。",
    })
    update_by_id(rows, "input_id", "E304-VIDEO-INTERCOM-DOORBELL", {
        "responsible_party": "现场负责人/物业门禁维护方/弱电设计",
        "required_evidence": "DNAKE 280M-S3 背部实际端子与线缆照片；实际供电为 802.3af PoE 或 DC 12 V 的记录；物业确认原位保留、装修期间拆下及复装条件；最终 I-503/E-304 安装定位",
        "notes": "准确型号已通过现场照片与官方产品页匹配关闭，不再要求拆旧箱或重复询问型号；网络地址不外发。",
    })
    update_by_id(rows, "input_id", "INT1-BATHG-TISSUE-WASTE-SCHEME", {
        "responsible_party": "室内设计",
        "required_evidence": "I-502 按台盆左侧全高服务塔、优选 400 mm 面宽（350–450 mm 可调）、250–280 mm 深（上限 300 mm）完成平立剖并证明不超过台盆前缘",
        "automatic_close_allowed": "yes",
        "notes": "业主设计决定已关闭；不再让全屋定制或业主二选一。",
    })
    update_by_id(rows, "input_id", "WFIN-DRY-BASEBOARD-SHADOW-GAP-SCOPE", {
        "responsible_party": "室内设计",
        "required_evidence": "业主已确认保留原参考：干区墙面同色齐平宽踢脚，上下各约 10 mm 阴影缝，门套／隐形门／柜体连续通缝",
        "automatic_close_allowed": "yes",
        "notes": "方向已关闭；中厨 VVD 60 mm 与湿区另行处理。",
    })
    missing_existing_rules = [
        ("E304-CABINET-IDENTITY", "owner_and_site_evidence", "业主/项目内部", "三江电气正面照片和业主确认：箱内无铭牌，当前不为读取背面铭牌拆箱；现状旧箱型号不阻塞按有效安装包络选择新箱", "yes", "现状品牌与旧型号取证边界已关闭。"),
        ("PLUM-PURCHASED-FAUCETS-20260820", "owner_purchase_evidence", "业主/项目内部", "Notion 交易成功订单截图及三款所选外观；各 1 件", "yes", "只关闭已购状态、数量和订单所示款式，不证明原厂身份或安装接口。"),
        ("PLUM-PURCHASED-FAUCETS-INSTALL-20260820", "seller_and_project_installation_evidence", "供货方/给排水设计/全屋定制", "MF287／CZ356／CZ028 准确 SKU、单位、阀体包络、埋深、完成面基准、接口、孔径／中心距和房间映射；P-202/I-502 项目图后供货方按图复核", "no", "MF287 与 CZ028 现有卖家图只关闭系统关系；CZ356 安装图仍缺。"),
        ("KITCHEN-H70FT-STORAGE", "site_measurement_and_project_detail", "室内设计/橱柜深化/电气设计/现场", "H70FT 到货实测 W×D×H、拆件高度、重量、电源线和附件包络；I-501/E-303 取放、承重、通风、防污和使用插座节点", "no", "采购、型号和岛台下存放位置已关闭。"),
        ("KITCHEN-TOOLS-DRAWERS", "site_measurement_and_project_detail", "室内设计/橱柜深化/五金/现场", "八项实物 SKU 和包络；K-D01～K-D04 摆样、抽屉内净尺寸、分隔、导轨、承重和取放空间", "no", "四类收纳分组已关闭，不由橱柜方重新分配。"),
        ("INT1-WINDOW-TREATMENTS", "site_measurement_product_quote_and_shopdrawing", "室内设计/窗饰供应安装/电气设计", "逐窗完成面复尺、分幅、操作与收拢、转角、头轨／侧轨、实体色卡、电源控制、分项报价和厂家安装图", "no", "日帘／百叶方向已确认；夜帘 A／B 两案尚未终选或下单。"),
        ("WIN-CRANK-RETROFIT", "product_interface_and_mockup_evidence", "门窗／开窗器供应方/工头/物业", "准确手摇器、推拉力、开启角、底座孔距、下轨基层与排水腔避让、防水保修；一樘非厨房窗样板及开关锁闭淋水检查", "no", "机械手摇方向、无电源控制和现场既有下轨实例已关闭。"),
        ("E303-STRONGBOX-IDENTITY", "owner_and_site_evidence", "项目内部", "现场铭牌和正面照片；P230／PX1 保持手写识读边界", "yes", "只关闭既有箱体可见身份，不证明回路能力。"),
        ("E303-STRONGBOX-CAPACITY", "site_and_electrical_design_evidence", "电气设计/现场电工", "总开关与分路断路器型号整定、回路标签、导线截面、剩余模数和新增负荷计算", "no", "63 A 箱体铭牌值不得冒充真实容量或可增容结论。"),
    ]
    for input_id, kind, party, evidence, automatic, notes in missing_existing_rules:
        item = {field: "" for field in fields}
        item.update({
            "input_id": input_id,
            "closeout_kind": kind,
            "responsible_party": party,
            "required_evidence": evidence,
            "automatic_close_allowed": automatic,
            "notes": notes,
        })
        upsert_full(rows, "input_id", item)
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
    legacy_rules = [
        ("LEGACY-RCP-CEILING-ACCESS", "project_detail_and_mockup", "室内设计/吊顶深化/机电安装/材料供应方", "RCP1/A-106 天花分区与 1:5 节点；材料防火防潮证据；龙骨／基层／固定；板块模数；整板检修和液压五金；设备拆出路线；1:1 样板", "旧稿方向保留，但 600／300 龙骨文字和网络图片不得直接下料。"),
        ("LEGACY-RCP-CURTAIN-COVE-WINERACK", "project_detail_and_structural_evidence", "室内设计/灯光/窗饰/全屋定制/结构或吊顶深化", "A-106/RCP1/I-501 联合节点；窗帘盒、灯槽、风口、电源和检修协调；酒架自重＋满载、顶部固定和加固计算／节点", "旧稿 200 mm 与约 1300 mm 为设计输入，最终以项目图和实测关闭。"),
        ("LEGACY-INT1-BAYWINDOW-WOOD", "project_detail_and_material_evidence", "室内设计/全屋定制/门窗防水施工", "I-504/D-601 飘窗平立剖和 1:5 节点；材料样板；基层、防潮、日晒、伸缩、检修及窗边防水证明", "不从旧稿图片反推加工尺寸。"),
        ("LEGACY-WFIN-BALCONY-DECK", "product_and_project_detail", "室内设计/材料供应/防水排水/现场", "准确塑木产品资料与样板；A-105/D-601 完成面、龙骨／支座、排水、检修、收边、防火和耐候节点", "旧稿只证明材料方向，不证明产品已选。"),
        ("LEGACY-M401-L-SUPPLY-AIR", "mep_design_and_external_review", "机电设计/日立安装方/风口供应方", "M-401/RCP1 风量、静压、噪声、有效开口、转角压损、风管、型材和检修图；安装方按图复核和调试记录", "造型方向不能替代性能核算。"),
    ]
    for input_id, kind, party, evidence, notes in legacy_rules:
        item = {field: "" for field in fields}
        item.update({
            "input_id": input_id,
            "closeout_kind": kind,
            "responsible_party": party,
            "required_evidence": evidence,
            "automatic_close_allowed": "no",
            "notes": notes,
        })
        upsert_full(rows, "input_id", item)
    write_csv(path, fields, rows)


def reconcile_equipment() -> None:
    path = DECISIONS / "equipment-register.csv"
    fields, rows = read_csv(path)
    updates = {
        "APP-016": {
            "manufacturer": "勒科斯（基准候选） / Franke（比较候选）",
            "model": "F50（基准候选） / Slim 50 CN LD370-E01B 134.0721.210（比较候选）",
            "variant": "F50 retained; Franke total height 365 mm; final product not selected",
            "procurement_status": "candidate",
            "decision_status": "partial",
            "identity_basis": "业主纠正 F50 未否决；Franke 中国官方资料确认 Slim 50 CN 准确 SKU、365 mm 总高和部分接口，只作比选",
            "source_ids": "OWNER-INBOX-20260814-001;OUTBOUND-FORM-PLUM-20260817;APP-016-FRANKE-SLIM50-CN-RESEARCH-20260818;OUTBOUND-FORM-KITCHEN-20260818",
            "notes": "F50 保留为基准候选，继续比较品质、资料、机身高度、接口、噪声、保修和可接受价格；两者均未最终选定／采购。项目先画 Foster 200 mm 深水槽下联合剖面，未给的法兰、开孔和接口中心继续不猜。",
        },
        "HVAC-001": {
            "manufacturer": "日立/海信日立",
            "model": "RPIZ-22FSLN5QDF/P 700x447x192",
            "variant": "18-50; QDF mirror orientation confirmed in formal IFC",
            "procurement_status": "existing",
            "decision_status": "confirmed",
            "identity_basis": "开发商随精装交付的既有设备；业主确认 A02/A03 QDF 与 QD 镜像关系和正式 IFC 建模一致；P02010Q 覆盖型号族",
            "notes": "不再做产品选型；项目负责 M-401/RCP1 路线，回风口兼检修口。说明书未给的接口中心不猜，旧 IFC 紫色管线不代表带保温完成外径。",
        },
        "HVAC-002": {
            "manufacturer": "日立/海信日立",
            "model": "RPIZ-22FSLN5QD/P 700x447x192",
            "variant": "18-50; QD orientation confirmed in formal IFC",
            "procurement_status": "existing",
            "decision_status": "confirmed",
            "identity_basis": "开发商随精装交付的既有设备；业主确认 A01/A04 QD 接管方向和正式 IFC 建模一致；P02010Q 覆盖型号族",
            "notes": "不再做产品选型；项目负责 M-401/RCP1 路线，回风口兼检修口。说明书未给的接口中心不猜。",
        },
        "HVAC-003": {
            "manufacturer": "日立/海信日立",
            "model": "RPIZ-71FSLN5QD/P 1180x447x192",
            "variant": "63-71; gas 15.88 mm; liquid 9.53 mm; VP25 drain OD 32 mm",
            "procurement_status": "existing",
            "decision_status": "confirmed",
            "identity_basis": "开发商既有设备；2026-08-19 现场铭牌关闭 A05 的 RPIZ-71FSLN5QD/P 身份；P02010Q 覆盖型号并给裸管规格",
            "source_ids": "RCP1-HVAC-IFACE-003;RCP1-HVAC-IFACE-005;OUTBOUND-FORM-HVAC-20260817;OWNER-HVAC-INSULATION-PHOTOS-20260819",
            "notes": "A05 不再是候选；项目按 P02010Q 画接口方向、风口、路线与检修净空。接口中心未提供时继续不猜。",
        },
        "HVAC-004": {
            "manufacturer": "日立/海信日立",
            "model": "existing developer-delivered unit; exact A06 mapping pending construction access",
            "variant": "site inventory includes RPIZ-32FSLN5QD/P and RPIZ-28FSLN5QD/P but location not bound",
            "procurement_status": "existing",
            "decision_status": "partial",
            "identity_basis": "A06 机位与开发商既有系统身份已确认；现场另拍到 28／32 型铭牌但无房间／A 编号，不能猜绑定",
            "source_ids": "RCP1-HVAC-IFACE-004;RCP1-A06-LEGACY-001;OUTBOUND-FORM-HVAC-20260817",
            "notes": "不要求业主拆遮挡补拍；施工正常开放遮挡后由安装方同框拍房间／机位环境与铭牌，再由项目绑定 A06。",
        },
        "HVAC-CTRL-001": {
            "manufacturer": "日立／海信日立",
            "model": "PC-P1HEQ",
            "variant": "5 个开发商既有控制点；一只最多控制 6 台；两芯≥0.75 mm²",
            "procurement_status": "existing",
            "decision_status": "confirmed",
            "identity_basis": "开发商随精装交付且正在既有日立系统使用；官方 P02037Q 确认产品和通用接线条件",
            "notes": "不再问采购或兼容性。五点控制分组由业主／项目设计画入 E-302；现场只追线记录实际机组、端子和底盒。",
        },
        "NET-AP-R09": {
            "manufacturer": "TP-Link",
            "model": "TL-XAP1500GE-PoE/DC 易展版（候选）",
            "variant": "Wi-Fi 6 AX1500; Ø184x40; cutout Ø155; 802.3at PoE; max 11.1 W",
            "procurement_status": "candidate",
            "decision_status": "candidate",
            "source_ids": "OWNER-INBOX-20260814-001;TPLINK-XAP1500GE-OFFICIAL-20260819;OUTBOUND-FORM-SITE-20260817",
            "identity_basis": "业主选择 TP-Link 嵌入式款式作为优先方向；官方页关闭外形、开孔和 PoE 参数",
            "notes": "Cat6 星型回弱电箱、AP 点不设 220 V。准确版本、控制器／网关、覆盖、通断、PoE 预算和吊顶机械适配待深化；未采购。",
        },
        "NET-AP-R14": {
            "manufacturer": "TP-Link",
            "model": "TL-XAP1500GE-PoE/DC 易展版（候选）",
            "variant": "Wi-Fi 6 AX1500; Ø184x40; cutout Ø155; 802.3at PoE; max 11.1 W",
            "procurement_status": "candidate",
            "decision_status": "candidate",
            "source_ids": "OWNER-INBOX-20260814-001;TPLINK-XAP1500GE-OFFICIAL-20260819;OUTBOUND-FORM-SITE-20260817",
            "identity_basis": "业主选择 TP-Link 嵌入式款式作为优先方向；官方页关闭外形、开孔和 PoE 参数",
            "notes": "Cat6 星型回弱电箱、AP 点不设 220 V。准确版本、控制器／网关、覆盖、通断、PoE 预算和吊顶机械适配待深化；未采购。",
        },
        "APP-019": {
            "manufacturer": "Musso（推荐候选） / Ninja（备选）",
            "model": "Mini Lussino 4080 / CREAMi 220 V（候选）",
            "variant": "project envelope 500W x 450D x 450H mm; ≥20 kg; 10A; no water/drain",
            "procurement_status": "candidate",
            "decision_status": "partial",
            "storage_location_candidate": "厨房台面／可取放柜格，待 I-501",
            "use_location_candidate": "厨房稳定台面，待 I-501",
            "identity_basis": "业主新增冰淇淋机；项目按两款家用台面机不利包络形成装修预留，不锁品牌",
            "source_ids": "OWNER-INPUT-APP019-20260815;OUTBOUND-FORM-KITCHEN-20260818",
            "notes": "按 500×450×450 mm、承重≥20 kg、一个可拔插 10A 接地插座和四周散热预留；不预留给排水，不按商用软冰机预留；未选／未购。",
        },
        "ACCESS-INTERCOM-001": {
            "manufacturer": "DNAKE／狄耐克",
            "model": "280M-S3",
            "variant": "10.1-inch Linux; approx 270x168x15 mm; 802.3af PoE or DC 12 V/2 A",
            "procurement_status": "existing",
            "decision_status": "confirmed",
            "source_ids": "OWNER-E304-DNAKE-PHOTO-20260817;DNAKE-280M-S3-OFFICIAL-20260819;OUTBOUND-FORM-SITE-20260817",
            "identity_basis": "业主现场正面照片与 DNAKE 官方产品资料匹配识别为 280M-S3；现场软件版本日期 2021-06-15",
            "notes": "型号不再重复询问。实际供电、背部端子、物业系统兼容、保留／拆下复装和最终坐标仍待现场；网络地址不投影。",
        },
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
    purchased_faucets = [
        {
            "equipment_id": "SAN-023",
            "domain": "SANITARY",
            "category": "concealed_basin_faucet",
            "item_name": "已购暗装面盆龙头",
            "manufacturer": "GESSI／gessi（订单商品标题；原厂身份待证）",
            "model": "MF287（卖家分享标识；非已验证厂家料号）",
            "variant": "拉丝金／双把",
            "quantity": "1",
            "procurement_status": "purchased_delivery_unverified",
            "decision_status": "partial",
            "use_location_candidate": "BATHG 客卫干区／Falper Sorgente 协调候选",
            "schedule_included": "yes",
            "selector_kind": "logical_input",
            "selector_value": "SAN-023",
            "source_ids": "OWNER-PLUM-PURCHASED-FAUCETS-20260820;OWNER-PLUM-FAUCETS-ORDER-IMAGE-20260820;OWNER-PLUM-MF287-INSTALL-IMAGE-20260820",
            "identity_basis": "订单截图确认交易成功、数量 1 和拉丝金／双把款式；MF287 只作卖家分享标识，尺寸图只关闭 G1/2 与部件关系",
            "confidence": "1.00",
            "human_review_required": "yes",
            "legacy_kind": "owner_product",
            "legacy_id": "MF287",
            "notes": "已购事实不是候选；准确厂家料号、204／200 图面单位、暗装阀体包络、允许埋深、完成面基准、中心关系和最终房间映射仍待供货方资料与项目图关闭。",
        },
        {
            "equipment_id": "SAN-024",
            "domain": "SANITARY",
            "category": "deck_mounted_basin_faucet",
            "item_name": "已购双孔台面面盆龙头",
            "manufacturer": "GESSI／gessi（订单商品标题；原厂身份待证）",
            "model": "CZ356（卖家分享标识；非已验证厂家料号）",
            "variant": "拉丝金半圆款／双孔分体",
            "quantity": "1",
            "procurement_status": "purchased_delivery_unverified",
            "decision_status": "partial",
            "use_location_candidate": "BATHM 主卫／antoniolupi Street 协调候选",
            "schedule_included": "yes",
            "selector_kind": "logical_input",
            "selector_value": "SAN-024",
            "source_ids": "OWNER-PLUM-PURCHASED-FAUCETS-20260820;OWNER-PLUM-FAUCETS-ORDER-IMAGE-20260820",
            "identity_basis": "订单截图确认交易成功、数量 1 和拉丝金半圆款；CZ356 只作卖家分享标识，当前没有单独安装图",
            "confidence": "1.00",
            "human_review_required": "yes",
            "legacy_kind": "owner_product",
            "legacy_id": "CZ356",
            "notes": "已购事实不是候选；准确厂家料号、两孔孔径／中心距、台面厚度、紧固、软管接口和最终房间映射保持 unknown。不得套用 Street 官方单 Ø38 mm 龙头孔。",
        },
        {
            "equipment_id": "SAN-025",
            "domain": "SANITARY",
            "category": "concealed_shower_system",
            "item_name": "已购暗装淋浴花洒",
            "manufacturer": "GESSI／gessi（订单商品标题；原厂身份待证）",
            "model": "CZ028（卖家分享标识；非已验证厂家料号）",
            "variant": "拉丝金色墙出式／顶喷＋手持关系",
            "quantity": "1",
            "procurement_status": "purchased_delivery_unverified",
            "decision_status": "partial",
            "use_location_candidate": "两个淋浴区之一；准确房间待 I-502",
            "schedule_included": "yes",
            "selector_kind": "logical_input",
            "selector_value": "SAN-025",
            "source_ids": "OWNER-PLUM-PURCHASED-FAUCETS-20260820;OWNER-PLUM-FAUCETS-ORDER-IMAGE-20260820;OWNER-PLUM-CZ028-INSTALL-IMAGE-20260820",
            "identity_basis": "订单截图确认交易成功、数量 1 和拉丝金色墙出式；CZ028 只作卖家分享标识，卖家图只关闭冷热进水、阀体、顶喷和手持的系统关系",
            "confidence": "1.00",
            "human_review_required": "yes",
            "legacy_kind": "owner_product",
            "legacy_id": "CZ028",
            "notes": "已购事实不是候选；项目有两个淋浴区但只购 1 套。准确房间、第二淋浴区产品、厂家料号、阀体包络、开槽单位／基准、接口、流量和检修仍待关闭。",
        },
    ]
    for values in purchased_faucets:
        item = {field: "" for field in fields}
        item.update(values)
        upsert_full(rows, "equipment_id", item)
    newly_absorbed_equipment = [
        ("APP-020", "APPLIANCE", "juicer", "已购惠人 H70FT 原汁机", "HUROM／惠人", "H70FT", "岛台下存放；台面使用", "purchased_delivery_unverified", "partial", "", "岛台下", "厨房台面", "OWNER-KITCHEN-H70FT-20260820;OUTBOUND-FORM-KITCHEN-20260818", "appliance", "型号、已购和存放位置已确认；到货实测、承重、取放、湿附件通风和使用插座待 I-501／E-303。"),
        ("KTOOL-001", "KITCHENWARE", "long_utensil", "Dreamfarm Clongs 食品夹", "Dreamfarm", "Clongs", "K-D01 长工具", "purchased_delivery_unverified", "partial", "", "K-D01", "厨房操作台", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "owner_kitchen_tool", "已购与收纳分组已关闭；只缺实物包络和抽屉摆样。"),
        ("KTOOL-002", "KITCHENWARE", "bowl_clip", "Kuhn Rikon 夹碗夹", "Kuhn Rikon", "订单 SKU 待到货复核", "K-D01 长工具", "purchased_delivery_unverified", "partial", "", "K-D01", "厨房操作台", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "owner_kitchen_tool", "已购与收纳分组已关闭；准确 SKU 和包络待实物复核。"),
        ("KTOOL-003", "KITCHENWARE", "wine_preservation", "Coravin Model 6+ 酒具套装", "Coravin", "Model 6+", "K-D02 独立酒具抽／柜格", "purchased_delivery_unverified", "partial", "", "K-D02", "餐厨区", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "owner_kitchen_tool", "连原收纳盒独立存放；只缺实物盒尺寸和取放验证。"),
        ("KTOOL-004", "KITCHENWARE", "wood_utensil_set", "Alessi Pots&Pans 木制厨具三件套", "Alessi", "Pots&Pans 3-piece", "K-D01 长工具", "purchased_delivery_unverified", "partial", "", "K-D01", "厨房操作台", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "owner_kitchen_tool", "已购与收纳分组已关闭；只缺实物摆样。"),
        ("KTOOL-005", "KITCHENWARE", "garlic_press", "Dreamfarm Garject 压蒜器", "Dreamfarm", "Garject", "K-D01 长工具", "purchased_delivery_unverified", "partial", "", "K-D01", "厨房操作台", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "owner_kitchen_tool", "已购；设独立短格，只缺实物摆样。"),
        ("KTOOL-006", "KITCHENWARE", "kitchen_shears", "Dreamfarm Bishears 厨房剪", "Dreamfarm", "Bishears", "K-D01 长工具", "purchased_delivery_unverified", "partial", "", "K-D01", "厨房操作台", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "owner_kitchen_tool", "已购；带护套独立短格，只缺实物摆样。"),
        ("KTOOL-007", "KITCHENWARE", "saucepan", "不锈钢单柄奶锅与玻璃沥水盖", "订单品牌待实物复核", "订单 SKU 待实物复核", "K-D03 深锅具抽", "purchased_delivery_unverified", "partial", "", "K-D03", "厨房操作台", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "owner_kitchen_tool", "锅体平放，玻璃盖软垫竖放；只缺包含锅柄的实测包络。"),
        ("KTOOL-008", "KITCHENWARE", "coffee_server_set", "HARIO 胡桃木 V60 玻璃套装", "HARIO", "V60 walnut glass set", "K-D04 易碎咖啡器具抽", "purchased_delivery_unverified", "partial", "", "K-D04", "厨房操作台", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "owner_kitchen_tool", "已购；防滑软垫和可替换定位块，只缺实物摆样。"),
        ("WINTR-001", "WINDOW_TREATMENT", "day_shade", "公共区日帘", "Hunter Douglas／亨特道格拉斯", "Luminette 萝美雅 K5-501", "产品家族与颜色方向确认", "not_selected", "partial", "", "", "公共区窗组", "OWNER-INT1-WINDOW-TREATMENTS-20260818", "owner_window_treatment", "未采购；只缺逐窗复尺、分幅、轨道、操作、色卡、报价和安装图。"),
        ("WINTR-002", "WINDOW_TREATMENT", "day_shade", "卧室日帘", "Hunter Douglas／亨特道格拉斯", "Silhouette 丝络雅 N31-204", "产品家族与颜色方向确认", "not_selected", "partial", "", "", "卧室窗组", "OWNER-INT1-WINDOW-TREATMENTS-20260818", "owner_window_treatment", "未采购；只缺逐窗复尺、分幅、轨道、操作、色卡、报价和安装图。"),
        ("WINTR-003", "WINDOW_TREATMENT", "wet_room_blind", "厨卫防潮百叶", "品牌待报价", "25 mm 铝百叶", "灰黑色方向确认", "not_selected", "partial", "", "", "厨卫窗组", "OWNER-INT1-WINDOW-TREATMENTS-20260818", "owner_window_treatment", "未采购；只缺逐窗复尺、耐湿样品、操作方式和安装图。"),
        ("WINTR-004", "WINDOW_TREATMENT", "bedroom_blackout", "卧室夜间遮光系统", "Hunter Douglas／亨特道格拉斯", "Duette LightLock（优先） / Silhouette Duolite（对照）", "两案比价，尚未终选", "not_selected", "candidate", "", "", "卧室窗组", "OWNER-INT1-WINDOW-TREATMENTS-20260818", "owner_window_treatment", "候选不是已选或已购；只缺实体遮光比较、报价和安装节点后的终选。"),
        ("WIN-OPENER-001", "DOOR_WINDOW", "manual_window_opener", "既有窗手摇开窗器改造", "准确产品待样板", "mechanical crank opener", "纯机械，不加电源与智能联动", "candidate", "partial", "", "", "非厨房既有窗，先做一樘样板", "OWNER-WIN-CRANK-RETROFIT-20260820;OWNER-WIN-CRANK-PHOTO-01-20260820;OWNER-WIN-CRANK-PHOTO-02-20260820", "owner_window_hardware", "22 mm 仅为现场几何依据；准确产品、受力、孔距、排水腔避让、防水保修和样板待关闭。"),
        ("ELEC-DB-EXISTING-001", "ELECTRICAL", "distribution_box", "既有强电配电箱", "深圳市三江电气有限公司", "P230（手写识读）", "箱号 PX1；铭牌 380/220 V、63 A", "existing", "partial", "", "", "玄关高柜现状强电箱", "E303-PHOTO-001;E303-PHOTO-002;ELEC-BOX-EVIDENCE-20260820", "owner_existing_equipment", "63 A 是箱体铭牌观察值，不是总开关整定、导线容量或可增容结论。"),
        ("HVAC-INS-EXISTING-001", "HVAC", "pipe_insulation", "既有华美 Class 1 橡塑保温管套", "Huamei／华美", "Class 1 Rubber Foam", "现场观察 11 种 ID×TK；逐段对应待施工开放", "existing", "ifc_observed", "", "", "开发商既有空调管路", "OWNER-HVAC-INSULATION-PHOTOS-20260819;HUAMEI-CLASS1-OFFICIAL-20260819", "hvac_insulation", "现场规格不是本项目最终计算厚度；最终按珠海露点、管温、导热系数、防火和逐段管径计算。"),
        ("PLUM-LAUNDRY-TRAP-001", "DRAINAGE", "double_appliance_trap", "洗烘双独立存水弯组合候选", "良舍高端进口卫浴（非吉博力官方旗舰店）", "Taobao item 703917654224 / SKU 5125125938318", "两台设备各自独立存水弯后下游合流", "candidate", "partial", "", "", "洗衣区墙排", "OWNER-PLUM-GEBERIT-DOUBLE-TRAP-20260819", "owner_product_component", "优先比选 SKU，不写成吉博力原厂套装或已采购；准确零件、接口、水封、防串水、标高、检修和双机试验待按图复核。"),
    ]
    for equipment_id, domain, category, item_name, manufacturer, model, variant, procurement, decision, storage_candidate, storage_confirmed, use_confirmed, source_ids, legacy_kind, notes in newly_absorbed_equipment:
        item = {field: "" for field in fields}
        item.update({
            "equipment_id": equipment_id, "domain": domain, "category": category,
            "item_name": item_name, "manufacturer": manufacturer, "model": model,
            "variant": variant, "quantity": "1", "procurement_status": procurement,
            "decision_status": decision, "storage_location_candidate": storage_candidate,
            "storage_location_confirmed": storage_confirmed, "use_location_confirmed": use_confirmed,
            "schedule_included": "no" if equipment_id == "ELEC-DB-EXISTING-001" else "yes",
            "selector_kind": "logical_input", "selector_value": equipment_id,
            "source_ids": source_ids, "identity_basis": notes, "confidence": "1.00",
            "human_review_required": "yes", "legacy_kind": legacy_kind,
            "legacy_id": equipment_id, "notes": notes,
        })
        upsert_full(rows, "equipment_id", item)
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
    for req_id in [f"REQ-{index:04d}" for index in range(290, 312)]:
        if any(row["requirement_id"] == req_id for row in rows):
            old_notes = next(row["notes"] for row in rows if row["requirement_id"] == req_id)
            update_by_id(rows, "requirement_id", req_id, {
                "status": "confirmed",
                "blocks_release": "no",
                "notes": append_sentence_once(old_notes, "开发商既有设备与 QD/QDF 镜像关系已确认；本参数不再作为选型候选。"),
            })
    hvac_a05 = {
        "REQ-0316": ("", "15.88", "mm", "confirmed", "no", "P02010Q 63–71 容量段气管外径。"),
        "REQ-0317": ("", "9.53", "mm", "confirmed", "no", "P02010Q 63–71 容量段液管外径。"),
        "REQ-0318": ("", "32", "mm", "confirmed", "no", "P02010Q VP25 排水管外径。"),
        "REQ-0319": ("1/100", "", "", "confirmed", "no", "P02010Q 冷凝水连续坡度边界。"),
        "REQ-0320": ("1/25", "", "", "confirmed", "no", "P02010Q 冷凝水连续坡度边界。"),
        "REQ-0321": ("QD orientation; project drawing required", "", "", "confirmed", "no", "A05 铭牌 QD 身份已关闭；项目负责画路线。"),
        "REQ-0322": ("official diagram basis; exact project center not provided", "", "", "pending", "yes", "不得推测接口中心；由项目图与厂家按图复核关闭。"),
    }
    for req_id, (text_value, number, unit, status, blocks, notes) in hvac_a05.items():
        update_by_id(rows, "requirement_id", req_id, {
            "value_text": text_value,
            "value_number": number,
            "unit": unit,
            "value_origin": "official_exact_model",
            "status": status,
            "source_id": "RCP1-HVAC-IFACE-003",
            "source_locator": "P02010Q；A05 RPIZ-71FSLN5QD/P",
            "blocks_release": blocks,
            "notes": notes,
        })
    update_by_id(rows, "requirement_id", "REQ-HVACCTRL-006", {
        "value_text": "existing developer-delivered controller operating with existing Hitachi system",
        "value_origin": "site_observed",
        "status": "confirmed",
        "blocks_release": "no",
        "source_id": "OUTBOUND-FORM-HVAC-20260817",
        "source_locator": "现行 02 表；业主确认精装房原配",
        "notes": "不再询问兼容性；五点物理映射和最终控制策略分别处理。",
    })
    update_by_id(rows, "requirement_id", "REQ-INTERCOM-004", {
        "value_text": "280M-S3",
        "value_origin": "research_conclusion",
        "status": "confirmed",
        "source_id": "DNAKE-280M-S3-OFFICIAL-20260819",
        "source_locator": "业主现场照片与官方外观／规格匹配",
        "blocks_release": "no",
        "notes": "型号不再重复询问；实际端子和供电仍由后续要求关闭。",
    })
    for req_id in ("REQ-0169", "REQ-0600", "REQ-0601", "REQ-0602"):
        old_notes = next(row["notes"] for row in rows if row["requirement_id"] == req_id)
        update_by_id(rows, "requirement_id", req_id, {
            "notes": old_notes.replace("F50 已否决；", "F50 保留为基准候选；").replace("最终同型号", "最终所选同型号"),
        })
    additions = [
        ("REQ-APP014-SLEEVE-20260818-001", "replaceable_service_sleeve_direction", "continuous_pull-through sleeve; install PE tube after fit-out", "user_input", "candidate", "no", "业主候选方向；不是已批准施工规格。"),
        ("REQ-APP014-SLEEVE-20260818-002", "supplied_tube_material_and_od", "unknown", "pending", "pending", "yes", "由西门子准确型号资料或书面回复关闭。"),
        ("REQ-APP014-SLEEVE-20260818-003", "sleeve_compatibility_and_bend_radius", "unknown", "pending", "pending", "yes", "不得用网络案例尺寸代填。"),
        ("REQ-APP014-SLEEVE-20260818-004", "sleeve_route_end_waterproof_and_leak_visibility", "unknown", "pending", "pending", "yes", "由 P-201／I-501 联合节点关闭。"),
        ("REQ-APP014-SLEEVE-20260818-005", "full_route_pull_through_mockup", "required", "project_candidate", "pending", "yes", "施工前完成全程抽换样板并留照片。"),
        ("REQ-APP019-ENVELOPE-20260820-001", "project_clear_width", "500", "research_conclusion", "confirmed", "no", "家用台面机候选不利包络；不是实购设备宽度。"),
        ("REQ-APP019-ENVELOPE-20260820-002", "project_clear_depth", "450", "research_conclusion", "confirmed", "no", "家用台面机候选不利包络；不是实购设备深度。"),
        ("REQ-APP019-ENVELOPE-20260820-003", "project_clear_height", "450", "research_conclusion", "confirmed", "no", "家用台面机候选不利包络；不是实购设备高度。"),
        ("REQ-APP019-ENVELOPE-20260820-004", "minimum_support_load", "20", "research_conclusion", "confirmed", "no", "稳定台面或可取放柜格承重下限。"),
        ("REQ-APP019-ENVELOPE-20260820-005", "wall_socket_rating", "10 A earthed accessible socket", "research_conclusion", "confirmed", "no", "最终按实购插头复核；当前家用候选可共用普通厨房小家电回路。"),
        ("REQ-APP019-ENVELOPE-20260820-006", "water_and_drain_reservation", "not required for current domestic candidates", "research_conclusion", "confirmed", "no", "若改商用软冰机须重新设计。"),
        ("REQ-NETAP-TPLINK-20260820-001", "candidate_model", "TL-XAP1500GE-PoE/DC 易展版", "official_exact_model", "candidate", "no", "业主优先方向，非已购。"),
        ("REQ-NETAP-TPLINK-20260820-002", "body_and_cutout", "body Ø184×40 mm; ceiling cutout Ø155 mm", "official_exact_model", "candidate", "yes", "吊顶深化须保留水晶头弯曲、散热和拆换空间。"),
        ("REQ-NETAP-TPLINK-20260820-003", "poe_and_power", "802.3at PoE; maximum 11.1 W", "official_exact_model", "candidate", "yes", "用于交换机端口和总 PoE 预算。"),
        ("REQ-NETAP-R14-TPLINK-20260820-001", "candidate_model", "TL-XAP1500GE-PoE/DC 易展版", "official_exact_model", "candidate", "no", "业主优先方向，非已购。"),
        ("REQ-NETAP-R14-TPLINK-20260820-002", "body_and_cutout", "body Ø184×40 mm; ceiling cutout Ø155 mm", "official_exact_model", "candidate", "yes", "吊顶深化须保留水晶头弯曲、散热和拆换空间。"),
        ("REQ-NETAP-R14-TPLINK-20260820-003", "poe_and_power", "802.3at PoE; maximum 11.1 W", "official_exact_model", "candidate", "yes", "用于交换机端口和总 PoE 预算。"),
    ]
    for req_id, requirement, text_value, basis, status, blocks, notes in additions:
        if req_id.startswith("REQ-APP019"):
            equipment_id = "APP-019"
            discipline = "ELEC/INT1/PLUM"
            source_id = "OUTBOUND-FORM-KITCHEN-20260818"
            locator = "现行 09 表；家用台面机候选不利包络"
        elif req_id.startswith("REQ-NETAP"):
            equipment_id = "NET-AP-R14" if req_id.startswith("REQ-NETAP-R14") else "NET-AP-R09"
            discipline = "NETWORK/ELEC/RCP1"
            source_id = "TPLINK-XAP1500GE-OFFICIAL-20260819"
            locator = "TP-Link 官方产品页"
        else:
            equipment_id = "APP-014"
            discipline = "PLUM/INT1"
            source_id = "OWNER-APP014-REPLACEABLE-SLEEVE-20260817"
            locator = "owner candidate and evidence boundary"
        item = {field: "" for field in fields}
        item.update({
            "requirement_id": req_id,
            "equipment_id": equipment_id,
            "discipline": discipline,
            "parameter_key": requirement,
            "value_text": text_value,
            "value_origin": basis,
            "status": status,
            "source_id": source_id,
            "source_locator": locator,
            "blocks_release": blocks,
            "notes": notes,
        })
        numeric_units = {
            "REQ-APP019-ENVELOPE-20260820-001": "mm",
            "REQ-APP019-ENVELOPE-20260820-002": "mm",
            "REQ-APP019-ENVELOPE-20260820-003": "mm",
            "REQ-APP019-ENVELOPE-20260820-004": "kg",
        }
        if req_id in numeric_units:
            item["value_text"] = ""
            item["value_number"] = text_value
            item["unit"] = numeric_units[req_id]
        if any(row["requirement_id"] == req_id for row in rows):
            update_by_id(rows, "requirement_id", req_id, item)
        else:
            rows.append(item)
    faucet_requirements = [
        # SAN-023 / MF287: purchase and visible configuration are confirmed;
        # installation dimensions remain intentionally unknown.
        ("REQ-SAN023-001", "SAN-023", "PROCUREMENT", "purchase_status", "purchased; transaction successful; quantity 1; delivery not evidenced", "owner_order_evidence", "confirmed", "OWNER-PLUM-FAUCETS-ORDER-IMAGE-20260820", "订单截图；2024-12-04", "no", "关闭采购事实和数量，不证明已收货、原厂身份或安装接口。"),
        ("REQ-SAN023-002", "SAN-023", "PLUM/INT1", "installation_type_and_finish", "concealed wall-mounted basin faucet; brushed gold; twin handle", "owner_order_evidence", "confirmed", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "采购事实与深化边界", "no", "功能类型和订单所示款式已确认。"),
        ("REQ-SAN023-003", "SAN-023", "PLUM", "seller_diagram_visible_interface", "G1/2 visible; control and spout relationship shown", "seller_installation_image", "candidate", "OWNER-PLUM-MF287-INSTALL-IMAGE-20260820", "MF287 卖家尺寸图", "no", "只记录图上明确可读关系；不据此冻结墙内中心。"),
        ("REQ-SAN023-004", "SAN-023", "PROCUREMENT", "exact_manufacturer_article_and_sku", "unknown", "pending", "pending", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "订单缺准确原厂料号", "yes", "MF287 是分享标识，不冒充厂家型号。"),
        ("REQ-SAN023-005", "SAN-023", "PLUM/INT1", "dimension_units_centres_and_valve_envelope", "unknown", "pending", "pending", "OWNER-PLUM-MF287-INSTALL-IMAGE-20260820", "204cm／200cm 单位不可信", "yes", "须由卖家在原图确认单位，并补阀体包络、允许埋深、完成面基准和公差。"),
        ("REQ-SAN023-006", "SAN-023", "PLUM/INT1", "final_room_mapping_and_project_coordinates", "unknown", "project_candidate", "pending", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "优先与客卫 Sorgente 协调", "yes", "项目先完成 I-502／P-202，再由供货方按图复核；接口中心不得推测。"),
        # SAN-024 / CZ356.
        ("REQ-SAN024-001", "SAN-024", "PROCUREMENT", "purchase_status", "purchased; transaction successful; quantity 1; delivery not evidenced", "owner_order_evidence", "confirmed", "OWNER-PLUM-FAUCETS-ORDER-IMAGE-20260820", "订单截图；2024-12-04", "no", "关闭采购事实和数量，不证明已收货、原厂身份或安装接口。"),
        ("REQ-SAN024-002", "SAN-024", "PLUM/INT1", "installation_type_and_finish", "deck-mounted two-hole split basin faucet; brushed-gold half-round style", "owner_order_evidence", "confirmed", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "采购事实与深化边界", "no", "功能类型和订单所示款式已确认。"),
        ("REQ-SAN024-003", "SAN-024", "PROCUREMENT", "exact_manufacturer_article_and_sku", "unknown", "pending", "pending", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "订单缺准确原厂料号", "yes", "CZ356 是分享标识，不冒充厂家型号。"),
        ("REQ-SAN024-004", "SAN-024", "PLUM/INT1", "hole_diameters_centres_and_countertop_range", "unknown", "pending", "pending", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "当前无 CZ356 安装图", "yes", "须补两孔孔径／中心距、允许台面厚度、台下紧固和软管抽换空间。"),
        ("REQ-SAN024-005", "SAN-024", "PLUM/INT1", "supply_interfaces_and_valve_access", "unknown", "pending", "pending", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "当前无 CZ356 安装图", "yes", "冷热接口、软管、阀门与检修必须进入联合剖面。"),
        ("REQ-SAN024-006", "SAN-024", "PLUM/INT1", "final_room_mapping_and_project_coordinates", "unknown", "project_candidate", "pending", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "优先与主卫 Street 协调", "yes", "不得套用 Street 官方单 Ø38 mm 孔；由 I-502／P-202 与定制盆加工图关闭。"),
        # SAN-025 / CZ028.
        ("REQ-SAN025-001", "SAN-025", "PROCUREMENT", "purchase_status", "purchased; transaction successful; quantity 1; delivery not evidenced", "owner_order_evidence", "confirmed", "OWNER-PLUM-FAUCETS-ORDER-IMAGE-20260820", "订单截图；2024-12-04", "no", "只购 1 套；订单不证明已收货；不得复制为两个淋浴区。"),
        ("REQ-SAN025-002", "SAN-025", "PLUM/INT1", "installation_type_and_finish", "concealed shower system; brushed-gold wall-outlet style", "owner_order_evidence", "confirmed", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "采购事实与深化边界", "no", "功能类型和订单所示款式已确认。"),
        ("REQ-SAN025-003", "SAN-025", "PLUM", "seller_diagram_system_relationship", "hot/cold supplies to concealed valve; overhead and hand shower shown", "seller_installation_image", "candidate", "OWNER-PLUM-CZ028-INSTALL-IMAGE-20260820", "CZ028 卖家安装关系图", "no", "只记录系统组成关系，不作为开槽或接口尺寸。"),
        ("REQ-SAN025-004", "SAN-025", "PROCUREMENT", "exact_manufacturer_article_and_sku", "unknown", "pending", "pending", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "订单缺准确原厂料号", "yes", "CZ028 是分享标识，不冒充厂家型号。"),
        ("REQ-SAN025-005", "SAN-025", "PLUM/INT1", "valve_envelope_depth_opening_units_and_datum", "unknown", "pending", "pending", "OWNER-PLUM-CZ028-INSTALL-IMAGE-20260820", "图示 30×10×6 无单位和基准", "yes", "须补阀体包络、允许埋深、完成面基准、公差、防水和检修。"),
        ("REQ-SAN025-006", "SAN-025", "PLUM", "interfaces_outputs_pressure_and_flow", "unknown", "pending", "pending", "OWNER-PLUM-CZ028-INSTALL-IMAGE-20260820", "卖家图只证明关系", "yes", "须补冷热及各出水口规格、顶喷固定、实际组件、压力和设计流量。"),
        ("REQ-SAN025-007", "SAN-025", "PLUM/INT1", "final_room_mapping_and_second_shower_product", "unknown", "pending", "pending", "OWNER-PLUM-PURCHASED-FAUCETS-20260820", "项目有两个淋浴区但只购 1 套", "yes", "由 I-502 确认本套房间归属；第二淋浴区产品另行关闭。"),
    ]
    for req_id, equipment_id, discipline, key, value, origin, status, source_id, locator, blocks, notes in faucet_requirements:
        item = {field: "" for field in fields}
        item.update({
            "requirement_id": req_id,
            "equipment_id": equipment_id,
            "discipline": discipline,
            "parameter_key": key,
            "value_text": value,
            "value_origin": origin,
            "status": status,
            "source_id": source_id,
            "source_locator": locator,
            "blocks_release": blocks,
            "notes": notes,
        })
        if any(row["requirement_id"] == req_id for row in rows):
            update_by_id(rows, "requirement_id", req_id, item)
        else:
            rows.append(item)
    update_by_id(rows, "requirement_id", "REQ-APP017-OWNER-20260818-003", {
        "value_text": "candidate Taobao item 703917654224 / SKU 5125125938318; two independent traps then downstream merge",
        "value_origin": "seller_installation_image",
        "status": "candidate",
        "source_id": "OWNER-PLUM-GEBERIT-DOUBLE-TRAP-20260819",
        "source_locator": "商家双存水弯组合图",
        "blocks_release": "yes",
        "notes": "不再写成完全 unknown；该 SKU 是优先候选，含非吉博力原厂配件，准确零件、接口、水封、包络、标高和双机试验待按图复核。",
    })
    absorbed_requirements = []
    def add_absorbed(req_id, equipment_id, discipline, key, value, origin, status, source_id, blocks, notes, number="", unit=""):
        origin = {
            "owner_confirmation": "user_input", "owner_decision": "user_input",
            "project_decision": "project_candidate", "owner_candidate_seller_evidence": "seller_installation_image",
            "official_product_family": "official_model_family", "engineering_calculation": "project_candidate",
            "project_requirement": "project_candidate",
        }.get(origin, origin)
        status = "observed" if status == "ifc_observed" else status
        absorbed_requirements.append((req_id, equipment_id, discipline, key, value, number, unit, origin, status, source_id, blocks, notes))

    add_absorbed("REQ-APP020-001", "APP-020", "PROCUREMENT", "purchase_status", "purchased; delivery not evidenced", "owner_confirmation", "confirmed", "OWNER-KITCHEN-H70FT-20260820", "no", "已购不等于已到货。")
    add_absorbed("REQ-APP020-002", "APP-020", "INT1", "storage_location", "VVD island base; storage only", "owner_decision", "confirmed", "OWNER-KITCHEN-H70FT-20260820", "no", "岛台下存放已关闭，不在关闭柜格内运行。")
    add_absorbed("REQ-APP020-003", "APP-020", "INT1", "arrival_dimensions_weight_and_accessories", "unknown", "pending", "pending", "OWNER-KITCHEN-H70FT-20260820", "yes", "实测整机 W×D×H、拆投料筒后高度、重量、电源线和附件盒。")
    add_absorbed("REQ-APP020-004", "APP-020", "INT1", "storage_detail", "front-access upright bay or full-extension tray; load, anti-soil mat and wet-accessory ventilation required", "project_decision", "confirmed", "OUTBOUND-FORM-KITCHEN-20260818", "yes", "由 I-501 给出净尺寸、导轨和取放路线后可加工。")
    add_absorbed("REQ-APP020-005", "APP-020", "ELEC", "use_power_location", "accessible 10 A earthed socket at countertop use position; no default socket in storage-only bay", "project_decision", "confirmed", "OUTBOUND-FORM-KITCHEN-20260818", "yes", "实购插头与功率待到货复核，不预留给排水。")

    tool_names = {
        "KTOOL-001": ("K-D01", "Dreamfarm Clongs"), "KTOOL-002": ("K-D01", "Kuhn Rikon bowl clip"),
        "KTOOL-003": ("K-D02", "Coravin Model 6+"), "KTOOL-004": ("K-D01", "Alessi Pots&Pans 3-piece"),
        "KTOOL-005": ("K-D01", "Dreamfarm Garject"), "KTOOL-006": ("K-D01", "Dreamfarm Bishears"),
        "KTOOL-007": ("K-D03", "saucepan and glass lid"), "KTOOL-008": ("K-D04", "HARIO V60 walnut glass set"),
    }
    for index, (equipment_id, (drawer, identity)) in enumerate(tool_names.items(), 1):
        add_absorbed(f"REQ-KTOOL-{index:03d}-01", equipment_id, "PROCUREMENT", "purchase_status", "purchased; delivery not evidenced", "owner_confirmation", "confirmed", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "no", f"{identity} 已购；未证明已到货。")
        add_absorbed(f"REQ-KTOOL-{index:03d}-02", equipment_id, "INT1", "storage_group", drawer, "project_decision", "confirmed", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "no", "抽屉分组已关闭，不交给橱柜方重新分配。")
        add_absorbed(f"REQ-KTOOL-{index:03d}-03", equipment_id, "INT1", "arrival_envelope_and_drawer_mockup", "unknown", "pending", "pending", "OWNER-KITCHEN-TOOLS-PURCHASED-20260820", "yes", "到货后记录准确 SKU／包络，在 I-501 冻结前完成八件实物摆样和导轨承重复核。")

    window_values = {
        "WINTR-001": "Luminette K5-501 public-area day shade",
        "WINTR-002": "Silhouette N31-204 bedroom day shade",
        "WINTR-003": "25 mm grey-black aluminium blind for kitchen and bathrooms",
        "WINTR-004": "Duette LightLock preferred / Silhouette Duolite comparison; not finally selected",
    }
    for index, (equipment_id, value) in enumerate(window_values.items(), 1):
        status = "candidate" if equipment_id == "WINTR-004" else "confirmed"
        add_absorbed(f"REQ-WINTR-{index:03d}-01", equipment_id, "INT1", "product_family_colour_direction", value, "owner_decision", status, "OWNER-INT1-WINDOW-TREATMENTS-20260818", "no", "产品方向不等于已选 SKU、已报价或已下单。")
        add_absorbed(f"REQ-WINTR-{index:03d}-02", equipment_id, "INT1/ELEC", "finished_measure_split_control_track_shopdrawing", "unknown", "pending", "pending", "OWNER-INT1-WINDOW-TREATMENTS-20260818", "yes", "逐窗复尺、分幅、收拢／操作侧、头轨／侧轨、电源控制、色卡、报价和安装图统一关闭。")

    add_absorbed("REQ-WINOPENER-001", "WIN-OPENER-001", "ARCH", "operation_and_power", "manual mechanical crank; no power, controls or smart linkage", "owner_decision", "confirmed", "OWNER-WIN-CRANK-RETROFIT-20260820", "no", "机械方向已关闭。")
    add_absorbed("REQ-WINOPENER-002", "WIN-OPENER-001", "ARCH", "observed_closed_edge", "", "site_observed", "ifc_observed", "OWNER-WIN-CRANK-PHOTO-01-20260820", "no", "约 22 mm 只作研究依据，不直接下料。", "22", "mm")
    add_absorbed("REQ-WINOPENER-003", "WIN-OPENER-001", "ARCH", "exact_product_force_angle_hole_spacing_substrate_waterproof_mockup", "unknown", "pending", "pending", "OWNER-WIN-CRANK-RETROFIT-20260820", "yes", "一樘非厨房窗样板同时验证锁闭、排水腔避让、淋水和保修。")

    add_absorbed("REQ-ELECDB-001", "ELEC-DB-EXISTING-001", "ELEC", "nameplate_identity", "Sanjiang Electric; handwritten P230; box PX1; 380/220 V", "site_observed", "confirmed", "E303-PHOTO-001", "no", "P230／PX1 为手写识读，订货或更换前近距复核。")
    add_absorbed("REQ-ELECDB-002", "ELEC-DB-EXISTING-001", "ELEC", "box_nameplate_current", "", "site_observed", "ifc_observed", "E303-PHOTO-001", "no", "箱体铭牌值，不是断路器整定或可增容结论。", "63", "A")
    add_absorbed("REQ-ELECDB-003", "ELEC-DB-EXISTING-001", "ELEC", "actual_breakers_wiring_modules_and_capacity", "unknown", "pending", "pending", "ELEC-BOX-EVIDENCE-20260820", "yes", "打开箱门记录总开关、分路、导线截面、剩余模数并完成新增负荷计算。")

    for index, spec in enumerate(("6x15", "10x15", "13x9", "13x15", "16x9", "16x15", "16x20", "20x15", "25x15", "32x9", "32x15"), 1):
        add_absorbed(f"REQ-HVACINS-{index:03d}", "HVAC-INS-EXISTING-001", "HVAC", f"observed_id_x_tk_{index:02d}", spec, "site_observed", "ifc_observed", "OWNER-HVAC-INSULATION-PHOTOS-20260819", "no", "ID=内径，TK=单边壁厚；只证明现场出现过该规格。")
    add_absorbed("REQ-HVACINS-020", "HVAC-INS-EXISTING-001", "HVAC", "official_product_boundary", "Huamei Class 1; lambda at 0 C <=0.034 W/(m.K); B1; standards per official source", "official_product_family", "confirmed", "HUAMEI-CLASS1-OFFICIAL-20260819", "no", "官方产品族性能不代表本项目逐段选型。")
    add_absorbed("REQ-HVACINS-021", "HVAC-INS-EXISTING-001", "HVAC", "final_segment_pipe_and_insulation_schedule", "unknown", "engineering_calculation", "pending", "OWNER-HVAC-INSULATION-PHOTOS-20260819", "yes", "由 M-401 按珠海露点、管温、导热系数、防火与逐段管径计算，再由日立／安装方按图复核。")

    add_absorbed("REQ-LAUNDRYTRAP-001", "PLUM-LAUNDRY-TRAP-001", "PLUM", "candidate_sku_and_topology", "Taobao item 703917654224 / SKU 5125125938318; two independent traps then downstream merge", "owner_candidate_seller_evidence", "candidate", "OWNER-PLUM-GEBERIT-DOUBLE-TRAP-20260819", "no", "优先比选，未采购，不是吉博力原厂成套批准。")
    add_absorbed("REQ-LAUNDRYTRAP-002", "PLUM-LAUNDRY-TRAP-001", "PLUM", "exact_parts_interfaces_water_seal_envelope_and_height", "unknown", "pending", "pending", "OWNER-PLUM-GEBERIT-DOUBLE-TRAP-20260819", "yes", "项目先出 P-201／I-503 联合剖面，再由商家与西门子对图确认。")
    add_absorbed("REQ-LAUNDRYTRAP-003", "PLUM-LAUNDRY-TRAP-001", "PLUM", "anti_crossflow_access_and_two_machine_test", "required", "project_requirement", "pending", "OWNER-PLUM-GEBERIT-DOUBLE-TRAP-20260819", "yes", "施工后同时排水试验，检查不串水、不溢水、可清扫和可更换。")

    for req_id, equipment_id, discipline, key, value, number, unit, origin, status, source_id, blocks, notes in absorbed_requirements:
        item = {field: "" for field in fields}
        item.update({
            "requirement_id": req_id, "equipment_id": equipment_id, "discipline": discipline,
            "parameter_key": key, "value_text": value, "value_number": number, "unit": unit,
            "value_origin": origin, "status": status, "source_id": source_id,
            "source_locator": "2026-08-20 SSOT absorption", "blocks_release": blocks, "notes": notes,
        })
        upsert_full(rows, "requirement_id", item)
    write_csv(path, fields, rows)


def reconcile_secondary_projections() -> None:
    drawing_path = DECISIONS / "drawing-register.csv"
    drawing_fields, drawing_rows = read_csv(drawing_path)
    drawing_notes = {
        "P-201": "Foster 1014850 官方参数已机械关闭；洗烘采用两台设备各自独存水弯后下游合流的双存水弯方向，优先比选 Taobao item 703917654224 / SKU 5125125938318。该组合未采购且包含非吉博力原厂配件；项目先出 P-201／I-503 联合剖面，商家与西门子再对图复核零件、接口、水封、防串水、标高、检修和双机试验。600 mm 石材盆由项目先画联合剖面。APP-016 保留 F50 基准候选并与 Franke Slim 50 CN 比较；说明书未给的中心和路线继续不猜。",
        "E-303": "设备插头、墙面插座、支路保护和独立回路分开表达。烤箱为 16A 插头／插座＋独立 C16；洗烘为两个 10A 插座共用一路 C16；冰淇淋机按家用候选包络预留。已购 H70FT 仅在台面使用位设可触及 10A 接地插座，纯存放柜格不默认加插座。既有强电箱为三江电气，手写 P230／PX1，铭牌 380/220 V、63 A；63 A 不得写成断路器整定或可增容结论，须由现场取证箱内断路器、导线、模数并计算负荷。",
        "E-304": "两个吸顶 AP 采用 Cat6 星型回弱电箱 PoE 拓扑；优先候选更新为 TP-Link TL-XAP1500GE-PoE/DC 易展版，准确外形 Ø184×40 mm、开孔 Ø155 mm、802.3at PoE、最大 11.1 W。现有对讲按 DNAKE 280M-S3 识别；型号不再询问，只补实际供电／端子、物业保留复装条件。弱电箱按有效安装包络、被动风道和可加装 120 mm 级风扇节点深化。",
        "M-401": "A01/A04 为开发商既有 RPIZ-22FSLN5QD/P；A02/A03 为既有 RPIZ-22FSLN5QDF/P，正式 IFC 镜像关系正确；A05 现场铭牌为 RPIZ-71FSLN5QD/P；A06 施工开放遮挡后再绑定铭牌。华美 Class 1 现场观察到 6×15、10×15、13×9、13×15、16×9、16×15、16×20、20×15、25×15、32×9、32×15 mm 的 ID×TK 规格；这些是既有现场观察，不是最终设计厚度。项目须按珠海露点、管温、导热系数、防火和逐段管径计算，分列裸管外径、保温壁厚和完成外径，再由日立／安装方按图复核。回风口兼检修口，冷凝水路线和管箍标高由项目出图、竣工验收。说明书未给的接口中心继续不猜。",
        "D-601": "干区墙脚采用业主已确认的原参考方向：墙面同色齐平宽踢脚，上下各约 10 mm 阴影缝，门套／隐形门／柜体连续通缝；不再作为 A/B 候选。D-601 仍须输出 1:5 截面、基层型材、公差、清洁背衬和转角展开并做 1:1 样板。中厨 VVD Pewter 60 mm 与湿区保持独立。",
        "S-701": "排程表同步登记 APP-020 惠人 H70FT（已购、到货未证）、八项已购厨房工具及 K-D01～K-D04 收纳分组、全屋窗饰方向、手摇开窗器改造、既有华美保温观察和洗烘双存水弯候选。夜帘、开窗器和双存水弯仍是候选，不得写成已选、已购或已批准；现场保温管套不得冒充最终设计厚度。",
    }
    for sheet, notes in drawing_notes.items():
        update_by_id(drawing_rows, "sheet_number", sheet, {"notes": notes})
    legacy_drawing_notes = {
        "A-105": "旧稿阳台塑木地板方向已登记；准确产品、完成面高度、排水、检修、收边、防火和耐候未关闭。公共区旧稿灰白洞石视觉与现行银白洞石岩板属于同向材料深化，不视为冲突。",
        "A-106": "旧稿厨房木饰面拼板、卫浴／阳台木条拼缝隐藏检修、反灯槽型材、200 mm 客厅窗帘盒和酒架上方加固方向已吸收。窗饰同步采用公共区 Luminette K5-501、卧室 Silhouette N31-204、厨卫 25 mm 灰黑铝百叶方向；夜帘为 LightLock 优先与 Duolite 对照，未终选或下单。完成面复尺、分幅、轨道、收拢、电源、色卡、报价和安装图待逐窗关闭。",
        "I-501": "旧稿酒架上方约 1300 mm 加固意图已登记；最终按自重、满载、顶部固定和吊顶体系出联合节点。已购惠人 H70FT 确认存放于岛台下，应画正面取放直立格／全拉出托盘、承重、防污垫和湿附件通风；八项已购工具按 K-D01 长工具、K-D02 Coravin、K-D03 锅具／锅盖、K-D04 HARIO 易碎器具画抽屉及摆样，不交给橱柜方重新分配。",
        "A-104": "全屋窗饰方向已确认但未下单，逐窗完成面、分幅和轨道等待厂家图关闭。既有窗增设纯机械手摇开窗器，不设电源和智能联动；约 22 mm 闭合边只作研究依据，须以准确产品和一樘非厨房窗样板关闭受力、孔距、排水腔避让、防水与保修。",
        "I-504": "旧稿次卧飘窗落地浅色木、客厅飘窗木盒与大白墙脱缝齐平的意图已登记；防潮、日晒、基层、伸缩、检修和窗边防水待平立剖与节点。",
    }
    for sheet, sentence in legacy_drawing_notes.items():
        row = next(item for item in drawing_rows if item["sheet_number"] == sheet)
        row["notes"] = append_sentence_once(row["notes"], sentence)
    write_csv(drawing_path, drawing_fields, drawing_rows)

    switch_path = DECISIONS / "e302-switch-product-review.csv"
    switch_fields, switch_rows = read_csv(switch_path)
    panel_candidates = {
        "CTRL-ENTRY-A": "JINK 候选产品族；至少 4 个基础功能，按两只相邻 86 模块预留；准确 SKU 未定",
        "CTRL-MASTER-A": "JINK 候选产品族；至少 2 个基础功能；准确 SKU 未定",
        "CTRL-MASTER-B": "JINK 候选产品族；至少 3 个基础功能并与 Entry A 形成实体有线双控；准确 SKU 未定",
    }
    for row in switch_rows:
        panel_id = row["panel_id"]
        if panel_id in panel_candidates:
            row["product_candidate"] = panel_candidates[panel_id]
            row["product_evidence_scope"] = "owner_function_count_confirmed_product_family_candidate_exact_sku_and_electrical_interface_pending"
            if panel_id == "CTRL-ENTRY-A":
                row["installation_context"] = "入户连续 2×86 模块预留区；最终可装一块或两块面板"
    write_csv(switch_path, switch_fields, switch_rows)

    detail_path = DECISIONS / "det1-detail-review.csv"
    detail_fields, detail_rows = read_csv(detail_path)
    update_by_id(detail_rows, "node_id", "D601-N04", {
        "confirmed_evidence": "中厨 VVD Pewter 60 mm 踢脚独立；干区已确认保留原参考：墙面同色齐平宽踢脚＋上下各约 10 mm 阴影缝，门套、隐形门和柜体连续通缝",
        "unresolved_material_or_product": "干区准确踢脚高度、厚度、材料、颜色、型材和清洁背衬待 D-601／样板；中厨只待项目样板与加工节点",
        "review_status": "owner_direction_confirmed_detail_pending_review",
    })
    write_csv(detail_path, detail_fields, detail_rows)

    wfin_path = DECISIONS / "wfin-open-issues.csv"
    wfin_fields, wfin_rows = read_csv(wfin_path)
    update_by_id(wfin_rows, "issue_id", "WFIN-R06", {
        "current_evidence": "业主已确认保留原参考：干区墙面同色齐平宽踢脚，上下各约 10 mm 阴影缝，墙下口型材和门墙柜连续通缝。中厨 VVD Pewter 60 mm 踢脚保持独立，湿区不套用。",
        "required_action_or_decision": "由 D-601 输出 1:5 墙脚节点、分房间适用表和门套／隐形门／柜体／转角展开，并以 1:1 实物样板关闭准确材料、尺寸、公差、耐撞、拖地水和积灰清洁。",
        "basis": "OWNER-WFIN-BASEBOARD-SHADOW-GAP-20260817；业主 2026-08-19 确认保留原参考做法。",
        "status": "owner_direction_confirmed_detail_and_mockup_pending",
        "stop_condition": "准确节点和样板批准前不得把约 10 mm 参考值直接作为加工公差、不得扩展到中厨或湿区、不得写正式 IFC 几何。",
    })
    write_csv(wfin_path, wfin_fields, wfin_rows)

    schedule_path = DECISIONS / "s701-schedule-review.csv"
    schedule_fields, schedule_rows = read_csv(schedule_path)
    update_by_id(schedule_rows, "schedule_id", "S701-APP-016", {
        "confirmed_scope": "使用位置=西厨 Foster 水槽柜；F50 未否决",
        "candidate_or_observed_scope": "F50=基准候选；Franke Slim 50 CN=365 mm 高比较候选；均未采购",
        "unresolved_for_release": "项目盆下联合剖面、准确铭牌功率、法兰、排水／洗碗机支管、拆换净空、噪声、保修和含安装报价",
        "review_status": "candidate_comparison_pending_project_section",
    })
    update_by_id(schedule_rows, "schedule_id", "S701-APP-019", {
        "confirmed_scope": "新增一台家用冰淇淋机；装修预留 500×450×450 mm、≥20 kg、10A 插座、散热，无给排水",
        "candidate_or_observed_scope": "Musso Mini Lussino 4080 推荐候选；Ninja CREAMi 220 V 备选；未选未购",
        "unresolved_for_release": "I-501 准确位置、取放路线和最终实购型号／插头复核；若改商用机须重做水电散热",
        "evidence_reference": "OWNER-INPUT-APP019-20260815;OUTBOUND-FORM-KITCHEN-20260818",
        "review_status": "project_envelope_confirmed_product_pending",
    })
    write_csv(schedule_path, schedule_fields, schedule_rows)


def main() -> None:
    reconcile_source_register()
    reconcile_owner_inputs()
    reconcile_closeout_rules()
    reconcile_equipment()
    reconcile_requirements()
    reconcile_secondary_projections()
    print(json.dumps({"status": "ok", "forms": len(FORM_SOURCES)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
