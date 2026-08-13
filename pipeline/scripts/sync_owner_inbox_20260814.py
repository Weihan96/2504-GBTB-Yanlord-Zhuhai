#!/usr/bin/env python3
"""Absorb the 2026-08-14 owner inbox snapshot into canonical project data.

The ignored inbox remains untouched.  This script only consumes the promoted,
hash-pinned evidence snapshot under drawings/evidence/.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "pipeline/decisions"
SNAPSHOT = ROOT / "drawings/evidence/OWNER-INPUT-20260814.md"
SNAPSHOT_SHA256 = "0dde421f925a8f16a0a52706db29b544d115cebd596b850cb9f18dea26da27ec"
IMAGE_DIR = ROOT / "drawings/evidence/owner-input-20260814"
BOM_FILES = {
    "appliance-input-register.csv",
    "elec-source-evidence.csv",
    "equipment-installation-requirements.csv",
    "equipment-register.csv",
    "source-evidence-register.csv",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    encoding = "utf-8-sig" if path.name in BOM_FILES else "utf-8"
    with path.open("w", encoding=encoding, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def by_key(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in rows}


def merge_ids(*values: str) -> str:
    result: list[str] = []
    for value in values:
        for item in (value or "").split(";"):
            item = item.strip()
            if item and item not in result:
                result.append(item)
    return ";".join(result)


def legacy_elec_projection(row: dict[str, str]) -> str:
    payload = {
        "evidence_id": row["source_id"],
        "discipline": row["discipline"],
        "sheet_id": row["sheet_id"],
        "decision_scope": row["decision_scope"],
        "source_kind": row["source_kind"],
        "source_document": row["source_url"] or row["local_path"] or row["source_document"],
        "source_sha256": row["sha256"],
        "source_locator": row["locator"],
        "evidence": row["evidence"],
        "proves": row["proves"],
        "does_not_prove": row["does_not_prove"],
        "status": row["status"],
        "confidence": row["confidence"],
        "review_required": row["review_required"],
        "formal_ifc_write_allowed": row["formal_ifc_write_allowed"],
        "notes": row["notes"],
    }
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))


def source_row(
    source_id: str,
    *,
    discipline: str,
    sheet_id: str,
    scope: str,
    kind: str,
    document: str,
    evidence: str,
    proves: str,
    does_not_prove: str,
    status: str,
    confidence: str = "1.00",
    review: str = "yes",
    url: str = "",
    local: str = "",
    digest: str = "not_applicable_live_reference",
    locator: str = "",
    manufacturer: str = "",
    model_scope: str = "",
    notes: str = "",
    project_to_elec: bool = False,
) -> dict[str, str]:
    row = {
        "source_id": source_id,
        "discipline": discipline,
        "sheet_id": sheet_id,
        "decision_scope": scope,
        "source_kind": kind,
        "source_document": document,
        "source_url": url,
        "local_path": local,
        "sha256": digest,
        "locator": locator,
        "evidence": evidence,
        "proves": proves,
        "does_not_prove": does_not_prove,
        "status": status,
        "confidence": confidence,
        "review_required": review,
        "formal_ifc_write_allowed": "no",
        "manufacturer": manufacturer,
        "model_scope": model_scope,
        "revision": "",
        "publication_date": "",
        "legacy_targets": "",
        "legacy_projection_json": "",
        "notes": notes,
    }
    if project_to_elec:
        row["legacy_targets"] = "elec-source-evidence.csv"
        row["legacy_projection_json"] = legacy_elec_projection(row)
    return row


def upsert_source(rows: list[dict[str, str]], row: dict[str, str]) -> None:
    index = by_key(rows, "source_id")
    if row["source_id"] in index:
        index[row["source_id"]].update(row)
    else:
        rows.append(row)


def next_requirement_id(rows: list[dict[str, str]]) -> str:
    maximum = max(int(re.search(r"(\d+)$", row["requirement_id"]).group(1)) for row in rows)
    return f"REQ-{maximum + 1:04d}"


def upsert_requirement(
    rows: list[dict[str, str]],
    equipment_id: str,
    parameter_key: str,
    *,
    discipline: str,
    value: str = "",
    unit: str = "",
    origin: str,
    status: str,
    source_id: str = "",
    blocks: str,
    notes: str = "",
) -> None:
    matches = [row for row in rows if row["equipment_id"] == equipment_id and row["parameter_key"] == parameter_key]
    row = matches[-1] if matches else {"requirement_id": next_requirement_id(rows)}
    numeric = value if re.fullmatch(r"-?\d+(?:\.\d+)?", value or "") else ""
    row.update({
        "equipment_id": equipment_id,
        "discipline": discipline,
        "parameter_key": parameter_key,
        "value_text": "" if numeric else value,
        "value_number": numeric,
        "unit": unit,
        "datum": "",
        "value_origin": origin,
        "status": status,
        "source_id": source_id,
        "source_locator": "",
        "blocks_release": blocks,
        "notes": notes,
    })
    if not matches:
        rows.append(row)


def append_equipment(rows: list[dict[str, str]], values: dict[str, str]) -> None:
    index = by_key(rows, "equipment_id")
    if values["equipment_id"] in index:
        index[values["equipment_id"]].update(values)
        return
    template = {key: "" for key in rows[0]}
    template.update(values)
    rows.append(template)


def append_named(rows: list[dict[str, str]], key: str, values: dict[str, str]) -> None:
    index = by_key(rows, key)
    if values[key] in index:
        index[values[key]].update(values)
    else:
        template = {field: "" for field in rows[0]}
        template.update(values)
        rows.append(template)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="verify evidence and report the planned scope without writing")
    args = parser.parse_args()

    if sha256(SNAPSHOT) != SNAPSHOT_SHA256:
        raise RuntimeError("promoted owner snapshot hash changed")
    image_paths = sorted(IMAGE_DIR.glob("*.webp"))
    if len(image_paths) != 18:
        raise RuntimeError(f"expected 18 promoted evidence images, found {len(image_paths)}")

    files = {
        name: read_csv(DECISIONS / name)
        for name in (
            "equipment-register.csv",
            "equipment-installation-requirements.csv",
            "source-evidence-register.csv",
            "owner-input-register.csv",
            "owner-input-closeout-rules.csv",
            "elec-design-rules.csv",
            "e302-switch-product-review.csv",
        )
    }
    equipment = files["equipment-register.csv"][1]
    requirements = files["equipment-installation-requirements.csv"][1]
    sources = files["source-evidence-register.csv"][1]
    owner_inputs = files["owner-input-register.csv"][1]
    closeout_rules = files["owner-input-closeout-rules.csv"][1]
    rules = files["elec-design-rules.csv"][1]
    switch_review = files["e302-switch-product-review.csv"][1]

    owner_source = "OWNER-INBOX-20260814-001"
    upsert_source(sources, source_row(
        owner_source,
        discipline="MULTI",
        sheet_id="OWNER-INPUT",
        scope="2026-08-14 业主补充输入完整快照",
        kind="owner_input_markdown_snapshot",
        document="业主输入-待同步.md",
        local="drawings/evidence/OWNER-INPUT-20260814.md",
        digest=SNAPSHOT_SHA256,
        locator="完整文档 01–09 节；原入口 tmp/owner-input-inbox/业主输入-待同步.md 保持不变",
        evidence="业主直接输入、明确确认、候选产品、研究结论与待确认事项的完整同步快照",
        proves="业主直接表达的功能、数量、使用场景、候选范围与确认程度",
        does_not_prove="厂家技术参数、商品在售状态、燃气公司准入、现场施工完成状态或未提供的型号接口",
        status="verified_owner_input_snapshot",
        review="no",
        notes="confirmed/tentative/unknown/research conclusion 必须按正文分别解释，不得把整份文档统一升级为确认。",
    ))

    image_groups = {
        "panasurface-switch-inset": ("E-302/INT1", "Panasurface 开关平嵌装饰件图片", "只证明装饰收口外观及与 86 型面板配合的商品展示", "不证明尺寸公差、JINK EGG 通用适配、底盒深度、阻燃、电气接线或施工批准"),
        "taobao-1062885293798": ("E-302", "JINK 2.5D Neo 淘宝商品主图", "只证明商品页展示 1–4 键、16A/20A/40A 等营销选项及 Matter 文案", "不证明准确 SKU、协议一致性、端子图、零线要求、证书映射或项目选定"),
        "taobao-964347781373": ("E-302", "JINK EGG 淘宝商品主图", "只证明商品页展示 EGG、Matter/Apple Home/场景控制等营销内容", "不证明准确 SKU、每路负载、端子图、底盒净深、传统有线双控或项目选定"),
        "taobao-992656775956": ("A-106", "Panasurface 报警器平嵌装饰件商品主图", "只证明装饰件存在多种设备名称和齐平外观展示且不含报警器", "不证明兼容尺寸、报警器型号、消防/燃气准入、安装位置或施工批准"),
    }
    for image_path in image_paths:
        stem = image_path.stem
        prefix = next(prefix for prefix in image_groups if stem.startswith(prefix))
        sheet, scope, proves, does_not = image_groups[prefix]
        sid = "OWNER-EVIDENCE-" + re.sub(r"[^A-Z0-9]+", "-", stem.upper()).strip("-")
        upsert_source(sources, source_row(
            sid,
            discipline="ELEC/INT1",
            sheet_id=sheet,
            scope=scope,
            kind="owner_supplied_product_image",
            document=str(image_path.relative_to(ROOT)),
            local=str(image_path.relative_to(ROOT)),
            digest=sha256(image_path),
            locator="业主 inbox evidence 原图按字节复制",
            evidence=f"{image_path.name} 本地视觉证据",
            proves=proves,
            does_not_prove=does_not,
            status="verified_local_owner_evidence_scope_limited",
            confidence="1.00",
            review="yes",
            notes="不得从无比例图片量取施工尺寸。",
            project_to_elec=True,
        ))

    web_sources = [
        source_row("TB-JINK-EGG-001", discipline="ELEC", sheet_id="E-302", scope="JINK EGG 淘宝商品页", kind="taobao_product_page", document="淘宝商品 964347781373", url="https://detail.tmall.com/item.htm?id=964347781373", evidence="商品页标题、属性型号 EGG 与零火三键/六键、锂电无线六键在售选项", proves="商品候选与在售 SKU 文案", does_not_prove="业主已选准确 SKU、官方端子图、每路负载、底盒净深或证书对应", status="commerce_page_research_candidate_only", confidence="0.90", manufacturer="JINK", model_scope="EGG", project_to_elec=True),
        source_row("JINK-EGG-OFFICIAL-001", discipline="ELEC", sheet_id="E-302", scope="JINK EGG 官方产品族", kind="official_product_web", document="JINK Switch EGG official product page", url="https://www.jinkhome.com/products/jink-switch-egg", evidence="官方产品族列出 3-Key/3-Relay、6-Key/3-Relay、Wireless Scene 6-Key、Matter over Thread、100–275V AC 与 85.8×86×37.55 mm 86 系列外形", proves="EGG 产品族功能、协议、供电和外形研究结论", does_not_prove="淘宝准确 SKU、每路/整机 2000W 口径、端子图、最低底盒净深、传统有线双控或项目选定", status="verified_official_family_research", confidence="0.95", manufacturer="JINK", model_scope="EGG family", project_to_elec=True),
        source_row("JINK-8E-CSA-001", discipline="ELEC", sheet_id="E-302", scope="JINK 8e Matter 证书边界", kind="official_certification_web", document="CSA JINK 8e Switch Series", url="https://csa-iot.org/csa_product/jink-8e-switch-series-3/", evidence="CSA 页面覆盖 JINK 8e 系列指定 Family SKU、Matter 1.1、Thread + Bluetooth", proves="8e 系列证书存在", does_not_prove="EGG 或 2.5D Neo 商品 SKU 与 8e Family SKU 相同或可借用该证书", status="verified_official_certificate_family_only", manufacturer="Longan Link/JINK", model_scope="JINK 8e family", project_to_elec=True),
        source_row("TB-JINK-2P5D-001", discipline="ELEC", sheet_id="E-302", scope="JINK 2.5D Neo 淘宝商品页", kind="taobao_product_page", document="淘宝商品 1062885293798", url="https://detail.tmall.com/item.htm?id=1062885293798", evidence="页面展示 1–3 GANG 16A 与单路 20A/40A 候选，属性与商品图对 Matter 传输协议表述冲突", proves="2.5D Neo 候选商品及页面内部协议冲突", does_not_prove="准确 SKU、最终传输协议、项目选定、接线图或证书", status="commerce_page_conflicting_research_candidate_only", confidence="0.80", manufacturer="JINK", model_scope="2.5D Neo", project_to_elec=True),
        source_row("TB-PANASURFACE-SWITCH-INSET-001", discipline="ELEC/INT1", sheet_id="E-302/I-501", scope="86 型开关平嵌装饰框", kind="taobao_product_page", document="淘宝商品 958552490234", url="https://item.taobao.com/item.htm?id=958552490234", evidence="商品标注 86 型、ABS、定制；一开/二开/三开表示相邻面板数量", proves="装饰框为饰面收口候选而非电气底盒", does_not_prove="未经商家按 EGG 实物确认即可通配、底盒深度、阻燃或接线", status="commerce_page_research_candidate_only", confidence="0.90", manufacturer="Panasurface", model_scope="86-type custom inset trim", project_to_elec=True),
        source_row("TB-PANASURFACE-AP-INSET-001", discipline="ELEC/INT1", sheet_id="E-304/RCP-1", scope="吸顶 AP 平嵌预埋件", kind="taobao_product_page", document="淘宝商品 1060053801442", url="https://item.taobao.com/item.htm?id=1060053801442", evidence="商品 SKU 仅明确 AP362E 开孔 215 与 RG-EAP262E 开孔 250", proves="两款指定 AP 的装饰预埋件候选", does_not_prove="兼容其他 AP、最终网络架构、安装深度、散热、承重或项目选定", status="commerce_page_research_candidate_only", confidence="0.90", manufacturer="Panasurface", model_scope="AP362E/RG-EAP262E inset", project_to_elec=True),
        source_row("HUAWEI-AP362E-OFFICIAL-001", discipline="ELEC", sheet_id="E-304", scope="Huawei eKit AP362E 官方参数", kind="official_product_web", document="Huawei eKit AP362E", url="https://ekit.huawei.com/ekit/front/ssr/en/product/1480064311765295872", evidence="约 Ø180×35 mm、1×GE、802.3af PoE、最大功耗约 9.4W", proves="AP362E 候选的机体包络、端口、PoE 与最大功耗", does_not_prove="业主已选、预埋件机械适配、PoE 交换机预算、无线覆盖或现场线缆通断", status="verified_official_exact_model_research", manufacturer="Huawei", model_scope="AP362E", project_to_elec=True),
        source_row("RUIJIE-EAP262E-OFFICIAL-001", discipline="ELEC", sheet_id="E-304", scope="Ruijie RG-EAP262(E) 官方供电方式", kind="official_product_web", document="Ruijie ceiling AP installation material", url="https://www.ruijie.com.cn/fw/wd/89870", evidence="RG-EAP262(E) 支持 802.3at PoE 或本地 12V DC", proves="RG-EAP262E 候选的可用供电方式", does_not_prove="业主已选、预埋件适配、最终采用本地电源或现场覆盖", status="verified_official_family_research", manufacturer="Ruijie", model_scope="RG-EAP262E", project_to_elec=True),
        source_row("APP-004-GS3-OFFICIAL-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-501", scope="La Marzocco GS3 候选", kind="official_product_pdf", document="La Marzocco GS3 official manual", url="https://lamarzocco.com/uk/en/wp-content/uploads/2023/06/MAN.2.1.01_GS3_EN_V1.3-1.pdf", evidence="约 410W×530D×355H mm、2120W，支持水箱与直连配置", proves="GS3 候选的产品包络、额定功率与供排水配置边界", does_not_prove="业主最终选择、实购版本、柜体接口或现场回路", status="verified_official_candidate_model", manufacturer="La Marzocco", model_scope="GS3", project_to_elec=True),
        source_row("APP-004-E1PRIMA-OFFICIAL-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-501", scope="Victoria Arduino E1 Prima EXP 候选", kind="official_product_pdf", document="Victoria Arduino E1 Prima EXP official quick guide", url="https://victoriaarduino.com/E1PrimaEXP_QuickGuide.pdf", evidence="约 411W×510D×379H mm；系列版本 1600/2600W；支持水箱与直连配置", proves="E1 Prima EXP 候选包络、版本功率范围与供排水配置边界", does_not_prove="中国 220V 实购版本功率、业主最终选择、柜体接口或现场回路", status="verified_official_candidate_model", manufacturer="Victoria Arduino", model_scope="E1 Prima EXP", project_to_elec=True),
        source_row("APP-014-SIEMENS-OFFICIAL-001", discipline="ELEC/PLUM", sheet_id="E-303/P-201/S-701", scope="Siemens WS7060BC1C/01 产品身份", kind="official_product_web", document="Siemens WS7060BC1C/01 official support page", url="https://www.siemens-home.bsh-group.com.hk/en/supportdetail/product/WS7060BC1C/01", evidence="西门子官方支持页存在准确 E-Nr. WS7060BC1C/01", proves="APP-014 准确厂家和 E-Nr. 身份", does_not_prove="业主已到货、额定功率、柜孔、给排水、电源和检修尺寸；这些仍须准确安装说明书", status="verified_official_exact_model_identity_only", manufacturer="Siemens", model_scope="WS7060BC1C/01", project_to_elec=True),
        source_row("APP-017-LG-OFFICIAL-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-501/S-701", scope="LG WashTower FN23BQH 官方身份", kind="official_product_web", document="LG FN23BQH official product page", url="https://www.lg.com/cn/washing-machines/lg-fn23bqh", evidence="13kg 洗+10kg 烘、600W×660D×1655H mm、曜岩黑", proves="FN23BQH 准确型号、容量和产品包络", does_not_prove="业主已购买、额定功率、进排水/插座坐标、柜门和搬运净空", status="verified_official_exact_model_research", manufacturer="LG", model_scope="FN23BQH", project_to_elec=True),
        source_row("TB-PANASURFACE-ALARM-INSET-001", discipline="ELEC/INT1", sheet_id="A-106", scope="烟感/燃气报警器平嵌装饰件", kind="taobao_product_page", document="淘宝商品 992656775956", url="https://item.taobao.com/item.htm?id=992656775956", evidence="商品列出多种设备名称且各 SKU 不含设备", proves="装饰收口候选存在", does_not_prove="任何报警器准确型号、兼容尺寸、消防/燃气准入或项目选定", status="commerce_page_research_candidate_only", confidence="0.90", manufacturer="Panasurface", model_scope="alarm inset trim", project_to_elec=True),
        source_row("GAS-HANWEI-KEA01-OFFICIAL-001", discipline="ELEC/GAS", sheet_id="A-106", scope="燃气公司咨询候选 Hanwei JT-KEA01", kind="official_product_web", document="Hanwei JT-KEA01 official page", url="https://hanwei.cn/pro_detail/JT-KEA01.html", evidence="家用甲烷探测、AC220V，页面列出阀门/输出选项和 GB15322.2-2019", proves="JT-KEA01 候选产品族公开参数", does_not_prove="珠海燃气公司批准、项目选定、准确输出配置或安装位置", status="official_research_candidate_for_authority_consultation_only", manufacturer="Hanwei", model_scope="JT-KEA01", project_to_elec=True),
        source_row("GAS-HANWEI-KWE-OFFICIAL-001", discipline="ELEC/GAS", sheet_id="A-106", scope="燃气公司咨询候选 Hanwei JT-KWE 系列", kind="official_product_web", document="Hanwei JT-KWE official page", url="https://hanwei.cn/pro_detail/JT-KWE.html", evidence="官方页面存在激光甲烷 JT-KWE 产品系列", proves="JT-KWE 可列为咨询候选系列", does_not_prove="准确子型号、珠海燃气公司批准、项目选定或施工接口", status="official_family_research_candidate_for_authority_consultation_only", manufacturer="Hanwei", model_scope="JT-KWE family", project_to_elec=True),
        source_row("GAS-CUBIC-AM5301-OFFICIAL-001", discipline="ELEC/GAS", sheet_id="A-106", scope="燃气公司咨询候选 Cubic JT-AM5301-JG", kind="official_product_web", document="Cubic JT-AM5301-JG official page", url="https://www.gassensor.com.cn/HomeAlarm/info_itemid_2007.html", evidence="家用甲烷探测，页面列出切断阀联动、NB-IoT 和 GB15322.2-2019", proves="JT-AM5301-JG 候选公开参数", does_not_prove="珠海燃气公司批准、项目选定、准确供电配置或安装位置", status="official_research_candidate_for_authority_consultation_only", manufacturer="Cubic", model_scope="JT-AM5301-JG", project_to_elec=True),
    ]
    for row in web_sources:
        upsert_source(sources, row)

    eq = by_key(equipment, "equipment_id")
    common_owner_ids = owner_source
    updates = {
        "APP-001": {"source_ids": merge_ids(eq["APP-001"]["source_ids"], common_owner_ids), "identity_basis": "业主确认 NS-01 使用场景；准确型号和功率未知", "notes": "NS-01 三面各一隐藏盖板插座；火锅可与 APP-002/APP-003 中一件同时使用。未提供型号和铭牌功率，不计算总负荷。"},
        "APP-002": {"source_ids": merge_ids(eq["APP-002"]["source_ids"], common_owner_ids), "identity_basis": "业主确认 NS-01 使用场景；准确型号和功率未知", "notes": "可与火锅同时使用；型号和铭牌功率保持 unknown。"},
        "APP-003": {"source_ids": merge_ids(eq["APP-003"]["source_ids"], common_owner_ids), "identity_basis": "业主确认 NS-01 使用场景；准确型号和功率未知", "notes": "可与火锅同时使用的备选小厨电；型号和铭牌功率保持 unknown。"},
        "APP-004": {"manufacturer": "La Marzocco / Victoria Arduino", "model": "GS3 / E1 Prima EXP（候选二选一）", "procurement_status": "candidate", "decision_status": "partial", "source_ids": merge_ids(common_owner_ids, "APP-004-GS3-OFFICIAL-001;APP-004-E1PRIMA-OFFICIAL-001"), "identity_basis": "业主给出两款候选；厂家资料仅作为候选参数，最终机型未定", "notes": "GS3 与 E1 Prima EXP 均为研究候选，不得作为已选产品。可预留可关闭净水支路和排水，但最终接口按实购机型关闭。"},
        "APP-005": {"manufacturer": "Hario", "model": "", "source_ids": merge_ids(common_owner_ids), "identity_basis": "业主提供 Hario 品牌和 1200W；准确型号待定", "notes": "1200W 为业主提供的候选功率，仍须按实购型号和铭牌复核；不与咖啡机同时使用。"},
        "APP-006": {"source_ids": merge_ids(eq["APP-006"]["source_ids"], common_owner_ids), "identity_basis": "业主确认 NS-02 使用场景；准确型号和功率未知", "notes": "与咖啡机同时使用；型号和铭牌功率保持 unknown。"},
        "APP-014": {"item_name": "双出水反渗透净饮一体机（直饮+净水）", "manufacturer": "Siemens", "model": "WS7060BC1C/01", "quantity": "1", "procurement_status": "selected", "decision_status": "partial", "source_ids": merge_ids(common_owner_ids, "APP-014-SIEMENS-OFFICIAL-001"), "identity_basis": "业主明确型号 WS7060BC1C；西门子官方支持页核实 E-Nr. WS7060BC1C/01；APP-015 为同一实物功能别名", "confidence": "1.00", "human_review_required": "yes", "notes": "一台实物同时承担直饮与净水功能；功率、给排水、电源、柜孔和检修条件仍须准确安装说明书。不得与 APP-015 重复统计。"},
        "APP-015": {"item_name": "净水功能别名（同 APP-014）", "manufacturer": "Siemens", "model": "WS7060BC1C/01", "quantity": "0", "procurement_status": "not_applicable", "decision_status": "superseded", "schedule_included": "no", "source_ids": merge_ids(common_owner_ids, "APP-014-SIEMENS-OFFICIAL-001"), "identity_basis": "业主明确 APP-015 与 APP-014 是同一台实物；本行仅为功能/历史编号别名", "confidence": "1.00", "human_review_required": "no", "notes": "Alias of APP-014。不得生成第二台设备，不重复统计电源、给水、排水、数量或回路。"},
        "APP-016": {"manufacturer": "勒科斯", "model": "F50（业主提供候选）", "procurement_status": "candidate", "decision_status": "partial", "source_ids": merge_ids(common_owner_ids), "identity_basis": "业主补充候选名称“勒科斯 F50”；尚无厂家安装图、订单或铭牌", "notes": "候选身份不等于已选/已购；额定功率、开关方式、法兰尺寸和排水接口保持 unknown。"},
        "APP-017": {"manufacturer": "LG", "model": "FN23BQH", "procurement_status": "candidate", "decision_status": "partial", "source_ids": merge_ids(eq["APP-017"]["source_ids"], common_owner_ids, "APP-017-LG-OFFICIAL-001"), "identity_basis": "正式 IFC 描述识别型号，LG 中国官网核实准确型号与包络；采购状态未确认", "notes": "13kg 洗+10kg 烘，600×660×1655mm 为官方产品包络；额定功率、进排水、插座、柜门和搬运净空仍待准确手册与现场复核。"},
        "SENSOR-001": {"source_ids": merge_ids(eq["SENSOR-001"]["source_ids"], common_owner_ids, "TB-PANASURFACE-ALARM-INSET-001"), "notes": "正式 IFC 只确认 A106-FIRE-R04 的 IfcSensor/FIRESENSOR 点位；最终感温/感烟/复合类型、准确产品、供电通信及装饰件适配仍待消防/设备方确认。"},
    }
    for equipment_id, values in updates.items():
        eq[equipment_id].update(values)

    new_equipment = [
        {"equipment_id":"CTRL-ENTRY-A","domain":"ASSEMBLY","category":"照明控制面板","item_name":"Entry A 智能照明控制面板","manufacturer":"JINK","model":"EGG 3-Key/3-Relay（研究优选，准确 SKU 未定）","variant":"wired_candidate","quantity":"1","procurement_status":"candidate","decision_status":"candidate","use_location_candidate":"入户门口墙面","use_location_confirmed":"","schedule_included":"yes","selector_kind":"logical_input","selector_value":"CTRL-ENTRY-A","source_ids":merge_ids(common_owner_ids,"TB-JINK-EGG-001;JINK-EGG-OFFICIAL-001;JINK-8E-CSA-001;TB-PANASURFACE-SWITCH-INSET-001"),"identity_basis":"业主确认两类照明控制用途；EGG 3-Key/3-Relay 仅为研究优选","confidence":"0.85","human_review_required":"yes","legacy_kind":"","legacy_id":"","notes":"准确 SKU、颜色、按键数、端子图、每路负载、传统有线双控、底盒净深和平嵌适配均未关闭。"},
        {"equipment_id":"CTRL-MASTER-A","domain":"ASSEMBLY","category":"照明控制面板","item_name":"Master A 主卧智能照明控制面板","manufacturer":"JINK","model":"EGG 3-Key/3-Relay 或 Wireless Scene 6-Key（研究候选）","variant":"direct_load_or_scene_only_candidate","quantity":"1","procurement_status":"candidate","decision_status":"candidate","use_location_candidate":"主卧内衣柜见光板","use_location_confirmed":"","schedule_included":"yes","selector_kind":"logical_input","selector_value":"CTRL-MASTER-A","source_ids":merge_ids(common_owner_ids,"TB-JINK-EGG-001;JINK-EGG-OFFICIAL-001;TB-PANASURFACE-SWITCH-INSET-001"),"identity_basis":"业主确认 Master A 独立角色；直接负载与纯场景方案尚未选择","confidence":"0.85","human_review_required":"yes","legacy_kind":"","legacy_id":"","notes":"见光板内直接承载负载前须核实阻燃背盒、固定、散热、检修、板厚开孔、各回路功率和启动浪涌。"},
        {"equipment_id":"CTRL-MASTER-B","domain":"ASSEMBLY","category":"照明控制面板","item_name":"Master B 公区双控面板","manufacturer":"JINK","model":"EGG 3-Key/3-Relay（研究优选，准确 SKU 未定）","variant":"wired_candidate","quantity":"1","procurement_status":"candidate","decision_status":"candidate","use_location_candidate":"主卧入口外墙面","use_location_confirmed":"","schedule_included":"yes","selector_kind":"logical_input","selector_value":"CTRL-MASTER-B","source_ids":merge_ids(common_owner_ids,"TB-JINK-EGG-001;JINK-EGG-OFFICIAL-001"),"identity_basis":"业主确认 Master B 与 Master A 为不同位置和角色；产品仅研究候选","confidence":"0.85","human_review_required":"yes","legacy_kind":"","legacy_id":"","notes":"客厅/书房/餐厅双控角色已确认；准确 SKU、传统有线双控拓扑、端子图和墙侧施工条件未关闭。"},
        {"equipment_id":"NET-AP-R09","domain":"NETWORK","category":"吸顶无线接入点","item_name":"主卧吸顶 AP","manufacturer":"Huawei / Ruijie","model":"AP362E / RG-EAP262E（候选）","variant":"PoE_recess_candidate","quantity":"1","procurement_status":"candidate","decision_status":"candidate","use_location_candidate":"R09 主卧天花候选点","use_location_confirmed":"","schedule_included":"yes","selector_kind":"logical_input","selector_value":"NET-AP-R09","source_ids":merge_ids(common_owner_ids,"TB-PANASURFACE-AP-INSET-001;HUAWEI-AP362E-OFFICIAL-001;RUIJIE-EAP262E-OFFICIAL-001"),"identity_basis":"预埋件只明确适配两款候选；若无锐捷生态暂优先研究 AP362E","confidence":"0.85","human_review_required":"yes","legacy_kind":"","legacy_id":"","notes":"最终型号/生态、CAT6 通断、PoE 交换机预算、覆盖、散热、检修和预埋件机械适配待关闭。"},
        {"equipment_id":"NET-AP-R14","domain":"NETWORK","category":"吸顶无线接入点","item_name":"次卧吸顶 AP","manufacturer":"Huawei / Ruijie","model":"AP362E / RG-EAP262E（候选）","variant":"PoE_recess_candidate","quantity":"1","procurement_status":"candidate","decision_status":"candidate","use_location_candidate":"R14 次卧天花候选点","use_location_confirmed":"","schedule_included":"yes","selector_kind":"logical_input","selector_value":"NET-AP-R14","source_ids":merge_ids(common_owner_ids,"TB-PANASURFACE-AP-INSET-001;HUAWEI-AP362E-OFFICIAL-001;RUIJIE-EAP262E-OFFICIAL-001"),"identity_basis":"预埋件只明确适配两款候选；若无锐捷生态暂优先研究 AP362E","confidence":"0.85","human_review_required":"yes","legacy_kind":"","legacy_id":"","notes":"最终型号/生态、CAT6 通断、PoE 交换机预算、覆盖、散热、检修和预埋件机械适配待关闭。"},
        {"equipment_id":"SENSOR-GAS-R04","domain":"SAFETY","category":"家用可燃气体探测器","item_name":"中厨燃气报警器","manufacturer":"","model":"","variant":"authority_approved_model_pending","quantity":"1","procurement_status":"not_selected","decision_status":"pending","use_location_candidate":"R04 中厨；准确位置待燃气公司","use_location_confirmed":"","schedule_included":"yes","selector_kind":"logical_input","selector_value":"SENSOR-GAS-R04","source_ids":merge_ids(common_owner_ids,"A106-GAS-OFFICIAL-001;GAS-HANWEI-KEA01-OFFICIAL-001;GAS-HANWEI-KWE-OFFICIAL-001;GAS-CUBIC-AM5301-OFFICIAL-001;TB-PANASURFACE-ALARM-INSET-001"),"identity_basis":"业主要求设置燃气报警器；准确型号须由当地燃气公司/主管要求确认","confidence":"1.00","human_review_required":"yes","legacy_kind":"","legacy_id":"","notes":"JT-GS838C-NBAC-H05、JT-KEA01/KEA31、JT-KWE 系列、JT-AM5301-JG 均只作为咨询候选，未获批准；不得锁定开孔、供电、联动或安装高度。"},
    ]
    for row in new_equipment:
        append_equipment(equipment, row)

    # Existing appliance requirements: preserve unknowns, add candidate-specific facts,
    # and mechanically remove APP-015 from all service counts.
    upsert_requirement(requirements,"APP-004","candidate_power_gs3",discipline="ELEC",value="2120",unit="W",origin="official_exact_model",status="candidate",source_id="APP-004-GS3-OFFICIAL-001",blocks="yes",notes="候选功率；最终未选 GS3 时不参与回路统计。")
    upsert_requirement(requirements,"APP-004","candidate_power_e1_prima_exp",discipline="ELEC",value="1600_or_2600_by_version",unit="W",origin="official_model_family",status="candidate",source_id="APP-004-E1PRIMA-OFFICIAL-001",blocks="yes",notes="须按中国 220V 实购版本锁值。")
    requirements[:] = [row for row in requirements if not (row["equipment_id"] == "APP-004" and row["parameter_key"] in {"candidate_water_mode", "candidate_drain_mode"})]
    upsert_requirement(requirements,"APP-004","candidate_gs3_water_mode",discipline="PLUM",value="reservoir_or_direct_plumb",origin="official_exact_model",status="candidate",source_id="APP-004-GS3-OFFICIAL-001",blocks="yes",notes="最终未选 GS3 时不参与施工接口。")
    upsert_requirement(requirements,"APP-004","candidate_gs3_drain_mode",discipline="PLUM",value="model_and_mode_dependent",origin="official_exact_model",status="candidate",source_id="APP-004-GS3-OFFICIAL-001",blocks="yes",notes="仅预留可关闭接口，不据此施工定版。")
    upsert_requirement(requirements,"APP-004","candidate_e1_prima_water_mode",discipline="PLUM",value="reservoir_or_direct_plumb",origin="official_exact_model",status="candidate",source_id="APP-004-E1PRIMA-OFFICIAL-001",blocks="yes",notes="最终未选 E1 Prima EXP 时不参与施工接口。")
    upsert_requirement(requirements,"APP-004","candidate_e1_prima_drain_mode",discipline="PLUM",value="model_and_mode_dependent",origin="official_exact_model",status="candidate",source_id="APP-004-E1PRIMA-OFFICIAL-001",blocks="yes",notes="仅预留可关闭接口，不据此施工定版。")
    upsert_requirement(requirements,"APP-005","rated_power",discipline="ELEC",value="1200",unit="W",origin="user_input",status="candidate",source_id=owner_source,blocks="yes",notes="业主提供；准确 Hario 型号和铭牌待补。")
    upsert_requirement(requirements,"APP-014","functional_alias",discipline="MULTI",value="APP-015",origin="user_input",status="confirmed",source_id=owner_source,blocks="no",notes="同一台实物的净水功能历史编号。")
    upsert_requirement(requirements,"APP-014","physical_unit_count",discipline="MULTI",value="1",unit="unit",origin="user_input",status="confirmed",source_id=owner_source,blocks="no",notes="APP-014/015 合计只统计一台。")
    upsert_requirement(requirements,"APP-014","exact_enr",discipline="MULTI",value="WS7060BC1C/01",origin="official_exact_model",status="confirmed",source_id="APP-014-SIEMENS-OFFICIAL-001",blocks="no")
    upsert_requirement(requirements,"APP-014","installation_manual",discipline="MULTI",value="unknown",origin="pending",status="pending",source_id="",blocks="yes",notes="需准确 E-Nr. 安装说明书关闭功率、柜孔、给排水、电源与检修。")
    for key, value in (("rated_power","0"),("simultaneous_group","APP-014_alias_no_count"),("water_required","not_applicable_alias"),("drain_required","not_applicable_alias"),("gas_required","not_applicable_alias"),("ventilation_required","not_applicable_alias")):
        discipline = "ELEC" if key in {"rated_power","simultaneous_group"} else "PLUM" if key in {"water_required","drain_required"} else "GAS" if key == "gas_required" else "HVAC"
        upsert_requirement(requirements,"APP-015",key,discipline=discipline,value=value,origin="user_input",status="not_applicable",source_id=owner_source,blocks="no",notes="功能别名，不新增实物或接口。")
    upsert_requirement(requirements,"APP-015","alias_of",discipline="MULTI",value="APP-014",origin="user_input",status="confirmed",source_id=owner_source,blocks="no",notes="唯一物理设备为 APP-014。")
    for key, discipline, unit in (("rated_power","ELEC","W"),("flange_size","PLUM","mm"),("drain_connection","PLUM","mm"),("switch_method","ELEC","")):
        upsert_requirement(requirements,"APP-016",key,discipline=discipline,value="",unit=unit,origin="pending",status="pending",blocks="yes",notes="勒科斯 F50 候选缺厂家安装图/铭牌，不按同类产品猜测。")
    for key, value, unit in (("product_width","600","mm"),("product_depth","660","mm"),("product_height","1655","mm"),("wash_capacity","13","kg"),("dry_capacity","10","kg")):
        upsert_requirement(requirements,"APP-017",key,discipline="INT1",value=value,unit=unit,origin="official_exact_model",status="confirmed",source_id="APP-017-LG-OFFICIAL-001",blocks="no")
    for equipment_id in ("CTRL-ENTRY-A","CTRL-MASTER-A","CTRL-MASTER-B"):
        upsert_requirement(requirements,equipment_id,"controlled_categories",discipline="ELEC",value="ambient_LED;accent_spotlights",origin="user_input",status="confirmed",source_id=owner_source,blocks="no")
        upsert_requirement(requirements,equipment_id,"exact_sku",discipline="ELEC",value="",origin="pending",status="pending",blocks="yes",notes="颜色、键数和有线/无线版本未锁。")
        upsert_requirement(requirements,equipment_id,"terminal_diagram",discipline="ELEC",value="",origin="pending",status="pending",blocks="yes",notes="零火文案不能替代端子定义和接线图。")
        upsert_requirement(requirements,equipment_id,"minimum_backbox_clear_depth",discipline="ELEC/INT1",value="",unit="mm",origin="pending",status="pending",blocks="yes",notes="须包含端子及导线弯曲空间。")
        upsert_requirement(requirements,equipment_id,"candidate_face_width",discipline="INT1",value="85.8",unit="mm",origin="official_model_family",status="candidate",source_id="JINK-EGG-OFFICIAL-001",blocks="yes")
        upsert_requirement(requirements,equipment_id,"candidate_face_height",discipline="INT1",value="86",unit="mm",origin="official_model_family",status="candidate",source_id="JINK-EGG-OFFICIAL-001",blocks="yes")
        upsert_requirement(requirements,equipment_id,"candidate_product_depth",discipline="INT1",value="37.55",unit="mm",origin="official_model_family",status="candidate",source_id="JINK-EGG-OFFICIAL-001",blocks="yes",notes="产品深度不等于底盒最小净深。")
        upsert_requirement(requirements,equipment_id,"inset_trim_fit",discipline="INT1",value="vendor_sample_confirmation_required",origin="project_candidate",status="candidate",source_id="TB-PANASURFACE-SWITCH-INSET-001",blocks="yes",notes="装饰框不替代电气底盒；一开/二开/三开指相邻面板数量。")
    for equipment_id in ("NET-AP-R09","NET-AP-R14"):
        upsert_requirement(requirements,equipment_id,"final_model",discipline="ELEC",value="",origin="pending",status="pending",blocks="yes")
        upsert_requirement(requirements,equipment_id,"home_run_cable",discipline="ELEC",value="CAT6_or_higher_to_weak_current_box",origin="project_candidate",status="candidate",source_id=owner_source,blocks="yes")
        upsert_requirement(requirements,equipment_id,"cable_continuity_test",discipline="ELEC",value="",origin="pending",status="pending",blocks="yes")
        upsert_requirement(requirements,equipment_id,"local_220v_for_current_candidates",discipline="ELEC",value="not_applicable_if_PoE",origin="project_candidate",status="candidate",source_id=owner_source,blocks="yes",notes="若改 FTTR/本地供电设备则重新设计，不能继承。")
        upsert_requirement(requirements,equipment_id,"candidate_ap362e_max_power",discipline="ELEC",value="9.4",unit="W",origin="official_exact_model",status="candidate",source_id="HUAWEI-AP362E-OFFICIAL-001",blocks="yes")
        upsert_requirement(requirements,equipment_id,"candidate_ap362e_power",discipline="ELEC",value="802.3af_PoE",origin="official_exact_model",status="candidate",source_id="HUAWEI-AP362E-OFFICIAL-001",blocks="yes")
        upsert_requirement(requirements,equipment_id,"candidate_rgeap262e_power",discipline="ELEC",value="802.3at_PoE_or_12V_DC",origin="official_exact_model",status="candidate",source_id="RUIJIE-EAP262E-OFFICIAL-001",blocks="yes")
    for key in ("exact_model","gas_company_approval","certification","supply","valve_interface","installation_height","clearance_to_appliances_and_openings","inset_trim_written_approval"):
        upsert_requirement(requirements,"SENSOR-GAS-R04",key,discipline="ELEC/GAS",value="",origin="pending",status="pending",blocks="yes",notes="须由珠海燃气公司/主管要求和最终设备资料关闭。")
    upsert_requirement(requirements,"SENSOR-GAS-R04","consultation_candidates",discipline="ELEC/GAS",value="JT-GS838C-NBAC-H05;JT-KEA01;JT-KEA31;JT-KWE_family;JT-AM5301-JG",origin="project_candidate",status="candidate",source_id=owner_source,blocks="yes",notes="仅供咨询，不代表批准或采购。")
    upsert_requirement(requirements,"SENSOR-001","final_detector_type_and_model",discipline="ELEC/FIRE",value="",origin="pending",status="pending",blocks="yes",notes="位置已确认；感温/感烟/复合和产品仍待消防/设备方。")
    upsert_requirement(requirements,"SENSOR-001","inset_trim_fit",discipline="INT1/FIRE",value="vendor_and_fire_authority_confirmation_required",origin="project_candidate",status="candidate",source_id="TB-PANASURFACE-ALARM-INSET-001",blocks="yes",notes="先定设备再按准确型号确认装饰件。")

    owner_updates = {
        "E302-ENTRY-PANEL": {"candidate_value":"JINK EGG 3-Key/3-Relay 作为两类照明直接控制的研究优选；准确 SKU 未定","user_value":"每个空间分别控制氛围 LED 与重点射灯；EGG 商品 ID 964347781373，颜色和键数未锁","status":"需证据","evidence_reference":merge_ids(owner_source,"TB-JINK-EGG-001;JINK-EGG-OFFICIAL-001;JINK-8E-CSA-001;TB-PANASURFACE-SWITCH-INSET-001"),"source_basis":"业主确认功能；官方产品族与淘宝页仅支持研究候选，CSA 8e 证书不得借给 EGG","notes":"先逐回路计算实际长度×W/m及驱动输入/数量并核对启动浪涌；2000W max 的整机/单路口径、端子图、底盒净深、传统有线双控和平嵌适配均未关闭。"},
        "E302-MASTER-PANEL": {"candidate_value":"Master A 直接负载时评估 EGG 3-Key/3-Relay，纯场景时评估 Wireless Scene 6-Key；Master B 研究优选 EGG 有线版","user_value":"Master A 与 Master B 为不同位置；各自控制已定义的两类照明/公区双控","status":"需证据","evidence_reference":merge_ids(owner_source,"TB-JINK-EGG-001;JINK-EGG-OFFICIAL-001;TB-PANASURFACE-SWITCH-INSET-001"),"source_basis":"业主确认角色；产品准确 SKU、接线和见光板施工条件未确认","notes":"Panasurface 套件只做收口，不替代 86 底盒；SKU 一开/二开/三开是相邻面板数。"},
        "E303-NS01-FORM": {"source_basis":"业主确认三个面各一隐藏盖板插座；功率必须逐台取证后再做回路统筹","notes":"APP-001～003 功率未齐，不计算总负荷、不锁回路数量、线径、保护器或 RCD。"},
        "E303-NS01-CIRCUIT": {"candidate_value":"先逐台记录 APP-001～003 额定功率；设计工况为火锅+APP-002/003 中一件","user_value":"不要假设，把每个功率都列出来，最后统一设计统筹","status":"待填写","evidence_reference":owner_source,"source_basis":"业主明确纠正旧聚合估值与两回路预设","notes":"旧 1.8–3.7kW 包络与 2 个独立回路候选失效；功率齐全前保持 pending。"},
        "E303-NS02-CIRCUIT": {"candidate_value":"咖啡机+磨豆机或手冲壶单独使用；先完成咖啡机二选一和磨豆机功率","user_value":"咖啡机与手冲壶不同时使用；Hario 1200W，准确型号待定","status":"待填写","evidence_reference":merge_ids(owner_source,"APP-004-GS3-OFFICIAL-001;APP-004-E1PRIMA-OFFICIAL-001"),"source_basis":"业主确认使用工况；候选咖啡机厂家功率不同，磨豆机功率未知","notes":"旧单回路容量候选不再关闭；按最终机型逐台功率统筹。"},
        "E304-AP-POWER": {"candidate_value":"AP362E 或 RG-EAP262E 均优先 PoE；CAT6 回弱电箱 PoE 交换机","user_value":"预埋件适配 AP362E 与 RG-EAP262E；若无锐捷生态暂优先评估 AP362E","status":"需证据","evidence_reference":merge_ids(owner_source,"TB-PANASURFACE-AP-INSET-001;HUAWEI-AP362E-OFFICIAL-001;RUIJIE-EAP262E-OFFICIAL-001"),"source_basis":"厂家资料证明两候选 PoE 能力；商品页只证明对应开孔候选","notes":"最终型号、生态、测线、PoE 预算、覆盖和预埋件施工条件未关闭；FTTR 需另选准确设备和预埋件。"},
        "E304-CABINET-DIMENSIONS": {"candidate_value":"弱电箱净深约 110mm；距入户墙 525mm","user_value":"已纠正：525mm，不是 525m；弱电箱净深约 110mm","unit":"mm","status":"需证据","evidence_reference":owner_source,"source_basis":"业主直接纠错和现场约测","notes":"525mm 已确认；110mm 为约值。设备厚度、插头突出、线缆弯曲、散热和检修净空仍待实测。"},
        "A106-FIRE-TYPE": {"candidate_value":"先由消防/设备方确定感温/感烟/复合及准确产品，再按实物确认装饰件","user_value":"准确型号 unknown；不因预埋件反向锁定报警器","status":"需证据","evidence_reference":merge_ids(owner_source,"A106-FIRE-OFFICIAL-001;A106-FIRE-OFFICIAL-002;A106-FIRE-OFFICIAL-003;TB-PANASURFACE-ALARM-INSET-001"),"source_basis":"业主最新输入撤回拼接型号作为选型依据；既有厂家核验只作负面匹配","notes":"正式 IFC 仅保留已确认点位；最终类型、产品、供电通信和厂家安装条件仍待确认。"},
        "A106-GAS-ALARM": {"candidate_value":"向珠海燃气公司书面确认准入型号、认证、联动阀、供电、安装位置与平嵌收口许可","user_value":"准确型号 unknown；列出的 JT-GS838C-NBAC-H05、汉威和四方光电型号只供燃气公司咨询","status":"需证据","evidence_reference":merge_ids(owner_source,"A106-GAS-OFFICIAL-001;GAS-HANWEI-KEA01-OFFICIAL-001;GAS-HANWEI-KWE-OFFICIAL-001;GAS-CUBIC-AM5301-OFFICIAL-001;TB-PANASURFACE-ALARM-INSET-001"),"source_basis":"业主要求燃气公司先确认；候选厂家页面不等于当地准入","notes":"不得把任何候选标记为批准；装饰件不含设备且不承担准入证明。"},
    }
    for input_id, values in owner_updates.items():
        by_key(owner_inputs,"input_id")[input_id].update(values)
    append_named(owner_inputs,"input_id",{"input_id":"APP014-APP015-ALIAS","workstream":"E-303/P-201/S-701","priority":"P0","blocks_release":"yes","question":"APP-014 西门子净饮机安装接口和 APP-015 别名去重","candidate_value":"APP-014 为唯一实物；APP-015 仅保留功能别名且数量为 0","user_value":"APP-014 与 APP-015 是同一台 Siemens WS7060BC1C/01","status":"需证据","evidence_reference":merge_ids(owner_source,"APP-014-SIEMENS-OFFICIAL-001"),"source_basis":"业主确认同一实物；官网确认 E-Nr. 身份","sync_target":"equipment SSOT","notes":"已机械去重电源、给水、排水和设备数量；仍缺准确安装说明书。"})
    append_named(owner_inputs,"input_id",{"input_id":"APP004-COFFEE-FINAL","workstream":"E-303/P-201/I-501","priority":"P0","blocks_release":"yes","question":"咖啡机最终选择 GS3 还是 E1 Prima EXP，并确认实购 220V 版本","candidate_value":"两款均保留水箱/直连兼容预留","user_value":"GS3 与 E1 Prima EXP 二选一，尚未决定","status":"待填写","evidence_reference":merge_ids(owner_source,"APP-004-GS3-OFFICIAL-001;APP-004-E1PRIMA-OFFICIAL-001"),"source_basis":"厂家候选资料","sync_target":"APP-004","notes":"不得把候选功率写成最终额定功率。"})
    append_named(owner_inputs,"input_id",{"input_id":"APP016-DISPOSER-DATA","workstream":"E-303/P-201/S-701","priority":"P0","blocks_release":"yes","question":"补勒科斯 F50 订单/铭牌/厂家安装图","candidate_value":"保持候选，不按同类产品估算","user_value":"勒科斯 F50","status":"需证据","evidence_reference":owner_source,"source_basis":"仅有业主候选名称","sync_target":"APP-016","notes":"额定功率、开关方式、法兰尺寸和排水接口均 unknown。"})
    append_named(owner_inputs,"input_id",{"input_id":"E302-BEDSIDE-220V","workstream":"E-302","priority":"P0","blocks_release":"yes","question":"各床位结合床/床头柜/软包深化准确点位、高度与面板组合","candidate_value":"每个可用床侧至少 2 个可同时使用的 220V 插座位","user_value":"confirmed；USB/Type-C 不计入 2 个 220V 插座位","status":"自定义确认","evidence_reference":owner_source,"source_basis":"业主直接确认","sync_target":"E-302 bedside sockets","notes":"数量规则已确认，准确点位/高度仍须结合双人床、床头柜、软包与窗帘深化。"})

    closeout_rows = [
        {"input_id":"APP014-APP015-ALIAS","closeout_kind":"product_interface_evidence","responsible_party":"设备方/给排水设计/电气设计","required_evidence":"WS7060BC1C/01 准确安装说明书，含功率、柜孔、给水、排水、电源与检修","automatic_close_allowed":"no","notes":"一台实物和别名去重已确认；施工接口仍须厂家资料。"},
        {"input_id":"APP004-COFFEE-FINAL","closeout_kind":"human_product_selection","responsible_party":"业主/设备方","required_evidence":"GS3 或 E1 Prima EXP 最终订单、准确 220V 版本铭牌及安装说明书","automatic_close_allowed":"no","notes":"候选功率不得代替最终铭牌。"},
        {"input_id":"APP016-DISPOSER-DATA","closeout_kind":"product_interface_evidence","responsible_party":"设备方/橱柜设计/给排水设计","required_evidence":"勒科斯 F50 订单或铭牌及厂家安装图，含功率、开关、法兰和排水接口","automatic_close_allowed":"no","notes":"同类产品常见值不能关闭。"},
        {"input_id":"E302-BEDSIDE-220V","closeout_kind":"site_or_drawing_evidence","responsible_party":"建筑设计/家具设计/电气设计","required_evidence":"逐房床、床头柜、软包和窗帘完成尺寸下的插座点位、高度与面板组合图","automatic_close_allowed":"no","notes":"每侧至少两个 220V 数量规则已确认；本规则只关闭准确定位。"},
    ]
    for row in closeout_rows:
        append_named(closeout_rules,"input_id",row)

    rule_index = by_key(rules,"rule_id")
    rule_index["ELEC-DES-001"].update({"value":"physical_wall_controls_required;smart_connected_features_candidate","basis":"业主仍要求实体墙面控制，同时补充 JINK Matter 智能面板候选；智能协议和场景尚未定版","review_required":"yes","status":"physical_interface_confirmed_smart_features_pending","notes":"实体按键必须可用；原“不采用智能场景”表述被最新业主输入部分取代。Matter 不是 KNX，准确 SKU 与行为验收未关闭。"})
    rule_index["ELEC-DES-020"].update({"value":"no_model_or_recess_locked","basis":"业主最新明确先定报警器准确型号再选同名预埋件，不得由装饰件反向选设备","status":"superseded_by_device_first_selection","notes":"旧小米烟感及 Ø105/Ø140×2/30mm 候选包络不得进入正式施工实体或开孔。"})
    rule_index["ELEC-DES-022"].update({"value":"no_model_or_recess_locked","basis":"业主最新明确燃气报警器型号 unknown 且须由燃气公司确认；装饰件不承担准入证明","status":"superseded_by_gas_authority_gate","notes":"旧小米燃气报警器及 Ø94/Ø120×2/43mm 候选包络不得锁定；由 ELEC-DES-041 取代。"})
    rule_index["ELEC-DES-032"].update({"basis":"业主确认 NS-01 三个面各一隐藏盖板插座；要求逐台记录功率后再统筹","confidence":"1.00","status":"user_form_confirmed_power_data_pending","notes":"APP-001～003 功率未齐，不再使用旧 2.6–4.9kW 聚合包络。"})
    rule_index["ELEC-DES-033"].update({"value":"pending_per_device_power_and_circuit_coordination","basis":"业主纠正：先逐台记录 APP-001～003 功率，再按火锅+一件小厨电同时工况统一回路统筹","confidence":"1.00","status":"pending_power_data_no_circuit_count_locked","notes":"旧 2_independent_circuits 候选已撤销；不得在功率齐全前锁定回路数量、线径、保护器或 RCD。"})
    rule_index["ELEC-DES-034"].update({"basis":"业主确认连接需求与使用位置；咖啡机二选一、手冲壶准确型号和磨豆机功率尚未关闭","status":"user_form_confirmed_product_and_power_pending","notes":"不再使用旧 2.35–4.0kW 聚合包络；最终模块数、墙侧、标高和线缆遮挡由餐边柜立面关闭。"})
    rule_index["ELEC-DES-035"].update({"value":"pending_selected_coffee_machine_and_grinder_power","basis":"业主确认咖啡机+磨豆机或手冲壶单独使用；GS3/E1 Prima EXP 候选功率不同且磨豆机未知","confidence":"1.00","status":"pending_power_data_no_circuit_count_locked","notes":"旧 1_circuit_capacity_candidate 已撤销；逐台功率齐全后再统一回路统筹。"})
    rule_index["ELEC-DES-036"].update({"value":"egg_family_preferred_for_entry_master;exact_sku_pending;2_5d_neo_alternative_only","basis":"业主确认各面板控制角色；JINK 官方 EGG 产品族与淘宝页支持研究优选；2.5D Neo 页面协议信息冲突","confidence":"0.85","status":"research_conclusion_exact_sku_wiring_and_certificate_pending","notes":"EGG/2.5D Neo 均未选定准确 SKU；CSA 8e 证书不得借用；端子图、每路负载、传统有线双控、底盒净深和平嵌适配未关闭。"})
    new_rules = [
        {"rule_id":"ELEC-DES-037","rule_kind":"smart_switch_load_calculation","locations":"R01;R10;ALL_LIGHTING_CIRCUITS","controlled_or_served_scope":"JINK EGG 候选继电器所控 LED/射灯回路","device_or_datum":"per_circuit_connected_load_and_inrush","value":"actual_length_x_W_per_m_plus_driver_input_and_inrush","basis":"业主纠正负载算法；示例 5W/m×200m=1000W 仅作算术说明","confidence":"1.00","review_required":"yes","status":"confirmed_method_product_limits_pending","automatic_ifc_write_allowed":"no","notes":"不得把 1000W 当作项目已测负荷；官方 16A/2000W max 的整机或单路口径未明，须逐回路核实。"},
        {"rule_id":"ELEC-DES-038","rule_kind":"network_candidate_architecture","locations":"R09;R14","controlled_or_served_scope":"主卧与次卧吸顶 AP","device_or_datum":"PoE_AP_candidate","value":"AP362E_or_RG-EAP262E;CAT6_home_run;no_local_220V_for_current_candidates","basis":"预埋件只明确两款候选；厂家资料证明 PoE 能力；业主接受网线回弱电箱作为推荐方案","confidence":"0.90","review_required":"yes","status":"research_conclusion_final_ecosystem_and_testing_pending","automatic_ifc_write_allowed":"no","notes":"若无锐捷控制器/网关生态暂优先研究 AP362E；不是华为 FTTR 光从机方案。"},
        {"rule_id":"ELEC-DES-039","rule_kind":"network_existing_measurement","locations":"R01","controlled_or_served_scope":"玄关弱电箱","device_or_datum":"entry_wall_distance_and_box_depth","value":"entry_wall_distance_525mm;box_depth_approx_110mm","basis":"业主明确纠正 525m 为 525mm，并提供弱电箱约 110mm 深度","confidence":"1.00","review_required":"yes","status":"owner_observed_dimensions_device_clearance_pending","automatic_ifc_write_allowed":"no","notes":"110mm 为约值；插头突出、线缆弯曲、设备厚度、散热与检修净空仍待实测。"},
        {"rule_id":"ELEC-DES-040","rule_kind":"bedside_socket_minimum","locations":"ALL_BEDS","controlled_or_served_scope":"所有左右两侧均可使用的床位","device_or_datum":"simultaneously_usable_220V_socket_positions_per_side","value":"minimum_2_per_usable_side;USB_TypeC_excluded","basis":"业主 2026-08-14 直接确认","confidence":"1.00","review_required":"yes","status":"confirmed_quantity_exact_positions_pending","automatic_ifc_write_allowed":"no","notes":"双人床两侧合计至少 4 个；结合床、床头柜、软包和窗帘确定准确点位/高度，插头不得妨碍贴墙和抽屉。"},
        {"rule_id":"ELEC-DES-041","rule_kind":"gas_alarm_authority_gate","locations":"R04","controlled_or_served_scope":"中厨家用可燃气体报警器","device_or_datum":"authority_approved_exact_model_and_installation","value":"required_before_selection_and_recess_detail","basis":"业主要求先向当地燃气公司确认；候选厂家资料不等于珠海准入","confidence":"1.00","review_required":"yes","status":"confirmed_process_exact_model_unknown","automatic_ifc_write_allowed":"no","notes":"候选型号仅供咨询；须关闭燃气种类、认证、供货安装责任、切断阀、供电、净距、平嵌许可和验收。"},
    ]
    for row in new_rules:
        append_named(rules,"rule_id",row)

    for row in switch_review:
        panel = row["panel_id"]
        if panel == "CTRL-ENTRY-A":
            row["product_candidate"] = "JINK Switch EGG 3-Key/3-Relay（研究优选，准确 SKU 未定）"
        elif panel == "CTRL-MASTER-A":
            row["product_candidate"] = "JINK EGG 3-Key/3-Relay 或 Wireless Scene 6-Key（按直接负载/纯场景二选一）"
        elif panel == "CTRL-MASTER-B":
            row["product_candidate"] = "JINK Switch EGG 3-Key/3-Relay（研究优选，准确 SKU 未定）"
        else:
            continue
        row["product_evidence_scope"] = "official_family_and_taobao_listing_only_exact_sku_terminal_diagram_unclosed"
        row["passed"] = "false"
        if row["gate"] == "neutral_and_wiring_diagram_verified":
            row["failure_reason"] = "零火文案只证明有线候选需 L/N；准确 SKU 端子定义、官方接线图及现场盒内线况未核实"
            row["required_evidence"] = "准确 SKU 官方端子图、零线要求、端子定义、底盒内线况及传统有线双控原理"
        elif row["gate"] == "rated_load_schedule_verified":
            row["failure_reason"] = "须按各回路实际长度×W/m、驱动输入与数量计算并核对启动浪涌；2000W max 的整机/单路口径未明"
            row["required_evidence"] = "准确 SKU 每路/整机额定口径、逐回路灯带长度与 W/m、驱动器输入/数量/浪涌及兼容性确认"
        elif row["gate"] == "box_and_joinery_interface_verified":
            row["failure_reason"] = "Panasurface 86 型装饰件只作平嵌收口且不是底盒；实物适配、阻燃背盒、端子弯线净深、固定与检修未关闭"
            row["required_evidence"] = "EGG 实物/完整尺寸图、Panasurface 书面适配确认、JINK 底盒净深与端子空间、基层和收口节点"

    if args.check:
        print(json.dumps({"snapshot_sha256": SNAPSHOT_SHA256, "evidence_images": len(image_paths), "mode": "check", "planned_tables": list(files)}, ensure_ascii=False, indent=2))
        return

    for name, (fields, _original_rows) in files.items():
        rows = {
            "equipment-register.csv": equipment,
            "equipment-installation-requirements.csv": requirements,
            "source-evidence-register.csv": sources,
            "owner-input-register.csv": owner_inputs,
            "owner-input-closeout-rules.csv": closeout_rules,
            "elec-design-rules.csv": rules,
            "e302-switch-product-review.csv": switch_review,
        }[name]
        write_csv(DECISIONS / name, fields, rows)

    print(json.dumps({
        "snapshot_sha256": SNAPSHOT_SHA256,
        "evidence_images": len(image_paths),
        "equipment_rows": len(equipment),
        "requirement_rows": len(requirements),
        "source_rows": len(sources),
        "owner_input_rows": len(owner_inputs),
        "closeout_rule_rows": len(closeout_rules),
        "design_rule_rows": len(rules),
        "switch_review_rows": len(switch_review),
        "mode": "write",
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
