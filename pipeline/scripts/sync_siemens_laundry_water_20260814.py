#!/usr/bin/env python3
"""Synchronize the 2026-08-14 Siemens laundry and drinking-water decisions.

This updates the canonical equipment SSOT plus owner-input closeout records.
It does not write IFC; the IFC batch is deliberately separate so geometry QA
can freeze and compare the formal model around that write boundary.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from equipment_ssot import projections, validate


ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "pipeline/decisions"
EVIDENCE_DIR = ROOT / "drawings/evidence"


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def merge_ids(*values: str) -> str:
    result: list[str] = []
    for value in values:
        for item in value.replace("；", ";").split(";"):
            item = item.strip()
            if item and item not in result:
                result.append(item)
    return ";".join(result)


def replace_row(rows: list[dict[str, str]], key: str, value: str, updates: dict[str, str]) -> None:
    matches = [row for row in rows if row[key] == value]
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one {key}={value}, got {len(matches)}")
    matches[0].update(updates)


def upsert_row(rows: list[dict[str, str]], key: str, value: str, payload: dict[str, str]) -> None:
    matches = [row for row in rows if row[key] == value]
    if len(matches) > 1:
        raise RuntimeError(f"duplicate {key}={value}")
    if matches:
        matches[0].update(payload)
    else:
        rows.append({key: value, **payload})


def evidence_projection(source: dict[str, str]) -> str:
    fields = [
        "evidence_id", "discipline", "sheet_id", "decision_scope", "source_kind",
        "source_document", "source_sha256", "source_locator", "evidence", "proves",
        "does_not_prove", "status", "confidence", "review_required",
        "formal_ifc_write_allowed", "notes",
    ]
    payload = {
        "evidence_id": source["source_id"],
        "discipline": source["discipline"],
        "sheet_id": source["sheet_id"],
        "decision_scope": source["decision_scope"],
        "source_kind": source["source_kind"],
        "source_document": source["source_url"] or source["source_document"],
        "source_sha256": source["sha256"],
        "source_locator": source["locator"],
        "evidence": source["evidence"],
        "proves": source["proves"],
        "does_not_prove": source["does_not_prove"],
        "status": source["status"],
        "confidence": source["confidence"],
        "review_required": source["review_required"],
        "formal_ifc_write_allowed": source["formal_ifc_write_allowed"],
        "notes": source["notes"],
    }
    return json.dumps({field: payload[field] for field in fields}, ensure_ascii=False, separators=(",", ":"))


def source_row(
    source_id: str,
    *,
    discipline: str,
    sheet_id: str,
    scope: str,
    kind: str,
    document: str,
    url: str,
    locator: str,
    evidence: str,
    proves: str,
    does_not_prove: str,
    status: str,
    manufacturer: str,
    model_scope: str,
    review_required: str = "no",
    formal_ifc_write_allowed: str = "yes",
    notes: str = "",
) -> dict[str, str]:
    local = EVIDENCE_DIR / document
    if not local.is_file():
        raise RuntimeError(f"missing evidence file: {local}")
    row = {
        "source_id": source_id,
        "discipline": discipline,
        "sheet_id": sheet_id,
        "decision_scope": scope,
        "source_kind": kind,
        "source_document": document,
        "source_url": url,
        "local_path": str(local.relative_to(ROOT)),
        "sha256": sha256(local),
        "locator": locator,
        "evidence": evidence,
        "proves": proves,
        "does_not_prove": does_not_prove,
        "status": status,
        "confidence": "1.00",
        "review_required": review_required,
        "formal_ifc_write_allowed": formal_ifc_write_allowed,
        "manufacturer": manufacturer,
        "model_scope": model_scope,
        "revision": "",
        "publication_date": "",
        "legacy_targets": "elec-source-evidence.csv",
        "legacy_projection_json": "",
        "notes": notes,
    }
    row["legacy_projection_json"] = evidence_projection(row)
    return row


def requirement_value(row: dict[str, str], value: str) -> None:
    try:
        float(value)
    except ValueError:
        row["value_text"] = value
        row["value_number"] = ""
    else:
        row["value_text"] = ""
        row["value_number"] = value


def upsert_requirement(
    rows: list[dict[str, str]],
    equipment_id: str,
    parameter_key: str,
    value: str,
    *,
    discipline: str,
    unit: str = "",
    origin: str,
    status: str,
    source_id: str = "",
    locator: str = "",
    blocks: str = "no",
    notes: str = "",
    force_text: bool = False,
) -> None:
    matches = [row for row in rows if row["equipment_id"] == equipment_id and row["parameter_key"] == parameter_key]
    if len(matches) > 1:
        raise RuntimeError(f"duplicate requirement {equipment_id}:{parameter_key}")
    if matches:
        row = matches[0]
    else:
        numbers = [int(row["requirement_id"].split("-")[-1]) for row in rows]
        row = {"requirement_id": f"REQ-{max(numbers, default=0) + 1:04d}"}
        rows.append(row)
    row.update({
        "equipment_id": equipment_id,
        "discipline": discipline,
        "parameter_key": parameter_key,
        "unit": unit,
        "datum": "",
        "value_origin": origin,
        "status": status,
        "source_id": source_id,
        "source_locator": locator,
        "blocks_release": blocks,
        "notes": notes,
    })
    if force_text:
        row["value_text"] = value
        row["value_number"] = ""
    else:
        requirement_value(row, value)


def main() -> None:
    equipment_fields, equipment = read_csv(DECISIONS / "equipment-register.csv")
    requirement_fields, requirements = read_csv(DECISIONS / "equipment-installation-requirements.csv")
    source_fields, sources = read_csv(DECISIONS / "source-evidence-register.csv")
    owner_fields, owner_inputs = read_csv(DECISIONS / "owner-input-register.csv")
    closeout_fields, closeout = read_csv(DECISIONS / "owner-input-closeout-rules.csv")

    new_sources = [
        source_row(
            "APP-017-SIEMENS-WASHER-WEB-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-503/S-701",
            scope="Siemens WG54M7D20W 产品身份与尺寸", kind="official_product_web_snapshot",
            document="APP-017-Siemens-WG54M7D20W-product-api.json",
            url="https://www.siemens-home.bsh-group.cn/productlist/product-detail.html?vib=WG54M7D20W",
            locator="official API product record and line drawing", evidence="iQ500 10 kg；598×848×600 mm；门关闭深 639 mm；90° 开门总深 1111 mm",
            proves="洗衣机准确型号、厂家产品尺寸、关闭与开门外廓", does_not_prove="叠放连接件兼容清单、接口中心坐标或现场安装完成状态",
            status="verified_official_exact_model", manufacturer="Siemens", model_scope="WG54M7D20W",
        ),
        source_row(
            "APP-017-SIEMENS-WASHER-MANUAL-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-503/S-701",
            scope="Siemens WG54M7D20W 使用及安装条件", kind="official_product_pdf",
            document="APP-017-Siemens-WG54M7D20W-manual-9002062046_A.pdf",
            url="https://cms-pos.bshg.com.cn/product_media_resource/Documents/9002062046_A.pdf",
            locator="pp.3, 7-8, 13, 24-25", evidence="220 V/50 Hz、1900 W、最小 10 A；3/4 in 进水；排水最高 1000 mm；羊毛程序 2 kg；同厂同宽深并用原厂连接件方可叠放",
            proves="洗衣机水电、排水、调平、叠放规则与羊毛程序限制", does_not_prove="开门铰链侧、维修前抽距离或接口中心坐标",
            status="verified_official_exact_model", manufacturer="Siemens", model_scope="WG54M7D20W",
        ),
        source_row(
            "APP-017-SIEMENS-DRYER-WEB-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-503/S-701",
            scope="Siemens WQ55M7U20W 产品身份与尺寸", kind="official_product_web_snapshot",
            document="APP-017-Siemens-WQ55M7U20W-product-api.json",
            url="https://www.siemens-home.bsh-group.cn/productlist/product-detail.html?vib=WQ55M7U20W",
            locator="official API product record and line drawing", evidence="iQ500 10 kg 热泵干衣机；598×842×600 mm；门关闭深 639 mm；90° 开门总深 1111 mm",
            proves="干衣机准确型号、厂家产品尺寸、关闭与开门外廓", does_not_prove="逐型号叠放兼容表、接口中心坐标或现场安装完成状态",
            status="verified_official_exact_model", manufacturer="Siemens", model_scope="WQ55M7U20W",
        ),
        source_row(
            "APP-017-SIEMENS-DRYER-MANUAL-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-503/S-701",
            scope="Siemens WQ55M7U20W 使用及安装条件", kind="official_product_pdf",
            document="APP-017-Siemens-WQ55M7U20W-manual-9001938966_A.pdf",
            url="https://cms-pos.bshg.com.cn/product_media_resource/Documents/9001938966_A.pdf",
            locator="pp.10-14, 20, 22, 27-28, 40", evidence="220 V/50 Hz、800 W、最小 10 A；开放柜叠放最小 630×640×1730 mm；连接件订货号 17008829；冷凝直排/Y 分配器；羊毛挂架篮随机附带且羊毛单件静置烘干",
            proves="干衣机电源、通风、调平、冷凝排水、原厂连接件订货号、开放柜最小尺寸与羊毛程序限制", does_not_prove="17008829 是否带抽板、与 WTZ27510 的商业型号映射、接口中心坐标或维修前抽毫米数",
            status="verified_official_exact_model", manufacturer="Siemens", model_scope="WQ55M7U20W",
        ),
        source_row(
            "APP-017-SIEMENS-WTZ27510-WEB-001", discipline="INT1", sheet_id="I-503/S-701",
            scope="Siemens WTZ27510 带抽板连接组件候选", kind="official_product_web_snapshot",
            document="APP-017-Siemens-WTZ27510-product-api.json",
            url="https://www.siemens-home.bsh-group.cn/productlist/product-detail.html?vib=WTZ27510",
            locator="official API product record", evidence="官网商品名为连接组件，带抽板；593×563×37 mm",
            proves="WTZ27510 商品身份、带抽板描述与附件外廓", does_not_prove="WTZ27510 等同于说明书订货号 17008829，或对本两台 E-Nr. 的逐型号兼容",
            status="verified_official_candidate_not_compatibility_approved", manufacturer="Siemens", model_scope="WTZ27510",
            review_required="yes", formal_ifc_write_allowed="no", notes="仅作采购候选；购买前由西门子书面确认 17008829/WTZ27510 映射及两台 E-Nr./FD 兼容。",
        ),
        source_row(
            "APP-017-SIEMENS-STACKKIT-MANUAL-001", discipline="INT1", sheet_id="I-503/S-701",
            scope="Siemens WTZ27510 系列连接件安装页", kind="official_product_pdf",
            document="APP-017-Siemens-stacking-kit-manual-9001824701_B.pdf",
            url="https://cms-pos.bshg.com.cn/product_media_resource/Documents/9001824701_B.pdf",
            locator="pp.1-2", evidence="安装页列出 WTZ27510 系列并给出固定与叠放步骤",
            proves="WTZ27510 属于西门子原厂叠放连接件系列", does_not_prove="17008829 与 WTZ27510 的订货映射、抽板配置或本两台设备逐型号兼容",
            status="verified_official_candidate_not_compatibility_approved", manufacturer="Siemens", model_scope="WTZ27510",
            review_required="yes", formal_ifc_write_allowed="no",
        ),
        source_row(
            "APP-014-SIEMENS-NEW-WEB-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-501/S-701",
            scope="Siemens WS7FSB0C1C 新款净饮机身份", kind="official_product_web_snapshot",
            document="APP-014-Siemens-WS7FSB0C1C-product-api.json",
            url="https://www.siemens-home.bsh-group.cn/productlist/product-detail.html?vib=WS7FSB0C1C",
            locator="official API product record", evidence="iQ700 水玲珑 ProS 杯满即停新款；594×455×550 mm；220 V；最大电流 10 A",
            proves="APP-014 新款优选型号身份、外廓与杯满即停功能", does_not_prove="已采购、现场安装完成或接口中心坐标",
            status="verified_official_exact_model", manufacturer="Siemens", model_scope="WS7FSB0C1C",
        ),
        source_row(
            "APP-014-SIEMENS-OLD-WEB-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-501/S-701",
            scope="Siemens WS7060BC1C 旧款比价备选身份", kind="official_product_web_snapshot",
            document="APP-014-Siemens-WS7060BC1C-product-api.json",
            url="https://www.siemens-home.bsh-group.cn/productlist/product-detail.html?vib=WS7060BC1C",
            locator="official API product record", evidence="iQ700 水玲珑 Pro 旧款；594×455×550 mm；220 V；最大电流 10 A",
            proves="旧款比价备选的准确身份与同平台外廓", does_not_prove="当前优选、已采购或具备新款杯满即停功能",
            status="verified_official_price_comparison_alternative", manufacturer="Siemens", model_scope="WS7060BC1C",
            review_required="yes", formal_ifc_write_allowed="no",
        ),
        source_row(
            "APP-014-SIEMENS-MANUAL-001", discipline="ELEC/PLUM/INT1", sheet_id="E-303/P-201/I-501/S-701",
            scope="Siemens WS7FSB0C1C/WS7060BC1C 共用安装说明书", kind="official_product_pdf",
            document="APP-014-Siemens-WS7FSB0C1C-WS7060BC1C-manual-8001325304_D.pdf",
            url="https://cms-pos.bshg.com.cn/product_media_resource/Documents/8001325304_D.pdf",
            locator="pp.3, 5-6, 18-20", evidence="594×455×550 mm；220 V/50 Hz、2200 W；0.1-0.4 MPa、5-38 °C；上部通风口至少 100 cm²；侧板孔不小于 80 mm；水源与排水 2-3 m 范围；插座位于左右或上方邻柜且不得在设备下方",
            proves="新旧款共用的准确水电、通风、柜孔与检修安装条件", does_not_prove="已采购、具体接口中心坐标或柜体现场完成状态",
            status="verified_official_exact_model_family_manual", manufacturer="Siemens", model_scope="WS7FSB0C1C;WS7060BC1C",
        ),
    ]
    for source in new_sources:
        upsert_row(sources, "source_id", source["source_id"], source)

    # Preserve the superseded LG evidence as explicit decision history only.
    replace_row(sources, "source_id", "APP-017-LG-OFFICIAL-001", {
        "status": "superseded_candidate_reference",
        "formal_ifc_write_allowed": "no",
        "does_not_prove": "当前选型、当前协调净空或任何西门子设备安装条件",
        "notes": "2026-08-14 由业主改选西门子独立洗衣机＋热泵干衣机上下叠放；FN23BQH 仅保留历史追溯。",
    })
    lg_source = next(row for row in sources if row["source_id"] == "APP-017-LG-OFFICIAL-001")
    lg_source["legacy_projection_json"] = evidence_projection(lg_source)

    stack_source_ids = merge_ids(
        "IFC-FORMAL-001", "OWNER-INBOX-20260814-001",
        "APP-017-SIEMENS-WASHER-WEB-001", "APP-017-SIEMENS-WASHER-MANUAL-001",
        "APP-017-SIEMENS-DRYER-WEB-001", "APP-017-SIEMENS-DRYER-MANUAL-001",
        "APP-017-SIEMENS-WTZ27510-WEB-001", "APP-017-SIEMENS-STACKKIT-MANUAL-001",
    )
    replace_row(equipment, "equipment_id", "APP-017", {
        "item_name": "西门子洗衣机＋热泵干衣机叠放套装",
        "manufacturer": "Siemens",
        "model": "WG54M7D20W + WQ55M7U20W",
        "variant": "iQ500 10 kg washer + iQ500 10 kg heat-pump dryer",
        "procurement_status": "selected",
        "decision_status": "partial",
        "storage_location_candidate": "中厨家政位",
        "use_location_candidate": "中厨家政位",
        "storage_location_confirmed": "",
        "use_location_confirmed": "中厨家政位",
        "source_ids": stack_source_ids,
        "identity_basis": "业主明确改选两台独立西门子设备上下叠放；官方精确型号资料核实尺寸、负载、水电和叠放规则",
        "confidence": "1.00",
        "human_review_required": "yes",
        "notes": "不是洗干一体机或 WashTower。有效协调净空 650W×800D×1900H mm；两个可检修插座；给排水优先侧置；不得六面密封；原厂连接件必须使用。说明书订货号 17008829 已确认，但其抽板状态及与 WTZ27510 的商业型号映射仍待西门子书面核对。",
    })

    component_defaults = {
        "domain": "APPLIANCE", "quantity": "1", "procurement_status": "selected", "decision_status": "partial",
        "storage_location_candidate": "中厨家政位", "use_location_candidate": "中厨家政位",
        "storage_location_confirmed": "", "use_location_confirmed": "中厨家政位", "schedule_included": "no",
        "selector_kind": "global_id", "ifc_type_name": "", "ifc_type_global_id": "", "confidence": "1.00",
        "human_review_required": "yes", "legacy_kind": "", "legacy_id": "",
    }
    upsert_row(equipment, "equipment_id", "APP-017-WASHER", {
        **component_defaults, "category": "洗衣机", "item_name": "西门子独立滚筒洗衣机", "manufacturer": "Siemens",
        "model": "WG54M7D20W", "variant": "iQ500 10 kg", "selector_value": "1Uzcf8kU9J2w9K40VkzCgV",
        "ifc_class": "IfcElectricAppliance", "ifc_global_ids": "1Uzcf8kU9J2w9K40VkzCgV",
        "source_ids": "APP-017-SIEMENS-WASHER-WEB-001;APP-017-SIEMENS-WASHER-MANUAL-001",
        "identity_basis": "业主主选型号；西门子官网及精确型号说明书核验",
        "notes": "APP-017 叠放套装的下部独立设备；产品尺寸不是接口中心。",
    })
    upsert_row(equipment, "equipment_id", "APP-017-DRYER", {
        **component_defaults, "category": "热泵干衣机", "item_name": "西门子独立热泵干衣机", "manufacturer": "Siemens",
        "model": "WQ55M7U20W", "variant": "iQ500 10 kg", "selector_value": "3_ZSeTAOPNThu9wyl2sOxQ",
        "ifc_class": "IfcElectricAppliance", "ifc_global_ids": "3_ZSeTAOPNThu9wyl2sOxQ",
        "source_ids": "APP-017-SIEMENS-DRYER-WEB-001;APP-017-SIEMENS-DRYER-MANUAL-001",
        "identity_basis": "业主主选型号；西门子官网及精确型号说明书核验",
        "notes": "APP-017 叠放套装的上部独立设备；产品尺寸不是接口中心。",
    })

    base_requirements = [
        ("rated_power", "2700", "ELEC", "W", "project_candidate", "confirmed", "APP-017-SIEMENS-WASHER-MANUAL-001", "两台额定功率算术合计；不自动决定共用回路。"),
        ("independent_power_plug_count", "2", "ELEC", "unit", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-MANUAL-001", "两台各一根电源线和插头；禁止排插/延长线。"),
        ("shared_branch_circuit_permission", "unknown", "ELEC", "", "pending", "pending", "", "厂家资料未明确允许共回路；由配电计算和当地规范关闭。"),
        ("water_required", "yes", "PLUM", "", "official_exact_model", "confirmed", "APP-017-SIEMENS-WASHER-MANUAL-001", "仅洗衣机需要生活给水。"),
        ("drain_required", "yes", "PLUM", "", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-MANUAL-001", "洗衣排水；干衣机可冷凝水盒或直接排水。"),
        ("ventilation_required", "yes", "HVAC", "", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-MANUAL-001", "开放柜、通风口畅通；不得六面密封。"),
        ("equipment_form", "two_independent_appliances_stacked", "MULTI", "", "user_input", "confirmed", "", "不得标作洗干一体机或 WashTower。"),
        ("stack_connector_order_number", "17008829", "INT1", "", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-MANUAL-001", "说明书订货号；是否带抽板 unknown。"),
        ("stack_connector_drawer_status", "unknown", "INT1", "", "pending", "pending", "", "说明书没有抽板描述。"),
        ("stack_connector_commercial_candidate", "WTZ27510", "INT1", "", "official_exact_model", "candidate", "APP-017-SIEMENS-WTZ27510-WEB-001", "官网称带抽板；尚未证明等同 17008829 或逐型号兼容。"),
        ("stack_body_height_sum", "1690", "INT1", "mm", "project_candidate", "confirmed", "APP-017-SIEMENS-DRYER-WEB-001", "两台厂家机身高度算术和，不得单独用于柜体放样。"),
        ("manufacturer_open_cabinet_width_min", "630", "INT1", "mm", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-MANUAL-001", "干衣机说明书叠放开放柜最小值。"),
        ("manufacturer_open_cabinet_depth_min", "640", "INT1", "mm", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-MANUAL-001", "干衣机说明书叠放开放柜最小值。"),
        ("manufacturer_open_cabinet_height_min", "1730", "INT1", "mm", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-MANUAL-001", "已经高于 1690 mm 纯机身和。"),
        ("coordination_clear_width", "650", "INT1", "mm", "user_input", "confirmed", "", "项目有效净空，不是厂家产品尺寸。"),
        ("coordination_clear_depth", "800", "INT1", "mm", "user_input", "confirmed", "", "项目有效净空，不是厂家产品尺寸。"),
        ("coordination_clear_height", "1900", "INT1", "mm", "user_input", "confirmed", "", "项目有效净空，含连接件、调平及安装维护余量。"),
        ("fixed_front_threshold_allowed", "no", "INT1", "", "user_input", "confirmed", "", "整套设备必须可向前抽出维修。"),
        ("removable_side_top_panels_required", "yes", "INT1", "", "user_input", "confirmed", "", "两侧与顶部可拆卸收口。"),
        ("service_connections_preferred_zone", "side_accessible_not_directly_behind", "ELEC/PLUM/INT1", "", "user_input", "confirmed", "", "避免软管和插头在机器正后方受挤压。"),
        ("service_interface_center_coordinates", "unknown", "ELEC/PLUM", "", "pending", "pending", "", "待厂家安装图/现场深化，不得从 IFC proxy 包围盒推测。"),
        ("front_open_total_depth", "1111", "INT1", "mm", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-WEB-001", "两台 90° 开门总深均为 1111 mm；前方取衣空间按此复核。"),
        ("dryer_condensate_discharge", "tank_or_direct_drain", "PLUM", "", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-MANUAL-001", "直排附件随机附带。"),
        ("shared_siphon_y_distributor", "15000490", "PLUM", "", "official_exact_model", "candidate", "APP-017-SIEMENS-DRYER-MANUAL-001", "官方 Y 分配器可使洗衣机与干衣机共用虹吸口；现场深化确认。"),
        ("wool_selection_reason", "washer wool 2kg; dryer wool rack one garment stationary", "MULTI", "", "official_exact_model", "confirmed", "APP-017-SIEMENS-DRYER-MANUAL-001", "产品选择理由，不是建筑安装条件。"),
    ]
    for key, value, discipline, unit, origin, status, source_id, notes in base_requirements:
        upsert_requirement(requirements, "APP-017", key, value, discipline=discipline, unit=unit, origin=origin, status=status, source_id=source_id, blocks="yes" if status in {"pending", "candidate"} and key in {"shared_branch_circuit_permission", "stack_connector_drawer_status", "stack_connector_commercial_candidate", "service_interface_center_coordinates"} else "no", notes=notes, force_text=key == "stack_connector_order_number")

    # APP-017 is now the business-level stacked set. Remove obsolete single-product
    # LG dimensions/capacities so they cannot be mistaken for current aggregate data.
    obsolete_parent_keys = {"product_width", "product_depth", "product_height", "wash_capacity", "dry_capacity"}
    requirements[:] = [
        row for row in requirements
        if not (row["equipment_id"] == "APP-017" and row["parameter_key"] in obsolete_parent_keys)
    ]

    component_requirements = {
        "APP-017-WASHER": [
            ("product_width", "598", "mm", "APP-017-SIEMENS-WASHER-WEB-001"), ("product_height", "848", "mm", "APP-017-SIEMENS-WASHER-WEB-001"),
            ("product_depth", "600", "mm", "APP-017-SIEMENS-WASHER-WEB-001"), ("door_closed_depth", "639", "mm", "APP-017-SIEMENS-WASHER-WEB-001"),
            ("door_open_90_depth", "1111", "mm", "APP-017-SIEMENS-WASHER-WEB-001"), ("leveling_height_add_max", "12", "mm", "APP-017-SIEMENS-WASHER-WEB-001"),
            ("rated_power", "1900", "W", "APP-017-SIEMENS-WASHER-MANUAL-001"), ("rated_voltage", "220", "V", "APP-017-SIEMENS-WASHER-MANUAL-001"),
            ("frequency", "50", "Hz", "APP-017-SIEMENS-WASHER-MANUAL-001"), ("minimum_fuse", "10", "A", "APP-017-SIEMENS-WASHER-MANUAL-001"),
            ("power_cord_length", "2100", "mm", "APP-017-SIEMENS-WASHER-MANUAL-001"), ("water_connection", "G3/4 cold water", "", "APP-017-SIEMENS-WASHER-MANUAL-001"),
            ("water_pressure_min", "100", "kPa", "APP-017-SIEMENS-WASHER-MANUAL-001"), ("water_pressure_max", "1000", "kPa", "APP-017-SIEMENS-WASHER-MANUAL-001"),
            ("inlet_hose_length", "1200", "mm", "APP-017-SIEMENS-WASHER-MANUAL-001"), ("drain_hose_length", "1500", "mm", "APP-017-SIEMENS-WASHER-MANUAL-001"),
            ("drain_height_max", "1000", "mm", "APP-017-SIEMENS-WASHER-MANUAL-001"), ("wool_program_load_max", "2", "kg", "APP-017-SIEMENS-WASHER-MANUAL-001"),
            ("door_handedness", "unknown", "", "APP-017-SIEMENS-WASHER-MANUAL-001"),
        ],
        "APP-017-DRYER": [
            ("product_width", "598", "mm", "APP-017-SIEMENS-DRYER-WEB-001"), ("product_height", "842", "mm", "APP-017-SIEMENS-DRYER-WEB-001"),
            ("product_depth", "600", "mm", "APP-017-SIEMENS-DRYER-WEB-001"), ("door_closed_depth", "639", "mm", "APP-017-SIEMENS-DRYER-WEB-001"),
            ("door_open_90_depth", "1111", "mm", "APP-017-SIEMENS-DRYER-WEB-001"), ("leveling_height_add_max", "15", "mm", "APP-017-SIEMENS-DRYER-WEB-001"),
            ("rated_power", "800", "W", "APP-017-SIEMENS-DRYER-MANUAL-001"), ("rated_voltage", "220", "V", "APP-017-SIEMENS-DRYER-MANUAL-001"),
            ("frequency", "50", "Hz", "APP-017-SIEMENS-DRYER-MANUAL-001"), ("minimum_fuse", "10", "A", "APP-017-SIEMENS-DRYER-MANUAL-001"),
            ("power_cord_length", "1450", "mm", "APP-017-SIEMENS-DRYER-MANUAL-001"), ("wool_rack_included", "yes", "", "APP-017-SIEMENS-DRYER-MANUAL-001"),
            ("wool_rack_program_load", "one_garment", "", "APP-017-SIEMENS-DRYER-MANUAL-001"), ("door_handedness", "customer_service_reversible", "", "APP-017-SIEMENS-DRYER-MANUAL-001"),
        ],
    }
    for equipment_id, values in component_requirements.items():
        for key, value, unit, source_id in values:
            pending = value == "unknown"
            upsert_requirement(requirements, equipment_id, key, value, discipline="ELEC/PLUM/INT1", unit=unit, origin="pending" if pending else "official_exact_model", status="pending" if pending else "confirmed", source_id=source_id, blocks="yes" if pending else "no", notes="厂家参数；不作为接口中心坐标。")

    water_sources = "OWNER-INBOX-20260814-001;APP-014-SIEMENS-NEW-WEB-001;APP-014-SIEMENS-OLD-WEB-001;APP-014-SIEMENS-MANUAL-001"
    for equipment_id, alias in (("APP-014", False), ("APP-015", True)):
        replace_row(equipment, "equipment_id", equipment_id, {
            "item_name": "净水功能别名（同 APP-014）" if alias else "西门子嵌入式净水饮水机（直饮+净水）",
            "manufacturer": "Siemens", "model": "WS7FSB0C1C",
            "quantity": "0" if alias else "1",
            "procurement_status": "not_applicable" if alias else "candidate",
            "decision_status": "superseded" if alias else "partial",
            "schedule_included": "no" if alias else "yes",
            "human_review_required": "no" if alias else "yes",
            "source_ids": water_sources,
            "identity_basis": "业主明确新款 WS7FSB0C1C 为优选、旧款 WS7060BC1C 仅作比价备选；两者官方资料核验",
            "notes": "Alias of APP-014；不得重复统计设备或接口。" if alias else "APP-015 为同一实物功能别名。新款为当前优选、旧款仅比价备选；均未因本次研究标记为已采购。共用官方安装说明书。",
        })
    for key, value, discipline, unit, origin, status, source_id, notes in [
        ("exact_enr", "WS7FSB0C1C", "MULTI", "", "official_exact_model", "confirmed", "APP-014-SIEMENS-NEW-WEB-001", "当前业主优选。"),
        ("price_comparison_alternative", "WS7060BC1C", "MULTI", "", "user_input", "confirmed", "APP-014-SIEMENS-OLD-WEB-001", "旧款只作比价备选，不是当前优选或已采购。"),
        ("installation_manual", "8001325304_D", "MULTI", "", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", "新旧款共用官方说明书。"),
        ("rated_power", "2200", "ELEC", "W", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("rated_voltage", "220", "ELEC", "V", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("frequency", "50", "ELEC", "Hz", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("rated_current_max", "10", "ELEC", "A", "official_exact_model", "confirmed", "APP-014-SIEMENS-NEW-WEB-001", ""),
        ("product_width", "594", "INT1", "mm", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("product_height", "455", "INT1", "mm", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("product_depth", "550", "INT1", "mm", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("water_pressure_min", "0.1", "PLUM", "MPa", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("water_pressure_max", "0.4", "PLUM", "MPa", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("inlet_water_temperature_min", "5", "PLUM", "°C", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("inlet_water_temperature_max", "38", "PLUM", "°C", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("ventilation_opening_area_min", "100", "INT1", "cm2", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", "橱柜上部通风口，不得遮盖。"),
        ("cabinet_side_opening_diameter_min", "80", "INT1", "mm", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("faucet_hole_diameter", "30", "INT1", "mm", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("service_distance_max", "3000", "ELEC/PLUM", "mm", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", "随机自带 3 m 管；改装须咨询安装工程师。"),
        ("power_cord_length", "2000", "ELEC", "mm", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("power_outlet_location", "left_right_or_upper_adjacent_cabinet_not_below_unit", "ELEC/INT1", "", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", "插头必须可自由插拔。"),
        ("water_required", "yes", "PLUM", "", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("drain_required", "yes", "PLUM", "", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", ""),
        ("ventilation_required", "yes", "HVAC/INT1", "", "official_exact_model", "confirmed", "APP-014-SIEMENS-MANUAL-001", "嵌柜须保持通风。"),
    ]:
        upsert_requirement(requirements, "APP-014", key, value, discipline=discipline, unit=unit, origin=origin, status=status, source_id=source_id, blocks="no", notes=notes)

    replace_row(owner_inputs, "input_id", "APP014-APP015-ALIAS", {
        "question": "APP-014 新款净饮机优选、旧款比价备选与 APP-015 别名去重",
        "candidate_value": "WS7FSB0C1C 为优选；WS7060BC1C 只作比价备选；APP-015 数量为 0",
        "user_value": "APP-014 优选 Siemens WS7FSB0C1C；旧款 WS7060BC1C 保留为比价备选；APP-015 与 APP-014 为同一实物",
        "status": "自定义确认",
        "evidence_reference": water_sources,
        "source_basis": "业主当前明确选型；西门子官网和共用官方说明书",
        "notes": "型号决策与水电柜孔条件已同步；现场接口中心与最终柜体安装图仍由深化关闭。",
    })
    upsert_row(owner_inputs, "input_id", "APP017-SIEMENS-STACK", {
        "workstream": "E-303/P-201/I-503/S-701", "priority": "P0", "blocks_release": "yes",
        "question": "APP-017 西门子独立洗衣机＋热泵干衣机叠放套装完整协调",
        "candidate_value": "WG54M7D20W + WQ55M7U20W；有效净空 650×800×1900 mm；原厂连接件",
        "user_value": "主选 WG54M7D20W + WQ55M7U20W；按 650W×800D×1900H mm 有效净空协调；侧置可检修水电；无固定门槛",
        "unit": "mm", "status": "自定义确认", "evidence_reference": stack_source_ids,
        "source_basis": "业主明确主选与协调净空；西门子精确型号官网及官方说明书",
        "sync_target": "equipment SSOT + IFC + E-303/P-201/I-503/S-701",
        "notes": "说明书订货号 17008829 已确认；抽板状态、与 WTZ27510 映射、两台 E-Nr./FD 逐型号兼容及接口中心仍待厂家/深化确认。",
    })
    replace_row(closeout, "input_id", "APP014-APP015-ALIAS", {
        "closeout_kind": "final_shop_drawing_and_site_interface",
        "responsible_party": "设备方/橱柜方/给排水设计/电气设计",
        "required_evidence": "WS7FSB0C1C 最终柜体安装图、现场给排水/电源位置及检修复核",
        "automatic_close_allowed": "no",
        "notes": "官方共用说明书已取得；APP-015 去重已关闭；剩余为项目落位与现场接口。",
    })
    upsert_row(closeout, "input_id", "APP017-SIEMENS-STACK", {
        "closeout_kind": "manufacturer_compatibility_and_rough_in",
        "responsible_party": "西门子客服/设备方/橱柜方/给排水设计/电气设计",
        "required_evidence": "书面确认 WG54M7D20W 与 WQ55M7U20W 对应 E-Nr./FD 可叠放、17008829 的商业型号与抽板配置；最终侧置插座/水龙头/排水和维修抽出尺寸",
        "automatic_close_allowed": "no",
        "notes": "不得以 WTZ27510 官网候选反向宣称 17008829 带抽板或逐型号已兼容。",
    })

    write_csv(DECISIONS / "equipment-register.csv", equipment_fields, equipment)
    write_csv(DECISIONS / "equipment-installation-requirements.csv", requirement_fields, requirements)
    write_csv(DECISIONS / "source-evidence-register.csv", source_fields, sources)
    write_csv(DECISIONS / "owner-input-register.csv", owner_fields, owner_inputs)
    write_csv(DECISIONS / "owner-input-closeout-rules.csv", closeout_fields, closeout)
    result = {"canonical": validate(ROOT), "projections": projections(ROOT)}
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
