#!/usr/bin/env python3
"""Upsert the accepted first electrical decision batch into project SSOT tables."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DECISIONS = ROOT / "pipeline/decisions"
RESEARCH_SOURCE_ID = "ELEC-RESEARCH-20260815-001"
RESEARCH_PATH = ROOT / "drawings/evidence/ELEC-第一组设计结论-20260815.md"


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


def main() -> None:
    source_fields, sources = read_csv(DECISIONS / "source-evidence-register.csv")
    projection = {
        "evidence_id": RESEARCH_SOURCE_ID,
        "discipline": "ELEC/INT1",
        "sheet_id": "E-302/E-303/E-304/I-501/I-503/S-701",
        "decision_scope": "第一组电气事实、工程结论与未知边界",
        "source_kind": "project_research_conclusion",
        "source_document": "drawings/evidence/ELEC-第一组设计结论-20260815.md",
        "source_sha256": sha256(RESEARCH_PATH),
        "source_locator": "语义边界；官方产品事实；工程研究结论；未知与外部证据门",
        "evidence": "将产品官方值与插座、支路保护、独立回路的研究结论分开；记录 NS-01/NS-02、三块控制面板和 PoE 拓扑",
        "proves": "本项目已采用的电气设计结论及其计算、状态和证据边界",
        "does_not_prove": "候选产品已购、终选 SKU、说明书未给的接口中心、现场网线通断或 PoE 最终功率预算",
        "status": "verified_research_conclusion",
        "confidence": "1.00",
        "review_required": "no",
        "formal_ifc_write_allowed": "no",
        "notes": "引用现有官方证据 ID；不将研究结论冒充厂家原文。",
    }
    source_row = {
        "source_id": RESEARCH_SOURCE_ID,
        "discipline": projection["discipline"],
        "sheet_id": projection["sheet_id"],
        "decision_scope": projection["decision_scope"],
        "source_kind": projection["source_kind"],
        "source_document": projection["source_document"],
        "source_url": "",
        "local_path": projection["source_document"],
        "sha256": projection["source_sha256"],
        "locator": projection["source_locator"],
        "evidence": projection["evidence"],
        "proves": projection["proves"],
        "does_not_prove": projection["does_not_prove"],
        "status": projection["status"],
        "confidence": projection["confidence"],
        "review_required": projection["review_required"],
        "formal_ifc_write_allowed": projection["formal_ifc_write_allowed"],
        "manufacturer": "",
        "model_scope": "",
        "revision": "2026-08-15",
        "publication_date": "2026-08-15",
        "legacy_targets": "elec-source-evidence.csv",
        "legacy_projection_json": json.dumps(projection, ensure_ascii=False, separators=(",", ":")),
        "notes": projection["notes"],
    }
    gb_row = {
        "source_id": "GB1002-2024-OFFICIAL-001",
        "discipline": "ELEC",
        "sheet_id": "E-303",
        "decision_scope": "家用单相插头插座型式、基本参数和尺寸",
        "source_kind": "official_national_standard_record",
        "source_document": "GB 1002-2024",
        "source_url": "https://openstd.samr.gov.cn/bzgk/std/newGbInfo?hcno=F8C9E208891B7BB5AF1B3E64933693C2",
        "local_path": "",
        "sha256": "not_applicable_live_official_standard",
        "locator": "国家标准公开系统；现行；2025-08-01 实施",
        "evidence": "GB 1002-2024 为现行强制性国家标准记录",
        "proves": "最终插头插座型式应按现行 GB 1002-2024 执行",
        "does_not_prove": "任一候选设备的实际插头额定值或已购状态",
        "status": "verified_official_source",
        "confidence": "1.00",
        "review_required": "no",
        "formal_ifc_write_allowed": "no",
        "manufacturer": "国家市场监督管理总局/国家标准化管理委员会",
        "model_scope": "GB 1002-2024",
        "revision": "2024",
        "publication_date": "2024-07-24",
        "legacy_targets": "",
        "legacy_projection_json": "",
        "notes": "网页为官方动态记录，因此 SHA-256 使用明确的 not_applicable 标记。",
    }
    source_by_id = {row["source_id"]: row for row in sources}
    source_by_id[RESEARCH_SOURCE_ID] = source_row
    source_by_id[gb_row["source_id"]] = gb_row
    write_csv(DECISIONS / "source-evidence-register.csv", source_fields, list(source_by_id.values()))

    master_fields, masters = read_csv(DECISIONS / "equipment-register.csv")
    master_notes = {
        "APP-001": "准确型号和名牌功率保持 unknown；NS-01 已按两个同时 10A 移动负载端口的不利工况设计两路 C16，未把回路容量冒充为产品功率。",
        "APP-002": "准确型号和名牌功率保持 unknown；与火锅同时使用时必须分路。",
        "APP-003": "准确型号和名牌功率保持 unknown；与火锅同时使用时必须分路。",
        "APP-004": "GS3 与 E1 Prima EXP 均保持候选，未采购、未终选；粗装按 2600W 不利值预留独立 C16 和接线盒，最终面板按实购插头匹配。",
        "APP-009": "候选洗碗机未采购；按官方 2000W/10A 和两台同时工况，每台各配 10A 接地插座、独立 C16 RCBO 和 2.5mm² 铜线。",
        "APP-010": "候选洗碗机未采购；按官方 2000W/10A 和两台同时工况，每台各配 10A 接地插座、独立 C16 RCBO 和 2.5mm² 铜线。",
        "APP-011": "已到货；官方 3400W/16A 插头。16A 三孔插座、独立 C16 RCBO、4mm² 铜线，插座设相邻可检修柜格而非机身正后。",
        "APP-012": "已到货；DC1.5V 电池点火，设备不需电源插座；原灶位插座仅作 10A 检修备用点，接厨房固定辅助回路。",
        "APP-013": "已到货；官方 385W，设计电流约 1.75A；10A 接地插座接厨房固定辅助回路，不独立成回路。",
        "APP-014": "新款为优选候选但未采购；2200W，10A 三孔邻柜插座，独立 C16 RCBO、2.5mm² 铜线和剩余电流保护；APP-015 仅为别名。",
        "APP-015": "APP-014 的功能别名，quantity=0、不入设备表、不重复计算功率、插座、回路或给排水。",
        "APP-017": "候选精确双机叠放；1900W+800W=2700W，约 12.3A@220V；一路洗衣区独立 C16 RCBO、2.5mm² 铜线和两个侧面可检修 10A 接地插座，不使用排插/延长线。",
    }
    for row in masters:
        if row["equipment_id"] in master_notes:
            row["source_ids"] = joined(row["source_ids"], RESEARCH_SOURCE_ID)
            row["notes"] = master_notes[row["equipment_id"]]
        if row["equipment_id"] in {"CTRL-ENTRY-A", "CTRL-MASTER-A", "CTRL-MASTER-B", "NET-AP-R09", "NET-AP-R14"}:
            row["source_ids"] = joined(row["source_ids"], RESEARCH_SOURCE_ID)
    write_csv(DECISIONS / "equipment-register.csv", master_fields, masters)

    req_fields, requirements = read_csv(DECISIONS / "equipment-installation-requirements.csv")
    obsolete_generic_ranges = {
        (equipment_id, parameter_key)
        for equipment_id in ("APP-001", "APP-002", "APP-003", "APP-004")
        for parameter_key in ("candidate_power_min", "candidate_power_max")
    }
    requirements = [
        row
        for row in requirements
        if (row["equipment_id"], row["parameter_key"]) not in obsolete_generic_ranges
    ]
    by_pair = {(row["equipment_id"], row["parameter_key"]): row for row in requirements}

    def upsert(equipment_id: str, key: str, value: str, *, unit: str = "", origin: str = "research_conclusion", status: str = "confirmed", source: str = RESEARCH_SOURCE_ID, blocks: str = "no", locator: str = "工程研究结论", notes: str = "") -> None:
        pair = (equipment_id, key)
        row = by_pair.get(pair)
        if row is None:
            row = {field: "" for field in req_fields}
            row["requirement_id"] = f"REQ-ELEC-20260815-{len([r for r in requirements if r['requirement_id'].startswith('REQ-ELEC-20260815-')]) + 1:03d}"
            row["equipment_id"] = equipment_id
            requirements.append(row)
            by_pair[pair] = row
        numeric = False
        try:
            float(value)
            numeric = value not in {"yes", "no"}
        except ValueError:
            pass
        row.update({
            "discipline": "ELEC",
            "parameter_key": key,
            "value_text": "" if numeric else value,
            "value_number": value if numeric else "",
            "unit": unit,
            "datum": "",
            "value_origin": origin,
            "status": status,
            "source_id": source,
            "source_locator": locator,
            "blocks_release": blocks,
            "notes": notes,
        })

    for equipment_id in ("APP-001", "APP-002", "APP-003", "APP-004", "APP-006", "APP-007", "APP-008", "APP-018"):
        if (equipment_id, "rated_power") in by_pair:
            by_pair[(equipment_id, "rated_power")]["blocks_release"] = "no"
            by_pair[(equipment_id, "rated_power")]["notes"] = "产品额定功率保持 unknown；相关支路已用端口上限或候选最不利值闭合粗装，不将设计值写成产品名牌值。"
    for equipment_id in ("APP-005", "APP-006"):
        for key in ("candidate_power_min", "candidate_power_max"):
            if (equipment_id, key) in by_pair:
                by_pair[(equipment_id, key)]["blocks_release"] = "no"
                by_pair[(equipment_id, key)]["notes"] = "仅作小功率辅助回路情景校核；NS-02 已按 10A 端口与候选最不利工况闭合粗装，准确型号只在设备登记时复核。"

    official = {
        "APP-011": ("APP-011-OFFICIAL-001", "官方规格表"),
        "APP-014": ("APP-014-SIEMENS-MANUAL-001", "pp.18-20"),
        "APP-017-W": ("APP-017-SIEMENS-WASHER-MANUAL-001", "p.24"),
        "APP-017-D": ("APP-017-SIEMENS-DRYER-MANUAL-001", "p.40"),
        "APP-013": ("APP-013-OFFICIAL-001", "安装资料"),
        "APP-012": ("APP-012-OFFICIAL-001", "官方产品资料"),
        "APP-009": ("APP-DW-SPEC-001", "官方规格表"),
    }
    upsert("APP-011", "appliance_plug_rating_a", "16", unit="A", origin="official_exact_model", source=official["APP-011"][0], locator=official["APP-011"][1])
    for key, value, unit in [("wall_socket_rating_a", "16", "A"), ("branch_breaker_rating_a", "16", "A"), ("conductor_cross_section_mm2", "4", "mm2")]: upsert("APP-011", key, value, unit=unit)
    for key, value in [("branch_breaker_curve", "C"), ("rcbo_required", "yes"), ("dedicated_branch_circuit", "yes"), ("socket_service_location", "adjacent_accessible_cabinet_not_directly_behind"), ("track_socket_allowed", "no")]: upsert("APP-011", key, value)
    upsert("APP-011", "interface_center_coordinates", "unknown", origin="pending", status="pending", source="", locator="", blocks="no", notes="说明书未提供接口中心坐标，不推测；项目采用相邻可检修柜格插座，不以机背中心坐标作为粗装停止条件。")

    upsert("APP-017", "washer_rated_power_w", "1900", unit="W", origin="official_exact_model", source=official["APP-017-W"][0], locator=official["APP-017-W"][1])
    upsert("APP-017", "dryer_rated_power_w", "800", unit="W", origin="official_exact_model", source=official["APP-017-D"][0], locator=official["APP-017-D"][1])
    for key, value, unit in [("simultaneous_design_load_w", "2700", "W"), ("simultaneous_design_current_a", "12.3", "A"), ("wall_socket_rating_a", "10", "A"), ("wall_socket_quantity", "2", "count"), ("branch_breaker_rating_a", "16", "A"), ("conductor_cross_section_mm2", "2.5", "mm2")]: upsert("APP-017", key, value, unit=unit)
    for key, value in [("shared_branch_circuit_permission", "yes"), ("dedicated_branch_circuit", "laundry_area_shared_by_washer_and_dryer"), ("rcbo_required", "yes"), ("extension_or_power_strip_allowed", "no"), ("socket_service_location", "two_independent_accessible_side_cabinet_positions_not_behind_appliances")]: upsert("APP-017", key, value)
    upsert("APP-017", "appliance_plug_rating_a", "unknown", origin="pending", status="pending", source="", locator="", blocks="no", notes="说明书给出最小保险 10A，但未明示插头额定标识；墙面 10A 插座为研究结论。")
    upsert("APP-017", "interface_center_coordinates", "unknown", origin="pending", status="pending", source="", locator="", blocks="yes", notes="水、排水和电源接口中心坐标未由说明书提供。")

    upsert("APP-014", "wall_socket_rating_a", "10", unit="A", origin="official_exact_model", source=official["APP-014"][0], locator=official["APP-014"][1])
    for key, value, unit in [("branch_breaker_rating_a", "16", "A"), ("conductor_cross_section_mm2", "2.5", "mm2")]: upsert("APP-014", key, value, unit=unit)
    for key, value in [("dedicated_branch_circuit", "yes"), ("rcbo_required", "yes"), ("residual_current_protection_required", "yes")]: upsert("APP-014", key, value)
    upsert("APP-014", "interface_center_coordinates", "unknown", origin="pending", status="pending", source="", locator="", blocks="yes", notes="安装图仅给服务区和距离，没有可写入的接口中心坐标。")

    upsert("APP-013", "design_current_a", "1.75", unit="A")
    upsert("APP-013", "wall_socket_rating_a", "10", unit="A")
    upsert("APP-013", "dedicated_branch_circuit", "no")
    upsert("APP-013", "assigned_branch_circuit", "kitchen_fixed_auxiliary")
    upsert("APP-012", "ignition_power_source", "DC_1.5V_battery", origin="official_exact_model", source=official["APP-012"][0], locator=official["APP-012"][1])
    upsert("APP-012", "wall_socket_required_for_appliance", "no")
    upsert("APP-012", "service_spare_socket_rating_a", "10", unit="A")
    upsert("APP-012", "service_spare_socket_circuit", "kitchen_fixed_auxiliary")

    for equipment_id in ("APP-009", "APP-010"):
        upsert(equipment_id, "wall_socket_rating_a", "10", unit="A")
        upsert(equipment_id, "branch_breaker_rating_a", "16", unit="A")
        upsert(equipment_id, "conductor_cross_section_mm2", "2.5", unit="mm2")
        upsert(equipment_id, "dedicated_branch_circuit", "yes")
        upsert(equipment_id, "rcbo_required", "yes")
        upsert(equipment_id, "simultaneous_operation_with_peer", "yes")

    upsert("APP-004", "design_reserve_power_w", "2600", unit="W")
    upsert("APP-004", "design_reserve_current_a", "11.8", unit="A")
    upsert("APP-004", "branch_breaker_rating_a", "16", unit="A")
    upsert("APP-004", "conductor_cross_section_mm2", "2.5", unit="mm2")
    upsert("APP-004", "dedicated_branch_circuit", "yes")
    upsert("APP-004", "rcbo_required", "yes")
    upsert("APP-004", "connection_box_reserved", "yes")
    upsert("APP-004", "wall_socket_rating_a", "unknown", origin="pending", status="pending", source="", locator="", blocks="no", notes="最终按实购机器插头匹配 10A/16A 面板，当前不冻结。")
    upsert("APP-004", "procurement_selection_required_before_rough_in", "no")
    upsert("APP-004", "water_and_drain_provision", "optional_valved_and_capped")

    for equipment_id in ("APP-001", "APP-002", "APP-003"):
        upsert(equipment_id, "wall_socket_rating_a", "10", unit="A")
        upsert(equipment_id, "ns01_branch_allocation", "hotpot_separate_from_second_simultaneous_appliance")
    upsert("APP-001", "ns01_total_branch_circuits", "2", unit="count")
    upsert("APP-001", "branch_breaker_rating_a", "16", unit="A")
    upsert("APP-001", "conductor_cross_section_mm2", "2.5", unit="mm2")
    upsert("APP-001", "rcbo_required", "yes")
    upsert("APP-001", "adverse_design_current_a", "20", unit="A")

    panel_requirements = {
        "CTRL-ENTRY-A": [("logical_panel_quantity", "1"), ("minimum_key_count", "4"), ("key_order_door_to_far", "living_two_way;study_two_way;dining_two_way;lighting_master"), ("controlled_loads", "living;study;dining_lighting"), ("scene_control_relation", "lighting_master_command_only_no_continuous_load_cutoff")],
        "CTRL-MASTER-B": [("logical_panel_quantity", "1"), ("minimum_key_count", "3"), ("key_order_door_to_far", "living_two_way;study_two_way;dining_two_way"), ("controlled_loads", "living;study;dining_lighting"), ("two_way_relation", "true_wired_two_way_with_CTRL_ENTRY_A_not_scene_substitution")],
        "CTRL-MASTER-A": [("logical_panel_quantity", "1"), ("minimum_key_count", "2"), ("key_order_door_to_far", "master_ambient_LED;master_accent_spotlights"), ("controlled_loads", "master_ambient_LED;master_accent_spotlights"), ("scene_control_relation", "direct_load_or_scene_command_by_final_product_underlying_load_switching_required")],
    }
    for equipment_id, items in panel_requirements.items():
        for key, value in items: upsert(equipment_id, key, value, unit="count" if key.endswith("quantity") or key.endswith("count") else "")
        upsert(equipment_id, "interface_center_coordinates", "unknown", origin="pending", status="pending", source="", locator="", blocks="yes", notes="待 M07/门套/见光板完成面证据确认真实安装面和净距。")

    for equipment_id in ("NET-AP-R09", "NET-AP-R14"):
        for key in (
            "home_run_cable",
            "local_220v_for_current_candidates",
            "candidate_ap362e_max_power",
            "candidate_ap362e_power",
            "candidate_rgeap262e_power",
            "cable_continuity_test",
        ):
            if (equipment_id, key) in by_pair:
                by_pair[(equipment_id, key)]["blocks_release"] = "no"
        if (equipment_id, "cable_continuity_test") in by_pair:
            by_pair[(equipment_id, "cable_continuity_test")]["notes"] = "历史重复字段；现行停止条件为 cable_continuity，不重复计入阻断。"
        for key, value in [("poe_supply_required", "yes"), ("local_220v_power_required", "no"), ("network_topology", "CAT6_star_home_run_to_weak_current_cabinet_PoE_switch"), ("poe_switch_port_class", "IEEE_802.3at_backward_compatible_with_802.3af")]: upsert(equipment_id, key, value)
        upsert(equipment_id, "final_endpoint_power_w", "unknown", origin="pending", status="pending", source="", locator="", blocks="yes", notes="按实购 AP 官方功率和交换机总 PoE 预算闭合。")
        upsert(equipment_id, "cable_continuity", "unknown", origin="pending", status="pending", source="", locator="", blocks="yes", notes="现场网线通断和端接顺序尚未测试。")

    write_csv(DECISIONS / "equipment-installation-requirements.csv", req_fields, requirements)

    owner_fields, owners = read_csv(DECISIONS / "owner-input-register.csv")
    owner_updates = {
        "E303-NS01-CIRCUIT": ("两个同时 10A 端口的不利工况为 20A；三个物理插座位分配两路独立 C16 RCBO/2.5mm²，火锅与第二台同时设备分路", "自定义确认", "设备名牌仍 unknown，但不阻断粗装回路。"),
        "E303-NS02-CIRCUIT": ("咖啡机按候选上限 2600W 预留独立 C16 RCBO/2.5mm² 和接线盒；手冲壶/磨豆机使用独立的小功率辅助回路", "自定义确认", "咖啡机仍选装；最终插座面板按实购插头。"),
        "E302-ENTRY-PANEL": ("1 个逻辑面板，最少 4 键：客厅/书房/餐厅双控+照明总控；总控不切断连续负载", "自定义确认", "准确 SKU、端子、负载与底盒仍由厂家/电气方闭合。"),
        "E302-MASTER-PANEL": ("Master A 最少 2 键：主卧氛围/重点；Master B 最少 3 键：客厅/书房/餐厅真实双控；两块各 1 个逻辑面板", "自定义确认", "真实安装面和产品证据门仍未闭合。"),
        "E304-AP-POWER": ("主卧/次卧 AP 均由弱电箱 PoE 交换机供电，CAT6 星型回箱，AP 点不另设 220V；交换端口按 802.3at 向下兼容 802.3af 预留", "自定义确认", "准确型号、端口/总功率预算和网线通断仍待厂家/现场。"),
    }
    for row in owners:
        if row["input_id"] in owner_updates:
            value, status, note = owner_updates[row["input_id"]]
            row["candidate_value"] = value
            row["status"] = status
            row["evidence_reference"] = joined(row["evidence_reference"], RESEARCH_SOURCE_ID)
            row["notes"] = note
        if row["input_id"] == "E303-NS02-FORM":
            row["source_basis"] = "业主确认三类连接需求与使用位置；旧 2.35–4.0kW 聚合包络已作废，负荷按咖啡机 2600W 独立预留和小功率辅助支路分开表达"
            row["notes"] = "安装型式仍待线性轨道或自制翻盖偏好；不得把旧聚合包络写回回路计算。"
    write_csv(DECISIONS / "owner-input-register.csv", owner_fields, owners)

    close_fields, closeouts = read_csv(DECISIONS / "owner-input-closeout-rules.csv")
    close_updates = {
        "E303-NS01-CIRCUIT": ("external_product_evidence", "设备商家/业主", "APP-001～003 实购后的型号和名牌功率，仅用于设备表完整性，不重开已闭合的两路 C16 粗装"),
        "E303-NS02-CIRCUIT": ("external_product_evidence", "咖啡机商家/电气方", "实购咖啡机型号、名牌和插头照片，用于最终面板匹配；粗装已按 2600W/C16 闭合"),
        "E302-ENTRY-PANEL": ("external_product_evidence", "智能面板厂家/电气方", "准确 SKU、端子图、每路负载/浪涌、证书编号、底盒净深和平嵌适配"),
        "E302-MASTER-PANEL": ("external_product_evidence", "智能面板厂家/照明设计", "Master A/B 准确 SKU、直接负载/场景方案、双控端子图、负载和底盒/见光板安装图"),
        "E304-AP-POWER": ("site_and_product_evidence", "网络设备方/现场", "实购 AP 准确型号和功耗、PoE 交换机端口/总预算、CAT6 通断与端接顺序测试"),
    }
    for row in closeouts:
        if row["input_id"] in close_updates:
            row["closeout_kind"], row["responsible_party"], row["required_evidence"] = close_updates[row["input_id"]]
            row["automatic_close_allowed"] = "no"
            row["notes"] = "功能/粗装结论已闭合；此门只关闭最终产品、现场测试或面板安装。"
    write_csv(DECISIONS / "owner-input-closeout-rules.csv", close_fields, closeouts)

    rule_fields, rules = read_csv(DECISIONS / "elec-design-rules.csv")
    rule_updates = {
        "ELEC-DES-031": ("CTRL-ENTRY-A;CTRL-MASTER-A;CTRL-MASTER-B", "three_distinct_logical_panels_with_fixed_function_maps", "Entry 4 keys; Master A 2 keys; Master B 3 keys", "research_conclusion_product_and_geometry_pending", "功能数量、键序、负载与场景关系已闭合；安装面与准确 SKU 仍待证据。"),
        "ELEC-DES-033": ("R03", "NS01_two_simultaneous_10A_port_adverse_case", "20A design demand; two independent C16 RCBO; 2.5mm2 Cu", "research_conclusion", "不利工况用端口上限，不虚构 APP-001～003 功率。"),
        "ELEC-DES-035": ("R06", "NS02_optional_coffee_worst_candidate_reserve", "2600W; dedicated C16 RCBO; 2.5mm2 Cu; connection box; faceplate by actual plug", "research_conclusion", "GS3/E1 Prima EXP 不需现在二选一；给排水只作可选预留。"),
        "ELEC-DES-038": ("R09;R14", "PoE_AP_architecture", "CAT6 star home run; weak-current-cabinet PoE switch; 802.3at ports backward-compatible with 802.3af; no local 220V", "research_conclusion_final_model_budget_and_testing_pending", "供电架构已闭合；准确 AP 功耗、总预算、通断和温升仍待证据。"),
    }
    for row in rules:
        if row["rule_id"] in rule_updates:
            row["locations"], row["device_or_datum"], row["value"], row["status"], row["notes"] = rule_updates[row["rule_id"]]
            row["basis"] = f"{RESEARCH_SOURCE_ID};业主 2026-08-15 第一组电气决策"
            row["confidence"] = "1.00"
            row["review_required"] = "no" if row["rule_id"] in {"ELEC-DES-033", "ELEC-DES-035"} else "yes"
    write_csv(DECISIONS / "elec-design-rules.csv", rule_fields, rules)

    switch_fields, switch_rows = read_csv(DECISIONS / "e302-switch-product-review.csv")
    for row in switch_rows:
        if row["gate"] == "control_role_closed":
            row["passed"] = "true"
            row["failure_reason"] = ""
            row["required_evidence"] = "已由 ELEC-RESEARCH-20260815-001 闭合逻辑面板数量、最少键数、键序、负载与场景关系"
        if row["panel_id"] == "CTRL-MASTER-A" and row["gate"] == "box_and_joinery_interface_verified":
            row["failure_reason"] = "见光板内的阻燃背盒、固定基层、端子弯线净深、散热和可检修性均未闭合"
            row["required_evidence"] = "准确面板完整尺寸/底盒图、端子空间、阻燃背盒、固定基层、通风散热和下方拆装检修节点"
    write_csv(DECISIONS / "e302-switch-product-review.csv", switch_fields, switch_rows)

    topology_fields, topology = read_csv(DECISIONS / "e304-network-topology.csv")
    for row in topology:
        if row["link_id"] in {"NET-LINK-002", "NET-LINK-004"}:
            row.update({
                "poe_required": "yes",
                "local_power_required": "no",
                "poe_standard": "IEEE_802.3at_port_backward_compatible_with_802.3af",
                "max_endpoint_power_w": "",
                "physical_port": "TBD_after_cable_test",
                "source_evidence_id": RESEARCH_SOURCE_ID,
                "verification_status": "poe_architecture_confirmed_final_model_power_budget_and_cable_continuity_pending",
                "status": "research_conclusion",
                "notes": "CAT6 星型回弱电箱 PoE 交换机；AP 点不设 220V。最终端口号、单端口功率和总预算按实购型号与现场测线闭合。",
            })
    write_csv(DECISIONS / "e304-network-topology.csv", topology_fields, topology)

    drawing_fields, drawings = read_csv(DECISIONS / "drawing-register.csv")
    drawing_notes = {
        "E-302": "Entry A/Master A/Master B 已分别闭合为 4/2/3 键逻辑功能，含三组真实双控、主卧氛围/重点和只作用于照明的总控命令；准确 SKU、端子、负载/浪涌、底盒、平嵌适配及 M07/门套净距仍阻断最终发布。",
        "E-303": "设备插头、墙面插座、支路保护和独立回路已拆分语义。烤箱为明确 16A 插头/16A 插座+独立 C16；洗烘为两个 10A 插座共用一路洗衣区 C16，非两个 16A 插座/两回路。NS-01 按两个 10A 端口同时的 20A 不利工况分两路 C16；NS-02 咖啡机按 2600W 预留独立 C16 与接线盒，面板按实购插头。说明书未给的接口中心仍 unknown。",
        "E-304": "两个吸顶 AP 已采用 CAT6 星型回弱电箱 PoE 交换机的拓扑，AP 点不另设 220V；交换端口按 802.3at 向下兼容 802.3af 预留。实购 AP 功耗、PoE 总预算、弱电箱净尺寸/柜门温升、网线通断和端接顺序仍待厂家/现场证据。",
    }
    for row in drawings:
        if row["sheet_number"] in drawing_notes:
            row["notes"] = drawing_notes[row["sheet_number"]]
    write_csv(DECISIONS / "drawing-register.csv", drawing_fields, drawings)


if __name__ == "__main__":
    main()
