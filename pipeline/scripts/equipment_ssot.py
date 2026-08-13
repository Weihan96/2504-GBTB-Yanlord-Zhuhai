#!/usr/bin/env python3
"""Build and enforce the equipment installation single source of truth.

The formal IFC is read-only. Legacy decision CSVs are migration inputs once and
then become compatibility projections generated from the canonical three-table
model.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import ifcopenshell
import ifcopenshell.util.element


MASTER = "pipeline/decisions/equipment-register.csv"
REQUIREMENTS = "pipeline/decisions/equipment-installation-requirements.csv"
SOURCES = "pipeline/decisions/source-evidence-register.csv"
SCHEMA = "pipeline/schemas/equipment-ssot.schema.json"
FORMAL_IFC = "2504 GBTB Yanlord Zhuhai.ifc"

APPLIANCE_FIELDS = [
    "appliance_id", "appliance_name", "category", "storage_location_candidate",
    "use_location_candidate", "storage_location_confirmed", "use_location_confirmed",
    "quantity", "rated_power_w", "simultaneous_group", "water_required",
    "drain_required", "gas_required", "ventilation_required", "model",
    "evidence_reference", "status", "notes",
]
FURNITURE_FIELDS = [
    "selector_kind", "selector_value", "manufacturer", "product_name",
    "product_variant", "product_category", "intended_use", "source_url",
    "identity_basis", "confidence", "human_review_required", "status",
]
FURNITURE_ROLE_FIELDS = [
    "selector_kind", "selector_value", "installation_role", "basis", "confidence",
    "human_review_required", "status",
]
HVAC_FIELDS = [
    "evidence_id", "equipment_ids", "equipment_global_ids", "ifc_type_global_id",
    "ifc_type_name", "ifc_model_text", "nominal_body_width_mm", "nominal_body_depth_mm",
    "nominal_body_height_mm", "capacity_group", "gas_pipe_od_mm", "liquid_pipe_od_mm",
    "drain_pipe_od_mm", "drain_slope_min", "drain_slope_max", "connection_side_basis",
    "port_coordinate_status", "source_document", "source_url", "source_sha256", "pdf_page",
    "evidence", "status", "confidence", "review_required", "formal_ifc_write_allowed", "notes",
]
ELEC_EVIDENCE_FIELDS = [
    "evidence_id", "discipline", "sheet_id", "decision_scope", "source_kind",
    "source_document", "source_sha256", "source_locator", "evidence", "proves",
    "does_not_prove", "status", "confidence", "review_required",
    "formal_ifc_write_allowed", "notes",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def split_ids(value: str) -> list[str]:
    return [part.strip() for part in re.split(r"[;；]", value or "") if part.strip()]


def joined(values: Iterable[str]) -> str:
    return ";".join(dict.fromkeys(value for value in values if value))


def status_for_appliance(value: str) -> tuple[str, str]:
    if value == "已确认":
        return "confirmed", "purchased_arrived"
    if value == "部分确认":
        return "partial", "candidate"
    return "pending", "not_selected"


def evidence_local_path(root: Path, document: str) -> str:
    candidates = [
        root / document,
        root / "drawings/evidence" / Path(document).name,
        root / "tmp/dwg" / Path(document).name,
        root.parent / "图纸/矩阵纵横" / Path(document).name,
        root.parent / "图纸" / Path(document).name,
        root.parent / "Blender文件" / Path(document).name,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate.relative_to(root)) if candidate.is_relative_to(root) else str(candidate)
    for candidate in sorted((root / "tmp/pdfs").glob(f"**/{Path(document).stem}*.pdf")):
        if candidate.is_file():
            return str(candidate.relative_to(root))
    return ""


def source_from_elec(root: Path, row: dict[str, str]) -> dict[str, str]:
    document = row["source_document"]
    is_url = document.startswith("http://") or document.startswith("https://")
    return {
        "source_id": row["evidence_id"], "discipline": row["discipline"],
        "sheet_id": row["sheet_id"], "decision_scope": row["decision_scope"],
        "source_kind": row["source_kind"], "source_document": document,
        "source_url": document if is_url else "", "local_path": "" if is_url else evidence_local_path(root, document),
        "sha256": row["source_sha256"], "locator": row["source_locator"],
        "evidence": row["evidence"], "proves": row["proves"],
        "does_not_prove": row["does_not_prove"], "status": row["status"],
        "confidence": row["confidence"], "review_required": row["review_required"],
        "formal_ifc_write_allowed": row["formal_ifc_write_allowed"],
        "manufacturer": "", "model_scope": "", "revision": "", "publication_date": "",
        "legacy_targets": "elec-source-evidence.csv",
        "legacy_projection_json": json.dumps(row, ensure_ascii=False, separators=(",", ":")),
        "notes": row["notes"],
    }


def source_from_hvac(row: dict[str, str]) -> dict[str, str]:
    return {
        "source_id": row["evidence_id"], "discipline": "HVAC/PLUM", "sheet_id": "RCP-1/M-401",
        "decision_scope": f"{row['equipment_ids']} 厂家接口", "source_kind": "official_product_pdf" if row["source_url"] else "legacy_model_evidence",
        "source_document": row["source_document"], "source_url": row["source_url"],
        "local_path": "", "sha256": row["source_sha256"], "locator": row["pdf_page"],
        "evidence": row["evidence"], "proves": "设备族、机体尺寸及已列明的冷媒/排水接口条件",
        "does_not_prove": row["notes"], "status": row["status"], "confidence": row["confidence"],
        "review_required": row["review_required"], "formal_ifc_write_allowed": row["formal_ifc_write_allowed"],
        "manufacturer": "日立/海信日立" if row["source_url"] else "", "model_scope": row["ifc_model_text"],
        "revision": "", "publication_date": "", "legacy_targets": "rcp1-hvac-equipment-interface-evidence.csv",
        "legacy_projection_json": json.dumps(row, ensure_ascii=False, separators=(",", ":")), "notes": row["notes"],
    }


def add_requirement(rows: list[dict[str, str]], equipment_id: str, key: str, value: str,
                    *, discipline: str = "MULTI", unit: str = "", origin: str = "project_candidate",
                    status: str = "candidate", source_id: str = "", blocks: str = "yes",
                    locator: str = "", notes: str = "") -> None:
    value = str(value or "")
    numeric = value if re.fullmatch(r"-?\d+(?:\.\d+)?", value) else ""
    text = "" if numeric else value
    rows.append({
        "requirement_id": f"REQ-{len(rows)+1:04d}", "equipment_id": equipment_id,
        "discipline": discipline, "parameter_key": key, "value_text": text,
        "value_number": numeric, "unit": unit, "datum": "",
        "value_origin": origin, "status": status, "source_id": source_id,
        "source_locator": locator, "blocks_release": blocks, "notes": notes,
    })


OFFICIAL_REQUIREMENTS: dict[str, list[tuple[str, str, str, str]]] = {
    "APP-009": [
        ("rated_voltage", "220", "V", "APP-DW-SPEC-001"), ("rated_current", "10", "A", "APP-DW-SPEC-001"),
        ("product_height", "775", "mm", "APP-DW-SPEC-001"), ("product_width", "598", "mm", "APP-DW-SPEC-001"),
        ("product_depth", "550", "mm", "APP-DW-SPEC-001"), ("niche_height_min", "780", "mm", "APP-DW-INSTALL-001"),
        ("niche_height_max", "835", "mm", "APP-DW-INSTALL-001"), ("niche_width_min", "600", "mm", "APP-DW-INSTALL-001"),
        ("niche_width_max", "608", "mm", "APP-DW-INSTALL-001"), ("niche_depth_min", "550", "mm", "APP-DW-INSTALL-001"),
        ("service_opening_width_min", "100", "mm", "APP-DW-INSTALL-001"), ("service_opening_height_min", "50", "mm", "APP-DW-INSTALL-001"),
        ("water_connection", "G3/4 cold water", "", "APP-DW-INSTALL-001"), ("drain_connection_od", "38", "mm", "APP-DW-INSTALL-001"),
        ("water_pressure_min", "0.05", "MPa", "APP-DW-INSTALL-001"), ("water_pressure_max", "1", "MPa", "APP-DW-INSTALL-001"),
        ("water_flow_min", "10", "L/min", "APP-DW-INSTALL-001"), ("cold_water_temperature_max", "25", "°C", "APP-DW-INSTALL-001"),
        ("door_panel_width", "594", "mm", "APP-DW-INSTALL-001"), ("door_panel_thickness_max", "20", "mm", "APP-DW-INSTALL-001"),
        ("door_panel_weight_min", "3.5", "kg", "APP-DW-INSTALL-001"), ("door_panel_weight_max", "8.5", "kg", "APP-DW-INSTALL-001"),
    ],
    "APP-011": [
        ("rated_voltage", "220", "V", "APP-011-OFFICIAL-001"), ("frequency", "50", "Hz", "APP-011-OFFICIAL-001"),
        ("rated_current", "16", "A", "APP-011-OFFICIAL-001"), ("product_height", "595", "mm", "APP-011-OFFICIAL-001"),
        ("product_width", "594", "mm", "APP-011-OFFICIAL-001"), ("product_depth", "548", "mm", "APP-011-OFFICIAL-001"),
        ("niche_width_min", "560", "mm", "APP-011-OFFICIAL-001"), ("niche_width_max", "568", "mm", "APP-011-OFFICIAL-001"),
        ("niche_height_min", "585", "mm", "APP-011-OFFICIAL-001"), ("niche_height_max", "595", "mm", "APP-011-OFFICIAL-001"),
        ("niche_depth_min", "550", "mm", "APP-011-OFFICIAL-001"), ("ventilation_area_min", "200", "cm2", "APP-011-OFFICIAL-001"),
    ],
    "APP-012": [
        ("gas_type", "天然气12T", "", "APP-012-PHOTO-001"), ("gas_pressure", "2.0", "kPa", "APP-012-PHOTO-001"),
        ("product_width", "920", "mm", "APP-012-OFFICIAL-001"), ("product_depth", "500", "mm", "APP-012-OFFICIAL-001"),
        ("product_height", "150", "mm", "APP-012-OFFICIAL-001"), ("worktop_cutout_width", "880", "mm", "APP-012-OFFICIAL-001"),
        ("worktop_cutout_depth", "455", "mm", "APP-012-OFFICIAL-001"), ("worktop_cutout_corner_radius", "15", "mm", "APP-012-OFFICIAL-001"),
        ("cabinet_ventilation_area_min", "100", "cm2", "APP-012-OFFICIAL-001"), ("clearance_below_min", "100", "mm", "APP-012-OFFICIAL-001"),
        ("separator_to_worktop_min", "120", "mm", "APP-012-OFFICIAL-001"),
    ],
    "APP-013": [
        ("rated_voltage", "220", "V", "APP-013-OFFICIAL-001"), ("frequency", "50", "Hz", "APP-013-OFFICIAL-001"),
        ("product_width", "895", "mm", "APP-013-OFFICIAL-001"), ("product_depth", "345", "mm", "APP-013-OFFICIAL-001"),
        ("product_height", "873", "mm", "APP-013-OFFICIAL-001"), ("duct_inner_diameter_min", "160", "mm", "APP-013-OFFICIAL-001"),
        ("wall_opening_diameter_candidate", "180", "mm", "APP-013-OFFICIAL-001"), ("ceiling_opening_diameter", "190", "mm", "APP-013-OFFICIAL-001"),
        ("access_opening_width", "400", "mm", "APP-013-OFFICIAL-001"), ("access_opening_height", "400", "mm", "APP-013-OFFICIAL-001"),
        ("clearance_to_gas_hob_pan_support_min", "650", "mm", "APP-013-OFFICIAL-001"),
    ],
}
OFFICIAL_REQUIREMENTS["APP-010"] = OFFICIAL_REQUIREMENTS["APP-009"]
CANDIDATE_POWER_RANGES = {
    "APP-001": (1500, 2200), "APP-002": (300, 1200), "APP-003": (800, 1500),
    "APP-004": (1200, 1800), "APP-005": (1000, 1800), "APP-006": (150, 400),
}


def migrate(root: Path) -> dict[str, Any]:
    master_path, req_path, source_path = (root / MASTER, root / REQUIREMENTS, root / SOURCES)
    if any(path.exists() for path in (master_path, req_path, source_path)):
        raise RuntimeError("canonical equipment tables already exist; migration is intentionally one-shot")
    schema = json.loads((root / SCHEMA).read_text(encoding="utf-8"))
    master_fields = schema["tables"]["equipment-register.csv"]["columns"]
    req_fields = schema["tables"]["equipment-installation-requirements.csv"]["columns"]
    source_fields = schema["tables"]["source-evidence-register.csv"]["columns"]
    appliances = read_csv(root / "pipeline/decisions/appliance-input-register.csv")
    furniture = read_csv(root / "pipeline/decisions/furniture-product-register.csv")
    doors = read_csv(root / "pipeline/decisions/a104-door-window-review.csv")
    hvac = read_csv(root / "pipeline/decisions/rcp1-hvac-equipment-interface-evidence.csv")
    elec_sources = read_csv(root / "pipeline/decisions/elec-source-evidence.csv")
    model = ifcopenshell.open(root / FORMAL_IFC)
    ifc_hash = sha256(root / FORMAL_IFC)
    masters: list[dict[str, str]] = []
    reqs: list[dict[str, str]] = []
    sources = [source_from_elec(root, row) for row in elec_sources] + [source_from_hvac(row) for row in hvac]
    sources.append({
        "source_id": "IFC-FORMAL-001", "discipline": "MULTI", "sheet_id": "", "decision_scope": "正式 IFC 对象身份与既有几何",
        "source_kind": "formal_ifc", "source_document": FORMAL_IFC, "source_url": "", "local_path": FORMAL_IFC,
        "sha256": ifc_hash, "locator": "GlobalId / type assignment / spatial containment", "evidence": "正式 IFC 当前对象身份与类型",
        "proves": "IFC 中对象存在、GlobalId、类型与当前容器", "does_not_prove": "最终产品选型或厂家安装要求",
        "status": "verified_formal_source", "confidence": "1.00", "review_required": "no", "formal_ifc_write_allowed": "no",
        "manufacturer": "", "model_scope": "", "revision": "", "publication_date": "", "legacy_targets": "", "legacy_projection_json": "",
        "notes": "由设备真源审计器逐次校验哈希",
    })
    existing_source_ids = {row["source_id"] for row in sources}
    for item in furniture:
        sid = f"FUR-SRC-{item['selector_value']}"
        if sid not in existing_source_ids:
            sources.append({
                "source_id": sid, "discipline": "INT1", "sheet_id": "S-701", "decision_scope": f"{item['selector_value']} 产品身份",
                "source_kind": "official_product_web", "source_document": item["source_url"], "source_url": item["source_url"], "local_path": "",
                "sha256": "not_applicable_live_official_web", "locator": "official product page", "evidence": item["identity_basis"],
                "proves": "厂家与产品系列身份", "does_not_prove": "项目最终规格、五金配置与安装尺寸", "status": "verified_official_source",
                "confidence": item["confidence"], "review_required": item["human_review_required"], "formal_ifc_write_allowed": "no",
                "manufacturer": item["manufacturer"], "model_scope": item["product_name"], "revision": "", "publication_date": "",
                "legacy_targets": "", "legacy_projection_json": "", "notes": "由原家具产品登记迁移",
            })
            existing_source_ids.add(sid)

    occurrence_owner: dict[str, str] = {}
    known_ifc_map = {"APP-011": "3PQOXKxgj6IftqWXFXMQXG", "APP-013": "288GLY62v8kPPydA1lAK8W"}
    for item in appliances:
        decision, procurement = status_for_appliance(item["status"])
        gids = known_ifc_map.get(item["appliance_id"], "")
        known_appliance_sources = {
            "APP-009": "APP-DW-SPEC-001;APP-DW-INSTALL-001;APP-DW-MANUAL-001",
            "APP-010": "APP-DW-SPEC-001;APP-DW-INSTALL-001;APP-DW-MANUAL-001",
        }
        source_ids = known_appliance_sources.get(item["appliance_id"], joined(split_ids(item["evidence_reference"]) if ";" in item["evidence_reference"] and "http" not in item["evidence_reference"] else []))
        domain = "NETWORK" if item["category"] == "网络设备" else "APPLIANCE"
        masters.append({
            "equipment_id": item["appliance_id"], "domain": domain, "category": item["category"], "item_name": item["appliance_name"],
            "manufacturer": "Siemens" if item["model"].startswith("Siemens") else "Huawei" if "Huawei" in item["appliance_name"] else "",
            "model": item["model"], "variant": "", "quantity": item["quantity"], "procurement_status": procurement,
            "decision_status": decision, "storage_location_candidate": item["storage_location_candidate"],
            "use_location_candidate": item["use_location_candidate"], "storage_location_confirmed": item["storage_location_confirmed"],
            "use_location_confirmed": item["use_location_confirmed"], "schedule_included": "yes", "selector_kind": "logical_input",
            "selector_value": item["appliance_id"], "ifc_class": "IfcElectricAppliance" if gids else "", "ifc_type_name": "",
            "ifc_type_global_id": "", "ifc_global_ids": gids, "source_ids": source_ids, "identity_basis": item["evidence_reference"],
            "confidence": "1.00" if decision == "confirmed" else "0.80", "human_review_required": "no" if decision == "confirmed" else "yes",
            "legacy_kind": "appliance", "legacy_id": item["appliance_id"], "notes": item["notes"],
        })
        if gids:
            occurrence_owner[gids] = item["appliance_id"]
        for key, value, unit, discipline in (
            ("rated_power", item["rated_power_w"], "W", "ELEC"), ("simultaneous_group", item["simultaneous_group"], "", "ELEC"),
            ("water_required", item["water_required"], "", "PLUM"), ("drain_required", item["drain_required"], "", "PLUM"),
            ("gas_required", item["gas_required"], "", "GAS"), ("ventilation_required", item["ventilation_required"], "", "HVAC"),
        ):
            add_requirement(reqs, item["appliance_id"], key, value, discipline=discipline, unit=unit,
                            origin="user_input", status="confirmed" if value not in {"", "model_dependent"} else "pending",
                            blocks="yes" if value in {"", "model_dependent"} else "no")
        for key, value, unit, sid in OFFICIAL_REQUIREMENTS.get(item["appliance_id"], []):
            add_requirement(reqs, item["appliance_id"], key, value, discipline="MULTI", unit=unit,
                            origin="official_exact_model", status="confirmed", source_id=sid, blocks="no")
        if item["appliance_id"] in CANDIDATE_POWER_RANGES:
            minimum, maximum = CANDIDATE_POWER_RANGES[item["appliance_id"]]
            add_requirement(reqs, item["appliance_id"], "candidate_power_min", str(minimum), discipline="ELEC", unit="W",
                            origin="project_candidate", status="candidate", blocks="yes", notes="常见值仅用于回路容量情景，不代表已选产品")
            add_requirement(reqs, item["appliance_id"], "candidate_power_max", str(maximum), discipline="ELEC", unit="W",
                            origin="project_candidate", status="candidate", blocks="yes", notes="常见值仅用于回路容量情景，不代表已选产品")

    furniture_types: dict[str, Any] = {}
    for product in model.by_type("IfcFurniture"):
        assigned = ifcopenshell.util.element.get_type(product)
        if assigned is not None:
            furniture_types.setdefault(assigned.Name or "", assigned)
    for item in furniture:
        assigned = furniture_types.get(item["selector_value"])
        products = [p for p in model.by_type("IfcFurniture") if (ifcopenshell.util.element.get_type(p) or object()) is assigned] if assigned else []
        gids = joined(p.GlobalId for p in products)
        eid = f"FUR-{item['selector_value']}"
        masters.append({
            "equipment_id": eid, "domain": "FURNITURE", "category": item["product_category"], "item_name": item["product_name"],
            "manufacturer": item["manufacturer"], "model": item["product_name"], "variant": item["product_variant"], "quantity": str(len(products) or 1),
            "procurement_status": "selected", "decision_status": "confirmed", "storage_location_candidate": "", "use_location_candidate": item["intended_use"],
            "storage_location_confirmed": "", "use_location_confirmed": item["intended_use"], "schedule_included": "yes", "selector_kind": item["selector_kind"],
            "selector_value": item["selector_value"], "ifc_class": "IfcFurniture", "ifc_type_name": item["selector_value"],
            "ifc_type_global_id": assigned.GlobalId if assigned else "", "ifc_global_ids": gids, "source_ids": f"FUR-SRC-{item['selector_value']}",
            "identity_basis": item["identity_basis"], "confidence": item["confidence"], "human_review_required": item["human_review_required"],
            "legacy_kind": "furniture_product", "legacy_id": item["selector_value"], "notes": "",
        })
        for gid in split_ids(gids): occurrence_owner[gid] = eid
        if item["intended_use"]:
            add_requirement(reqs, eid, "intended_use", item["intended_use"], discipline="INT1", origin="user_input", status="confirmed", blocks="no")

    for item in doors:
        eid = f"DW-{item['candidate_id']}"
        masters.append({
            "equipment_id": eid, "domain": "DOOR_WINDOW", "category": item["ifc_class"], "item_name": f"{item['candidate_tag']} · {item['name']}",
            "manufacturer": "", "model": "", "variant": "", "quantity": "1", "procurement_status": "not_selected",
            "decision_status": "ifc_observed", "storage_location_candidate": "", "use_location_candidate": item["space_candidates"],
            "storage_location_confirmed": "", "use_location_confirmed": "", "schedule_included": "yes", "selector_kind": "global_id",
            "selector_value": item["global_id"], "ifc_class": item["ifc_class"], "ifc_type_name": "", "ifc_type_global_id": "",
            "ifc_global_ids": item["global_id"], "source_ids": "IFC-FORMAL-001", "identity_basis": item["basis"], "confidence": item["confidence"],
            "human_review_required": item["review_required"], "legacy_kind": "a104", "legacy_id": item["candidate_id"], "notes": item["review_question"],
        })
        occurrence_owner[item["global_id"]] = eid
        for key, value, unit in (("nominal_width", item["nominal_width_mm"], "mm"), ("nominal_height", item["nominal_height_mm"], "mm"),
                                 ("operation_type", item["operation_type"], ""), ("host_relation", item["host_relation"], ""),
                                 ("sill_or_threshold_z", item["sill_or_threshold_z_mm"], "mm")):
            add_requirement(reqs, eid, key, value, discipline="ARCH", unit=unit, origin="ifc_observed",
                            status="observed" if value else "pending", source_id="IFC-FORMAL-001", blocks="yes" if not value else "no")

    for item in hvac:
        eid = f"HVAC-{item['evidence_id'].split('-')[-1]}"
        masters.append({
            "equipment_id": eid, "domain": "HVAC", "category": "ducted_indoor_unit", "item_name": item["equipment_ids"],
            "manufacturer": "日立/海信日立" if item["source_url"] else "", "model": item["ifc_model_text"], "variant": item["capacity_group"],
            "quantity": str(len(split_ids(item["equipment_ids"]))), "procurement_status": "candidate", "decision_status": "candidate",
            "storage_location_candidate": "", "use_location_candidate": "吊顶内", "storage_location_confirmed": "", "use_location_confirmed": "",
            "schedule_included": "yes", "selector_kind": "ifc_type_global_id" if item["ifc_type_global_id"] else "logical_input",
            "selector_value": item["ifc_type_global_id"] or item["equipment_ids"], "ifc_class": "IfcElectricAppliance" if item["equipment_global_ids"] else "",
            "ifc_type_name": item["ifc_type_name"], "ifc_type_global_id": item["ifc_type_global_id"], "ifc_global_ids": item["equipment_global_ids"],
            "source_ids": item["evidence_id"], "identity_basis": item["evidence"], "confidence": item["confidence"],
            "human_review_required": item["review_required"], "legacy_kind": "hvac_interface", "legacy_id": item["evidence_id"], "notes": item["notes"],
        })
        for gid in split_ids(item["equipment_global_ids"]): occurrence_owner[gid] = eid
        for key, unit in (("nominal_body_width_mm", "mm"), ("nominal_body_depth_mm", "mm"), ("nominal_body_height_mm", "mm"),
                          ("capacity_group", ""), ("gas_pipe_od_mm", "mm"), ("liquid_pipe_od_mm", "mm"), ("drain_pipe_od_mm", "mm"),
                          ("drain_slope_min", ""), ("drain_slope_max", ""), ("connection_side_basis", ""), ("port_coordinate_status", "")):
            value = item[key]
            add_requirement(reqs, eid, key.removesuffix("_mm"), value, discipline="HVAC/PLUM", unit=unit,
                            origin="official_model_family" if item["source_url"] else "project_candidate",
                            status="candidate" if value else "pending", source_id=item["evidence_id"], blocks="yes" if not value else "no")

    def register_ifc_groups(cls: str, domain: str, prefix: str, selected_names: set[str] | None = None) -> None:
        groups: dict[tuple[str, str, str], list[Any]] = defaultdict(list)
        for product in model.by_type(cls):
            if product.GlobalId in occurrence_owner:
                continue
            if selected_names is not None and (product.Name or "") not in selected_names:
                continue
            assigned = ifcopenshell.util.element.get_type(product)
            key = (assigned.GlobalId if assigned else "", assigned.Name if assigned else "", product.Name or "")
            groups[key].append(product)
        for index, ((type_gid, type_name, occurrence_name), products) in enumerate(sorted(groups.items()), 1):
            eid = f"{prefix}-{index:03d}"
            gids = joined(p.GlobalId for p in products)
            masters.append({
                "equipment_id": eid, "domain": domain, "category": cls, "item_name": type_name or occurrence_name or cls,
                "manufacturer": "", "model": type_name, "variant": "", "quantity": str(len(products)), "procurement_status": "existing",
                "decision_status": "ifc_observed", "storage_location_candidate": "", "use_location_candidate": "", "storage_location_confirmed": "",
                "use_location_confirmed": "", "schedule_included": "no", "selector_kind": "ifc_type_global_id" if type_gid else "global_id",
                "selector_value": type_gid or products[0].GlobalId, "ifc_class": cls, "ifc_type_name": type_name, "ifc_type_global_id": type_gid,
                "ifc_global_ids": gids, "source_ids": "IFC-FORMAL-001", "identity_basis": "formal IFC inventory",
                "confidence": "1.00", "human_review_required": "yes", "legacy_kind": "ifc_inventory", "legacy_id": "", 
                "notes": "对象已机械纳入覆盖；最终产品与厂家安装条件待确认",
            })
            for gid in split_ids(gids): occurrence_owner[gid] = eid

    register_ifc_groups("IfcElectricAppliance", "APPLIANCE", "ELECIFC")
    register_ifc_groups("IfcFurniture", "FURNITURE", "FURIFC")
    register_ifc_groups("IfcSanitaryTerminal", "SANITARY", "SAN")
    register_ifc_groups("IfcWasteTerminal", "DRAINAGE", "WASTE")
    register_ifc_groups("IfcSensor", "SAFETY", "SENSOR")
    register_ifc_groups("IfcElementAssembly", "ASSEMBLY", "PLUMASM", {"WC", "WC_01", "Drain Center"})

    write_csv(master_path, master_fields, masters)
    write_csv(req_path, req_fields, reqs)
    write_csv(source_path, source_fields, sources)
    return {"master_count": len(masters), "requirement_count": len(reqs), "source_count": len(sources), "ifc_sha256": ifc_hash}


def load_canonical(root: Path) -> tuple[dict[str, Any], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    schema = json.loads((root / SCHEMA).read_text(encoding="utf-8"))
    return schema, read_csv(root / MASTER), read_csv(root / REQUIREMENTS), read_csv(root / SOURCES)


def validate(root: Path) -> dict[str, Any]:
    schema, masters, reqs, sources = load_canonical(root)
    errors: list[str] = []
    tables = [("equipment-register.csv", masters), ("equipment-installation-requirements.csv", reqs), ("source-evidence-register.csv", sources)]
    for name, rows in tables:
        spec = schema["tables"][name]
        path = root / "pipeline/decisions" / name
        with path.open(newline="", encoding="utf-8-sig") as stream:
            actual = next(csv.reader(stream))
        if actual != spec["columns"]:
            errors.append(f"{name}: header differs from schema")
        keys = [row[spec["primary_key"]] for row in rows]
        if len(keys) != len(set(keys)) or any(not key for key in keys):
            errors.append(f"{name}: primary key is blank or duplicated")
        for row in rows:
            for field in spec.get("required", []):
                if not row.get(field, "").strip(): errors.append(f"{name}:{row.get(spec['primary_key'])}: missing {field}")
            for field, allowed in spec.get("enums", {}).items():
                if row.get(field, "") not in allowed: errors.append(f"{name}:{row.get(spec['primary_key'])}: invalid {field}={row.get(field)}")
    master_ids = {row["equipment_id"] for row in masters}
    source_ids = {row["source_id"] for row in sources}
    for row in reqs:
        if row["equipment_id"] not in master_ids: errors.append(f"{row['requirement_id']}: missing equipment {row['equipment_id']}")
        if row["source_id"] and row["source_id"] not in source_ids: errors.append(f"{row['requirement_id']}: missing source {row['source_id']}")
        if bool(row["value_text"]) == bool(row["value_number"]) and row["status"] != "pending":
            errors.append(f"{row['requirement_id']}: exactly one value column is required")
        if row["value_origin"] in {"official_exact_model", "official_model_family"} and not row["source_id"]:
            errors.append(f"{row['requirement_id']}: official value lacks source")
    for row in masters:
        for sid in split_ids(row["source_ids"]):
            if sid not in source_ids: errors.append(f"{row['equipment_id']}: missing source {sid}")
        if row["procurement_status"] == "purchased_arrived" and (not row["model"] or not row["source_ids"]):
            errors.append(f"{row['equipment_id']}: arrived product lacks model/source")
    for row in sources:
        local = row["local_path"]
        expected = row["sha256"]
        exact_hash = bool(re.fullmatch(r"[0-9a-f]{64}", expected))
        not_applicable_hash = expected.startswith("not_applicable_")
        if not exact_hash and not not_applicable_hash:
            errors.append(
                f"{row['source_id']}: sha256 must be 64 lowercase hexadecimal characters or an explicit not_applicable_* marker"
            )
        if exact_hash:
            resolved = local or evidence_local_path(root, row["source_document"])
            path = Path(resolved) if Path(resolved).is_absolute() else root / resolved
            if not path.is_file():
                resolved = evidence_local_path(root, row["source_document"])
                path = Path(resolved) if Path(resolved).is_absolute() else root / resolved
            if not path.is_file(): errors.append(f"{row['source_id']}: local evidence missing: {local or row['source_document']}")
            elif sha256(path) != expected: errors.append(f"{row['source_id']}: local evidence hash mismatch")
        if row["legacy_projection_json"]:
            try: json.loads(row["legacy_projection_json"])
            except json.JSONDecodeError: errors.append(f"{row['source_id']}: invalid legacy projection JSON")
    if errors:
        raise RuntimeError("equipment SSOT validation failed:\n- " + "\n- ".join(errors))
    return {"master_count": len(masters), "requirement_count": len(reqs), "source_count": len(sources), "schema_version": schema["schema_version"]}


def requirement_map(reqs: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = defaultdict(dict)
    for row in reqs:
        result[row["equipment_id"]][row["parameter_key"]] = row["value_number"] or row["value_text"]
    return result


def appliance_projection_rows(root: Path) -> list[dict[str, str]]:
    """Return the owner-appliance compatibility view from canonical data."""
    _, masters, reqs, _ = load_canonical(root)
    reqmap = requirement_map(reqs)
    result: list[dict[str, str]] = []
    for row in masters:
        if row["legacy_kind"] != "appliance":
            continue
        values = reqmap[row["equipment_id"]]
        status = "已确认" if row["decision_status"] == "confirmed" else "部分确认" if row["decision_status"] == "partial" else "待填写"
        result.append({
            "appliance_id": row["legacy_id"], "appliance_name": row["item_name"], "category": row["category"],
            "storage_location_candidate": row["storage_location_candidate"], "use_location_candidate": row["use_location_candidate"],
            "storage_location_confirmed": row["storage_location_confirmed"], "use_location_confirmed": row["use_location_confirmed"],
            "quantity": row["quantity"], "rated_power_w": values.get("rated_power", ""), "simultaneous_group": values.get("simultaneous_group", ""),
            "water_required": values.get("water_required", ""), "drain_required": values.get("drain_required", ""),
            "gas_required": values.get("gas_required", ""), "ventilation_required": values.get("ventilation_required", ""),
            "model": row["model"], "evidence_reference": row["identity_basis"], "status": status, "notes": row["notes"],
            "_candidate_power_range_w": [values.get("candidate_power_min", ""), values.get("candidate_power_max", "")],
        })
    return result


def apply_owner_appliances(root: Path, appliances: list[dict[str, str]]) -> None:
    """Apply validated workbook appliance rows to canonical tables."""
    schema, masters, reqs, sources = load_canonical(root)
    master_by_legacy = {row["legacy_id"]: row for row in masters if row["legacy_kind"] == "appliance"}
    source_ids = {row["source_id"] for row in sources}
    req_by_key = {(row["equipment_id"], row["parameter_key"]): row for row in reqs}
    for item in appliances:
        master = master_by_legacy[item["appliance_id"]]
        decision, procurement = status_for_appliance(item["status"])
        prior_evidence = master["identity_basis"]
        prior_source_ids = master["source_ids"]
        master.update({
            "item_name": item["appliance_name"], "category": item["category"], "model": item["model"],
            "quantity": item["quantity"], "procurement_status": procurement, "decision_status": decision,
            "storage_location_candidate": item["storage_location_candidate"], "use_location_candidate": item["use_location_candidate"],
            "storage_location_confirmed": item["storage_location_confirmed"], "use_location_confirmed": item["use_location_confirmed"],
            "identity_basis": item["evidence_reference"], "human_review_required": "no" if decision == "confirmed" else "yes", "notes": item["notes"],
        })
        references = split_ids(item["evidence_reference"]) if ";" in item["evidence_reference"] and "http" not in item["evidence_reference"] else []
        canonical_refs = split_ids(prior_source_ids) if item["evidence_reference"] == prior_evidence else [sid for sid in references if sid in source_ids]
        if item["evidence_reference"] and not canonical_refs:
            sid = f"OWNER-SRC-{item['appliance_id']}-{hashlib.sha256(item['evidence_reference'].encode()).hexdigest()[:10]}"
            if sid not in source_ids:
                sources.append({
                    "source_id": sid, "discipline": "MULTI", "sheet_id": "OWNER-INPUT", "decision_scope": f"{item['appliance_id']} 业主设备输入",
                    "source_kind": "owner_supplied_reference", "source_document": item["evidence_reference"],
                    "source_url": item["evidence_reference"] if item["evidence_reference"].startswith("http") else "", "local_path": "",
                    "sha256": "not_applicable_owner_reference", "locator": "owner workbook", "evidence": item["evidence_reference"],
                    "proves": "业主提供的设备参考", "does_not_prove": "官方安装条件或现场安装完成状态", "status": "owner_supplied_unverified",
                    "confidence": "0.80", "review_required": "yes", "formal_ifc_write_allowed": "no", "manufacturer": "", "model_scope": item["model"],
                    "revision": "", "publication_date": "", "legacy_targets": "", "legacy_projection_json": "", "notes": "同步脚本自动登记；须补官方资料核验",
                })
                source_ids.add(sid)
            canonical_refs = [sid]
        master["source_ids"] = joined(canonical_refs)
        for key, value in (("rated_power", item["rated_power_w"]), ("simultaneous_group", item["simultaneous_group"]),
                           ("water_required", item["water_required"]), ("drain_required", item["drain_required"]),
                           ("gas_required", item["gas_required"]), ("ventilation_required", item["ventilation_required"])):
            requirement = req_by_key[(item["appliance_id"], key)]
            requirement["value_number"] = value if re.fullmatch(r"-?\d+(?:\.\d+)?", value or "") else ""
            requirement["value_text"] = "" if requirement["value_number"] else value
            requirement["status"] = "confirmed" if value not in {"", "model_dependent"} else "pending"
            requirement["blocks_release"] = "no" if requirement["status"] == "confirmed" else "yes"
    write_csv(root / MASTER, schema["tables"]["equipment-register.csv"]["columns"], masters)
    write_csv(root / REQUIREMENTS, schema["tables"]["equipment-installation-requirements.csv"]["columns"], reqs)
    write_csv(root / SOURCES, schema["tables"]["source-evidence-register.csv"]["columns"], sources)
    validate(root)


def projections(root: Path, *, check: bool = False) -> dict[str, Any]:
    validate(root)
    schema, masters, reqs, sources = load_canonical(root)
    reqmap = requirement_map(reqs)
    appliance_rows = [
        {field: row.get(field, "") for field in APPLIANCE_FIELDS}
        for row in appliance_projection_rows(root)
    ]
    furniture_rows = []
    source_by_id = {row["source_id"]: row for row in sources}
    for row in masters:
        if row["legacy_kind"] != "furniture_product": continue
        furniture_rows.append({
            "selector_kind": row["selector_kind"], "selector_value": row["selector_value"], "manufacturer": row["manufacturer"],
            "product_name": row["item_name"], "product_variant": row["variant"], "product_category": row["category"],
            "intended_use": row["use_location_confirmed"], "source_url": source_by_id[row["source_ids"]]["source_url"],
            "identity_basis": row["identity_basis"], "confidence": row["confidence"], "human_review_required": row["human_review_required"], "status": "confirmed",
        })
    source_projections: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in sources:
        if not row["legacy_projection_json"]: continue
        payload = json.loads(row["legacy_projection_json"])
        for target in split_ids(row["legacy_targets"]): source_projections[target].append(payload)
    outputs: list[tuple[Path, list[str], list[dict[str, str]]]] = [
        (root / "pipeline/decisions/appliance-input-register.csv", APPLIANCE_FIELDS, appliance_rows),
        (root / "pipeline/decisions/furniture-product-register.csv", FURNITURE_FIELDS, furniture_rows),
        (root / "pipeline/decisions/elec-source-evidence.csv", ELEC_EVIDENCE_FIELDS, source_projections["elec-source-evidence.csv"]),
        (root / "pipeline/decisions/rcp1-hvac-equipment-interface-evidence.csv", HVAC_FIELDS, source_projections["rcp1-hvac-equipment-interface-evidence.csv"]),
        (root / "pipeline/decisions/furniture-installation-role.csv", FURNITURE_ROLE_FIELDS, source_projections["furniture-installation-role.csv"]),
    ]
    drift: list[str] = []
    for path, fields, rows in outputs:
        if check:
            current = read_csv(path)
            if current != rows: drift.append(str(path.relative_to(root)))
        else:
            write_csv(path, fields, rows)
    if drift: raise RuntimeError("compatibility projection drift: " + ", ".join(drift))
    return {"appliance_rows": len(appliance_rows), "furniture_rows": len(furniture_rows),
            "elec_evidence_rows": len(source_projections["elec-source-evidence.csv"]), "hvac_evidence_rows": len(source_projections["rcp1-hvac-equipment-interface-evidence.csv"]),
            "furniture_role_rows": len(source_projections["furniture-installation-role.csv"]),
            "mode": "check" if check else "write"}


def audit_ifc(root: Path, report_path: Path) -> dict[str, Any]:
    validate(root)
    _, masters, _, _ = load_canonical(root)
    model = ifcopenshell.open(root / FORMAL_IFC)
    summary, coverage = assert_ifc_coverage(root, masters, model)
    report = {"mode": "read_only_equipment_ssot_ifc_coverage", "summary": summary,
              "coverage": coverage}
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def assert_ifc_coverage(root: Path, masters: list[dict[str, str]], model: Any) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Prove canonical equipment identities still match the formal IFC.

    This intentionally does not call ``validate`` first, so it can guard a
    controlled source-hash refresh after a legitimate, geometry-safe IFC save.
    """
    scoped = ["IfcElectricAppliance", "IfcFurniture", "IfcSanitaryTerminal", "IfcWasteTerminal", "IfcSensor", "IfcDoor", "IfcWindow"]
    expected: dict[str, str] = {}
    for cls in scoped:
        for product in model.by_type(cls): expected[product.GlobalId] = cls
    for product in model.by_type("IfcElementAssembly"):
        if (product.Name or "") in {"WC", "WC_01", "Drain Center"}: expected[product.GlobalId] = "IfcElementAssembly"
    for row in masters:
        if row["ifc_class"] != "IfcBuildingElementProxy":
            continue
        for gid in split_ids(row["ifc_global_ids"]):
            product = model.by_guid(gid)
            if not product or product.is_a() != "IfcBuildingElementProxy":
                raise RuntimeError(f"registered equipment proxy is missing or has the wrong IFC class: {gid}")
            expected[gid] = "IfcBuildingElementProxy"
    owners: dict[str, list[str]] = defaultdict(list)
    for row in masters:
        for gid in split_ids(row["ifc_global_ids"]): owners[gid].append(row["equipment_id"])
    missing = sorted(set(expected) - set(owners))
    duplicate = {gid: ids for gid, ids in owners.items() if gid in expected and len(ids) != 1}
    unexpected = sorted(set(owners) - set(expected))
    if missing or duplicate or unexpected:
        raise RuntimeError(f"IFC coverage failed: missing={missing}, duplicate={duplicate}, unexpected={unexpected}")
    summary = {"formal_ifc_sha256": sha256(root / FORMAL_IFC), "scope_count": len(expected),
               "covered_count": len(expected), "missing_count": 0, "duplicate_count": 0,
               "class_counts": dict(Counter(expected.values()))}
    coverage = [{"global_id": gid, "ifc_class": expected[gid], "equipment_id": owners[gid][0]} for gid in sorted(expected)]
    return summary, coverage


def refresh_ifc_source_hashes(root: Path, *, dry_run: bool = False) -> dict[str, Any]:
    """Refresh only evidence rows whose local source is the formal IFC.

    Product identity coverage is checked before any CSV write. This keeps an
    elevation-only or relationship-only IFC save from silently accepting lost,
    duplicated, or newly unregistered equipment identities.
    """
    schema, masters, _, sources = load_canonical(root)
    model = ifcopenshell.open(root / FORMAL_IFC)
    coverage, _ = assert_ifc_coverage(root, masters, model)
    current_hash = coverage["formal_ifc_sha256"]
    targets = [row for row in sources if row["source_document"] == FORMAL_IFC and row["local_path"] == FORMAL_IFC]
    if not targets:
        raise RuntimeError("no canonical evidence rows point explicitly to the formal IFC")
    prior_hashes = sorted({row["sha256"] for row in targets})
    changed = [row for row in targets if row["sha256"] != current_hash]
    if not dry_run and changed:
        for row in changed:
            row["sha256"] = current_hash
        write_csv(root / SOURCES, schema["tables"]["source-evidence-register.csv"]["columns"], sources)
        validate(root)
    return {
        "formal_ifc_sha256": current_hash,
        "scoped_ifc_object_count": coverage["scope_count"],
        "target_source_count": len(targets),
        "updated_source_count": len(changed),
        "prior_hashes": prior_hashes,
        "dry_run": dry_run,
    }


def reconcile(root: Path, report_path: Path) -> dict[str, Any]:
    validate(root)
    _, masters, requirements, sources = load_canonical(root)
    legacy_counts = {
        "appliance-input-register.csv": len(read_csv(root / "pipeline/decisions/appliance-input-register.csv")),
        "furniture-product-register.csv": len(read_csv(root / "pipeline/decisions/furniture-product-register.csv")),
        "elec-source-evidence.csv": len(read_csv(root / "pipeline/decisions/elec-source-evidence.csv")),
        "rcp1-hvac-equipment-interface-evidence.csv": len(read_csv(root / "pipeline/decisions/rcp1-hvac-equipment-interface-evidence.csv")),
        "furniture-installation-role.csv": len(read_csv(root / "pipeline/decisions/furniture-installation-role.csv")),
    }
    canonical_counts = {
        "appliance_records": sum(row["legacy_kind"] == "appliance" for row in masters),
        "furniture_product_records": sum(row["legacy_kind"] == "furniture_product" for row in masters),
        "hvac_interface_records": sum(row["legacy_kind"] == "hvac_interface" for row in masters),
        "elec_evidence_projection_records": sum("elec-source-evidence.csv" in split_ids(row["legacy_targets"]) for row in sources),
        "hvac_evidence_projection_records": sum("rcp1-hvac-equipment-interface-evidence.csv" in split_ids(row["legacy_targets"]) for row in sources),
        "furniture_role_projection_records": sum("furniture-installation-role.csv" in split_ids(row["legacy_targets"]) for row in sources),
        "furniture_role_requirements": sum(row["parameter_key"] == "installation_role" for row in requirements),
    }
    expected = {
        "appliance_records": legacy_counts["appliance-input-register.csv"],
        "furniture_product_records": legacy_counts["furniture-product-register.csv"],
        "hvac_interface_records": legacy_counts["rcp1-hvac-equipment-interface-evidence.csv"],
        "elec_evidence_projection_records": legacy_counts["elec-source-evidence.csv"],
        "hvac_evidence_projection_records": legacy_counts["rcp1-hvac-equipment-interface-evidence.csv"],
        "furniture_role_projection_records": legacy_counts["furniture-installation-role.csv"],
        "furniture_role_requirements": legacy_counts["furniture-installation-role.csv"],
    }
    if canonical_counts != expected:
        raise RuntimeError(f"migration reconciliation failed: canonical={canonical_counts}, expected={expected}")
    projections(root, check=True)
    summary = {"legacy_counts": legacy_counts, "canonical_migrated_counts": canonical_counts,
               "row_count_parity": True, "compatibility_projection_parity": True}
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps({"mode": "equipment_ssot_migration_reconciliation", "summary": summary}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def normalize_source_paths(root: Path) -> dict[str, Any]:
    schema, masters, _, sources = load_canonical(root)
    updated = 0
    for row in sources:
        resolved = evidence_local_path(root, row["source_document"])
        if resolved and row["local_path"] != resolved:
            row["local_path"] = resolved
            updated += 1
    write_csv(root / SOURCES, schema["tables"]["source-evidence-register.csv"]["columns"], sources)
    link_updates = 0
    for row in masters:
        if row["equipment_id"] in {"APP-009", "APP-010"}:
            expected = "APP-DW-SPEC-001;APP-DW-INSTALL-001;APP-DW-MANUAL-001"
            if row["source_ids"] != expected:
                row["source_ids"] = expected
                link_updates += 1
    write_csv(root / MASTER, schema["tables"]["equipment-register.csv"]["columns"], masters)
    validate(root)
    return {"updated_paths": updated, "updated_links": link_updates, "source_count": len(sources)}


def migrate_furniture_roles(root: Path) -> dict[str, Any]:
    schema, masters, reqs, sources = load_canonical(root)
    roles = read_csv(root / "pipeline/decisions/furniture-installation-role.csv")
    existing_source_ids = {row["source_id"] for row in sources}
    existing_keys = {(row["equipment_id"], row["parameter_key"]) for row in reqs}
    by_type = {row["ifc_type_name"]: row for row in masters if row["domain"] == "FURNITURE" and row["ifc_type_name"]}
    by_global = {gid: row for row in masters if row["domain"] == "FURNITURE" for gid in split_ids(row["ifc_global_ids"])}
    added_requirements = 0
    added_sources = 0
    for role in roles:
        master = by_type.get(role["selector_value"]) if role["selector_kind"] == "type_name" else by_global.get(role["selector_value"])
        if master is None:
            raise RuntimeError(f"furniture role selector is not represented in canonical IFC coverage: {role}")
        sid = f"FUR-ROLE-SRC-{hashlib.sha256((role['selector_kind'] + ':' + role['selector_value']).encode()).hexdigest()[:12]}"
        if sid not in existing_source_ids:
            sources.append({
                "source_id": sid, "discipline": "INT1", "sheet_id": "INT1/S-701", "decision_scope": f"{role['selector_value']} 安装角色",
                "source_kind": "project_ifc_classification", "source_document": FORMAL_IFC, "source_url": "", "local_path": FORMAL_IFC,
                "sha256": sha256(root / FORMAL_IFC), "locator": f"{role['selector_kind']}={role['selector_value']}", "evidence": role["basis"],
                "proves": f"installation_role={role['installation_role']}", "does_not_prove": "最终产品安装节点、连接件或现场完成状态",
                "status": role["status"], "confidence": role["confidence"], "review_required": role["human_review_required"],
                "formal_ifc_write_allowed": "no", "manufacturer": "", "model_scope": role["selector_value"], "revision": "", "publication_date": "",
                "legacy_targets": "furniture-installation-role.csv", "legacy_projection_json": json.dumps(role, ensure_ascii=False, separators=(",", ":")),
                "notes": "由原家具安装角色登记无损迁移",
            })
            existing_source_ids.add(sid)
            added_sources += 1
        key = (master["equipment_id"], "installation_role")
        if key not in existing_keys:
            add_requirement(reqs, master["equipment_id"], "installation_role", role["installation_role"], discipline="INT1",
                            origin="project_candidate", status="confirmed" if role["status"] == "implemented" else "candidate",
                            source_id=sid, blocks="yes" if role["human_review_required"] == "yes" else "no",
                            notes=f"{role['basis']} | confidence={role['confidence']}")
            existing_keys.add(key)
            added_requirements += 1
    write_csv(root / REQUIREMENTS, schema["tables"]["equipment-installation-requirements.csv"]["columns"], reqs)
    write_csv(root / SOURCES, schema["tables"]["source-evidence-register.csv"]["columns"], sources)
    validate(root)
    return {"role_rows": len(roles), "added_requirements": added_requirements, "added_sources": added_sources}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["migrate", "validate", "sync-projections", "check-projections", "audit-ifc", "refresh-ifc-source-hashes", "reconcile", "normalize-sources", "migrate-furniture-roles", "all"])
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--report", type=Path, default=Path("build/equipment-ssot/ifc-coverage.json"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    result: dict[str, Any] = {}
    if args.command == "migrate":
        result = {"canonical": migrate(root), "furniture_roles": migrate_furniture_roles(root),
                  "source_paths": normalize_source_paths(root)}
    elif args.command == "validate": result = validate(root)
    elif args.command == "sync-projections": result = projections(root)
    elif args.command == "check-projections": result = projections(root, check=True)
    elif args.command == "audit-ifc": result = audit_ifc(root, args.report if args.report.is_absolute() else root / args.report)
    elif args.command == "refresh-ifc-source-hashes": result = refresh_ifc_source_hashes(root, dry_run=args.dry_run)
    elif args.command == "reconcile": result = reconcile(root, root / "build/equipment-ssot/migration-reconciliation.json")
    elif args.command == "normalize-sources": result = normalize_source_paths(root)
    elif args.command == "migrate-furniture-roles": result = migrate_furniture_roles(root)
    else:
        result = {"validate": validate(root), "projections": projections(root, check=True),
                  "ifc_coverage": audit_ifc(root, args.report if args.report.is_absolute() else root / args.report),
                  "migration_reconciliation": reconcile(root, root / "build/equipment-ssot/migration-reconciliation.json")}
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
