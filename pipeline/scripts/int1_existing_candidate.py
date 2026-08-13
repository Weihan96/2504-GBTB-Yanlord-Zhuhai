#!/usr/bin/env python3
"""Generate the INT1 existing-object candidate register without writing IFC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element

from furniture_anchor_audit import (
    classify_furniture,
    find_role_decision,
    find_selector_record,
    product_identity,
    read_product_register,
    read_role_decisions,
)
from equipment_ssot import load_canonical, split_ids, validate as validate_equipment_ssot


EXPECTED_FURNITURE = 89
EXPECTED_FIXED = 70
EXPECTED_LOOSE = 19
EXPECTED_SHEETS = {"I-501", "I-502", "I-503", "I-504"}
KITCHEN_PRODUCT_IDS = ("APP-009", "APP-010", "APP-011", "APP-012", "APP-013")

KITCHEN_APPLIANCE_IDS = {
    "0UOnmuAdP1MPy6p3olwiEU",  # WD01 water dispenser
    "3PQOXKxgj6IftqWXFXMQXG",  # OV01 oven
    "288GLY62v8kPPydA1lAK8W",  # HD01 Siemens hood
}
KITCHEN_SANITARY_IDS = {
    "2C9JmXaZXCX9I$vMds4lYK",  # SIN01 island sink
    "1$3UG05MPFmumejbg1KJ$I",  # FAU01 island faucet
}
KITCHEN_CONTEXT_IDS = {
    "2Ca2tGerPBj8kvUTjlUtDl",  # kitchen sliding-door pier and storage
    "3NIZVJBir6rQlVwhKs84Ca",  # cabinet assembly
}
BATHROOM_CONTEXT_IDS = {
    "3ARl_CqPrCWQrA6_$W073W",  # guest bathroom depressed slab
    "1UixnpnQb97wuNAGrTmMJd",  # master bathroom depressed slab
    "04rs0EDjn2EvxytEQSxWRB",  # confirmed new WAL170 bathroom wall
    "3jha5L04zBjh0$pMl_tHLy",  # confirmed new WAL170 bathroom wall
}
WASHTOWER_ID = "3JAkt8PsX7vPfGKWLK5EKp"
GEBERIT_FLUSH_PLATE_IDS = {
    "2gFgcOYEXEaQWAzcKulTFt",
    "2lDPsdQevFSfeOThtjSlPG",
}

CSV_FIELDS = [
    "record_kind",
    "sheet_id",
    "global_id",
    "ifc_class",
    "object_name",
    "type_name",
    "type_description",
    "container",
    "installation_role",
    "bbox_min_mm",
    "bbox_max_mm",
    "dimensions_mm",
    "manufacturer",
    "product_name",
    "product_variant",
    "candidate_use",
    "dimension_status",
    "review_status",
    "human_review_required",
    "basis",
    "confidence",
    "equipment_id",
    "equipment_decision_status",
    "equipment_procurement_status",
    "equipment_source_ids",
    "installation_requirement_count",
    "blocking_requirement_count",
    "stop_condition",
    "source_ifc_sha256",
]

BLOCKERS = [
    {
        "issue_id": "INT1-BLOCK-001",
        "scope": "equipment_installation_drawings",
        "condition": (
            "Product installation evidence and project rough-in locations are incomplete; "
            "available exact-model constraints are listed separately."
        ),
        "stop_condition": (
            "Do not freeze cabinet openings, dedicated services, ventilation, "
            "heat-clearance or access dimensions from the existing model envelope."
        ),
    },
    {
        "issue_id": "INT1-BLOCK-002",
        "scope": "gas_meter_site_survey",
        "condition": (
            "Gas meter, valve and pipe dimensions and authority requirements "
            "have not been verified on site."
        ),
        "stop_condition": (
            "Do not freeze enclosure depth, openings, ventilation, observation "
            "or access provisions."
        ),
    },
    {
        "issue_id": "INT1-BLOCK-003",
        "scope": "sanitary_rough_in_drawings",
        "condition": "Manufacturer sanitary rough-in drawings are missing.",
        "stop_condition": (
            "Do not infer water, waste, fixing or wall/floor connection centres "
            "from complex product-shell geometry."
        ),
    },
    {
        "issue_id": "INT1-BLOCK-004",
        "scope": "hardware_motion_envelopes",
        "condition": (
            "Door, drawer, mirror-cabinet and kitchen-pier hardware motion "
            "envelopes are not confirmed."
        ),
        "stop_condition": (
            "Do not publish opening, collision or ergonomic clearances for "
            "unmodelled hardware."
        ),
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rounded(values: list[float]) -> list[float]:
    return [round(value, 6) for value in values]


def world_bbox(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
    bbox_registry: dict[str, tuple[list[float], list[float], list[float]]],
) -> tuple[list[float], list[float], list[float]]:
    if product.GlobalId in bbox_registry:
        return bbox_registry[product.GlobalId]
    try:
        shape = ifcopenshell.geom.create_shape(settings, product)
        values = list(shape.geometry.verts)
    except RuntimeError:
        values = []
    if not values:
        child_boxes = [
            world_bbox(settings, child, bbox_registry)
            for child in ifcopenshell.util.element.get_decomposition(product)
            if getattr(child, "ObjectPlacement", None)
        ]
        if not child_boxes:
            raise RuntimeError(f"product has no tessellated vertices: {product.GlobalId}")
        minimum = [min(box[0][axis] for box in child_boxes) for axis in range(3)]
        maximum = [max(box[1][axis] for box in child_boxes) for axis in range(3)]
        dimensions = [maximum[axis] - minimum[axis] for axis in range(3)]
        result = rounded(minimum), rounded(maximum), rounded(dimensions)
        bbox_registry[product.GlobalId] = result
        return result
    vertices = [
        [values[index] * 1000.0, values[index + 1] * 1000.0, values[index + 2] * 1000.0]
        for index in range(0, len(values), 3)
    ]
    minimum = [min(vertex[axis] for vertex in vertices) for axis in range(3)]
    maximum = [max(vertex[axis] for vertex in vertices) for axis in range(3)]
    dimensions = [maximum[axis] - minimum[axis] for axis in range(3)]
    result = rounded(minimum), rounded(maximum), rounded(dimensions)
    bbox_registry[product.GlobalId] = result
    return result


def read_bbox_registry(
    path: Path,
    expected_source_hash: str,
) -> dict[str, tuple[list[float], list[float], list[float]]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("source", {}).get("sha256") != expected_source_hash:
        raise RuntimeError("construction-surface audit does not match the formal IFC")
    result = {}
    for record in data.get("records", []):
        result[record["global_id"]] = (
            rounded(record["bbox_min_mm"]),
            rounded(record["bbox_max_mm"]),
            rounded(record["dimensions_mm"]),
        )
    return result


def container_name(product: ifcopenshell.entity_instance) -> str:
    container = ifcopenshell.util.element.get_container(product)
    return str(getattr(container, "Name", "") or "")


def assigned_type(product: ifcopenshell.entity_instance) -> Any:
    return ifcopenshell.util.element.get_type(product)


def sheet_for_furniture(container: str) -> str:
    if container in {"KITCHEN", "VVD"}:
        return "I-501"
    if container == "BATHM":
        return "I-502"
    return "I-504"


def object_record(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
    sheet_id: str,
    source_hash: str,
    bbox_registry: dict[str, tuple[list[float], list[float], list[float]]],
    installation_role: str,
    candidate_use: str,
    basis: str,
    confidence: float,
    identity: dict[str, Any] | None = None,
    equipment: dict[str, Any] | None = None,
    review_status: str = "existing_model_candidate",
) -> dict[str, Any]:
    type_object = assigned_type(product)
    minimum, maximum, dimensions = world_bbox(settings, product, bbox_registry)
    identity = identity or {}
    equipment = equipment or {}
    return {
        "record_kind": "existing_object",
        "sheet_id": sheet_id,
        "global_id": product.GlobalId,
        "ifc_class": product.is_a(),
        "object_name": product.Name or "",
        "type_name": getattr(type_object, "Name", "") or "",
        "type_description": (
            getattr(type_object, "Description", "")
            or getattr(product, "Description", "")
            or ""
        ),
        "container": container_name(product),
        "installation_role": installation_role,
        "bbox_min_mm": minimum,
        "bbox_max_mm": maximum,
        "dimensions_mm": dimensions,
        "manufacturer": identity.get("manufacturer", ""),
        "product_name": identity.get("product_name", ""),
        "product_variant": identity.get("product_variant", ""),
        "candidate_use": candidate_use,
        "dimension_status": "existing_world_bbox_not_fabrication_dimension",
        "review_status": review_status,
        "human_review_required": False,
        "basis": basis,
        "confidence": confidence,
        "equipment_id": equipment.get("equipment_id", ""),
        "equipment_decision_status": equipment.get("decision_status", ""),
        "equipment_procurement_status": equipment.get("procurement_status", ""),
        "equipment_source_ids": equipment.get("source_ids", []),
        "installation_requirement_count": len(equipment.get("requirements", [])),
        "blocking_requirement_count": sum(
            row.get("blocks_release") == "yes"
            for row in equipment.get("requirements", [])
        ),
        "stop_condition": "Do not infer fabrication openings or service connection centres from this envelope.",
        "source_ifc_sha256": source_hash,
    }


def blocker_csv_record(blocker: dict[str, str], source_hash: str) -> dict[str, Any]:
    return {
        "record_kind": "blocker",
        "sheet_id": "I-501/I-502/I-503/I-504",
        "global_id": "",
        "ifc_class": "",
        "object_name": blocker["issue_id"],
        "type_name": blocker["scope"],
        "type_description": blocker["condition"],
        "container": "",
        "installation_role": "",
        "bbox_min_mm": "",
        "bbox_max_mm": "",
        "dimensions_mm": "",
        "manufacturer": "",
        "product_name": "",
        "product_variant": "",
        "candidate_use": "human_or_site_input_required",
        "dimension_status": "not_available",
        "review_status": "BLOCK",
        "human_review_required": True,
        "basis": blocker["condition"],
        "confidence": 1.0,
        "equipment_id": "",
        "equipment_decision_status": "",
        "equipment_procurement_status": "",
        "equipment_source_ids": "",
        "installation_requirement_count": "",
        "blocking_requirement_count": "",
        "stop_condition": blocker["stop_condition"],
        "source_ifc_sha256": source_hash,
    }


def encode_csv_value(value: Any) -> Any:
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, bool):
        return "yes" if value else "no"
    return value


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({field: encode_csv_value(record.get(field, "")) for field in CSV_FIELDS})


def write_summary(path: Path, report: dict[str, Any]) -> None:
    summary = report["summary"]
    lines = [
        "# INT1 Existing-Object Candidate",
        "",
        f"- Source IFC SHA-256: `{report['source']['ifc_sha256']}`",
        f"- Furniture: {summary['furniture_total']} = {summary['fixed_furniture']} fixed + {summary['loose_furniture']} loose",
        f"- Existing-object records: {summary['existing_object_records']}",
        f"- Candidate sheets: {', '.join(summary['sheet_ids'])}",
        f"- Candidate generation: {'PASS' if report['gates']['candidate_generation_pass'] else 'FAIL'}",
        f"- INT1 completion gate: BLOCK ({summary['block_count']} disclosed blockers)",
        "",
        "## Blocking inputs",
        "",
    ]
    for blocker in report["blockers"]:
        lines.append(
            f"- **{blocker['issue_id']} / {blocker['scope']}**: {blocker['stop_condition']}"
        )
    lines.extend(
        [
            "",
            "> Bounding boxes are current world-coordinate model envelopes, not fabrication dimensions. This candidate does not write IFC or infer openings, service centres, rough-in points, or hardware motion.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--decision-csv", required=True, type=Path)
    parser.add_argument("--construction-audit", required=True, type=Path)
    return parser.parse_args()


def equipment_owner_index(
    equipment_rows: list[dict[str, str]],
    requirement_rows: list[dict[str, str]],
) -> dict[str, dict[str, Any]]:
    requirements: dict[str, list[dict[str, str]]] = {}
    for row in requirement_rows:
        requirements.setdefault(row["equipment_id"], []).append(row)
    owners: dict[str, dict[str, Any]] = {}
    for row in equipment_rows:
        snapshot = {
            **row,
            "source_ids": split_ids(row["source_ids"]),
            "requirements": requirements.get(row["equipment_id"], []),
        }
        for global_id in split_ids(row["ifc_global_ids"]):
            if global_id in owners:
                raise RuntimeError(f"equipment SSOT assigns {global_id} more than once")
            owners[global_id] = snapshot
    return owners


def requirement_snapshot(row: dict[str, str]) -> dict[str, Any]:
    return {
        "requirement_id": row["requirement_id"],
        "discipline": row["discipline"],
        "parameter_key": row["parameter_key"],
        "value": row["value_number"] or row["value_text"],
        "unit": row["unit"],
        "value_origin": row["value_origin"],
        "status": row["status"],
        "source_id": row["source_id"],
        "blocks_release": row["blocks_release"] == "yes",
    }


def kitchen_product_snapshot(
    equipment_id: str,
    equipment_by_id: dict[str, dict[str, str]],
    requirements_by_id: dict[str, list[dict[str, str]]],
) -> dict[str, Any]:
    equipment = equipment_by_id[equipment_id]
    requirements = [
        requirement_snapshot(row)
        for row in requirements_by_id.get(equipment_id, [])
        if row["status"] == "confirmed"
        and row["value_origin"] in {"official_exact_model", "user_input"}
    ]
    official_source_ids = sorted(
        {row["source_id"] for row in requirements if row["source_id"]}
    )
    return {
        "equipment_id": equipment_id,
        "sheet_id": "I-501",
        "item_name": equipment["item_name"],
        "manufacturer": equipment["manufacturer"],
        "model": equipment["model"],
        "procurement_status": equipment["procurement_status"],
        "decision_status": equipment["decision_status"],
        "use_location": (
            equipment["use_location_confirmed"]
            or equipment["use_location_candidate"]
        ),
        "official_source_ids": official_source_ids,
        "requirements": requirements,
        "product_constraints_available": bool(official_source_ids),
        "project_interface_status": "unlocated",
        "final_product_confirmed": (
            equipment["procurement_status"] == "purchased_arrived"
            and equipment["decision_status"] == "confirmed"
        ),
        "release_blockers": [
            "project rough-in XYZ, valve, drain adaptor, hose path and opening position are unlocated",
            *(
                []
                if equipment["procurement_status"] == "purchased_arrived"
                else ["final procurement and installed-unit identity are unconfirmed"]
            ),
        ],
    }


def main() -> None:
    args = parse_args()
    source_path = args.input.resolve()
    source_hash = sha256(source_path)
    root = Path(__file__).resolve().parents[2]
    validate_equipment_ssot(root)
    _, equipment_rows, requirement_rows, _ = load_canonical(root)
    owners = equipment_owner_index(equipment_rows, requirement_rows)
    equipment_by_id = {row["equipment_id"]: row for row in equipment_rows}
    requirements_by_id: dict[str, list[dict[str, str]]] = {}
    for row in requirement_rows:
        requirements_by_id.setdefault(row["equipment_id"], []).append(row)
    model = ifcopenshell.open(source_path)
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    canonical_register = root / "pipeline/decisions/equipment-register.csv"
    role_decisions = read_role_decisions(canonical_register)
    product_records = read_product_register(canonical_register)
    bbox_registry = read_bbox_registry(args.construction_audit, source_hash)

    records: list[dict[str, Any]] = []
    furniture_role_counts: Counter[str] = Counter()
    furniture_sheet_counts: Counter[str] = Counter()
    furniture_ids: set[str] = set()

    for furniture in model.by_type("IfcFurniture"):
        type_object = assigned_type(furniture)
        type_name = getattr(type_object, "Name", None)
        type_predefined = getattr(type_object, "PredefinedType", None)
        decision = find_role_decision(furniture.GlobalId, type_name, role_decisions)
        role = classify_furniture(type_predefined, decision)
        installation_role = role["installation_role"]
        if installation_role not in {"fixed_furniture", "loose_furniture"}:
            raise RuntimeError(f"unresolved furniture role: {furniture.GlobalId}")
        identity_record = find_selector_record(
            furniture.GlobalId, type_name, product_records
        )
        identity = product_identity(identity_record)
        sheet_id = sheet_for_furniture(container_name(furniture))
        records.append(
            object_record(
                settings,
                furniture,
                sheet_id,
                source_hash,
                bbox_registry,
                installation_role,
                "fixed_joinery_detail" if installation_role == "fixed_furniture" else "clearance_context_only",
                role["basis"],
                float(role["confidence"]),
                identity,
                owners.get(furniture.GlobalId),
            )
        )
        furniture_role_counts[installation_role] += 1
        furniture_sheet_counts[sheet_id] += 1
        furniture_ids.add(furniture.GlobalId)

    if len(furniture_ids) != EXPECTED_FURNITURE:
        raise RuntimeError(f"expected {EXPECTED_FURNITURE} unique furniture objects")
    if furniture_role_counts != Counter(
        {"fixed_furniture": EXPECTED_FIXED, "loose_furniture": EXPECTED_LOOSE}
    ):
        raise RuntimeError(f"unexpected furniture role split: {furniture_role_counts}")

    for global_id in sorted(KITCHEN_APPLIANCE_IDS):
        product = model.by_guid(global_id)
        records.append(
            object_record(
                settings, product, "I-501", source_hash, bbox_registry, "fixed_equipment",
                "equipment_envelope_coordination", "Exact IFC GlobalId and assigned type preserve the existing equipment identity.", 1.0,
                equipment=owners.get(global_id),
            )
        )
    for global_id in sorted(KITCHEN_SANITARY_IDS):
        product = model.by_guid(global_id)
        records.append(
            object_record(
                settings, product, "I-501", source_hash, bbox_registry, "fixed_sanitary_fixture",
                "kitchen_fixture_envelope_coordination", "Exact IFC GlobalId and assigned sanitary type preserve the existing fixture identity.", 1.0,
                equipment=owners.get(global_id),
            )
        )
    for global_id in sorted(KITCHEN_CONTEXT_IDS):
        product = model.by_guid(global_id)
        records.append(
            object_record(
                settings, product, "I-501", source_hash, bbox_registry, "fixed_context",
                "kitchen_joinery_context", "Previously confirmed kitchen pier or cabinet assembly identity; world geometry is read only.", 1.0,
            )
        )

    for sanitary in model.by_type("IfcSanitaryTerminal"):
        if sanitary.GlobalId in KITCHEN_SANITARY_IDS:
            continue
        review_status = (
            "semantic_correction_required_flush_plate_not_wcseat"
            if sanitary.GlobalId in GEBERIT_FLUSH_PLATE_IDS
            else "existing_model_candidate"
        )
        records.append(
            object_record(
                settings, sanitary, "I-502", source_hash, bbox_registry, "fixed_sanitary_fixture",
                "bathroom_fixture_envelope_coordination", "Existing assigned sanitary type is retained for grouping; it does not prove rough-in centres.", 0.8,
                equipment=owners.get(sanitary.GlobalId),
                review_status=review_status,
            )
        )
    for waste in model.by_type("IfcWasteTerminal"):
        records.append(
            object_record(
                settings, waste, "I-502", source_hash, bbox_registry, "fixed_waste_terminal",
                "bathroom_drain_envelope_coordination", "Existing waste-terminal identity and world geometry; connection centres remain unconfirmed.", 0.8,
                equipment=owners.get(waste.GlobalId),
            )
        )
    for global_id in sorted(BATHROOM_CONTEXT_IDS):
        product = model.by_guid(global_id)
        records.append(
            object_record(
                settings, product, "I-502", source_hash, bbox_registry, "fixed_context",
                "bathroom_wall_or_depressed_slab_context", "Previously confirmed bathroom wall or 50 mm depressed-slab context.", 1.0,
            )
        )

    washtower = model.by_guid(WASHTOWER_ID)
    records.append(
        object_record(
            settings, washtower, "I-503", source_hash, bbox_registry, "fixed_equipment",
            "laundry_equipment_envelope_coordination", "Exact IFC GlobalId, LG WashTower name and existing description preserve the current model identity.", 1.0,
            equipment=owners.get(WASHTOWER_ID),
        )
    )

    object_records = [record for record in records if record["record_kind"] == "existing_object"]
    ssot_scoped_records = [
        record for record in object_records
        if record["ifc_class"] in {
            "IfcFurniture", "IfcElectricAppliance", "IfcSanitaryTerminal", "IfcWasteTerminal"
        }
    ]
    missing_ssot = sorted(record["global_id"] for record in ssot_scoped_records if not record["equipment_id"])
    if missing_ssot:
        raise RuntimeError(f"INT1 scoped objects missing from equipment SSOT: {missing_ssot}")
    object_ids = [record["global_id"] for record in object_records]
    if len(object_ids) != len(set(object_ids)):
        duplicates = [key for key, count in Counter(object_ids).items() if count > 1]
        raise RuntimeError(f"duplicate INT1 object assignment: {duplicates}")
    if set(record["sheet_id"] for record in object_records) != EXPECTED_SHEETS:
        raise RuntimeError("INT1 candidate must cover I-501 through I-504")
    if furniture_sheet_counts != Counter({"I-501": 63, "I-502": 1, "I-504": 25}):
        raise RuntimeError(f"unexpected furniture sheet allocation: {furniture_sheet_counts}")

    csv_records = [*object_records, *(blocker_csv_record(row, source_hash) for row in BLOCKERS)]
    write_csv(args.decision_csv, csv_records)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sheet_counts = Counter(record["sheet_id"] for record in object_records)
    summary = {
        "furniture_total": sum(furniture_role_counts.values()),
        "fixed_furniture": furniture_role_counts["fixed_furniture"],
        "loose_furniture": furniture_role_counts["loose_furniture"],
        "furniture_by_sheet": dict(sorted(furniture_sheet_counts.items())),
        "existing_object_records": len(object_records),
        "records_by_sheet": dict(sorted(sheet_counts.items())),
        "sheet_ids": sorted(sheet_counts),
        "block_count": len(BLOCKERS),
        "geberit_flush_plate_semantic_corrections": sum(
            record["global_id"] in GEBERIT_FLUSH_PLATE_IDS for record in object_records
        ),
        "equipment_ssot_linked_object_records": sum(bool(record["equipment_id"]) for record in object_records),
        "equipment_ssot_unlinked_scoped_records": len(missing_ssot),
    }
    candidate_generation_pass = (
        summary["furniture_total"] == EXPECTED_FURNITURE
        and summary["fixed_furniture"] == EXPECTED_FIXED
        and summary["loose_furniture"] == EXPECTED_LOOSE
        and set(summary["sheet_ids"]) == EXPECTED_SHEETS
        and summary["block_count"] == 4
        and summary["geberit_flush_plate_semantic_corrections"] == 2
    )
    report = {
        "mode": "read_only_existing_object_candidate",
        "source": {
            "ifc": str(source_path), "ifc_sha256": source_hash,
            "equipment_ssot": {
                "equipment_register_sha256": sha256(root / "pipeline/decisions/equipment-register.csv"),
                "installation_requirements_sha256": sha256(root / "pipeline/decisions/equipment-installation-requirements.csv"),
                "source_evidence_sha256": sha256(root / "pipeline/decisions/source-evidence-register.csv"),
            },
        },
        "summary": summary,
        "gates": {
            "candidate_generation_pass": candidate_generation_pass,
            "ifc_write_allowed": False,
            "fabrication_dimensions_ready": False,
            "int1_completion_pass": False,
        },
        "blockers": BLOCKERS,
        "kitchen_product_installation_requirements": [
            kitchen_product_snapshot(
                equipment_id, equipment_by_id, requirements_by_id
            )
            for equipment_id in KITCHEN_PRODUCT_IDS
        ],
        "records": object_records,
    }
    if not candidate_generation_pass:
        raise RuntimeError("INT1 existing-object candidate gate failed")
    (args.output_dir / "int1-existing-report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    write_summary(args.output_dir / "summary.md", report)
    print(json.dumps({"summary": summary, "gates": report["gates"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
