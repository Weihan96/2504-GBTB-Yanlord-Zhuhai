#!/usr/bin/env python3
"""Classify furniture origins without moving complex furniture geometry."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.util.placement
import ifcopenshell.util.unit


EXPLICIT_LOOSE_PREDEFINED_TYPES = {"BED", "CHAIR", "SOFA", "TABLE"}
ROLE_DECISION_FIELDS = {
    "selector_kind",
    "selector_value",
    "installation_role",
    "basis",
    "confidence",
    "human_review_required",
    "status",
}
PRODUCT_REGISTER_FIELDS = {
    "selector_kind",
    "selector_value",
    "manufacturer",
    "product_name",
    "product_variant",
    "product_category",
    "intended_use",
    "source_url",
    "identity_basis",
    "confidence",
    "human_review_required",
    "status",
}


def read_role_decisions(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if rows and "equipment_id" in rows[0]:
        source_path = path.with_name("source-evidence-register.csv")
        with source_path.open(newline="", encoding="utf-8-sig") as handle:
            rows = [
                json.loads(row["legacy_projection_json"])
                for row in csv.DictReader(handle)
                if "furniture-installation-role.csv" in row["legacy_targets"].split(";")
            ]
    if not rows or set(rows[0]) != ROLE_DECISION_FIELDS:
        raise RuntimeError(f"invalid furniture role decision schema: {path}")
    seen: set[tuple[str, str]] = set()
    result: list[dict[str, Any]] = []
    for row in rows:
        selector = (row["selector_kind"], row["selector_value"])
        if row["selector_kind"] not in {"global_id", "type_name"}:
            raise RuntimeError(f"invalid furniture selector: {selector}")
        if selector in seen:
            raise RuntimeError(f"duplicate furniture selector: {selector}")
        seen.add(selector)
        if row["installation_role"] not in {"fixed_furniture", "loose_furniture"}:
            raise RuntimeError(f"invalid furniture role: {row}")
        if row["human_review_required"] not in {"yes", "no"}:
            raise RuntimeError(f"invalid furniture review flag: {row}")
        confidence = float(row["confidence"])
        if not 0.0 <= confidence <= 1.0:
            raise RuntimeError(f"invalid furniture confidence: {row}")
        result.append({**row, "confidence": confidence})
    return result


def read_product_register(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if rows and "equipment_id" in rows[0]:
        source_path = path.with_name("source-evidence-register.csv")
        with source_path.open(newline="", encoding="utf-8-sig") as handle:
            source_by_id = {row["source_id"]: row for row in csv.DictReader(handle)}
        rows = [
            {
                "selector_kind": row["selector_kind"], "selector_value": row["selector_value"],
                "manufacturer": row["manufacturer"], "product_name": row["item_name"],
                "product_variant": row["variant"], "product_category": row["category"],
                "intended_use": row["use_location_confirmed"],
                "source_url": source_by_id.get(row["source_ids"], {}).get("source_url", ""),
                "identity_basis": row["identity_basis"], "confidence": row["confidence"],
                "human_review_required": row["human_review_required"], "status": "confirmed",
            }
            for row in rows if row["legacy_kind"] == "furniture_product"
        ]
    if not rows or set(rows[0]) != PRODUCT_REGISTER_FIELDS:
        raise RuntimeError(f"invalid furniture product register schema: {path}")
    seen: set[tuple[str, str]] = set()
    result: list[dict[str, Any]] = []
    for row in rows:
        selector = (row["selector_kind"], row["selector_value"])
        if row["selector_kind"] not in {"global_id", "type_name"}:
            raise RuntimeError(f"invalid furniture product selector: {selector}")
        if selector in seen:
            raise RuntimeError(f"duplicate furniture product selector: {selector}")
        seen.add(selector)
        if row["human_review_required"] not in {"yes", "no"}:
            raise RuntimeError(f"invalid furniture product review flag: {row}")
        confidence = float(row["confidence"])
        if not 0.0 <= confidence <= 1.0:
            raise RuntimeError(f"invalid furniture product confidence: {row}")
        result.append({**row, "confidence": confidence})
    return result


def find_role_decision(
    global_id: str,
    type_name: str | None,
    decisions: list[dict[str, Any]],
) -> dict[str, Any] | None:
    global_matches = [
        row
        for row in decisions
        if row["selector_kind"] == "global_id"
        and row["selector_value"] == global_id
    ]
    if global_matches:
        return global_matches[0]
    type_matches = [
        row
        for row in decisions
        if row["selector_kind"] == "type_name"
        and row["selector_value"] == (type_name or "")
    ]
    return type_matches[0] if type_matches else None


def find_selector_record(
    global_id: str,
    type_name: str | None,
    records: list[dict[str, Any]],
) -> dict[str, Any] | None:
    global_matches = [
        row
        for row in records
        if row["selector_kind"] == "global_id"
        and row["selector_value"] == global_id
    ]
    if global_matches:
        return global_matches[0]
    type_matches = [
        row
        for row in records
        if row["selector_kind"] == "type_name"
        and row["selector_value"] == (type_name or "")
    ]
    return type_matches[0] if type_matches else None


def product_identity(record: dict[str, Any] | None) -> dict[str, Any] | None:
    if record is None or record["status"] != "confirmed":
        return None
    return {
        key: record[key]
        for key in (
            "manufacturer",
            "product_name",
            "product_variant",
            "product_category",
            "intended_use",
            "source_url",
            "identity_basis",
            "confidence",
            "human_review_required",
        )
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify_furniture(
    predefined_type: str | None,
    decision: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if decision is not None and decision["status"] == "implemented":
        shared = {
            "basis": decision["basis"],
            "confidence": decision["confidence"],
            "human_review_required": decision["human_review_required"] == "yes",
            "role_decision_selector": {
                "kind": decision["selector_kind"],
                "value": decision["selector_value"],
            },
            "review_required": decision["human_review_required"] == "yes",
            "automatic_write_allowed": False,
        }
        if decision["installation_role"] == "loose_furniture":
            return {
                "installation_role": "loose_furniture",
                "normalization_disposition": "controlled_exception_keep_complex_origin",
                **shared,
            }
        return {
            "installation_role": "fixed_furniture",
            "normalization_disposition": "fixed_installation_anchor_required",
            **shared,
        }
    normalized = (predefined_type or "").upper()
    if normalized in EXPLICIT_LOOSE_PREDEFINED_TYPES:
        return {
            "installation_role": "loose_furniture",
            "normalization_disposition": "controlled_exception_keep_complex_origin",
            "basis": (
                "assigned IfcFurnitureType.PredefinedType explicitly identifies "
                f"{normalized}; no fixed installation anchor is required"
            ),
            "review_required": False,
            "automatic_write_allowed": False,
            "confidence": 1.0,
            "human_review_required": False,
            "role_decision_selector": None,
        }
    return {
        "installation_role": "fixed_or_loose_requires_review",
        "normalization_disposition": "installation_anchor_review",
        "basis": (
            "USERDEFINED, SHELF or missing furniture type does not prove whether "
            "the object is fixed joinery or loose furniture"
        ),
        "review_required": True,
        "automatic_write_allowed": False,
        "confidence": 0.5,
        "human_review_required": True,
        "role_decision_selector": None,
    }


def origin_mm(model: ifcopenshell.file, product: Any) -> list[float]:
    matrix = ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement)
    factor = ifcopenshell.util.unit.calculate_unit_scale(model) * 1000.0
    return [float(matrix[index, 3]) * factor for index in range(3)]


def nearest_integer_delta(values: list[float]) -> list[float]:
    return [round(value) - value for value in values]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--role-decisions", type=Path)
    parser.add_argument("--product-register", type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= 1.0:
        raise SystemExit("--tolerance-mm must be greater than zero and at most 1 mm")
    source_path = args.input.resolve()
    model = ifcopenshell.open(source_path)
    role_decisions = read_role_decisions(args.role_decisions)
    product_records = read_product_register(args.product_register)
    records: list[dict[str, Any]] = []
    for furniture in model.by_type("IfcFurniture"):
        assigned_type = ifcopenshell.util.element.get_type(furniture)
        predefined_type = getattr(assigned_type, "PredefinedType", None)
        type_name = getattr(assigned_type, "Name", None)
        role_decision = find_role_decision(
            furniture.GlobalId,
            type_name,
            role_decisions,
        )
        product_record = find_selector_record(
            furniture.GlobalId,
            type_name,
            product_records,
        )
        classification = classify_furniture(predefined_type, role_decision)
        current = origin_mm(model, furniture)
        delta = nearest_integer_delta(current)
        container = ifcopenshell.util.element.get_container(furniture)
        records.append(
            {
                "global_id": furniture.GlobalId,
                "name": furniture.Name,
                "type_global_id": getattr(assigned_type, "GlobalId", None),
                "type_name": type_name,
                "type_predefined_type": predefined_type,
                "container": getattr(container, "Name", None),
                "origin_mm": current,
                "nearest_integer_delta_mm": delta,
                "nearest_integer_shift_mm": math.sqrt(sum(value * value for value in delta)),
                "origin_within_tolerance": max(abs(value) for value in delta)
                <= args.tolerance_mm,
                "product_identity": product_identity(product_record),
                **classification,
            }
        )
    records.sort(key=lambda record: record["global_id"])
    summary = {
        "furniture": len(records),
        "explicit_loose_controlled_exceptions": sum(
            record["type_predefined_type"] in EXPLICIT_LOOSE_PREDEFINED_TYPES
            for record in records
        ),
        "loose_controlled_exceptions": sum(
            record["normalization_disposition"]
            == "controlled_exception_keep_complex_origin"
            for record in records
        ),
        "fixed_installation_anchor_required": sum(
            record["normalization_disposition"]
            == "fixed_installation_anchor_required"
            for record in records
        ),
        "installation_role_review": sum(
            record["review_required"] for record in records
        ),
        "confirmed_product_identities": sum(
            record["product_identity"] is not None for record in records
        ),
        "product_identity_review": sum(
            bool(record["product_identity"])
            and record["product_identity"]["human_review_required"] == "yes"
            for record in records
        ),
        "automatic_write_allowed": sum(
            record["automatic_write_allowed"] for record in records
        ),
        "by_predefined_type": dict(
            sorted(
                Counter(
                    record["type_predefined_type"] or "MISSING" for record in records
                ).items()
            )
        ),
    }
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-furniture-anchor-classification",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": model.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "summary": summary,
        "records": records,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **summary}, ensure_ascii=False))


if __name__ == "__main__":
    main()
