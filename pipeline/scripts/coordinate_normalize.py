#!/usr/bin/env python3
"""Audit IFC coordinates and remove sub-0.01 mm floating-point noise.

Only IfcLengthMeasure values already within 0.01 mm of an integer are
normalized. Intentional decimals, directions, angles, ratios, and uncertain
placements remain unchanged for Blender review.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import ifcopenshell
import ifcopenshell.util.placement
import ifcopenshell.util.unit
FLOAT_NOISE_TOLERANCE_MM = 0.01
DEFAULT_ORIGIN_REVIEW_TOLERANCE_MM = 0.1
CARDINAL_TOLERANCE = 1e-6
FLOAT_NOISE_EPSILON_MM = 1e-9


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def residual_mm(value: float) -> float:
    return abs(value - round(value))


def is_integer_mm(value: float, tolerance_mm: float = FLOAT_NOISE_TOLERANCE_MM) -> bool:
    return residual_mm(value) <= tolerance_mm


def flatten_numbers(value: Any) -> Iterable[float]:
    if isinstance(value, bool):
        return
    if isinstance(value, (int, float)):
        yield float(value)
        return
    if isinstance(value, (tuple, list)):
        for item in value:
            yield from flatten_numbers(item)


def clean_noise_nested(value: Any) -> tuple[Any, int]:
    """Round only sub-0.01 mm tails, preserving intentional decimals."""
    if isinstance(value, bool):
        return value, 0
    if isinstance(value, (int, float)):
        number = float(value)
        residual = residual_mm(number)
        if FLOAT_NOISE_EPSILON_MM < residual <= FLOAT_NOISE_TOLERANCE_MM:
            return float(round(number)), 1
        return value, 0
    if isinstance(value, tuple):
        cleaned = [clean_noise_nested(item) for item in value]
        return tuple(item[0] for item in cleaned), sum(item[1] for item in cleaned)
    if isinstance(value, list):
        cleaned = [clean_noise_nested(item) for item in value]
        return [item[0] for item in cleaned], sum(item[1] for item in cleaned)
    return value, 0


def clean_floating_noise(ifc: ifcopenshell.file) -> dict[str, Any]:
    """Normalize near-integer IfcLengthMeasure values in project millimetres.

    Directions, angles, ratios, and intentional length decimals farther than
    0.01 mm from an integer are untouched. Explicitly unit-scoped select values
    are skipped because their raw number may not use the project's millimetres.
    """
    if abs(ifcopenshell.util.unit.calculate_unit_scale(ifc) - 0.001) > 1e-12:
        raise RuntimeError("Floating-noise cleanup currently requires an IFC model using millimetres.")

    changed_entities: set[int] = set()
    changed_scalars = 0
    by_attribute: Counter[str] = Counter()

    for entity in ifc:
        declaration = entity.wrapped_data.declaration().as_entity()
        for index in range(len(entity)):
            value = entity[index]
            if value is None:
                continue
            attribute = declaration.attribute_by_index(index)
            attribute_type = str(attribute.type_of_attribute())
            key = f"{entity.is_a()}.{attribute.name()}"

            if isinstance(value, ifcopenshell.entity_instance):
                is_length = value.id() == 0 and "LengthMeasure" in value.is_a()
                if not is_length:
                    continue
                if entity.is_a() in {"IfcMeasureWithUnit", "IfcPropertySingleValue"} and getattr(entity, "Unit", None):
                    continue
                cleaned, count = clean_noise_nested(value.wrappedValue)
                if count:
                    value.wrappedValue = cleaned
            else:
                if "IfcLengthMeasure" not in attribute_type:
                    continue
                cleaned, count = clean_noise_nested(value)
                if count:
                    entity[index] = cleaned

            if count:
                changed_entities.add(entity.id())
                changed_scalars += count
                by_attribute[key] += count

    return {
        "tolerance_mm": FLOAT_NOISE_TOLERANCE_MM,
        "changed_entities": len(changed_entities),
        "changed_scalars": changed_scalars,
        "by_attribute": dict(by_attribute.most_common()),
    }


def product_origin_records(
    ifc: ifcopenshell.file, origin_tolerance_mm: float
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for product in ifc.by_type("IfcProduct"):
        placement = getattr(product, "ObjectPlacement", None)
        if not placement or not placement.is_a("IfcLocalPlacement"):
            continue
        matrix = ifcopenshell.util.placement.get_local_placement(placement)
        xyz = [float(matrix[index][3]) for index in range(3)]
        target = [float(round(value)) for value in xyz]
        delta = [target[index] - xyz[index] for index in range(3)]
        records.append(
            {
                "step_id": product.id(),
                "global_id": getattr(product, "GlobalId", None),
                "ifc_class": product.is_a(),
                "name": getattr(product, "Name", None),
                "current_mm": xyz,
                "target_mm": target,
                "delta_mm": delta,
                "max_residual_mm": max(residual_mm(value) for value in xyz),
                "shift_mm": math.sqrt(sum(value * value for value in delta)),
                "within_review_tolerance": all(
                    is_integer_mm(value, origin_tolerance_mm) for value in xyz
                ),
            }
        )
    return records


def length_value_audit(ifc: ifcopenshell.file) -> dict[str, Any]:
    total = 0
    noninteger = 0
    by_attribute: Counter[str] = Counter()
    noninteger_by_attribute: Counter[str] = Counter()
    for entity in ifc:
        declaration = entity.wrapped_data.declaration().as_entity()
        for index in range(len(entity)):
            value = entity[index]
            if value is None:
                continue
            attribute = declaration.attribute_by_index(index)
            attribute_type = str(attribute.type_of_attribute())
            # Select-valued attributes list every permitted measure in their
            # schema text. Count their concrete wrapped type, not the entire
            # select declaration (for example, do not mistake an angle unit
            # conversion factor for a length).
            if isinstance(value, ifcopenshell.entity_instance):
                is_length = value.id() == 0 and "LengthMeasure" in value.is_a()
            else:
                is_length = "IfcLengthMeasure" in attribute_type
            if not is_length:
                continue
            key = f"{entity.is_a()}.{attribute.name()}"
            for number in flatten_numbers(value):
                total += 1
                by_attribute[key] += 1
                if not is_integer_mm(number):
                    noninteger += 1
                    noninteger_by_attribute[key] += 1
    return {
        "total": total,
        "integer": total - noninteger,
        "noninteger": noninteger,
        "top_noninteger_attributes": [
            {
                "attribute": key,
                "noninteger": count,
                "total": by_attribute[key],
            }
            for key, count in noninteger_by_attribute.most_common(20)
        ],
    }


def direction_audit(ifc: ifcopenshell.file) -> dict[str, int]:
    total = 0
    exact_cardinal = 0
    near_cardinal = 0
    intentional_noncardinal = 0
    for direction in ifc.by_type("IfcDirection"):
        values = tuple(float(value) for value in direction.DirectionRatios)
        target, distance = nearest_cardinal(values)
        total += 1
        if distance <= 1e-12:
            exact_cardinal += 1
        elif distance <= CARDINAL_TOLERANCE:
            near_cardinal += 1
        else:
            intentional_noncardinal += 1
    return {
        "total": total,
        "exact_cardinal": exact_cardinal,
        "near_cardinal": near_cardinal,
        "intentional_noncardinal": intentional_noncardinal,
    }


def nearest_cardinal(values: tuple[float, ...]) -> tuple[tuple[float, ...], float]:
    norm = math.sqrt(sum(value * value for value in values))
    if not norm:
        return values, math.inf
    unit = tuple(value / norm for value in values)
    candidates: list[tuple[float, ...]] = []
    for axis in range(len(unit)):
        for sign in (-1.0, 1.0):
            candidates.append(tuple(sign if index == axis else 0.0 for index in range(len(unit))))
    target = min(candidates, key=lambda item: math.sqrt(sum((unit[i] - item[i]) ** 2 for i in range(len(unit)))))
    distance = math.sqrt(sum((unit[index] - target[index]) ** 2 for index in range(len(unit))))
    return target, distance


def audit(
    ifc: ifcopenshell.file,
    source: Path,
    origin_tolerance_mm: float = DEFAULT_ORIGIN_REVIEW_TOLERANCE_MM,
) -> dict[str, Any]:
    unit_scale = ifcopenshell.util.unit.calculate_unit_scale(ifc)
    origins = product_origin_records(ifc, origin_tolerance_mm)
    by_class: Counter[str] = Counter(record["ifc_class"] for record in origins)
    noninteger_by_class: Counter[str] = Counter(
        record["ifc_class"]
        for record in origins
        if not record["within_review_tolerance"]
    )
    shifts = [record["shift_mm"] for record in origins]
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {"path": str(source), "sha256": sha256(source)},
        "ifc_schema": ifc.schema,
        "unit_scale_to_m": unit_scale,
        "float_noise_tolerance_mm": FLOAT_NOISE_TOLERANCE_MM,
        "origin_review_tolerance_mm": origin_tolerance_mm,
        "product_origins": {
            "total": len(origins),
            "within_review_tolerance": sum(
                record["within_review_tolerance"] for record in origins
            ),
            "over_review_tolerance": sum(
                not record["within_review_tolerance"] for record in origins
            ),
            "max_candidate_shift_mm": max(shifts, default=0.0),
            "by_class": dict(sorted(by_class.items())),
            "noninteger_by_class": dict(sorted(noninteger_by_class.items())),
            "records": origins,
        },
        "raw_length_values": length_value_audit(ifc),
        "directions": direction_audit(ifc),
    }


def entity_counts(ifc: ifcopenshell.file) -> Counter[str]:
    return Counter(entity.is_a() for entity in ifc)


def build_noise_candidate(source: Path, output: Path, report_path: Path) -> dict[str, Any]:
    ifc = ifcopenshell.open(source)
    before_counts = entity_counts(ifc)
    before_guids = sorted(
        entity.GlobalId for entity in ifc.by_type("IfcRoot") if getattr(entity, "GlobalId", None)
    )
    before_sha256 = sha256(source)
    changes = clean_floating_noise(ifc)
    output.parent.mkdir(parents=True, exist_ok=True)
    ifc.write(output)
    reopened = ifcopenshell.open(output)
    remaining = clean_floating_noise(reopened)
    after_guids = sorted(
        entity.GlobalId for entity in reopened.by_type("IfcRoot") if getattr(entity, "GlobalId", None)
    )
    qa = {
        "schema_unchanged": reopened.schema == ifc.schema,
        "entity_counts_unchanged": entity_counts(reopened) == before_counts,
        "global_ids_unchanged": after_guids == before_guids,
        "noise_cleanup_idempotent": remaining["changed_scalars"] == 0,
    }
    report = {
        "mode": "clean-noise",
        "source": {"path": str(source), "sha256": before_sha256},
        "candidate": {"path": str(output), "sha256": sha256(output)},
        "changes": changes,
        "qa": qa,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    if not all(qa.values()):
        raise RuntimeError(f"Floating-noise candidate failed QA. See {report_path}")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("audit", "clean-noise"))
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument(
        "--origin-tolerance-mm",
        type=float,
        default=DEFAULT_ORIGIN_REVIEW_TOLERANCE_MM,
        help="Review threshold for product origins in audit mode (default: 0.1 mm).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = args.input.resolve()
    if args.mode == "audit":
        report = audit(ifcopenshell.open(source), source, args.origin_tolerance_mm)
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    elif args.mode == "clean-noise":
        if args.output is None:
            raise SystemExit("--output is required in clean-noise mode")
        report = build_noise_candidate(source, args.output.resolve(), args.report.resolve())
    if "qa" in report:
        summary = report["qa"]
    else:
        summary = {
            key: value
            for key, value in report["product_origins"].items()
            if key != "records"
        }
    print(json.dumps({"report": str(args.report), "summary": summary}, ensure_ascii=False))


if __name__ == "__main__":
    main()
