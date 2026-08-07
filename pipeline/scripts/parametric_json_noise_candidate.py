#!/usr/bin/env python3
"""Remove sub-0.01 mm integer tails from Bonsai door/window JSON metadata."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell

from geometry_alignment_audit import geometry_difference_audit, sha256


TARGET_PSETS = {"BBIM_Door", "BBIM_Window"}


def snap_near_integer(
    value: Any,
    tolerance_mm: float,
    path: tuple[str, ...] = (),
) -> tuple[Any, list[dict[str, Any]]]:
    if isinstance(value, bool):
        return value, []
    if isinstance(value, dict):
        result = {}
        changes = []
        for key, item in value.items():
            normalized, item_changes = snap_near_integer(
                item, tolerance_mm, (*path, key)
            )
            result[key] = normalized
            changes.extend(item_changes)
        return result, changes
    if isinstance(value, list):
        result = []
        changes = []
        for index, item in enumerate(value):
            normalized, item_changes = snap_near_integer(
                item, tolerance_mm, (*path, str(index))
            )
            result.append(normalized)
            changes.extend(item_changes)
        return result, changes
    if isinstance(value, (int, float)):
        number = float(value)
        if not math.isfinite(number):
            raise RuntimeError(f"Non-finite JSON number at {'.'.join(path)}")
        target = round(number)
        delta = target - number
        if 0.0 < abs(delta) <= tolerance_mm:
            return target, [
                {
                    "path": ".".join(path),
                    "before": number,
                    "after": target,
                    "delta_mm": delta,
                }
            ]
    return value, []


def data_properties(model: ifcopenshell.file):
    for product_type in [
        *model.by_type("IfcDoorType"),
        *model.by_type("IfcWindowType"),
    ]:
        for pset in product_type.HasPropertySets or ():
            if not pset.is_a("IfcPropertySet") or pset.Name not in TARGET_PSETS:
                continue
            matches = [
                prop
                for prop in pset.HasProperties
                if prop.is_a("IfcPropertySingleValue") and prop.Name == "Data"
            ]
            if len(matches) != 1:
                raise RuntimeError(
                    f"{product_type.GlobalId} {pset.Name} must contain one Data property"
                )
            related_products = [
                product
                for relation in product_type.Types
                for product in relation.RelatedObjects
            ]
            yield product_type, related_products, pset, matches[0]


def apply_json_cleanup(
    model: ifcopenshell.file, tolerance_mm: float
) -> list[dict[str, Any]]:
    results = []
    for product_type, related_products, pset, prop in data_properties(model):
        raw = prop.NominalValue.wrappedValue
        parsed = json.loads(raw)
        normalized, changes = snap_near_integer(parsed, tolerance_mm)
        if not changes:
            continue
        prop.NominalValue = model.create_entity(
            "IfcText", json.dumps(normalized, ensure_ascii=False)
        )
        results.append(
            {
                "product_class": (
                    "IfcDoor"
                    if product_type.is_a("IfcDoorType")
                    else "IfcWindow"
                ),
                "type_class": product_type.is_a(),
                "type_global_id": product_type.GlobalId,
                "affected_product_global_ids": sorted(
                    product.GlobalId for product in related_products
                ),
                "pset": pset.Name,
                "property_id": prop.id(),
                "change_count": len(changes),
                "maximum_absolute_delta_mm": max(
                    abs(change["delta_mm"]) for change in changes
                ),
                "changes": changes,
            }
        )
    return results


def residual_count(model: ifcopenshell.file, tolerance_mm: float) -> int:
    result = 0
    for _, _, _, prop in data_properties(model):
        _, changes = snap_near_integer(
            json.loads(prop.NominalValue.wrappedValue), tolerance_mm
        )
        result += len(changes)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.01)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= 0.1:
        raise SystemExit("--tolerance-mm must be greater than 0 and at most 0.1")
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    results = apply_json_cleanup(candidate, args.tolerance_mm)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)
    difference = geometry_difference_audit(
        candidate,
        source,
        str(source_path),
        tolerance_mm=0.1,
        classes=["IfcDoor", "IfcWindow"],
        global_ids=[],
    )
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    residual = residual_count(candidate, args.tolerance_mm)
    affected_product_ids = {
        global_id
        for result in results
        for global_id in result["affected_product_global_ids"]
    }
    gates = {
        "changed_type_property_sets": len(results),
        "changed_products": len(affected_product_ids),
        "changed_values": sum(result["change_count"] for result in results),
        "residual_values": residual,
        "door_window_geometry_over_tolerance": difference["over_tolerance"],
        "schema_equal": source.schema == candidate.schema,
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
        "entity_count_equal": len(list(source)) == len(list(candidate)),
    }
    gates["pass"] = (
        gates["changed_products"] > 0
        and gates["changed_values"] > 0
        and residual == 0
        and difference["over_tolerance"] == 0
        and gates["schema_equal"]
        and gates["root_global_ids_equal"]
        and gates["entity_count_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-parametric-json-noise-candidate",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "candidate": {
            "path": str(args.output),
            "sha256": sha256(args.output),
            "schema": candidate.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "results": results,
        "geometry_difference": difference,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
