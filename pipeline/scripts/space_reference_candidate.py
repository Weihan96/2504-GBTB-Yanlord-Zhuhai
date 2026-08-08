#!/usr/bin/env python3
"""Generate a read-only R01-R22 Space reference candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
import ifcopenshell.geom
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--register", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bbox_centre_mm(settings: ifcopenshell.geom.settings, product: ifcopenshell.entity_instance) -> list[float]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3)) * 1000.0
    if len(vertices) == 0:
        raise RuntimeError(f"Space {product.GlobalId} has no Body geometry")
    return ((vertices.min(axis=0) + vertices.max(axis=0)) / 2.0).tolist()


def current_reference(space: ifcopenshell.entity_instance) -> str:
    for relationship in space.IsDefinedBy or ():
        if not relationship.is_a("IfcRelDefinesByProperties"):
            continue
        definition = relationship.RelatingPropertyDefinition
        if not definition.is_a("IfcPropertySet") or definition.Name != "Pset_SpaceCommon":
            continue
        for prop in definition.HasProperties or ():
            if prop.Name == "Reference" and prop.NominalValue is not None:
                return str(prop.NominalValue.wrappedValue or "")
    return ""


def main() -> None:
    args = parse_args()
    model = ifcopenshell.open(args.input)
    if model.schema != "IFC4":
        raise RuntimeError(f"expected IFC4, got {model.schema}")

    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    records = []
    for space in model.by_type("IfcSpace"):
        records.append(
            {
                "global_id": space.GlobalId,
                "long_name": str(space.LongName or space.Name or ""),
                "name": str(space.Name or ""),
                "current_reference": current_reference(space),
                "centre_mm": bbox_centre_mm(settings, space),
            }
        )
    if len(records) != 22:
        raise RuntimeError(f"expected 22 IfcSpace objects, got {len(records)}")
    plan_centre = np.mean([record["centre_mm"][:2] for record in records], axis=0)
    foyers = [record for record in records if record["long_name"] == "玄关"]
    if len(foyers) != 1:
        raise RuntimeError(f"expected exactly one 玄关 Space, got {len(foyers)}")
    foyer = foyers[0]
    start_angle = math.atan2(
        foyer["centre_mm"][1] - plan_centre[1],
        foyer["centre_mm"][0] - plan_centre[0],
    )
    for record in records:
        angle = math.atan2(
            record["centre_mm"][1] - plan_centre[1],
            record["centre_mm"][0] - plan_centre[0],
        )
        record["clockwise_angle_deg"] = math.degrees((start_angle - angle) % (2.0 * math.pi))
    records.sort(key=lambda record: (round(record["clockwise_angle_deg"], 9), record["global_id"]))
    for index, record in enumerate(records, 1):
        record["candidate_reference"] = f"R{index:02d}"

    existing_count = sum(bool(record["current_reference"]) for record in records)
    existing_references_match_candidate = all(
        record["current_reference"] == record["candidate_reference"] for record in records
    )
    if existing_count not in (0, 22):
        raise RuntimeError(f"formal IFC contains a partial Space Reference set: {existing_count}/22")
    if existing_count == 22 and not existing_references_match_candidate:
        raise RuntimeError("formal IFC Space References do not match the confirmed R01-R22 candidate")

    if records[0]["long_name"] != "玄关" or records[0]["candidate_reference"] != "R01":
        raise RuntimeError("R01 is not assigned to 玄关")
    references = [record["candidate_reference"] for record in records]
    if len(set(references)) != 22 or references != [f"R{i:02d}" for i in range(1, 23)]:
        raise RuntimeError("candidate References are not the complete unique R01-R22 sequence")

    args.register.parent.mkdir(parents=True, exist_ok=True)
    with args.register.open("w", encoding="utf-8-sig", newline="") as handle:
        fieldnames = [
            "candidate_reference",
            "space_global_id",
            "space_long_name",
            "space_name",
            "current_reference",
            "centre_x_mm",
            "centre_y_mm",
            "centre_z_mm",
            "clockwise_angle_deg",
            "basis",
            "confidence",
            "review_required",
            "status",
            "formal_ifc_write_allowed",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "candidate_reference": record["candidate_reference"],
                    "space_global_id": record["global_id"],
                    "space_long_name": record["long_name"],
                    "space_name": record["name"],
                    "current_reference": record["current_reference"],
                    "centre_x_mm": f'{record["centre_mm"][0]:.3f}',
                    "centre_y_mm": f'{record["centre_mm"][1]:.3f}',
                    "centre_z_mm": f'{record["centre_mm"][2]:.3f}',
                    "clockwise_angle_deg": f'{record["clockwise_angle_deg"]:.6f}',
                    "basis": "用户确认玄关为 R01；按 22 个 Space 世界包围盒中心相对总体中心的顺时针极角排序",
                    "confidence": "1.00",
                    "review_required": "no",
                    "status": "implemented" if existing_count == 22 else "confirmed_candidate",
                    "formal_ifc_write_allowed": "no",
                }
            )

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {"path": str(args.input), "sha256": sha256(args.input), "schema": model.schema},
        "automatic_ifc_write_allowed": False,
        "ordering_rule": "玄关 R01；Space bbox centre around plan centroid, clockwise",
        "space_count": len(records),
        "unique_reference_count": len(set(references)),
        "r01_long_name": records[0]["long_name"],
        "formal_ifc_reference_count": existing_count,
        "existing_references_match_candidate": existing_count == 0 or existing_references_match_candidate,
        "long_names_preserved": True,
        "records": records,
        "qa": {
            "space_count_22": len(records) == 22,
            "references_complete_unique": len(set(references)) == 22,
            "foyer_is_r01": records[0]["long_name"] == "玄关",
            "formal_ifc_unchanged": True,
            "existing_reference_state_valid": existing_count == 0 or existing_references_match_candidate,
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        f"Space Reference candidate: {len(records)} spaces, R01={records[0]['long_name']}, "
        f"formal references {existing_count}/22, IFC unchanged"
    )


if __name__ == "__main__":
    main()
