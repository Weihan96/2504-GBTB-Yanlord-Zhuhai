#!/usr/bin/env python3
"""Read-only audit of the A-105 50 mm finish zone and depressed slabs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import numpy as np


REFERENCE_FLOORS = {
    "F11": {"global_id": "1ogoq1VJP4vgBTqBCMFyCn", "confirmed_material": "地板"},
    "F20": {"global_id": "11c$NwzxL9hASgCgwDxM$w", "confirmed_material": "银白洞石岩板"},
    "F21": {"global_id": "3BwA5Rnkf4Bf3v3vTmeWTw", "confirmed_material": "地板"},
}
DEPRESSED_SLABS = {
    "3ARl_CqPrCWQrA6_$W073W": "1Ro3zU6lXBjuo7VDM3mTa7",
    "1UixnpnQb97wuNAGrTmMJd": "2F667JcWnAm9Hchuxrf8vg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--floor-register", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--expected-build-up-mm", type=float, default=50.0)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def world_bbox(settings: ifcopenshell.geom.settings, product: Any) -> dict[str, list[float]]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3)) * 1000.0
    minimum = vertices.min(axis=0)
    maximum = vertices.max(axis=0)
    return {
        "min_mm": minimum.tolist(),
        "max_mm": maximum.tolist(),
        "dimensions_mm": (maximum - minimum).tolist(),
    }


def bbox_overlap_area_mm2(first: dict[str, list[float]], second: dict[str, list[float]]) -> float:
    width = max(0.0, min(first["max_mm"][0], second["max_mm"][0]) - max(first["min_mm"][0], second["min_mm"][0]))
    height = max(0.0, min(first["max_mm"][1], second["max_mm"][1]) - max(first["min_mm"][1], second["min_mm"][1]))
    return width * height


def opening_inventory(model: ifcopenshell.file, settings: ifcopenshell.geom.settings) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for slab in model.by_type("IfcSlab"):
        if not slab.HasOpenings:
            continue
        slab_bbox = world_bbox(settings, slab)
        for relation in slab.HasOpenings:
            opening = relation.RelatedOpeningElement
            opening_bbox = world_bbox(settings, opening)
            records.append(
                {
                    "slab_global_id": slab.GlobalId,
                    "slab_name": str(slab.Name or ""),
                    "slab_bbox": slab_bbox,
                    "opening_global_id": opening.GlobalId,
                    "opening_name": str(opening.Name or ""),
                    "opening_bbox": opening_bbox,
                    "recess_below_slab_top_mm": slab_bbox["max_mm"][2] - opening_bbox["min_mm"][2],
                }
            )
    return records


def read_reference_floors(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = {row["candidate_id"]: row for row in csv.DictReader(handle)}
    records: list[dict[str, Any]] = []
    for candidate_id, expected in REFERENCE_FLOORS.items():
        row = rows.get(candidate_id)
        if row is None or row["global_id"] != expected["global_id"]:
            raise RuntimeError(f"{candidate_id} identity drift in A-105 floor register")
        records.append(
            {
                "candidate_id": candidate_id,
                "global_id": row["global_id"],
                "primary_space": row["primary_space"],
                "confirmed_material": expected["confirmed_material"],
                "bbox": json.loads(row["bbox"]),
            }
        )
    return records


def audit(args: argparse.Namespace) -> dict[str, Any]:
    model = ifcopenshell.open(args.input)
    if model.schema != "IFC4":
        raise RuntimeError(f"expected IFC4, found {model.schema}")
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    openings = opening_inventory(model, settings)
    reference_floors = read_reference_floors(args.floor_register)

    for floor in reference_floors:
        finish_z = float(floor["bbox"]["min_mm"][2])
        evidence = []
        for opening in openings:
            overlap = bbox_overlap_area_mm2(floor["bbox"], opening["opening_bbox"])
            if overlap <= 0.0:
                continue
            build_up = finish_z - float(opening["opening_bbox"]["min_mm"][2])
            evidence.append(
                {
                    "slab_global_id": opening["slab_global_id"],
                    "opening_global_id": opening["opening_global_id"],
                    "xy_bbox_overlap_area_mm2": overlap,
                    "opening_min_z_mm": opening["opening_bbox"]["min_mm"][2],
                    "finish_zone_depth_mm": build_up,
                    "matches_expected_build_up": abs(build_up - args.expected_build_up_mm) <= args.tolerance_mm,
                }
            )
        floor["finish_reference_z_mm"] = finish_z
        floor["opening_evidence"] = evidence
        floor["matching_opening_count"] = sum(item["matches_expected_build_up"] for item in evidence)
        floor["passes"] = (
            abs(finish_z) <= args.tolerance_mm
            and floor["matching_opening_count"] >= 1
        )

    depressed = []
    opening_by_pair = {(row["slab_global_id"], row["opening_global_id"]): row for row in openings}
    for slab_id, opening_id in DEPRESSED_SLABS.items():
        evidence = opening_by_pair.get((slab_id, opening_id))
        if evidence is None:
            raise RuntimeError(f"missing depressed slab/opening pair {slab_id}/{opening_id}")
        depth = float(evidence["recess_below_slab_top_mm"])
        depressed.append(
            {
                "slab_global_id": slab_id,
                "opening_global_id": opening_id,
                "slab_top_z_mm": evidence["slab_bbox"]["max_mm"][2],
                "opening_min_z_mm": evidence["opening_bbox"]["min_mm"][2],
                "recess_below_slab_top_mm": depth,
                "passes": abs(depth - args.expected_build_up_mm) <= args.tolerance_mm,
            }
        )

    gates = {
        "ifc_schema_is_ifc4": model.schema == "IFC4",
        "reference_floor_count": len(reference_floors),
        "reference_floor_at_shared_ffl_count": sum(abs(float(row["finish_reference_z_mm"])) <= args.tolerance_mm for row in reference_floors),
        "reference_floor_with_50mm_opening_evidence_count": sum(row["matching_opening_count"] >= 1 for row in reference_floors),
        "depressed_slab_count": len(depressed),
        "depressed_slab_50mm_recess_count": sum(row["passes"] for row in depressed),
        "automatic_ifc_write_allowed": False,
    }
    gates["mechanical_pass"] = (
        gates["reference_floor_count"] == 3
        and gates["reference_floor_at_shared_ffl_count"] == 3
        and gates["reference_floor_with_50mm_opening_evidence_count"] == 3
        and gates["depressed_slab_count"] == 2
        and gates["depressed_slab_50mm_recess_count"] == 2
    )
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a105-floor-build-up-audit",
        "source": {"ifc": str(args.input.resolve()), "ifc_sha256": sha256(args.input), "schema": model.schema},
        "expected_build_up_mm": args.expected_build_up_mm,
        "tolerance_mm": args.tolerance_mm,
        "method": "World-space AABB overlap identifies evidence only; finish depth is shared FFL Z minus overlapping slab Opening minimum Z. No IFC geometry or semantics are written.",
        "reference_floors": reference_floors,
        "depressed_slabs": depressed,
        "gates": gates,
    }


def main() -> None:
    args = parse_args()
    report = audit(args)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(args.report), **report["gates"]}, ensure_ascii=False))
    if not report["gates"]["mechanical_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
