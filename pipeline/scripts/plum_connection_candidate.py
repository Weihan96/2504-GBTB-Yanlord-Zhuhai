#!/usr/bin/env python3
"""Generate read-only sanitary-to-drainage proximity candidates."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom


EXPECTED_IFC_SHA256 = "7521c09991f3d0c7b7d91ca2324fd55ad961d8e32e9e3e9a9777a4cc19b06e81"
OFF_PLAN_WASTE_TERMINAL_ID = "3IVqCnhGr51hY4LrOq_5G_"


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument("--p201", type=Path, default=root / "build/plum/p201-demand-endpoints.json")
    parser.add_argument("--centerlines", type=Path, default=root / "build/coordinate-normalization/flow-segment-centerline-audit.json")
    parser.add_argument("--output", type=Path, default=root / "build/plum/plum-connection-candidate.json")
    parser.add_argument("--contact-mm", type=float, default=10.0)
    parser.add_argument("--near-mm", type=float, default=150.0)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def world_bbox_mm(settings: ifcopenshell.geom.settings, product: Any) -> tuple[list[float], list[float]]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    values = list(shape.geometry.verts)
    points = [
        [values[index] * 1000.0, values[index + 1] * 1000.0, values[index + 2] * 1000.0]
        for index in range(0, len(values), 3)
    ]
    return (
        [min(point[axis] for point in points) for axis in range(3)],
        [max(point[axis] for point in points) for axis in range(3)],
    )


def point_bbox_distance(point: list[float], minimum: list[float], maximum: list[float]) -> tuple[float, float]:
    offsets = [max(minimum[axis] - point[axis], 0.0, point[axis] - maximum[axis]) for axis in range(3)]
    return math.hypot(offsets[0], offsets[1]), math.sqrt(sum(value * value for value in offsets))


def branch_records(centerlines: dict[str, Any]) -> list[dict[str, Any]]:
    branches = []
    for product in centerlines["records"]:
        for index, component in enumerate(product["components"], 1):
            branches.append({
                "branch_id": f"{product['global_id']}:B{index:02d}",
                "product_global_id": product["global_id"],
                "product_name": product["name"],
                "classification": product["classification"],
                "component_index": index - 1,
                "end_centres_mm": component["derived_end_centres_mm"],
                "derived_centreline_is_write_authority": False,
            })
    return branches


def proximity_class(distance_mm: float, contact_mm: float, near_mm: float) -> tuple[str, float]:
    if distance_mm <= contact_mm:
        return "geometric_contact_candidate", 0.90
    if distance_mm <= near_mm:
        return "nearby_connection_candidate", 0.70
    return "remote_nearest_branch_only", 0.35


def main() -> int:
    args = parse_args()
    source_hash = sha256(args.ifc)
    if source_hash != EXPECTED_IFC_SHA256:
        raise RuntimeError(f"formal IFC hash changed: {source_hash}")
    p201 = read_json(args.p201)
    centerlines = read_json(args.centerlines)
    if p201["source_ifc_sha256"] != source_hash or centerlines["source"]["sha256"] != source_hash:
        raise RuntimeError("P-201 or centerline report is stale")

    model = ifcopenshell.open(args.ifc)
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    branches = branch_records(centerlines)
    records = []
    for endpoint in p201["demand_endpoints"]:
        if not endpoint["service_demand_candidate"]:
            continue
        product = model.by_guid(endpoint["global_id"])
        minimum, maximum = world_bbox_mm(settings, product)
        candidates = []
        for branch in branches:
            distances = [point_bbox_distance(point, minimum, maximum) for point in branch["end_centres_mm"]]
            plan_distance, distance = min(distances, key=lambda pair: pair[1])
            candidates.append({
                "branch_id": branch["branch_id"],
                "product_global_id": branch["product_global_id"],
                "minimum_plan_distance_to_fixture_bbox_mm": round(plan_distance, 6),
                "minimum_3d_distance_to_fixture_bbox_mm": round(distance, 6),
            })
        candidates.sort(key=lambda row: (row["minimum_3d_distance_to_fixture_bbox_mm"], row["branch_id"]))
        nearest = candidates[0]
        classification, confidence = proximity_class(
            nearest["minimum_3d_distance_to_fixture_bbox_mm"], args.contact_mm, args.near_mm
        )
        records.append({
            "global_id": endpoint["global_id"],
            "type_name": endpoint["type_name"],
            "type_predefined_type": endpoint["type_predefined_type"],
            "candidate_role": endpoint["candidate_role"],
            "fixture_bbox_min_mm": [round(value, 6) for value in minimum],
            "fixture_bbox_max_mm": [round(value, 6) for value in maximum],
            "nearest_drainage_branch": nearest,
            "proximity_classification": classification,
            "confidence": confidence,
            "basis": "minimum distance from derived drainage component end centres to the current fixture world bbox",
            "review_required": True,
            "automatic_ifc_write_allowed": False,
        })

    counts = Counter(row["proximity_classification"] for row in records)
    report = {
        "source_ifc_sha256": source_hash,
        "tolerances_mm": {"contact": args.contact_mm, "near": args.near_mm},
        "summary": {
            "registered_sanitary_terminals": len(p201["demand_endpoints"]),
            "service_demand_candidates": len(records),
            "non_service_components": len(p201["demand_endpoints"]) - len(records),
            "derived_drainage_components": len(branches),
            "controlled_pvc110_branches": sum(branch["classification"] == "controlled_disconnected_bundle" for branch in branches),
            "proximity_counts": dict(sorted(counts.items())),
            "off_plan_waste_terminal_id": OFF_PLAN_WASTE_TERMINAL_ID,
        },
        "branches": branches,
        "records": records,
        "gates": {
            "all_service_candidates_have_nearest_branch_evidence": len(records) == 24,
            "derived_points_are_formal_connectors": False,
            "automatic_ifc_write_allowed": False,
            "construction_release_ready": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
