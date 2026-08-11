#!/usr/bin/env python3
"""Generate read-only kitchen electrical role candidates from current IFC and source drawings."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


SOCKET_ROLE_CANDIDATES = {
    "P001": ("c_face_purifier_or_spare_socket_group", 0.60, "one of three adjacent low sockets matching page 7 purifier plus two spare outlets"),
    "P002": ("c_face_purifier_or_spare_socket_group", 0.60, "one of three adjacent low sockets matching page 7 purifier plus two spare outlets"),
    "P003": ("c_face_purifier_or_spare_socket_group", 0.60, "one of three adjacent low sockets matching page 7 purifier plus two spare outlets"),
    "P004": ("refrigerator_socket_candidate", 0.85, "upper B-face point beside the first W570D660H1780 refrigerator volume"),
    "P005": ("direct_drinking_machine_socket_candidate", 0.75, "B-face point adjacent to WD01 and within the refrigerator/direct-drinking appliance bank"),
    "P006": ("oven_socket_candidate", 0.70, "B-face point in the oven/direct-drinking appliance bank; individual assignment needs elevation review"),
    "P007": ("refrigerator_socket_candidate", 0.85, "upper B-face point beside the second W570D660H1780 refrigerator volume"),
    "P008": ("existing_or_spare_service_socket_candidate", 0.80, "low B-face point below the cooking zone; the arrived ER9EPA33MP/01 hob is battery-ignited with 0 W mains load, so this point is not assigned to that product"),
    "P009": ("hood_socket_candidate", 0.95, "high B-face point 215 mm from HD01; latest page 6 explicitly circles the hood outlet"),
    "P010": ("h1200_wall_or_track_power_candidate", 0.65, "one of two adjacent H1173 points matching page 9 H1200 wall/track power intent"),
    "P011": ("h1200_wall_or_track_power_candidate", 0.65, "one of two adjacent H1173 points matching page 9 H1200 wall/track power intent"),
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument("--elec", type=Path, default=root / "build/elec/elec-existing-candidate.json")
    parser.add_argument("--source-audit", type=Path, default=root / "build/mep-positioning/source-audit.json")
    parser.add_argument("--output", type=Path, default=root / "build/elec/elec-positioning-candidate.json")
    parser.add_argument("--expected-ifc-sha256", help="Optional caller-frozen source hash")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def point_bbox_distance(point: list[float], bbox: dict[str, list[float]]) -> float:
    offsets = [
        max(bbox["min_mm"][axis] - point[axis], 0.0, point[axis] - bbox["max_mm"][axis])
        for axis in range(3)
    ]
    return math.sqrt(sum(value * value for value in offsets))


def proxy_identity(record: dict[str, Any]) -> tuple[str, float, str]:
    name = record["name"]
    if record["global_id"] in {"1SNXhKeZb7r9Y0MmzcQnc3", "2d2Vw3ZSn0exH1seMBMiVf"}:
        return (
            "island_dishwasher",
            1.00,
            "user confirmed the two SJ45ZB24MC proxies are two separate island dishwashers; power, water, drainage and service interfaces remain pending",
        )
    if name == "W570D660H1780":
        return "refrigerator_volume_candidate", 0.90, "two matching equipment volumes align with the two refrigerators and two refrigerator outlets in the source drawing"
    if name == "SJ45ZB24MC":
        return "dishwasher_model_candidate", 0.75, "unconfirmed SJ45ZB24MC proxy; exact appliance identity and interface remain unverified"
    if name == "LECOASE F50":
        return "food_waste_disposer", 0.95, "existing ObjectType is FOOD_WASTE_DISPOSER"
    if name.startswith("Hole W"):
        return "appliance_service_hole_candidate", 0.80, "existing name and low horizontal envelope identify a service-hole candidate without formal IfcOpeningElement semantics"
    if name == "Electric flue check valve":
        return "electric_flue_check_valve", 0.95, "existing name identifies the service component; power and maintenance interface remain unknown"
    return "unresolved_equipment_proxy", 0.30, "no supported identity rule"


def main() -> int:
    args = parse_args()
    source_hash = sha256(args.ifc)
    if args.expected_ifc_sha256 and source_hash != args.expected_ifc_sha256:
        raise RuntimeError(
            f"formal IFC hash changed: expected {args.expected_ifc_sha256}, got {source_hash}"
        )
    elec = read_json(args.elec)
    source = read_json(args.source_audit)
    if elec["source"]["sha256"] != source_hash or source["source_ifc_sha256"] != source_hash:
        raise RuntimeError("ELEC or MEP source report is stale")

    e303 = elec["sheets"]["E-303"]
    equipment = e303["typed_equipment"] + e303["proxy_handoffs"]
    socket_records = []
    for socket in e303["sockets"]:
        role, confidence, basis = SOCKET_ROLE_CANDIDATES[socket["candidate_id"]]
        centre = socket["bbox"]["centre_mm"]
        nearest = sorted(
            (
                {
                    "candidate_id": item["candidate_id"],
                    "global_id": item["global_id"],
                    "name": item["assigned_type"]["name"] or item["name"],
                    "distance_to_bbox_mm": round(point_bbox_distance(centre, item["bbox"]), 6),
                }
                for item in equipment
                if not (item["assigned_type"]["name"] or "").startswith("AC")
            ),
            key=lambda row: (row["distance_to_bbox_mm"], row["global_id"]),
        )
        socket_records.append({
            "candidate_id": socket["candidate_id"],
            "global_id": socket["global_id"],
            "space": socket["candidate_space"]["long_name"],
            "centre_mm": centre,
            "candidate_role": role,
            "nearest_non_hvac_equipment": nearest[0],
            "basis": basis,
            "confidence": confidence,
            "review_required": True,
            "automatic_ifc_write_allowed": False,
        })

    proxy_records = []
    for proxy in e303["proxy_handoffs"]:
        role, confidence, basis = proxy_identity(proxy)
        proxy_records.append({
            "candidate_id": proxy["candidate_id"],
            "global_id": proxy["global_id"],
            "current_name": proxy["name"],
            "space": proxy["candidate_space"]["long_name"],
            "candidate_role": role,
            "basis": basis,
            "confidence": confidence,
            "review_required": True,
            "automatic_ifc_write_allowed": False,
        })

    role_counts = Counter(row["candidate_role"] for row in socket_records)
    proxy_role_counts = Counter(row["candidate_role"] for row in proxy_records)
    report = {
        "source_ifc_sha256": source_hash,
        "source_pdf_sha256": source["source_drawings"]["latest"]["sha256"],
        "summary": {
            "existing_lights": len(elec["sheets"]["E-301"]["lights"]),
            "existing_sockets": len(socket_records),
            "typed_equipment": len(e303["typed_equipment"]),
            "proxy_handoffs": len(proxy_records),
            "socket_role_counts": dict(sorted(role_counts.items())),
            "proxy_role_counts": dict(sorted(proxy_role_counts.items())),
            "switch_instances": len(elec["sheets"]["E-302"]["instances"]),
            "network_instances": len(elec["sheets"]["E-304"]["instances"]),
        },
        "socket_candidates": socket_records,
        "proxy_identity_candidates": proxy_records,
        "gates": {
            "all_existing_kitchen_sockets_have_source_role_candidates": len(socket_records) == 11,
            "all_proxy_handoffs_have_identity_candidates": all(row["candidate_role"] != "unresolved_equipment_proxy" for row in proxy_records),
            "whole_home_socket_positioning_complete": False,
            "switch_positioning_complete": False,
            "network_positioning_complete": False,
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
