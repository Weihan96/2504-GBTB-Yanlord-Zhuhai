#!/usr/bin/env python3
"""Compare developer handover MEP points with the current renovation model."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import ifcopenshell
from ifcopenshell.util.element import get_psets
from shapely.geometry import Point, box


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ifc", type=Path, default=root / "2504 GBTB Yanlord Zhuhai.ifc")
    parser.add_argument(
        "--developer",
        type=Path,
        default=root / "build/mep-positioning/developer-handover-mep-candidate.json",
    )
    parser.add_argument("--elec", type=Path, default=root / "build/elec/elec-existing-candidate.json")
    parser.add_argument("--plum", type=Path, default=root / "build/plum/p202-existing-location-register.json")
    parser.add_argument(
        "--walls",
        type=Path,
        default=root / "pipeline/decisions/a103-wall-status-review.csv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "build/mep-positioning/mep-renovation-delta-candidate.json",
    )
    parser.add_argument("--wall-near-mm", type=float, default=200.0)
    parser.add_argument("--current-point-near-mm", type=float, default=150.0)
    parser.add_argument("--plum-coordination-mm", type=float, default=600.0)
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


def read_walls(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    records = []
    for row in rows:
        records.append(
            {
                "global_id": row["global_id"],
                "status": row["current_status"],
                "candidate_status": row["candidate_status"],
                "footprint": box(
                    float(row["min_x_mm"]),
                    float(row["min_y_mm"]),
                    float(row["max_x_mm"]),
                    float(row["max_y_mm"]),
                ),
            }
        )
    return records


def nearest_wall(position: list[float], walls: list[dict[str, Any]]) -> dict[str, Any]:
    point = Point(position)
    distance, wall = min(
        ((record["footprint"].distance(point), record) for record in walls),
        key=lambda item: (item[0], item[1]["global_id"]),
    )
    return {
        "global_id": wall["global_id"],
        "current_status": wall["status"],
        "candidate_status": wall["candidate_status"],
        "plan_distance_mm": round(distance, 6),
        "basis": "minimum 2D distance to the current confirmed wall footprint bbox",
    }


def nearest_record(position: list[float], records: list[dict[str, Any]], point_getter) -> dict[str, Any]:
    distance, record = min(
        ((math.dist(position, point_getter(row)), row) for row in records),
        key=lambda item: (item[0], item[1]["global_id"]),
    )
    return {"distance_mm": round(distance, 6), "record": record}


def current_spaces(model: ifcopenshell.file) -> list[dict[str, str]]:
    records = []
    for space in model.by_type("IfcSpace"):
        reference = str(get_psets(space).get("Pset_SpaceCommon", {}).get("Reference", ""))
        records.append({"reference": reference, "global_id": space.GlobalId, "long_name": str(space.LongName or "")})
    return sorted(records, key=lambda row: row["reference"])


def main() -> int:
    args = parse_args()
    ifc_hash = sha256(args.ifc)
    if args.expected_ifc_sha256 and ifc_hash != args.expected_ifc_sha256:
        raise RuntimeError(
            f"formal IFC hash changed: expected {args.expected_ifc_sha256}, got {ifc_hash}"
        )
    developer = read_json(args.developer)
    elec = read_json(args.elec)
    plum = read_json(args.plum)
    if {developer["source_ifc_sha256"], elec["source"]["sha256"], plum["source_ifc_sha256"]} != {ifc_hash}:
        raise RuntimeError("MEP input report is stale")

    walls = read_walls(args.walls)
    sockets = elec["sheets"]["E-303"]["sockets"]
    plum_objects = plum["objects"]
    electrical = []
    plumbing = []
    for source in developer["records"]:
        if source["kind"] in {"developer_electrical_point", "developer_switch_or_control"}:
            wall = nearest_wall(source["ifc_position_mm"], walls)
            nearest_socket = None
            if source["kind"] == "developer_electrical_point":
                match = nearest_record(
                    source["ifc_position_mm"],
                    sockets,
                    lambda row: row["bbox"]["centre_mm"][:2],
                )
                nearest_socket = {
                    "global_id": match["record"]["global_id"],
                    "candidate_id": match["record"]["candidate_id"],
                    "distance_mm": match["distance_mm"],
                }
            if nearest_socket and nearest_socket["distance_mm"] <= args.current_point_near_mm:
                disposition = "current_socket_near_handover_reference"
                confidence = 0.80
            elif wall["current_status"] == "EXISTING" and wall["plan_distance_mm"] <= args.wall_near_mm:
                disposition = "existing_wall_reference_candidate"
                confidence = 0.65
            else:
                disposition = "relocation_or_replacement_review"
                confidence = 0.85
            electrical.append(
                {
                    "candidate_id": source["candidate_id"],
                    "kind": source["kind"],
                    "candidate_role": source["candidate_role"],
                    "ifc_position_mm": source["ifc_position_mm"],
                    "installation_height_mm": source.get("installation_height_mm"),
                    "candidate_space": source["candidate_space"],
                    "nearest_current_wall": wall,
                    "nearest_current_socket": nearest_socket,
                    "disposition_candidate": disposition,
                    "disposition_confidence": confidence,
                    "review_required": True,
                    "automatic_ifc_write_allowed": False,
                }
            )
        else:
            match = nearest_record(
                source["ifc_position_mm"],
                plum_objects,
                lambda row: row["geometry"]["location_candidate_mm"][:2],
            )
            if match["distance_mm"] <= args.current_point_near_mm:
                disposition = "current_plum_geometry_near_handover_reference"
                confidence = 0.75
            elif match["distance_mm"] <= args.plum_coordination_mm:
                disposition = "plum_coordination_zone_candidate"
                confidence = 0.55
            else:
                disposition = "source_only_or_relocated_plum_reference"
                confidence = 0.80
            plumbing.append(
                {
                    "candidate_id": source["candidate_id"],
                    "kind": source["kind"],
                    "candidate_role": source["candidate_role"],
                    "ifc_position_mm": source["ifc_position_mm"],
                    "candidate_space": source["candidate_space"],
                    "nearest_current_plum_object": {
                        "global_id": match["record"]["global_id"],
                        "ifc_class": match["record"]["ifc_class"],
                        "type_or_name": match["record"].get("type_name") or match["record"].get("name"),
                        "distance_mm": match["distance_mm"],
                        "is_formal_connector": False,
                    },
                    "disposition_candidate": disposition,
                    "disposition_confidence": confidence,
                    "review_required": True,
                    "automatic_ifc_write_allowed": False,
                }
            )

    model = ifcopenshell.open(args.ifc)
    spaces = current_spaces(model)
    current_light_counts = Counter(row["candidate_space"]["long_name"] for row in elec["sheets"]["E-301"]["lights"])
    current_socket_counts = Counter(row["candidate_space"]["long_name"] for row in sockets)
    room_program = []
    for space in spaces:
        source_electrical = [row for row in electrical if row["candidate_space"]["candidate_reference"] == space["reference"]]
        source_plumbing = [row for row in plumbing if row["candidate_space"]["candidate_reference"] == space["reference"]]
        room_program.append(
            {
                **space,
                "current_light_instances": current_light_counts[space["long_name"]],
                "current_socket_instances": current_socket_counts[space["long_name"]],
                "developer_power_or_other_electrical_references": sum(row["kind"] == "developer_electrical_point" for row in source_electrical),
                "developer_switch_or_control_references": sum(row["kind"] == "developer_switch_or_control" for row in source_electrical),
                "developer_weak_point_references": sum(
                    row["candidate_role"] in {
                        "wireless_ap_wall_point_reference",
                        "tv_point_reference",
                        "tv_telephone_box_reference",
                    }
                    for row in source_electrical
                ),
                "developer_water_or_drain_references": len(source_plumbing),
                "whole_home_design_status": "candidate_program_only",
                "automatic_ifc_write_allowed": False,
            }
        )

    electrical_counts = Counter(row["disposition_candidate"] for row in electrical)
    plumbing_counts = Counter(row["disposition_candidate"] for row in plumbing)
    report = {
        "source_ifc_sha256": ifc_hash,
        "tolerances_mm": {
            "current_point_near": args.current_point_near_mm,
            "current_wall_near": args.wall_near_mm,
            "plum_coordination": args.plum_coordination_mm,
        },
        "summary": {
            "developer_electrical_and_control_points": len(electrical),
            "developer_plumbing_points": len(plumbing),
            "electrical_dispositions": dict(sorted(electrical_counts.items())),
            "plumbing_dispositions": dict(sorted(plumbing_counts.items())),
            "rooms": len(room_program),
            "rooms_without_current_socket_instance": sum(row["current_socket_instances"] == 0 for row in room_program),
            "current_switch_instances": len(elec["sheets"]["E-302"]["instances"]),
            "current_network_instances": len(elec["sheets"]["E-304"]["instances"]),
        },
        "electrical": electrical,
        "plumbing": plumbing,
        "room_program": room_program,
        "gates": {
            "all_22_spaces_programmed": len(room_program) == 22 and all(row["reference"] for row in room_program),
            "developer_points_are_current_design": False,
            "rough_in_connectors_inferred": False,
            "whole_home_positioning_complete": False,
            "automatic_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
