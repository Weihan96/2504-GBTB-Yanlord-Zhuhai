"""Audit an exported Blender HVAC constraint preview against IFC-backed anchors."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--route-report", type=Path, required=True)
    parser.add_argument("--preview", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    args = parser.parse_args()

    formal_sha = sha256(args.input)
    route_report = json.loads(args.route_report.read_text(encoding="utf-8"))
    preview = json.loads(args.preview.read_text(encoding="utf-8"))
    if route_report["source"]["ifc_sha256"] != formal_sha:
        raise RuntimeError("route report does not match the formal IFC")
    if preview["source_ifc_sha256"] != formal_sha:
        raise RuntimeError("Blender preview candidate does not match the formal IFC")
    if preview.get("formal_ifc_write_allowed") is not False:
        raise RuntimeError("Blender preview must not authorize formal IFC writes")

    expected_routes = {
        route["route_id"]: route
        for route in route_report["confirmed_route_graph"]
        if len(route["waypoints"]) >= 2
    }
    preview_routes = {route["route_id"]: route for route in preview["routes"]}
    if set(preview_routes) != set(expected_routes):
        raise RuntimeError("Blender preview route set drift")

    anchor_checks = []
    route_checks = []
    maximum_anchor_deviation = 0.0
    bend_count = 0
    for route_id, expected in expected_routes.items():
        actual = preview_routes[route_id]
        expected_by_id = {row["anchor_id"]: row for row in expected["waypoints"]}
        expected_fixed_order = [row["anchor_id"] for row in expected["waypoints"]]
        actual_fixed_order = [
            row["anchor_id"]
            for row in actual["anchors"]
            if row["anchor_kind"] != "blender_bend"
        ]
        if actual_fixed_order != expected_fixed_order:
            raise RuntimeError(f"fixed waypoint order drift: {route_id}")
        for row in actual["anchors"]:
            if row["anchor_kind"] == "blender_bend":
                bend_count += 1
                anchor_checks.append({
                    "route_id": route_id,
                    "anchor_id": row["anchor_id"],
                    "anchor_kind": "blender_bend",
                    "world_mm": row["world_mm"],
                    "status": "human_review_required_before_register_update",
                })
                continue
            expected_row = expected_by_id[row["anchor_id"]]
            deviation = math.dist(row["world_mm"], expected_row["centre_mm"])
            maximum_anchor_deviation = max(maximum_anchor_deviation, deviation)
            anchor_checks.append({
                "route_id": route_id,
                "anchor_id": row["anchor_id"],
                "anchor_kind": row["anchor_kind"],
                "global_id": row["global_id"],
                "deviation_from_current_ifc_anchor_mm": deviation,
                "within_tolerance": deviation <= args.tolerance_mm,
            })
        segment_lengths = [
            math.dist(first["world_mm"], second["world_mm"])
            for first, second in zip(actual["anchors"], actual["anchors"][1:])
        ]
        route_checks.append({
            "route_id": route_id,
            "anchor_order": [row["anchor_id"] for row in actual["anchors"]],
            "segment_lengths_mm": segment_lengths,
            "total_length_mm": sum(segment_lengths),
            "status": "constraint_skeleton_not_fabrication_geometry",
        })

    fixed_anchor_checks = [row for row in anchor_checks if row["anchor_kind"] != "blender_bend"]
    fixed_anchors_pass = all(row["within_tolerance"] for row in fixed_anchor_checks)
    output = {
        "mode": "read_only_rcp1_hvac_preview_audit",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source": {
            "ifc": str(args.input.resolve()),
            "ifc_sha256": formal_sha,
            "route_report": str(args.route_report.resolve()),
            "preview_candidate": str(args.preview.resolve()),
        },
        "tolerance_mm": args.tolerance_mm,
        "route_checks": route_checks,
        "anchor_checks": anchor_checks,
        "summary": {
            "route_count": len(route_checks),
            "fixed_anchor_check_count": len(fixed_anchor_checks),
            "temporary_bend_count": bend_count,
            "maximum_fixed_anchor_deviation_mm": maximum_anchor_deviation,
        },
        "gates": {
            "source_hash_matches": True,
            "route_set_matches": True,
            "fixed_anchor_order_matches": True,
            "fixed_anchors_within_tolerance": fixed_anchors_pass,
            "temporary_bends_require_human_review": bend_count > 0,
            "fabrication_geometry_ready": False,
            "formal_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output["summary"], ensure_ascii=False))
    return 0 if fixed_anchors_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
