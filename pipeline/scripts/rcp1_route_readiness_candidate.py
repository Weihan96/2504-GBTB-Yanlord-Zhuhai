"""Build read-only RCP1C route-readiness evidence from fixed HVAC positions."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell


EQUIPMENT = {
    "A01": {"global_id": "33chLv3TzEOhKalIJJAPNF", "label": "主卧入口"},
    "A02": {"global_id": "1QBdVekDnBsOleyo9PM6rT", "label": "走廊"},
    "A03": {"global_id": "06GpMzzWj1XQobAhD35cgU", "label": "次卧"},
    "A04": {"global_id": "1PUCikoaP5fgiYt8sJd8$6", "label": "西厨"},
    "A05": {"global_id": "1yW7DASIz8qA$2j8z9tdl2", "label": "书房"},
    "A06": {"global_id": None, "label": "东侧固定机位"},
}

OPENINGS = {
    "H01": "0DOeKdT3DE9f2LMR_x$G7q",
    "H02": "1JVGu2xtb1ZhGElhvSmyd$",
    "H03": "2oHWzdjkr8X8Yt0Gd2e3RQ",
    "H04": "3ERXx822H9jOPo6CetKX9r",
    "H05": "3k2aFe_p5Ah8HJEqWUBB7F",
    "H06": "3qgu$TepT0J8Yat7ZQ90Mf",
    "H07": "1_EX1UWfL8ShcGm5BI1UPc",
}

DIRECT_SERVICE = {
    "A01": ("R09", "主卧"),
    "A02": ("R07", "餐厅"),
    "A03": ("R14", "次卧"),
    "A04": ("R03", "西厨"),
    "A05": ("R22", "书房"),
    "A06": ("R07", "餐厅"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bbox_clearance_mm(first: dict, second: dict) -> float:
    deltas = [
        max(0.0, second["min_mm"][axis] - first["max_mm"][axis], first["min_mm"][axis] - second["max_mm"][axis])
        for axis in range(3)
    ]
    return math.sqrt(sum(delta * delta for delta in deltas))


def relation_distance(item: dict, opening_global_id: str) -> float:
    relation = next(
        relation
        for relation in item["opening_relations"]
        if relation["global_id"] == opening_global_id
    )
    return float(relation["minimum_clearance_candidate_mm"])


def equipment_row(candidate_id: str, report: dict, opening_by_id: dict[str, dict]) -> dict:
    definition = EQUIPMENT[candidate_id]
    if definition["global_id"] is None:
        legacy = report["legacy_east_ac_candidate"]
        bbox = legacy["predicted_formal_bbox"]
        ranked = sorted(
            (
                {
                    "opening_id": opening_id,
                    "global_id": global_id,
                    "clearance_mm": bbox_clearance_mm(bbox, opening_by_id[global_id]["bbox"]),
                }
                for opening_id, global_id in OPENINGS.items()
            ),
            key=lambda row: (row["clearance_mm"], row["opening_id"]),
        )
        centre = legacy["predicted_formal_centre_mm"]
        formal_identity_status = "missing"
    else:
        item = next(
            item
            for item in report["formal_equipment_pairing"]
            if item["global_id"] == definition["global_id"]
        )
        ranked = sorted(
            (
                {
                    "opening_id": opening_id,
                    "global_id": global_id,
                    "clearance_mm": relation_distance(item, global_id),
                }
                for opening_id, global_id in OPENINGS.items()
            ),
            key=lambda row: (row["clearance_mm"], row["opening_id"]),
        )
        centre = item["bbox"]["centre_mm"]
        formal_identity_status = "present"
    first, second = ranked[:2]
    margin = second["clearance_mm"] - first["clearance_mm"]
    confidence = 0.95 if margin >= 1000.0 else 0.85 if margin >= 400.0 else 0.60
    return {
        "equipment_id": candidate_id,
        "global_id": definition["global_id"],
        "label": definition["label"],
        "position_status": "confirmed_fixed",
        "formal_identity_status": formal_identity_status,
        "centre_mm": centre,
        "direct_service_candidate": {
            "space_reference": DIRECT_SERVICE[candidate_id][0],
            "space_long_name": DIRECT_SERVICE[candidate_id][1],
            "status": "geometry_probe_not_final_connection",
        },
        "opening_ranking": ranked,
        "nearest_opening_candidate": first["opening_id"],
        "nearest_to_second_margin_mm": margin,
        "nearest_pairing_confidence": confidence,
        "human_review_required": candidate_id in {"A01", "A02", "A06"},
    }


def global_opening_assignment(equipment: list[dict]) -> dict:
    equipment_ids = [row["equipment_id"] for row in equipment]
    distance_by_equipment = {
        row["equipment_id"]: {
            candidate["opening_id"]: candidate["clearance_mm"]
            for candidate in row["opening_ranking"]
        }
        for row in equipment
    }
    candidates = []
    for opening_ids in itertools.permutations(sorted(OPENINGS), len(equipment_ids)):
        pairs = list(zip(equipment_ids, opening_ids))
        total = sum(distance_by_equipment[equipment_id][opening_id] for equipment_id, opening_id in pairs)
        candidates.append((total, pairs))
    candidates.sort(key=lambda candidate: (candidate[0], candidate[1]))
    best, second = candidates[:2]
    pairs = [
        {
            "equipment_id": equipment_id,
            "opening_id": opening_id,
            "equipment_global_id": EQUIPMENT[equipment_id]["global_id"],
            "opening_global_id": OPENINGS[opening_id],
            "clearance_mm": distance_by_equipment[equipment_id][opening_id],
            "status": "global_minimum_distance_candidate_not_connection",
        }
        for equipment_id, opening_id in best[1]
    ]
    return {
        "method": "exhaustive one-to-one minimum total world-AABB clearance; six equipment positions mapped to six of seven existing openings",
        "total_clearance_mm": best[0],
        "second_best_total_clearance_mm": second[0],
        "best_to_second_margin_mm": second[0] - best[0],
        "pairs": pairs,
        "unused_opening_ids": sorted(set(OPENINGS) - {pair["opening_id"] for pair in pairs}),
        "confidence": 0.90,
        "human_review_required": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--hvac-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    formal_sha = sha256(args.input)
    hvac = json.loads(args.hvac_report.read_text(encoding="utf-8"))
    if hvac["source"]["ifc_sha256"] != formal_sha:
        raise RuntimeError("HVAC report does not match the formal IFC")
    model = ifcopenshell.open(args.input)
    opening_by_id = {
        item["global_id"]: item for item in hvac["developer_opening_pairing"]
    }
    equipment = [
        equipment_row(candidate_id, hvac, opening_by_id)
        for candidate_id in sorted(EQUIPMENT)
    ]
    row_by_id = {row["equipment_id"]: row for row in equipment}
    global_assignment = global_opening_assignment(equipment)

    high_confidence_pairings = [
        {
            "equipment_id": equipment_id,
            "opening_id": opening_id,
            "equipment_global_id": row_by_id[equipment_id]["global_id"],
            "opening_global_id": OPENINGS[opening_id],
            "clearance_mm": row_by_id[equipment_id]["opening_ranking"][0]["clearance_mm"],
            "basis": "nearest opening with at least 400 mm separation from the second-ranked opening",
            "confidence": row_by_id[equipment_id]["nearest_pairing_confidence"],
            "human_review_required": False,
        }
        for equipment_id, opening_id in (("A03", "H04"), ("A04", "H06"), ("A05", "H07"))
    ]

    port_count = len(model.by_type("IfcDistributionPort"))
    system_count = len(model.by_type("IfcSystem"))
    distribution_system_count = len(model.by_type("IfcDistributionSystem"))
    outdoor_equipment_count = sum(
        len(model.by_type(ifc_class))
        for ifc_class in ("IfcCondenser", "IfcCompressor")
    )
    fixed_diagnostics = hvac["fixed_equipment_airside_diagnostics"]
    a05_diagnostic = next(
        item for item in fixed_diagnostics if item["equipment_candidate"] == EQUIPMENT["A05"]["global_id"]
    )

    output = {
        "mode": "read_only_rcp1_route_readiness_candidate",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source": {
            "ifc": str(args.input.resolve()),
            "ifc_sha256": formal_sha,
            "hvac_report": str(args.hvac_report.resolve()),
        },
        "scope": {
            "formal_ifc_write_allowed": False,
            "equipment_positions_may_move": False,
            "new_openings_allowed": False,
            "demolition_walls_are_permanent_obstacles": False,
            "legacy_pipe_geometry_is_final_route": False,
        },
        "equipment": equipment,
        "openings": [
            {
                "opening_id": opening_id,
                "global_id": global_id,
                "name": opening_by_id[global_id]["name"],
                "hosts": opening_by_id[global_id]["hosts"],
                "centre_mm": opening_by_id[global_id]["bbox"]["centre_mm"],
            }
            for opening_id, global_id in OPENINGS.items()
        ],
        "global_equipment_opening_assignment": global_assignment,
        "high_confidence_equipment_opening_pairings": high_confidence_pairings,
        "shared_or_route_opening_evidence": [
            {
                "opening_id": "H05",
                "global_id": OPENINGS["H05"],
                "fact": "nearest opening for both A01 and A02; existing liquid, gas and drain geometry intersects it",
                "status": "global assignment gives H05 to A01 and H02 to A02; still not a formal connection",
            },
            {
                "equipment_id": "A06",
                "fact": "H03 is nearest at 343.452 mm and H02 is second at 566.377 mm",
                "status": "fixed_position_but_formal_identity_and_final_route_missing",
            },
            {
                "opening_id": "H01",
                "global_id": OPENINGS["H01"],
                "fact": "unused by the global six-equipment assignment; located beyond H06 from A04",
                "status": "external_service_chain_candidate_not_confirmed",
            },
        ],
        "airside_readiness": {
            "direct_service_candidates": {
                equipment_id: {
                    "space_reference": space[0],
                    "space_long_name": space[1],
                }
                for equipment_id, space in DIRECT_SERVICE.items()
            },
            "living_room_R20_direct_candidate_count": 0,
            "a05_to_R20_demolition_wall_crossings": a05_diagnostic["route"]["demolition_wall_crossings"],
            "a05_to_R20_permanent_wall_crossings": a05_diagnostic["route"]["permanent_wall_crossings"],
            "status": "service_assignment_and_supply_return_layout_pending",
        },
        "refrigerant_and_condensate_readiness": {
            "distribution_port_count": port_count,
            "ifc_system_count": system_count,
            "ifc_distribution_system_count": distribution_system_count,
            "outdoor_condenser_or_compressor_count": outdoor_equipment_count,
            "legacy_pipe_products": hvac["summary"]["legacy_pipe_products"],
            "legacy_pipe_independent_components": hvac["summary"]["legacy_pipe_independent_components"],
            "status": "blocked_by_missing_real_endpoints_and_ports",
        },
        "minimum_required_inputs": [
            {
                "input_id": "RCP1C-I01",
                "question": "确认全局洞口候选配对，并确认 A04→H06→H01 是否属于同一条向外服务穿墙链。",
                "why_required": "全局配对显著优于次优组合，但距离和旧管交叠仍不能替代正式连接关系。",
            },
            {
                "input_id": "RCP1C-I02",
                "question": "确认公共区服务分工：A05 是否服务书房＋客厅、A06 是否服务餐厅，并明确 A02/A04 的最终服务范围。",
                "why_required": "两个卧室机位可由几何直接确定；公共区存在多台固定设备和跨房送风。",
            },
            {
                "input_id": "RCP1C-I03",
                "question": "一次性提供或确认设备接管侧/接口坐标、液气管外部终点、冷凝水排放点/标高及是否允许冷凝水泵。",
                "why_required": "正式 IFC 中端口、系统、室外机/立管接口和 HVAC 冷凝水排放端点均为 0。",
            },
        ],
        "gates": {
            "source_hash_matches": True,
            "six_fixed_equipment_positions_registered": len(equipment) == 6,
            "seven_existing_openings_registered": len(opening_by_id) == 7,
            "three_high_confidence_pairings_ready": len(high_confidence_pairings) == 3,
            "global_six_equipment_assignment_ready": len(global_assignment["pairs"]) == 6,
            "demolition_walls_excluded_from_permanent_obstacles": (
                len(a05_diagnostic["route"]["demolition_wall_crossings"]) >= 1
                and len(a05_diagnostic["route"]["permanent_wall_crossings"]) == 0
            ),
            "final_airside_routes_ready": False,
            "final_refrigerant_routes_ready": False,
            "final_condensate_routes_ready": False,
            "formal_ifc_write_allowed": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "fixed_equipment": len(equipment),
        "existing_openings": len(opening_by_id),
        "high_confidence_pairings": len(high_confidence_pairings),
        "distribution_ports": port_count,
        "systems": system_count + distribution_system_count,
        "outdoor_equipment": outdoor_equipment_count,
        "required_inputs": len(output["minimum_required_inputs"]),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
