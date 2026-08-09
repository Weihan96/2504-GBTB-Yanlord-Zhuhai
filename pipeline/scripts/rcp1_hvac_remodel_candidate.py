#!/usr/bin/env python3
"""Build a read-only remodel HVAC pairing candidate from verified IFC evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import numpy as np


EXPECTED_IFC_SHA256 = "7521c09991f3d0c7b7d91ca2324fd55ad961d8e32e9e3e9a9777a4cc19b06e81"
EXPECTED_AC_IDS = {
    "1yW7DASIz8qA$2j8z9tdl2",
    "1QBdVekDnBsOleyo9PM6rT",
    "33chLv3TzEOhKalIJJAPNF",
    "06GpMzzWj1XQobAhD35cgU",
    "1PUCikoaP5fgiYt8sJd8$6",
}
EXPECTED_OPENING_IDS = {
    "3qgu$TepT0J8Yat7ZQ90Mf",
    "2oHWzdjkr8X8Yt0Gd2e3RQ",
    "1JVGu2xtb1ZhGElhvSmyd$",
    "3k2aFe_p5Ah8HJEqWUBB7F",
    "3ERXx822H9jOPo6CetKX9r",
    "0DOeKdT3DE9f2LMR_x$G7q",
    "1_EX1UWfL8ShcGm5BI1UPc",
}
EXPECTED_PIPE_IDS = {
    "0hHnbLj0X4jPDz4o3QJo1l",
    "10Wm8ivdX7dAVfz4cV8l5Q",
    "0Ik2RcgGbFOhdYTJPgh5AQ",
    "1hZRB0eOX8OA8rcjke67P0",
    "0f2ZLauDH8lRnYj6oervDm",
}
EXPECTED_DIFFUSER_IDS = {
    "16Ey9Flj9BK9VRun$ozzjH",
    "3Bv_Kl3jDC5RvaUMdyge1U",
}
EXPECTED_REVIEW_IDS = {
    "RCP1-HVAC-EQUIPMENT-001",
    "RCP1-HVAC-EAST-001",
    "RCP1-HVAC-OPENING-001",
    "RCP1-HVAC-PIPE-001",
    "RCP1-HVAC-AIRSIDE-001",
    "RCP1-HVAC-ROUTE-001",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument(
        "--rcp1-report",
        type=Path,
        default=Path("build/rcp1/rcp1-existing-candidate.json"),
    )
    parser.add_argument(
        "--coordination-report",
        type=Path,
        default=Path("build/rcp1/coordination-report.json"),
    )
    parser.add_argument(
        "--m401-report",
        type=Path,
        default=Path("build/rcp1/m401-existing-report.json"),
    )
    parser.add_argument(
        "--legacy-report",
        type=Path,
        default=Path("build/rcp1/legacy-base-audit.json"),
    )
    parser.add_argument(
        "--review",
        type=Path,
        default=Path("pipeline/decisions/rcp1-hvac-remodel-review.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("build/rcp1/hvac-remodel-candidate.json"),
    )
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_review(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ids = {row["review_id"] for row in rows}
    if ids != EXPECTED_REVIEW_IDS:
        raise RuntimeError(f"HVAC review boundary drift: {sorted(ids)}")
    return rows


def geometry_settings() -> ifcopenshell.geom.settings:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    return settings


def world_mesh_mm(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> tuple[np.ndarray, np.ndarray]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3)) * 1000.0
    faces = np.asarray(shape.geometry.faces, dtype=int).reshape((-1, 3))
    if not len(vertices) or not len(faces):
        raise RuntimeError(f"missing Body mesh: {product.GlobalId}")
    return vertices, faces


def bounds(vertices: np.ndarray) -> dict[str, list[float]]:
    minimum = vertices.min(axis=0)
    maximum = vertices.max(axis=0)
    return {
        "min_mm": minimum.tolist(),
        "max_mm": maximum.tolist(),
        "centre_mm": ((minimum + maximum) / 2.0).tolist(),
        "dimensions_mm": (maximum - minimum).tolist(),
    }


def connected_components(vertices: np.ndarray, faces: np.ndarray) -> list[dict[str, Any]]:
    parents = list(range(len(vertices)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(first: int, second: int) -> None:
        first_root = find(first)
        second_root = find(second)
        if first_root != second_root:
            parents[second_root] = first_root

    used: set[int] = set()
    for face in faces:
        indices = [int(value) for value in face]
        used.update(indices)
        union(indices[0], indices[1])
        union(indices[0], indices[2])

    groups: dict[int, list[int]] = {}
    for index in used:
        groups.setdefault(find(index), []).append(index)

    records = []
    for indices in groups.values():
        component_vertices = vertices[np.asarray(indices, dtype=int)]
        records.append(
            {
                "vertex_count": len(indices),
                "bbox": bounds(component_vertices),
            }
        )
    records.sort(key=lambda item: tuple(item["bbox"]["centre_mm"]))
    for index, record in enumerate(records, start=1):
        record["component_id"] = f"C{index:02d}"
    return records


def relation_index(coordination: dict[str, Any]) -> dict[tuple[str, str], dict[str, Any]]:
    return {tuple(sorted(pair["pair"])): pair for pair in coordination["pairs"]}


def compact_relation(
    index: dict[tuple[str, str], dict[str, Any]],
    first: str,
    second: str,
) -> dict[str, Any]:
    relation = index.get(tuple(sorted((first, second))))
    if relation is None:
        raise RuntimeError(f"missing coordination relation: {first} / {second}")
    return {
        "global_id": second,
        "geometry_state": relation["geometry_state"],
        "minimum_clearance_candidate_mm": relation["minimum_clearance_candidate_mm"],
        "minimum_clearance_is_lower_bound": relation["minimum_clearance_is_lower_bound"],
        "clearance_method": relation["clearance_method"],
    }


def relation_rank(relation: dict[str, Any]) -> tuple[int, float, str]:
    state_rank = {"intersecting": 0, "touching": 1, "near": 2, "separated": 3}
    distance = relation["minimum_clearance_candidate_mm"]
    return (state_rank.get(relation["geometry_state"], 9), float(distance), relation["global_id"])


def inventory_map(report: dict[str, Any], key: str) -> dict[str, dict[str, Any]]:
    return {item["global_id"]: item for item in report["inventory"][key]}


def main() -> int:
    args = parse_args()
    if args.tolerance_mm <= 0:
        raise RuntimeError("tolerance must be positive")

    formal_sha = sha256(args.input)
    if formal_sha != EXPECTED_IFC_SHA256:
        raise RuntimeError(f"formal IFC SHA drift: {formal_sha}")

    rcp1 = read_json(args.rcp1_report)
    coordination = read_json(args.coordination_report)
    m401 = read_json(args.m401_report)
    legacy = read_json(args.legacy_report)
    review = read_review(args.review)

    report_hashes = {
        rcp1["source"]["sha256"],
        coordination["source"]["sha256"],
        m401["source"]["ifc_sha256"],
        legacy["source"]["formal_ifc_sha256"],
    }
    if report_hashes != {formal_sha}:
        raise RuntimeError(f"stale HVAC evidence: {sorted(report_hashes)}")
    if not legacy["gates"]["legacy_base_consistency_pass"]:
        raise RuntimeError("legacy HVAC base audit does not pass")

    ac_inventory = inventory_map(rcp1, "typed_high_equipment")
    opening_inventory = inventory_map(rcp1, "named_high_openings")
    pipe_inventory = inventory_map(rcp1, "high_flow_segments")
    dcl_inventory = inventory_map(rcp1, "dcl_proxies")
    if set(ac_inventory) != EXPECTED_AC_IDS:
        raise RuntimeError("formal AC set drift")
    if set(opening_inventory) != EXPECTED_OPENING_IDS:
        raise RuntimeError("developer AC opening set drift")
    if set(pipe_inventory) != EXPECTED_PIPE_IDS:
        raise RuntimeError("legacy pipe set drift")
    if not EXPECTED_DIFFUSER_IDS <= set(dcl_inventory):
        raise RuntimeError("confirmed diffuser set drift")

    model = ifcopenshell.open(args.input)
    settings = geometry_settings()
    pipe_components = []
    for global_id in sorted(EXPECTED_PIPE_IDS):
        product = model.by_guid(global_id)
        vertices, faces = world_mesh_mm(settings, product)
        components = connected_components(vertices, faces)
        expected_count = next(
            item["formal_component_count"]
            for item in legacy["pipe_comparisons"]
            if item["global_id"] == global_id
        )
        if len(components) != expected_count:
            raise RuntimeError(f"pipe component count drift: {global_id}")
        pipe_components.append(
            {
                "global_id": global_id,
                "name": pipe_inventory[global_id]["name"],
                "component_count": len(components),
                "components": components,
                "role": "legacy_design_base_only",
            }
        )

    relations = relation_index(coordination)
    equipment_pairing = []
    for global_id in sorted(EXPECTED_AC_IDS):
        item = ac_inventory[global_id]
        pipe_relations = sorted(
            (compact_relation(relations, global_id, target) for target in EXPECTED_PIPE_IDS),
            key=relation_rank,
        )
        opening_relations = sorted(
            (compact_relation(relations, global_id, target) for target in EXPECTED_OPENING_IDS),
            key=relation_rank,
        )
        equipment_pairing.append(
            {
                "global_id": global_id,
                "assigned_type": item["assigned_type"],
                "bbox": item["bbox"],
                "primary_space_candidate": item["primary_space_candidate"],
                "pipe_relations": pipe_relations,
                "opening_relations": opening_relations,
                "nearest_opening_candidate": opening_relations[0]["global_id"],
                "basis": "formal world meshes; proximity ranks a review candidate and does not create a connection",
                "confidence": 0.85,
                "human_review_required": True,
            }
        )

    opening_pairing = []
    for global_id in sorted(EXPECTED_OPENING_IDS):
        item = opening_inventory[global_id]
        equipment_relations = sorted(
            (compact_relation(relations, global_id, target) for target in EXPECTED_AC_IDS),
            key=relation_rank,
        )
        pipe_relations = sorted(
            (compact_relation(relations, global_id, target) for target in EXPECTED_PIPE_IDS),
            key=relation_rank,
        )
        opening_pairing.append(
            {
                "global_id": global_id,
                "name": item["name"],
                "bbox": item["bbox"],
                "hosts": item["hosts"],
                "primary_space_candidate": item["primary_space_candidate"],
                "equipment_relations": equipment_relations,
                "pipe_relations": pipe_relations,
                "nearest_equipment_candidate": equipment_relations[0]["global_id"],
                "nearest_pipe_candidate": pipe_relations[0]["global_id"],
                "basis": "formal world meshes and existing opening host; no port or connection inferred",
                "confidence": 0.80,
                "human_review_required": True,
            }
        )

    airside_pairing = []
    for global_id in sorted(EXPECTED_DIFFUSER_IDS):
        item = dcl_inventory[global_id]
        equipment_relations = sorted(
            (compact_relation(relations, global_id, target) for target in EXPECTED_AC_IDS),
            key=relation_rank,
        )
        airside_pairing.append(
            {
                "global_id": global_id,
                "name": item["name"],
                "bbox": item["bbox"],
                "nearest_equipment_candidate": equipment_relations[0]["global_id"],
                "equipment_relations": equipment_relations,
                "confirmed_identity": (
                    "corner_ac_outlet"
                    if global_id == "16Ey9Flj9BK9VRun$ozzjH"
                    else "legacy_defective_ac_outlet"
                ),
                "airside_role": "pending_supply_return_exhaust_design",
                "basis": "confirmed object identity plus formal world-mesh proximity; no airflow inferred",
                "confidence": 0.90,
                "human_review_required": True,
            }
        )

    existing_intersections = [
        relation
        for relation in coordination["pairs"]
        if relation["geometry_state"] == "intersecting"
        and set(relation["pair"]) <= (EXPECTED_AC_IDS | EXPECTED_PIPE_IDS | EXPECTED_OPENING_IDS)
        and (set(relation["pair"]) & EXPECTED_PIPE_IDS)
    ]

    output = {
        "mode": "read_only_rcp1_hvac_remodel_pairing_candidate",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source": {
            "ifc": str(args.input.resolve()),
            "ifc_sha256": formal_sha,
            "legacy_blend_sha256": legacy["source"]["legacy_blend_sha256"],
        },
        "scope": {
            "formal_ifc_write_allowed": False,
            "existing_geometry_moved": False,
            "new_openings_allowed": False,
            "proximity_is_connection": False,
            "legacy_pipe_geometry_is_final_route": False,
        },
        "summary": {
            "formal_ac_candidates": len(EXPECTED_AC_IDS),
            "legacy_east_ac_candidates_without_global_id": 1,
            "developer_openings": len(EXPECTED_OPENING_IDS),
            "legacy_pipe_products": len(EXPECTED_PIPE_IDS),
            "legacy_pipe_independent_components": sum(
                item["component_count"] for item in pipe_components
            ),
            "confirmed_outlet_identity_candidates": len(EXPECTED_DIFFUSER_IDS),
            "existing_pipe_intersections_in_scope": len(existing_intersections),
            "decision_rows": len(review),
        },
        "formal_equipment_pairing": equipment_pairing,
        "legacy_east_ac_candidate": legacy["legacy_east_ac_candidate"],
        "developer_opening_pairing": opening_pairing,
        "airside_pairing": airside_pairing,
        "legacy_pipe_components": pipe_components,
        "existing_intersections": existing_intersections,
        "human_review_bundle": [
            {
                "question_id": "RCP1B-Q01",
                "question": "最终采用正式 IFC 的 5 台设备，还是加入旧 blend 的东侧第 6 机位？",
                "required_evidence": "服务房间、回风可达性、检修包络和既有洞口配对",
            },
            {
                "question_id": "RCP1B-Q02",
                "question": "逐台确认服务房间及送风、回风形式；转角风口不得仅凭几何自动判定送风或回风。",
                "required_evidence": "设备能力、房间负荷、风量、风口尺寸和回风路径",
            },
            {
                "question_id": "RCP1B-Q03",
                "question": "确认室外机或立管接口、冷凝水排放点和设备厂家接管位置后，才能生成装修后管线路由。",
                "required_evidence": "液管/气管管径与保温、冷凝水坡度、吊架、检修和穿洞顺序",
            },
        ],
        "gates": {
            "all_source_hashes_current": report_hashes == {formal_sha},
            "formal_object_sets_exact": True,
            "legacy_pipe_topology_verified": legacy["gates"][
                "five_pipe_topology_and_component_bounds_within_tolerance"
            ],
            "developer_opening_constraint_recorded": any(
                row["review_id"] == "RCP1-HVAC-OPENING-001"
                and row["review_status"] == "confirmed_constraint"
                for row in review
            ),
            "candidate_ready_for_blender_review": True,
            "formal_ifc_write_allowed": False,
            "hvac_design_ready": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output["summary"], ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
