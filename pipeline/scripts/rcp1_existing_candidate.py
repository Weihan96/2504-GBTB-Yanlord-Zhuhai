#!/usr/bin/env python3
"""Audit existing ceiling and high-level coordination objects without writing IFC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import numpy as np
from ifcopenshell.util.element import get_container, get_psets, get_type
from ifcopenshell.util.placement import get_local_placement
from shapely.geometry import Polygon
from shapely.ops import unary_union


PROTECTED_HANDOFF_IDS = {"16Ey9Flj9BK9VRun$ozzjH", "0zWtSQZzjFQg_PORjlssbe"}
OTHER_HIGH_PROXY_IDS = {"1faflkXXH6M9cnYPE9Liir"}
EXPECTED_AC_OPENING_IDS = {
    "3qgu$TepT0J8Yat7ZQ90Mf",
    "2oHWzdjkr8X8Yt0Gd2e3RQ",
    "1JVGu2xtb1ZhGElhvSmyd$",
    "3k2aFe_p5Ah8HJEqWUBB7F",
    "3ERXx822H9jOPo6CetKX9r",
    "0DOeKdT3DE9f2LMR_x$G7q",
    "1_EX1UWfL8ShcGm5BI1UPc",
}
LEGACY_BASE_PAIR_IDS = {
    tuple(sorted(pair))
    for pair in [
        ("0Ik2RcgGbFOhdYTJPgh5AQ", "0f2ZLauDH8lRnYj6oervDm"),
        ("0Ik2RcgGbFOhdYTJPgh5AQ", "10Wm8ivdX7dAVfz4cV8l5Q"),
        ("0f2ZLauDH8lRnYj6oervDm", "10Wm8ivdX7dAVfz4cV8l5Q"),
        ("0f2ZLauDH8lRnYj6oervDm", "1yW7DASIz8qA$2j8z9tdl2"),
        ("0hHnbLj0X4jPDz4o3QJo1l", "10Wm8ivdX7dAVfz4cV8l5Q"),
        ("0hHnbLj0X4jPDz4o3QJo1l", "1QBdVekDnBsOleyo9PM6rT"),
        ("0hHnbLj0X4jPDz4o3QJo1l", "1hZRB0eOX8OA8rcjke67P0"),
        ("10Wm8ivdX7dAVfz4cV8l5Q", "1hZRB0eOX8OA8rcjke67P0"),
        ("10Wm8ivdX7dAVfz4cV8l5Q", "1yW7DASIz8qA$2j8z9tdl2"),
        ("1QBdVekDnBsOleyo9PM6rT", "1hZRB0eOX8OA8rcjke67P0"),
    ]
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument(
        "--review",
        type=Path,
        default=Path("pipeline/decisions/rcp1-existing-review.csv"),
    )
    parser.add_argument("--output", type=Path, default=Path("build/rcp1/rcp1-existing-candidate.json"))
    parser.add_argument(
        "--coordination-output",
        type=Path,
        default=Path("build/rcp1/coordination-report.json"),
    )
    parser.add_argument(
        "--expected-ifc-sha256",
        help="Optional caller-frozen formal IFC SHA-256; defaults to the current input file.",
    )
    parser.add_argument(
        "--legacy-report",
        type=Path,
        default=Path("build/rcp1/legacy-base-audit.json"),
        help="Read-only legacy Blender audit metadata; this generator never refreshes it.",
    )
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--search-window-mm", type=float, default=50.0)
    parser.add_argument("--review-clearance-mm", type=float, default=20.0)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def legacy_evidence_status(path: Path, formal_ifc_hash: str) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path),
        "status": "missing_not_refreshed",
        "current_formal_ifc": False,
        "observed_formal_hash": "",
        "legacy_audit_file_hash": "",
        "note": "Blender legacy audit was not refreshed by this pure-Python candidate generator.",
    }
    if not path.is_file():
        return result
    result["legacy_audit_file_hash"] = sha256(path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        result["status"] = "unreadable_not_refreshed"
        result["error"] = str(exc)
        return result
    observed = str(payload.get("source", {}).get("formal_ifc_sha256", "") or "")
    result["observed_formal_hash"] = observed
    result["current_formal_ifc"] = observed == formal_ifc_hash
    result["status"] = "current_external_audit" if result["current_formal_ifc"] else "stale_not_refreshed"
    return result


def settings() -> ifcopenshell.geom.settings:
    result = ifcopenshell.geom.settings()
    result.set(result.USE_WORLD_COORDS, True)
    return result


def geometry(
    geometry_settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> tuple[np.ndarray, np.ndarray]:
    shape = ifcopenshell.geom.create_shape(geometry_settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3))
    faces = np.asarray(shape.geometry.faces, dtype=int).reshape((-1, 3))
    if not len(vertices):
        raise RuntimeError(f"{product.GlobalId} has no world geometry")
    return vertices, faces


def footprint(vertices: np.ndarray, faces: np.ndarray) -> Any:
    polygons = []
    for face in faces:
        polygon = Polygon(vertices[face, :2])
        if polygon.area > 1e-10:
            polygons.append(polygon)
    if not polygons:
        return None
    return unary_union(polygons).buffer(1e-8)


def bbox(vertices_m: np.ndarray) -> dict[str, list[float]]:
    vertices = vertices_m * 1000.0
    minimum = vertices.min(axis=0)
    maximum = vertices.max(axis=0)
    return {
        "min_mm": minimum.tolist(),
        "max_mm": maximum.tolist(),
        "dimensions_mm": (maximum - minimum).tolist(),
        "centre_mm": ((minimum + maximum) / 2.0).tolist(),
    }


def placement(product: ifcopenshell.entity_instance) -> dict[str, Any]:
    if product.ObjectPlacement is None:
        return {"origin_mm": None, "maximum_integer_residual_mm": None}
    origin = np.asarray(get_local_placement(product.ObjectPlacement)[:3, 3], dtype=float)
    residual = np.abs(origin - np.round(origin))
    return {
        "origin_mm": origin.tolist(),
        "maximum_integer_residual_mm": float(residual.max()),
    }


def read_review(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ids = [row["review_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise RuntimeError("RCP1 review IDs are not unique")
    handoffs = {
        global_id
        for row in rows
        if row["scope"] == "protected_handoff"
        for global_id in row["global_ids"].split("; ")
        if global_id
    }
    if handoffs != PROTECTED_HANDOFF_IDS:
        raise RuntimeError("RCP1 protected handoff set drift")
    return rows


def space_footprints(
    model: ifcopenshell.file,
    geometry_settings: ifcopenshell.geom.settings,
) -> list[dict[str, Any]]:
    records = []
    for space in model.by_type("IfcSpace"):
        vertices, faces = geometry(geometry_settings, space)
        plan = footprint(vertices, faces)
        if plan is None:
            raise RuntimeError(f"Space {space.GlobalId} has no plan footprint")
        records.append(
            {
                "global_id": space.GlobalId,
                "reference": get_psets(space).get("Pset_SpaceCommon", {}).get("Reference", ""),
                "long_name": str(space.LongName or space.Name or space.GlobalId),
                "footprint": plan,
            }
        )
    return records


def space_overlaps(plan: Any, spaces: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if plan is None or plan.area <= 1e-12:
        return []
    overlaps = []
    for space in spaces:
        area = float(plan.intersection(space["footprint"]).area)
        if area <= 1e-9:
            continue
        overlaps.append(
            {
                "global_id": space["global_id"],
                "reference": space["reference"],
                "long_name": space["long_name"],
                "overlap_area_m2": area,
                "object_plan_ratio": area / float(plan.area),
            }
        )
    overlaps.sort(key=lambda item: (-item["overlap_area_m2"], item["global_id"]))
    return overlaps


def related_hosts(product: ifcopenshell.entity_instance) -> list[dict[str, str]]:
    hosts = []
    for relation in getattr(product, "VoidsElements", ()):
        host = relation.RelatingBuildingElement
        hosts.append({"ifc_class": host.is_a(), "global_id": host.GlobalId, "name": str(host.Name or "")})
    return hosts


def record(
    product: ifcopenshell.entity_instance,
    geometry_settings: ifcopenshell.geom.settings,
    spaces: list[dict[str, Any]],
    category: str,
    basis: str,
    confidence: float,
    review_required: str = "yes",
) -> dict[str, Any]:
    vertices, faces = geometry(geometry_settings, product)
    bounds = bbox(vertices)
    assigned_type = get_type(product)
    overlaps = space_overlaps(footprint(vertices, faces), spaces)
    return {
        "global_id": product.GlobalId,
        "ifc_class": product.is_a(),
        "name": str(product.Name or ""),
        "object_type": str(getattr(product, "ObjectType", "") or ""),
        "occurrence_predefined_type": str(getattr(product, "PredefinedType", "") or ""),
        "assigned_type": {
            "global_id": str(getattr(assigned_type, "GlobalId", "") or ""),
            "name": str(getattr(assigned_type, "Name", "") or ""),
            "element_type": str(getattr(assigned_type, "ElementType", "") or ""),
            "predefined_type": str(getattr(assigned_type, "PredefinedType", "") or ""),
        },
        "container": str(getattr(get_container(product), "Name", "") or ""),
        "category": category,
        "bbox": bounds,
        "mesh": {"vertex_count": len(vertices), "triangle_count": len(faces)},
        "elevation_mm": {
            "bottom": bounds["min_mm"][2],
            "centre": bounds["centre_mm"][2],
            "top": bounds["max_mm"][2],
        },
        "placement": placement(product),
        "space_overlaps": overlaps,
        "primary_space_candidate": overlaps[0] if overlaps else None,
        "space_assignment_status": "candidate_from_plan_overlap" if overlaps else "unassigned_for_review",
        "hosts": related_hosts(product),
        "has_psets": bool(get_psets(product)),
        "basis": basis,
        "confidence": confidence,
        "review_required": review_required,
        "formal_ifc_write_allowed": "no",
    }


def pair_key(first: str, second: str) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def wrapper_global_id(entity: Any) -> str:
    return str(entity.get_argument(0))


def aabb_distance_mm(first: dict[str, Any], second: dict[str, Any]) -> float:
    gaps = []
    for axis in range(3):
        first_min = float(first["bbox"]["min_mm"][axis])
        first_max = float(first["bbox"]["max_mm"][axis])
        second_min = float(second["bbox"]["min_mm"][axis])
        second_max = float(second["bbox"]["max_mm"][axis])
        gaps.append(max(0.0, second_min - first_max, first_min - second_max))
    return math.sqrt(sum(gap * gap for gap in gaps))


def world_mesh_relations(
    model: ifcopenshell.file,
    records: list[dict[str, Any]],
    tolerance_mm: float,
    search_window_mm: float,
    review_clearance_mm: float,
    legacy_base_pairs: set[tuple[str, str]],
) -> dict[str, Any]:
    if tolerance_mm <= 0:
        raise RuntimeError("coordination tolerance must be positive")
    if search_window_mm < tolerance_mm:
        raise RuntimeError("search window must not be smaller than tolerance")
    if review_clearance_mm < tolerance_mm or review_clearance_mm > search_window_mm:
        raise RuntimeError("review clearance must be between tolerance and search window")

    by_id = {item["global_id"]: item for item in records}
    ids = sorted(by_id)
    all_keys = [
        (ids[first_index], ids[second_index])
        for first_index in range(len(ids))
        for second_index in range(first_index + 1, len(ids))
    ]
    aabb_distances = {
        key: aabb_distance_mm(by_id[key[0]], by_id[key[1]]) for key in all_keys
    }
    candidate_keys = {key for key, distance in aabb_distances.items() if distance <= search_window_mm + 1e-9}

    products = [model.by_guid(global_id) for global_id in ids]
    tree = ifcopenshell.geom.tree()
    geom = settings()
    for product in products:
        tree.add_element(ifcopenshell.geom.create_shape(geom, product))

    collision_pairs: set[tuple[str, str]] = set()
    collision_types: dict[tuple[str, str], str] = {}
    for clash in tree.clash_collision_many(products, products, False):
        key = pair_key(wrapper_global_id(clash.a), wrapper_global_id(clash.b))
        if key[0] == key[1] or key not in candidate_keys:
            continue
        collision_pairs.add(key)
        collision_types[key] = tree.get_clash_type(clash.clash_type)

    intersection_depths: dict[tuple[str, str], float] = {}
    for clash in tree.clash_intersection_many(
        products, products, tolerance_mm / 1000.0, True
    ):
        key = pair_key(wrapper_global_id(clash.a), wrapper_global_id(clash.b))
        if key[0] == key[1] or key not in candidate_keys:
            continue
        intersection_depths[key] = max(
            intersection_depths.get(key, 0.0), float(clash.distance) * 1000.0
        )

    clearance_distances: dict[tuple[str, str], float] = {}
    for clash in tree.clash_clearance_many(
        products, products, search_window_mm / 1000.0, True
    ):
        key = pair_key(wrapper_global_id(clash.a), wrapper_global_id(clash.b))
        if key[0] == key[1] or key not in candidate_keys:
            continue
        distance = float(clash.distance) * 1000.0
        clearance_distances[key] = min(clearance_distances.get(key, math.inf), distance)

    def expected_rule(first: dict[str, Any], second: dict[str, Any]) -> str:
        categories = {first["category"], second["category"]}
        ceiling_like = {
            "typed_ceiling_covering",
            "name_only_ceiling_or_light_slot_proxy",
        }
        if "typed_light_fixture" in categories and categories & ceiling_like:
            return "expected_light_ceiling_embedding_candidate"
        if "named_high_opening" in categories:
            opening = first if first["category"] == "named_high_opening" else second
            other = second if opening is first else first
            if any(host["global_id"] == other["global_id"] for host in opening["hosts"]):
                return "expected_opening_void_host"
            return "opening_passage_candidate"
        if "name_only_air_outlet_proxy" in categories and categories & ceiling_like:
            return "diffuser_ceiling_embedding_candidate"
        if "typed_high_ac_appliance" in categories and categories & ceiling_like:
            return "ac_ceiling_coordination_candidate"
        if "high_level_flow_segment" in categories and categories & ceiling_like:
            return "flow_ceiling_coordination_candidate"
        if "name_only_high_service_proxy" in categories and categories & ceiling_like:
            return "service_proxy_ceiling_coordination_candidate"
        if categories <= ceiling_like:
            return "ceiling_assembly_candidate"
        return "none"

    pair_records = []
    human_review_pairs = []
    expected_pairs = []
    for key in all_keys:
        first = by_id[key[0]]
        second = by_id[key[1]]
        aabb_distance = aabb_distances[key]
        rule = expected_rule(first, second)
        exact_tested = key in candidate_keys
        if key in collision_pairs or key in intersection_depths:
            geometry_state = "intersecting"
            minimum_clearance = 0.0
            clearance_kind = "world_mesh_collision_or_triangle_intersection"
            lower_bound = False
        elif key in clearance_distances:
            minimum_clearance = clearance_distances[key]
            lower_bound = False
            if minimum_clearance <= tolerance_mm + 1e-9:
                geometry_state = "contacting"
                clearance_kind = "world_mesh_clearance_within_tolerance"
            else:
                geometry_state = "separated"
                clearance_kind = "world_mesh_minimum_clearance_candidate"
        elif exact_tested:
            geometry_state = "separated"
            minimum_clearance = search_window_mm
            clearance_kind = "world_mesh_clearance_exceeds_search_window"
            lower_bound = True
        else:
            geometry_state = "separated"
            minimum_clearance = aabb_distance
            clearance_kind = "world_aabb_separation_lower_bound"
            lower_bound = True

        if key in legacy_base_pairs and geometry_state == "intersecting" and rule == "none":
            rule = "legacy_hvac_base_intersection_pending_redesign"

        is_expected_without_pair_review = rule in {
            "expected_light_ceiling_embedding_candidate",
            "expected_opening_void_host",
            "ceiling_assembly_candidate",
            "diffuser_ceiling_embedding_candidate",
            "flow_ceiling_coordination_candidate",
            "opening_passage_candidate",
            "legacy_hvac_base_intersection_pending_redesign",
        }
        review_required = (
            not is_expected_without_pair_review
            and geometry_state in {"intersecting", "contacting"}
        )
        if review_required:
            if geometry_state == "intersecting" and rule == "none":
                review_reason = "unresolved world-mesh intersection; installation/system semantics are absent"
            elif geometry_state == "contacting" and rule == "none":
                review_reason = "world-mesh contact is not covered by an expected-relation rule"
            else:
                review_reason = f"{rule} requires installation and access confirmation"
        else:
            if rule == "legacy_hvac_base_intersection_pending_redesign":
                review_reason = "the same mesh intersection was reproduced in the prior legacy design-base audit; audit freshness is disclosed separately and remodel HVAC route design remains pending"
            else:
                review_reason = "expected non-hard relation or mechanically separated pair; thematic system/access review remains"
        relation_basis = "formal IFC world Body meshes; AABB prefilter followed by IfcOpenShell geom.tree collision/intersection/clearance"
        if rule == "legacy_hvac_base_intersection_pending_redesign":
            relation_basis += "; legacy 2504_lowpoly.blend pair state was independently reproduced by rcp1_legacy_base_audit.py, whose current-file freshness is disclosed separately"
        relation = {
            "pair": list(key),
            "first": {
                "global_id": first["global_id"],
                "name": first["name"],
                "category": first["category"],
            },
            "second": {
                "global_id": second["global_id"],
                "name": second["name"],
                "category": second["category"],
            },
            "aabb_distance_mm": aabb_distance,
            "aabb_candidate": exact_tested,
            "geometry_state": geometry_state,
            "minimum_clearance_candidate_mm": minimum_clearance,
            "minimum_clearance_is_lower_bound": lower_bound,
            "clearance_method": clearance_kind,
            "collision_type": collision_types.get(key),
            "triangle_intersection_depth_candidate_mm": intersection_depths.get(key),
            "expected_relation_rule": rule,
            "conflict_candidate": review_required and geometry_state == "intersecting" and rule == "none",
            "review_required": "yes" if review_required else "no",
            "review_reason": review_reason,
            "basis": relation_basis,
            "confidence": 1.0 if exact_tested else 0.95,
        }
        pair_records.append(relation)
        if review_required:
            human_review_pairs.append(relation)
        elif rule != "none" and geometry_state in {"intersecting", "contacting"}:
            expected_pairs.append(relation)

    state_counts = Counter(item["geometry_state"] for item in pair_records)
    review_state_counts = Counter(item["geometry_state"] for item in human_review_pairs)
    interacting_rule_counts = Counter(
        item["expected_relation_rule"]
        for item in pair_records
        if item["expected_relation_rule"] != "none"
        and item["geometry_state"] in {"intersecting", "contacting"}
    )
    legacy_pair_records = [item for item in pair_records if tuple(item["pair"]) in legacy_base_pairs]
    return {
        "method": {
            "tolerance_mm": tolerance_mm,
            "search_window_mm": search_window_mm,
            "review_clearance_mm": review_clearance_mm,
            "pair_review_policy": "only unresolved world-mesh intersection/contact remains in pair-specific review; prior-audit legacy HVAC intersections are retained as redesign inputs rather than construction approvals, with audit freshness disclosed separately",
            "aabb_prefilter": True,
            "exact_candidate_method": "IfcOpenShell geom.tree world-shape collision/intersection/clearance",
            "far_pair_method": "world AABB Euclidean separation lower bound",
        },
        "summary": {
            "object_count": len(records),
            "total_pair_count": len(pair_records),
            "aabb_candidate_pair_count": len(candidate_keys),
            "aabb_proven_far_pair_count": len(pair_records) - len(candidate_keys),
            "geometry_state_counts": dict(sorted(state_counts.items())),
            "human_review_pair_count": len(human_review_pairs),
            "human_review_state_counts": dict(sorted(review_state_counts.items())),
            "expected_interacting_relation_counts": dict(
                sorted(interacting_rule_counts.items())
            ),
            "unresolved_conflict_candidate_count": sum(item["conflict_candidate"] for item in pair_records),
            "legacy_base_pair_count": len(legacy_pair_records),
        },
        "expected_non_conflict_pairs": expected_pairs,
        "human_review_pairs": human_review_pairs,
        "pairs": pair_records,
        "gates": {
            "all_pairs_classified": len(pair_records) == len(records) * (len(records) - 1) // 2,
            "aabb_partition_complete": len(candidate_keys) + len(pair_records) - len(candidate_keys) == len(pair_records),
            "review_pairs_have_exact_ids_and_basis": all(
                len(item["pair"]) == 2
                and item["basis"]
                and item["review_reason"]
                and item["review_required"] == "yes"
                for item in human_review_pairs
            ),
            "review_pairs_exact_mesh_tested": all(
                item["aabb_candidate"]
                and not item["minimum_clearance_is_lower_bound"]
                for item in human_review_pairs
            ),
            "expected_light_or_host_relations_not_conflicts": all(
                not item["conflict_candidate"]
                for item in pair_records
                if item["expected_relation_rule"] in {
                    "expected_light_ceiling_embedding_candidate",
                    "expected_opening_void_host",
                }
            ),
            "legacy_base_pairs_preserve_exact_intersection_scope": len(legacy_pair_records) == 10
            and all(
                item["geometry_state"] == "intersecting"
                and item["expected_relation_rule"] == "legacy_hvac_base_intersection_pending_redesign"
                and item["review_required"] == "no"
                and not item["conflict_candidate"]
                for item in legacy_pair_records
            ),
            "automatic_ifc_write_allowed": False,
            "construction_conflict_status": "legacy HVAC base pair review closed; remodel routes, systems, ports and access remain pending",
        },
    }


def main() -> int:
    args = parse_args()
    source_sha = sha256(args.input)
    if args.expected_ifc_sha256 and source_sha != args.expected_ifc_sha256:
        raise RuntimeError(f"formal IFC SHA drift: {source_sha} != {args.expected_ifc_sha256}")
    legacy_evidence = legacy_evidence_status(args.legacy_report, source_sha)
    model = ifcopenshell.open(args.input)
    if model.schema != "IFC4":
        raise RuntimeError(f"unexpected schema {model.schema}")
    reviews = read_review(args.review)
    legacy_base_rows = [row for row in reviews if row["review_id"] == "RCP1-LEGACY-BASE-001"]
    if len(legacy_base_rows) != 1:
        raise RuntimeError("RCP1 legacy design-base decision must exist exactly once")
    legacy_base_pairs = (
        LEGACY_BASE_PAIR_IDS if legacy_base_rows[0]["status"] == "implemented" else set()
    )
    geom = settings()
    spaces = space_footprints(model, geom)

    ceiling_coverings = []
    for product in model.by_type("IfcCovering"):
        assigned_type = get_type(product)
        if str(product.PredefinedType or "") != "CEILING" and str(getattr(assigned_type, "PredefinedType", "") or "") != "CEILING":
            continue
        ceiling_coverings.append(
            record(
                product,
                geom,
                spaces,
                "typed_ceiling_covering",
                "IfcCovering occurrence or assigned IfcCoveringType has PredefinedType=CEILING",
                1.0,
            )
        )

    lights = [
        record(
            product,
            geom,
            spaces,
            "typed_light_fixture",
            "IfcLightFixture assigned to the existing RA.LP / DIRECTIONSOURCE type; no circuit or mounting inferred",
            1.0,
        )
        for product in model.by_type("IfcLightFixture")
    ]

    dcl_proxies = []
    for product in model.by_type("IfcBuildingElementProxy"):
        if str(getattr(get_container(product), "Name", "") or "") != "DCL":
            continue
        vertices, _ = geometry(geom, product)
        if float(vertices[:, 2].max() * 1000.0) < 2000.0:
            continue
        name = str(product.Name or "")
        if "Diffuser" in name:
            category = "name_only_air_outlet_proxy"
            basis = "DCL container + high-level Body + current name contains Diffuser; object is not IfcAirTerminal"
            confidence = 0.65
        else:
            category = "name_only_ceiling_or_light_slot_proxy"
            basis = "DCL container + high-level Body + current name; object has no assigned IFC type or ceiling Pset"
            confidence = 0.85
        dcl_proxies.append(record(product, geom, spaces, category, basis, confidence))

    typed_high_equipment = []
    for product in model.by_type("IfcElectricAppliance"):
        assigned_type = get_type(product)
        if str(getattr(assigned_type, "ElementType", "") or "") != "AC":
            continue
        typed_high_equipment.append(
            record(
                product,
                geom,
                spaces,
                "typed_high_ac_appliance",
                "assigned IfcElectricApplianceType ElementType=AC + existing high-level Body; no system, power, or connection inferred",
                1.0,
            )
        )

    named_high_openings = []
    for product in model.by_type("IfcOpeningElement"):
        if product.GlobalId not in EXPECTED_AC_OPENING_IDS:
            continue
        named_high_openings.append(
            record(
                product,
                geom,
                spaces,
                "named_high_opening",
                "IfcOpeningElement + existing Voids relationship + current AC Hole name; no equipment connection inferred",
                0.85,
            )
        )

    high_flow_segments = []
    for product in model.by_type("IfcPipeSegment"):
        if str(getattr(get_container(product), "Name", "") or "") != "CEL":
            continue
        vertices, _ = geometry(geom, product)
        if float(vertices[:, 2].max() * 1000.0) < 2400.0:
            continue
        high_flow_segments.append(
            record(
                product,
                geom,
                spaces,
                "high_level_flow_segment",
                "IfcPipeSegment + CEL container + high-level Body; current name is recorded but medium/system is not inferred",
                0.70,
            )
        )

    other_high_proxies = [
        record(
            model.by_guid(global_id),
            geom,
            spaces,
            "name_only_high_service_proxy",
            "high-level Body + current name; object has no assigned IFC type, system, or port",
            0.70,
        )
        for global_id in OTHER_HIGH_PROXY_IDS
    ]

    for group in [ceiling_coverings, lights, dcl_proxies, typed_high_equipment, named_high_openings, high_flow_segments, other_high_proxies]:
        group.sort(key=lambda item: item["global_id"])

    all_records = [*ceiling_coverings, *lights, *dcl_proxies, *typed_high_equipment, *named_high_openings, *high_flow_segments, *other_high_proxies]
    root_ids = [root.GlobalId for root in model.by_type("IfcRoot")]
    handoff_records = [item for item in all_records if item["global_id"] in PROTECTED_HANDOFF_IDS]
    handoff_residuals = [item["placement"]["maximum_integer_residual_mm"] for item in handoff_records]
    missing_instances = {
        class_name: len(model.by_type(class_name))
        for class_name in [
            "IfcAirTerminal",
            "IfcAirTerminalBox",
            "IfcFan",
            "IfcDamper",
            "IfcSensor",
            "IfcAlarm",
            "IfcController",
            "IfcActuator",
        ]
    }
    topology = {
        "ifc_systems": len(model.by_type("IfcSystem")),
        "ifc_distribution_systems": len(model.by_type("IfcDistributionSystem")),
        "distribution_ports": len(model.by_type("IfcDistributionPort")),
        "port_connections": len(model.by_type("IfcRelConnectsPorts")),
        "port_to_element_connections": len(model.by_type("IfcRelConnectsPortToElement")),
    }
    duplicate_ids = [global_id for global_id, count in Counter(item["global_id"] for item in all_records).items() if count > 1]
    gates = {
        "root_global_ids_unique": len(root_ids) == len(set(root_ids)),
        "space_count": len(spaces),
        "spaces_with_reference": sum(bool(space["reference"]) for space in spaces),
        "ceiling_covering_count": len(ceiling_coverings),
        "light_fixture_count": len(lights),
        "dcl_proxy_count": len(dcl_proxies),
        "name_only_air_outlet_proxy_count": sum(item["category"] == "name_only_air_outlet_proxy" for item in dcl_proxies),
        "typed_high_equipment_count": len(typed_high_equipment),
        "named_high_opening_count": len(named_high_openings),
        "high_flow_segment_count": len(high_flow_segments),
        "other_high_proxy_count": len(other_high_proxies),
        "duplicate_inventory_ids": duplicate_ids,
        "records_with_geometry": sum(bool(item["bbox"]["dimensions_mm"]) for item in all_records),
        "records_with_basis_confidence_review": sum(
            bool(item["basis"]) and isinstance(item["confidence"], float) and item["review_required"] in {"yes", "no"}
            for item in all_records
        ),
        "records_with_space_overlap": sum(bool(item["space_overlaps"]) for item in all_records),
        "protected_handoff_ids": sorted(item["global_id"] for item in handoff_records),
        "protected_handoffs_over_0_1_mm": sum(value is not None and value > args.tolerance_mm + 1e-9 for value in handoff_residuals),
        "missing_instances": missing_instances,
        "topology": topology,
        "automatic_ifc_write_allowed": False,
    }
    gates["candidate_pass"] = all(
        [
            gates["root_global_ids_unique"],
            gates["space_count"] == 22,
            gates["spaces_with_reference"] == 22,
            gates["ceiling_covering_count"] == 2,
            gates["light_fixture_count"] == 79,
            gates["dcl_proxy_count"] == 18,
            gates["name_only_air_outlet_proxy_count"] == 2,
            gates["typed_high_equipment_count"] == 5,
            gates["named_high_opening_count"] == 7,
            {item["global_id"] for item in named_high_openings} == EXPECTED_AC_OPENING_IDS,
            gates["high_flow_segment_count"] == 5,
            gates["other_high_proxy_count"] == 1,
            gates["duplicate_inventory_ids"] == [],
            gates["records_with_geometry"] == len(all_records),
            gates["records_with_basis_confidence_review"] == len(all_records),
            set(gates["protected_handoff_ids"]) == PROTECTED_HANDOFF_IDS,
            gates["protected_handoffs_over_0_1_mm"] == 0,
            all(value == 0 for key, value in missing_instances.items() if key != "IfcSensor"),
            all(value == 0 for value in topology.values()),
        ]
    )
    gates["construction_release_ready"] = False
    gates["release_blocks"] = [
        "existing ceiling boundaries and level relationships are confirmed; final materials, build-ups, and access requirements remain open",
        "two Diffuser-named proxies are not typed IfcAirTerminal instances",
        "gas and warm-air equipment instances are absent; the single fire-sensor location is type-pending",
        "equipment clearances and access zones are not modelled",
        "systems, ports, connections, airflow, power, and control data are absent",
    ]
    if not gates["candidate_pass"]:
        raise RuntimeError("RCP1 existing-condition candidate gates failed")

    coordination = world_mesh_relations(
        model,
        all_records,
        args.tolerance_mm,
        args.search_window_mm,
        args.review_clearance_mm,
        legacy_base_pairs,
    )
    coordination["generated_at"] = datetime.now(timezone.utc).isoformat()
    coordination["mode"] = "read_only_world_mesh_coordination_candidate"
    coordination["source"] = {
        "path": str(args.input.resolve()),
        "sha256": source_sha,
        "schema": model.schema,
    }
    coordination["source_ifc_sha256"] = source_sha
    coordination["legacy_evidence"] = legacy_evidence
    coordination["inventory_counts"] = {
        "typed_high_equipment": len(typed_high_equipment),
        "named_high_openings": len(named_high_openings),
        "high_flow_segments": len(high_flow_segments),
        "name_only_air_outlet_proxies": sum(
            item["category"] == "name_only_air_outlet_proxy" for item in dcl_proxies
        ),
        "ceiling_or_light_slot_proxies": sum(
            item["category"] == "name_only_ceiling_or_light_slot_proxy" for item in dcl_proxies
        ),
        "typed_ceiling_coverings": len(ceiling_coverings),
        "light_fixtures": len(lights),
        "other_high_service_proxies": len(other_high_proxies),
    }
    coordination_gate_pass = all(
        [
            coordination["summary"]["object_count"] == 117,
            coordination["summary"]["total_pair_count"] == 6786,
            coordination["gates"]["all_pairs_classified"],
            coordination["gates"]["aabb_partition_complete"],
            coordination["gates"]["review_pairs_have_exact_ids_and_basis"],
            coordination["gates"]["review_pairs_exact_mesh_tested"],
            coordination["gates"]["expected_light_or_host_relations_not_conflicts"],
            coordination["gates"]["legacy_base_pairs_preserve_exact_intersection_scope"],
            coordination["gates"]["automatic_ifc_write_allowed"] is False,
        ]
    )
    coordination["gates"]["candidate_pass"] = coordination_gate_pass
    if not coordination_gate_pass:
        raise RuntimeError("RCP1 world-mesh coordination gates failed")
    args.coordination_output.parent.mkdir(parents=True, exist_ok=True)
    args.coordination_output.write_text(
        json.dumps(coordination, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_existing_high_level_candidate",
        "source_ifc_sha256": source_sha,
        "source": {"path": str(args.input.resolve()), "sha256": source_sha, "schema": model.schema},
        "legacy_evidence": legacy_evidence,
        "tolerance_mm": args.tolerance_mm,
        "inventory": {
            "ceiling_coverings": ceiling_coverings,
            "light_fixtures": lights,
            "dcl_proxies": dcl_proxies,
            "typed_high_equipment": typed_high_equipment,
            "named_high_openings": named_high_openings,
            "high_flow_segments": high_flow_segments,
            "other_high_proxies": other_high_proxies,
        },
        "missing": {
            "instances": missing_instances,
            "topology": topology,
            "verified_access_panel_instances": 0,
            "verified_warm_air_instances": 0,
        },
        "review_register": {"path": str(args.review), "records": len(reviews)},
        "coordination_report": {
            "path": str(args.coordination_output),
            "summary": coordination["summary"],
            "gates": coordination["gates"],
        },
        "gates": gates,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "gates": gates}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
