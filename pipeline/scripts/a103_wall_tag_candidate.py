#!/usr/bin/env python3
"""Build a geometry-neutral A-103 wall-tag candidate for one-rule review."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.util.placement

import p0_ids_metadata_candidate as p0


ROOT = Path(__file__).resolve().parents[2]
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
REGISTER = ROOT / "pipeline/decisions/a103-wall-status-review.csv"
EXPECTED_SOURCE_IDS_PASS = 499
EXPECTED_CANDIDATE_IDS_PASS = 587


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--output", type=Path, default=ROOT / "build/candidates/2504-GBTB-a103-wall-tags.ifc")
    parser.add_argument("--report", type=Path, default=ROOT / "build/a103/a103-wall-tag-candidate.json")
    parser.add_argument("--expected-ifc-sha256")
    return parser.parse_args()


def centre(row: dict[str, str]) -> tuple[float, float]:
    return (
        (float(row["min_x_mm"]) + float(row["max_x_mm"])) / 2.0,
        (float(row["min_y_mm"]) + float(row["max_y_mm"])) / 2.0,
    )


def tag_records(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for status, prefix, expected_count in (("EXISTING", "EW", 84), ("NEW", "NW", 4)):
        selected = [row for row in rows if row["current_status"] == status]
        if len(selected) != expected_count:
            raise RuntimeError(f"expected {expected_count} {status} walls, found {len(selected)}")
        selected.sort(key=lambda row: (-centre(row)[1], centre(row)[0], row["global_id"]))
        for index, row in enumerate(selected, 1):
            x, y = centre(row)
            records.append(
                {
                    "global_id": row["global_id"],
                    "status": status,
                    "candidate_tag": f"{prefix}{index:02d}",
                    "centre_x_mm": x,
                    "centre_y_mm": y,
                    "ordering": "north_to_south_then_west_to_east_then_global_id",
                }
            )
    return records


def main() -> int:
    args = parse_args()
    source_path = args.input.resolve()
    output_path = args.output.resolve()
    report_path = args.report.resolve()
    if source_path == output_path:
        raise RuntimeError("candidate output must not overwrite the formal IFC")
    source_hash = p0.sha256(source_path)
    if args.expected_ifc_sha256 and source_hash != args.expected_ifc_sha256:
        raise RuntimeError("formal IFC hash changed from the caller-frozen value")

    source = ifcopenshell.open(source_path)
    rows = p0.read_csv(REGISTER)
    if len(rows) != 88 or len({row["global_id"] for row in rows}) != 88:
        raise RuntimeError("A-103 register must contain 88 unique walls")
    if {row["source_ifc_sha256"] for row in rows} != {source_hash}:
        raise RuntimeError("A-103 register is stale against the formal IFC")
    records = tag_records(rows)
    tags = [record["candidate_tag"] for record in records]
    if len(tags) != len(set(tags)):
        raise RuntimeError("wall-tag candidate contains duplicates")

    target_ids = {record["global_id"] for record in records}
    final_wall_ids = {
        wall.GlobalId
        for wall in source.by_type("IfcWall")
        if str(
            ifcopenshell.util.element.get_psets(wall)
            .get("Pset_WallCommon", {})
            .get("Status")
            or ""
        ).strip().upper()
        not in {"DEMOLISH", "DEMOLISHED"}
    }
    if target_ids != final_wall_ids:
        raise RuntimeError("tag targets do not exactly equal the final-built wall set")
    occupied_tags = {
        str(element.Tag)
        for element in source.by_type("IfcElement")
        if element.GlobalId not in target_ids and getattr(element, "Tag", None)
    }
    collisions = sorted(set(tags) & occupied_tags)
    if collisions:
        raise RuntimeError(f"candidate wall tags collide with existing element tags: {collisions}")

    source_ids = p0.run_ids(source_path)
    source_relationships = p0.relationship_fingerprint(source)
    source_products = {
        product.GlobalId: product
        for product in source.by_type("IfcProduct")
        if not (product.is_a("IfcAnnotation") and product.ObjectType == "DRAWING")
    }
    source_graphs = {global_id: p0.physical_graph(product) for global_id, product in source_products.items()}
    source_matrices = {
        global_id: ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement)
        for global_id, product in source_products.items()
    }

    candidate = ifcopenshell.open(source_path)
    for record in records:
        candidate.by_guid(record["global_id"]).Tag = record["candidate_tag"]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(output_path)
    if p0.sha256(source_path) != source_hash:
        raise RuntimeError("formal IFC changed while building the wall-tag candidate")

    candidate = ifcopenshell.open(output_path)
    candidate_products = {product.GlobalId: product for product in candidate.by_type("IfcProduct")}
    graph_bad: list[str] = []
    placement_bad: list[dict[str, Any]] = []
    for global_id, source_product in sorted(source_products.items()):
        product = candidate_products.get(global_id)
        if product is None or source_graphs[global_id] != p0.physical_graph(product):
            graph_bad.append(global_id)
            continue
        delta = p0.matrix_max_delta(
            source_matrices[global_id],
            ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
        )
        if delta != 0.0:
            placement_bad.append({"global_id": global_id, "maximum_delta": delta})

    candidate_ids = p0.run_ids(output_path)
    candidate_tags_correct = all(
        candidate.by_guid(record["global_id"]).Tag == record["candidate_tag"]
        for record in records
    )
    gates = {
        "source_ifc_hash_preserved": p0.sha256(source_path) == source_hash,
        "existing_wall_tag_count": sum(record["status"] == "EXISTING" for record in records),
        "new_wall_tag_count": sum(record["status"] == "NEW" for record in records),
        "candidate_tags_unique": len(tags) == len(set(tags)),
        "candidate_tags_do_not_collide": not collisions,
        "candidate_tags_correct": candidate_tags_correct,
        "protected_product_count": len(source_products),
        "protected_products_geometry_exact": not graph_bad,
        "maximum_world_geometry_delta_mm": 0.0 if not graph_bad and not placement_bad else math.inf,
        "maximum_placement_matrix_delta": max((record["maximum_delta"] for record in placement_bad), default=0.0),
        "fills_voids_relationships_unchanged": p0.relationship_fingerprint(candidate) == source_relationships,
        "source_ids_pass": source_ids["total_checks_pass"],
        "candidate_ids_pass": candidate_ids["total_checks_pass"],
        "candidate_ids_fail": candidate_ids["total_checks_fail"],
        "formal_ifc_write_allowed": False,
    }
    passed = (
        gates["source_ifc_hash_preserved"]
        and gates["existing_wall_tag_count"] == 84
        and gates["new_wall_tag_count"] == 4
        and gates["candidate_tags_unique"]
        and gates["candidate_tags_do_not_collide"]
        and gates["candidate_tags_correct"]
        and gates["protected_products_geometry_exact"]
        and gates["maximum_placement_matrix_delta"] == 0.0
        and gates["fills_voids_relationships_unchanged"]
        and gates["source_ids_pass"] == EXPECTED_SOURCE_IDS_PASS
        and gates["candidate_ids_pass"] == EXPECTED_CANDIDATE_IDS_PASS
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_a103_wall_tag_candidate",
        "source_ifc_sha256": source_hash,
        "scheme": {
            "existing": "EW01-EW84",
            "new": "NW01-NW04",
            "ordering": "north_to_south_then_west_to_east_then_global_id",
            "review_required": True,
        },
        "records": records,
        "gates": gates,
        "pass": passed,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"candidate": str(output_path), "report": str(report_path), "gates": gates, "pass": passed}, ensure_ascii=False))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
