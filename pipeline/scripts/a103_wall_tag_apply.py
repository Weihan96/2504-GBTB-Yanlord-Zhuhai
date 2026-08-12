#!/usr/bin/env python3
"""Apply the reviewed A-103 wall-tag scheme through an atomic, fail-closed write."""

from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.util.placement

import a103_wall_tag_candidate as candidate_builder
import p0_ids_metadata_candidate as p0


ROOT = Path(__file__).resolve().parents[2]
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
CANDIDATE_REPORT = ROOT / "build/a103/a103-wall-tag-candidate.json"
DEFAULT_REPORT = ROOT / "build/a103/a103-wall-tag-postwrite.json"
APPROVAL_TOKEN = "APPROVE-A103-EW-NW"
EXPECTED_SOURCE_IDS_PASS = 499
EXPECTED_POSTWRITE_IDS_PASS = 587


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--candidate-report", type=Path, default=CANDIDATE_REPORT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--expected-ifc-sha256", required=True)
    parser.add_argument("--approval-token", required=True)
    parser.add_argument(
        "--replace-source",
        action="store_true",
        help="Required when output is the formal input path; replacement remains atomic.",
    )
    return parser.parse_args()


def encoded(value: Any) -> Any:
    if hasattr(value, "id"):
        return ["entity", value.is_a(), value.id()]
    if isinstance(value, (tuple, list)):
        return [encoded(item) for item in value]
    return value


def direct_signature(entity: Any, ignored: set[str] | None = None) -> dict[str, Any]:
    ignored = ignored or set()
    return {
        entity.attribute_name(index): encoded(entity[index])
        for index in range(len(entity))
        if entity.attribute_name(index) not in ignored
    }


def root_signatures(model: ifcopenshell.file, target_ids: set[str]) -> dict[str, dict[str, Any]]:
    signatures: dict[str, dict[str, Any]] = {}
    for root in model.by_type("IfcRoot"):
        ignored = {"Tag"} if root.GlobalId in target_ids and root.is_a("IfcElement") else set()
        signatures[root.GlobalId] = {
            "class": root.is_a(),
            "attributes": direct_signature(root, ignored),
        }
    return signatures


def entity_counts(model: ifcopenshell.file) -> dict[str, int]:
    return {
        "entities": sum(1 for _ in model),
        "roots": len(model.by_type("IfcRoot")),
        "products": len(model.by_type("IfcProduct")),
        "walls": len(model.by_type("IfcWall")),
    }


def final_wall_ids(model: ifcopenshell.file) -> set[str]:
    return {
        wall.GlobalId
        for wall in model.by_type("IfcWall")
        if str(
            ifcopenshell.util.element.get_psets(wall)
            .get("Pset_WallCommon", {})
            .get("Status")
            or ""
        ).strip().upper()
        not in {"DEMOLISH", "DEMOLISHED"}
    }


def write_atomic(model: ifcopenshell.file, output: Path) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.stem}-", suffix=".ifc", dir=output.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        model.write(temporary)
        return temporary
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    args = parse_args()
    if args.approval_token != APPROVAL_TOKEN:
        raise RuntimeError("A-103 wall-tag write requires the exact explicit approval token")

    source_path = args.input.resolve()
    output_path = args.output.resolve()
    report_path = args.report.resolve()
    replacing_source = source_path == output_path
    if replacing_source != args.replace_source:
        raise RuntimeError("--replace-source must be supplied if and only if output equals input")

    source_hash = p0.sha256(source_path)
    if source_hash != args.expected_ifc_sha256:
        raise RuntimeError(
            f"formal IFC hash changed: expected {args.expected_ifc_sha256}, got {source_hash}"
        )
    candidate = json.loads(args.candidate_report.read_text(encoding="utf-8"))
    if not candidate.get("pass") or candidate.get("source_ifc_sha256") != source_hash:
        raise RuntimeError("A-103 wall-tag candidate is failed or stale")
    if candidate.get("scheme") != {
        "existing": "EW01-EW84",
        "new": "NW01-NW04",
        "ordering": "north_to_south_then_west_to_east_then_global_id",
        "review_required": True,
    }:
        raise RuntimeError("A-103 wall-tag candidate scheme changed after review")

    records = candidate.get("records") or []
    if len(records) != 88:
        raise RuntimeError("A-103 wall-tag candidate must contain exactly 88 records")
    rows = p0.read_csv(candidate_builder.REGISTER)
    regenerated = candidate_builder.tag_records(rows)
    if records != regenerated:
        raise RuntimeError("A-103 wall-tag candidate no longer matches the decision register")
    target_ids = {record["global_id"] for record in records}
    tags = [record["candidate_tag"] for record in records]
    if len(target_ids) != 88 or len(tags) != len(set(tags)):
        raise RuntimeError("A-103 wall-tag candidate IDs or tags are not unique")

    source = ifcopenshell.open(source_path)
    if source.schema != "IFC4" or final_wall_ids(source) != target_ids:
        raise RuntimeError("A-103 targets do not equal the formal final-built wall set")
    source_counts = entity_counts(source)
    source_relationships = p0.relationship_fingerprint(source)
    source_roots = root_signatures(source, target_ids)
    source_products = {
        product.GlobalId: product
        for product in source.by_type("IfcProduct")
        if not (product.is_a("IfcAnnotation") and product.ObjectType == "DRAWING")
    }
    source_graphs = {
        global_id: p0.physical_graph(product)
        for global_id, product in source_products.items()
    }
    source_placements = {
        global_id: ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement)
        for global_id, product in source_products.items()
    }
    source_ids = p0.run_ids(source_path)
    if source_ids["total_checks_pass"] != EXPECTED_SOURCE_IDS_PASS:
        raise RuntimeError("formal IFC IDS baseline changed before A-103 write")

    working = ifcopenshell.open(source_path)
    for record in records:
        working.by_guid(record["global_id"]).Tag = record["candidate_tag"]
    temporary = write_atomic(working, output_path)
    try:
        saved = ifcopenshell.open(temporary)
        saved_counts = entity_counts(saved)
        saved_products = {product.GlobalId: product for product in saved.by_type("IfcProduct")}
        graph_bad: list[str] = []
        placement_bad: list[dict[str, Any]] = []
        for global_id in sorted(source_products):
            product = saved_products.get(global_id)
            if product is None or source_graphs[global_id] != p0.physical_graph(product):
                graph_bad.append(global_id)
                continue
            delta = p0.matrix_max_delta(
                source_placements[global_id],
                ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
            )
            if delta != 0.0:
                placement_bad.append({"global_id": global_id, "maximum_delta": delta})

        saved_roots = root_signatures(saved, target_ids)
        tags_correct = all(
            saved.by_guid(record["global_id"]).Tag == record["candidate_tag"]
            for record in records
        )
        saved_ids = p0.run_ids(temporary)
        gates = {
            "approval_token_exact": True,
            "source_hash_frozen": p0.sha256(source_path) == source_hash,
            "candidate_report_current_and_passed": True,
            "target_count": len(target_ids),
            "tags_unique": len(tags) == len(set(tags)),
            "tags_correct_after_reload": tags_correct,
            "entity_counts_unchanged": saved_counts == source_counts,
            "all_root_attributes_except_target_tags_unchanged": saved_roots == source_roots,
            "all_product_physical_graphs_unchanged": not graph_bad,
            "maximum_placement_matrix_delta": max(
                (item["maximum_delta"] for item in placement_bad), default=0.0
            ),
            "maximum_world_geometry_delta_mm": (
                0.0 if not graph_bad and not placement_bad else math.inf
            ),
            "fills_voids_relationships_unchanged": (
                p0.relationship_fingerprint(saved) == source_relationships
            ),
            "source_ids_pass": source_ids["total_checks_pass"],
            "postwrite_ids_pass": saved_ids["total_checks_pass"],
            "postwrite_ids_fail": saved_ids["total_checks_fail"],
        }
        passed = (
            gates["source_hash_frozen"]
            and gates["tags_correct_after_reload"]
            and gates["entity_counts_unchanged"]
            and gates["all_root_attributes_except_target_tags_unchanged"]
            and gates["all_product_physical_graphs_unchanged"]
            and gates["maximum_placement_matrix_delta"] == 0.0
            and gates["fills_voids_relationships_unchanged"]
            and gates["postwrite_ids_pass"] == EXPECTED_POSTWRITE_IDS_PASS
        )
        if not passed:
            raise RuntimeError(f"A-103 postwrite gates failed: {json.dumps(gates, ensure_ascii=False)}")

        output_hash = p0.sha256(temporary)
        os.replace(temporary, output_path)
        if p0.sha256(output_path) != output_hash:
            raise RuntimeError("atomic A-103 output failed disk hash verification")
        report = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "mode": "formal_replace" if replacing_source else "validated_output",
            "source_ifc": str(source_path),
            "output_ifc": str(output_path),
            "source_ifc_sha256": source_hash,
            "output_ifc_sha256": output_hash,
            "candidate_report": str(args.candidate_report.resolve()),
            "scheme": candidate["scheme"],
            "counts_before": source_counts,
            "counts_after": saved_counts,
            "gates": gates,
            "formal_ifc_write_performed": replacing_source,
            "pass": True,
        }
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(report, ensure_ascii=False))
        return 0
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
