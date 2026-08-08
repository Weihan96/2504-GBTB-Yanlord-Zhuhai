#!/usr/bin/env python3
"""Generate read-only PLUM existing-object candidate registries.

The script records observable IFC objects. It deliberately does not infer
cold/hot-water connections, drainage connectivity, rough-in points, pipe
systems, or installation semantics from tessellated shells.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.util.placement
import numpy as np

from flow_segment_centerline_audit import mesh_components
from geometry_alignment_audit import geometry_settings, world_mesh_mm


EXPECTED_IFC_SHA256 = "c7295688003f3f36775a25f6adc2c9878e203c52c980f8e46faed66a3537c4a8"
PVC110_IDS = ("178mqyyzzFowLcbXcH6prO", "0bfVg4Ys1CevZs$qxhkXTo")
ASSEMBLY_IDS = ("1FICW5lCjEIQPvN1oi4cE2", "0a3r2aWdbFi9OlDGcQ00CE", "1fEFh83259DRb927ZP9yW7")
DIRECT_DRAIN_IDS = (
    "3AehGtDgv4cwDJtCqrtkrC", "1J6buM4aD3LO4vfb0AtPRz",
    "1vFwjs$QnFcOPrrmIAjNrG", "178mqyyzzFowLcbXcH6prO",
    "2hIKyAC2X3RhKEThQtEyDv", "1NwqlTj$T8DAd8wveTVU0J",
    "0bfVg4Ys1CevZs$qxhkXTo",
)
DRAIN_PIPE_TYPE_NAMES = {"TEE01", "ELB01", "PIP01", "DRA01"}
EXPECTED_REVIEW_IDS = {
    "PLUM-EXISTING-001", "PLUM-DATA-001", "PLUM-P201-001",
    "PLUM-P202-001", "PLUM-PVC110-001", "PLUM-ANCHOR-001",
    "PLUM-SITE-001", "PLUM-HOTWATER-001", "PLUM-KITCHEN-DRAIN-001",
}
GEOMETRY_SETTINGS = geometry_settings()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument("--review", type=Path, default=Path("pipeline/decisions/plum-existing-review.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("build/plum"))
    parser.add_argument("--geometry-baseline-ref", default="4700894")
    parser.add_argument("--expected-ifc-sha256", default=EXPECTED_IFC_SHA256)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assigned_type(product: Any) -> Any | None:
    return ifcopenshell.util.element.get_type(product)


def container_record(product: Any) -> dict[str, Any] | None:
    container = ifcopenshell.util.element.get_container(product)
    if not container:
        return None
    return {
        "global_id": getattr(container, "GlobalId", None),
        "ifc_class": container.is_a(),
        "name": container.Name,
        "long_name": getattr(container, "LongName", None),
    }


def origin_mm(product: Any) -> list[float] | None:
    if not product.ObjectPlacement:
        return None
    matrix = ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement)
    return [float(value * 1000.0) for value in matrix[:3, 3]]


def descendant_leaves(product: Any) -> list[Any]:
    children = [
        child
        for relation in getattr(product, "IsDecomposedBy", ()) or ()
        for child in relation.RelatedObjects
    ]
    if not children:
        return [product]
    leaves: list[Any] = []
    for child in children:
        leaves.extend(descendant_leaves(child))
    return leaves


def geometry_record(product: Any) -> dict[str, Any]:
    leaves = descendant_leaves(product)
    represented_leaf_ids = [
        leaf.GlobalId for leaf in leaves if getattr(leaf, "Representation", None)
    ]
    return {
        "has_world_geometry": bool(represented_leaf_ids),
        "represented_leaf_global_ids": represented_leaf_ids,
        "bbox_min_mm": None,
        "bbox_max_mm": None,
        "bbox_size_mm": None,
        "location_candidate_mm": origin_mm(product),
        "location_basis": "world ObjectPlacement for existing-object registration only; not a connector or rough-in point",
    }


def product_record(product: Any) -> dict[str, Any]:
    product_type = assigned_type(product)
    return {
        "global_id": product.GlobalId,
        "ifc_class": product.is_a(),
        "name": product.Name,
        "predefined_type": getattr(product, "PredefinedType", None),
        "type_global_id": getattr(product_type, "GlobalId", None),
        "type_name": getattr(product_type, "Name", None),
        "container": container_record(product),
        "object_origin_mm": origin_mm(product),
        "geometry": geometry_record(product),
    }


def drainage_products(model: Any) -> list[Any]:
    products: list[Any] = []
    for product in model.by_type("IfcFlowSegment"):
        product_type = assigned_type(product)
        type_name = getattr(product_type, "Name", None)
        if product.GlobalId in DIRECT_DRAIN_IDS or type_name in DRAIN_PIPE_TYPE_NAMES:
            products.append(product)
    return sorted(products, key=lambda product: product.GlobalId)


def load_baseline(root: Path, ref: str, ifc_relative: Path) -> Any:
    if ifc_relative.is_absolute():
        raise RuntimeError("git geometry baseline requires a repository-relative IFC path")
    payload = subprocess.check_output(
        ["git", "show", f"{ref}:{ifc_relative.as_posix()}"], cwd=root
    )
    return ifcopenshell.file.from_string(payload.decode("utf-8"))


def pvc110_record(current: Any, baseline: Any, global_id: str) -> dict[str, Any]:
    current_product = current.by_guid(global_id)
    baseline_product = baseline.by_guid(global_id)
    current_vertices, current_faces = world_mesh_mm(GEOMETRY_SETTINGS, current_product)
    baseline_vertices, _ = world_mesh_mm(GEOMETRY_SETTINGS, baseline_product)
    current_points = np.asarray(current_vertices, dtype=float)
    baseline_points = np.asarray(baseline_vertices, dtype=float)
    components = mesh_components(current_vertices, current_faces)
    compatible = current_points.shape == baseline_points.shape
    maximum_change = (
        float(np.max(np.linalg.norm(current_points - baseline_points, axis=1)))
        if compatible else None
    )
    return {
        "global_id": global_id,
        "ifc_product_preserved": current_product.is_a() == "IfcFlowSegment",
        "mesh_component_count": len(components),
        "component_vertex_counts": [len(component) for component in components],
        "baseline_vertex_count": len(baseline_points),
        "current_vertex_count": len(current_points),
        "vertex_order_compatible": compatible,
        "maximum_world_vertex_change_mm": maximum_change,
        "world_geometry_unchanged": compatible and maximum_change == 0.0,
        "formal_ifc_write_authorized": False,
    }


def read_review(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    ids = {row["issue_id"] for row in rows}
    if ids != EXPECTED_REVIEW_IDS or len(ids) != len(rows):
        raise RuntimeError("PLUM review register is incomplete or contains duplicate issue IDs")
    required = ("basis", "confidence", "review_required", "status", "stop_condition")
    if any(not row[field] for row in rows for field in required):
        raise RuntimeError("PLUM review rows require basis, confidence, review status, and stop condition")
    return rows


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    resolve = lambda path: path if path.is_absolute() else root / path
    ifc_path = resolve(args.ifc)
    review_path = resolve(args.review)
    output_dir = resolve(args.output_dir)
    ifc_sha = sha256(ifc_path)
    if ifc_sha != args.expected_ifc_sha256:
        raise RuntimeError(
            f"formal IFC hash changed: expected {args.expected_ifc_sha256}, got {ifc_sha}"
        )

    model = ifcopenshell.open(ifc_path)
    review_rows = read_review(review_path)
    sanitary = sorted(model.by_type("IfcSanitaryTerminal"), key=lambda product: product.GlobalId)
    waste = sorted(model.by_type("IfcWasteTerminal"), key=lambda product: product.GlobalId)
    assemblies = [model.by_guid(global_id) for global_id in ASSEMBLY_IDS]
    drainage = drainage_products(model)
    ports = model.by_type("IfcDistributionPort")
    systems = model.by_type("IfcSystem")
    distribution_systems = model.by_type("IfcDistributionSystem")
    connects_ports = model.by_type("IfcRelConnectsPorts")
    connects_elements = model.by_type("IfcRelConnectsPortToElement")

    if len(sanitary) != 27 or len(waste) != 3 or len(assemblies) != 3 or len(drainage) != 16:
        raise RuntimeError("PLUM observable object counts changed; review classification before regeneration")

    record_cache: dict[str, dict[str, Any]] = {}

    def record(product: Any) -> dict[str, Any]:
        if product.GlobalId not in record_cache:
            record_cache[product.GlobalId] = product_record(product)
        return record_cache[product.GlobalId]

    p201 = {
        "candidate": "P-201 service-demand endpoint register",
        "source_ifc_sha256": ifc_sha,
        "scope": "observable sanitary product locations only",
        "not_in_scope": [
            "cold-water connection inference", "hot-water connection inference",
            "pipe routing", "pipe sizing", "system connectivity", "rough-in anchor writes",
        ],
        "demand_endpoints": [
            {
                **record(product),
                "connection_requirement": "unknown",
                "is_ifc_distribution_port": False,
                "candidate_is_write_authority": False,
            }
            for product in sanitary
        ],
    }
    p202_products = sanitary + waste + assemblies + drainage
    p202 = {
        "candidate": "P-202 existing drainage and sanitary location register",
        "source_ifc_sha256": ifc_sha,
        "scope": "existing product identity, type, container, and world geometry",
        "not_in_scope": [
            "connectivity inference", "rough-in point inference", "drainage redesign",
            "PVC110 product splitting", "formal IFC writes",
        ],
        "objects": [record(product) for product in p202_products],
    }

    baseline = load_baseline(root, args.geometry_baseline_ref, args.ifc)
    pvc110 = [pvc110_record(model, baseline, global_id) for global_id in PVC110_IDS]
    missing_connectivity_declared = not any(
        (ports, systems, distribution_systems, connects_ports, connects_elements)
    )
    qa = {
        "expected_counts_pass": len(sanitary) == 27 and len(waste) == 3
        and len(assemblies) == 3 and len(drainage) == 16,
        "unique_candidate_global_ids": len({p.GlobalId for p in p202_products}) == len(p202_products),
        "distribution_data": {
            "IfcDistributionPort": len(ports),
            "IfcSystem": len(systems),
            "IfcDistributionSystem": len(distribution_systems),
            "IfcRelConnectsPorts": len(connects_ports),
            "IfcRelConnectsPortToElement": len(connects_elements),
            "status": "data_missing" if missing_connectivity_declared else "present",
            "connectivity_qa_passed": False,
        },
        "pvc110_product_count": len(pvc110),
        "pvc110_branch_count": sum(row["mesh_component_count"] for row in pvc110),
        "pvc110_world_geometry_unchanged": all(row["world_geometry_unchanged"] for row in pvc110),
        "review_register_complete": len(review_rows) == len(EXPECTED_REVIEW_IDS),
    }
    qa["candidate_registry_pass"] = all((
        qa["expected_counts_pass"], qa["unique_candidate_global_ids"],
        missing_connectivity_declared, qa["pvc110_product_count"] == 2,
        qa["pvc110_branch_count"] == 6, qa["pvc110_world_geometry_unchanged"],
        qa["review_register_complete"],
    ))
    qa["construction_release_pass"] = False

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {"path": str(ifc_path), "ifc_sha256": ifc_sha, "schema": model.schema},
        "geometry_baseline": {"git_ref": args.geometry_baseline_ref},
        "inventory": {
            "sanitary_terminal_count": len(sanitary),
            "waste_terminal_count": len(waste),
            "sanitary_assembly_count": len(assemblies),
            "recognizable_drainage_product_count": len(drainage),
            "p201_demand_endpoint_count": len(p201["demand_endpoints"]),
            "p202_existing_object_count": len(p202["objects"]),
        },
        "pvc110": pvc110,
        "qa": qa,
        "stop_conditions": [
            "water pressure, valve-room equipment space, and filtration/boosting conditions are unconfirmed",
            "hot-water route, pipe size, insulation, and waiting-time strategy are unconfirmed",
            "manufacturer rough-in points and IFC ports are absent; no placement write is authorized",
            "kitchen drain adapter location and maintainable access require site verification",
        ],
        "outputs": {
            "p201": "build/plum/p201-demand-endpoints.json",
            "p202": "build/plum/p202-existing-location-register.json",
            "review": "pipeline/decisions/plum-existing-review.csv",
        },
    }
    write_json(output_dir / "p201-demand-endpoints.json", p201)
    write_json(output_dir / "p202-existing-location-register.json", p202)
    write_json(output_dir / "plum-report.json", report)
    print(json.dumps(report["inventory"], ensure_ascii=False))
    print(json.dumps(qa, ensure_ascii=False))
    if not qa["candidate_registry_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
