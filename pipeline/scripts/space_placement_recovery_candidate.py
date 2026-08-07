#!/usr/bin/env python3
"""Recover collapsed IfcSpace placements from a verified IFC baseline."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api.geometry
import ifcopenshell.util.element
import ifcopenshell.util.placement
import numpy as np

from geometry_alignment_audit import geometry_difference_audit, sha256


def space_matrix(space: ifcopenshell.entity_instance) -> np.ndarray:
    if not space.ObjectPlacement or not space.ObjectPlacement.is_a("IfcLocalPlacement"):
        raise RuntimeError(f"{space.GlobalId} does not have an IfcLocalPlacement")
    return np.array(
        ifcopenshell.util.placement.get_local_placement(space.ObjectPlacement),
        dtype=float,
    )


def placement_evidence(
    baseline_matrix: np.ndarray, tolerance_mm: float
) -> tuple[bool, str]:
    translation = baseline_matrix[:3, 3]
    integer_delta = np.round(translation) - translation
    rotation = baseline_matrix[:3, :3]
    rotation_delta = float(np.max(np.abs(rotation - np.eye(3))))
    accepted = bool(
        np.max(np.abs(integer_delta)) <= tolerance_mm
        and rotation_delta <= 1e-9
    )
    return (
        accepted,
        "baseline placement has an integer-millimetre translation and identity rotation"
        if accepted
        else "baseline placement is not an integer-millimetre identity-rotation target",
    )


def space_semantic_signature(space: ifcopenshell.entity_instance) -> dict[str, Any]:
    return {
        "name": space.Name,
        "long_name": space.LongName,
        "object_type": space.ObjectType,
        "composition_type": space.CompositionType,
        "psets": ifcopenshell.util.element.get_psets(space),
    }


def apply_recovery(
    candidate: ifcopenshell.file,
    baseline: ifcopenshell.file,
    tolerance_mm: float,
) -> list[dict[str, Any]]:
    candidate_spaces = {space.GlobalId: space for space in candidate.by_type("IfcSpace")}
    baseline_spaces = {space.GlobalId: space for space in baseline.by_type("IfcSpace")}
    if set(candidate_spaces) != set(baseline_spaces):
        raise RuntimeError("Current and baseline IfcSpace GlobalId sets differ")
    results = []
    for global_id in sorted(candidate_spaces):
        space = candidate_spaces[global_id]
        baseline_space = baseline_spaces[global_id]
        before = space_matrix(space)
        target = space_matrix(baseline_space)
        accepted, basis = placement_evidence(target, tolerance_mm)
        if not accepted:
            raise RuntimeError(f"Unsafe baseline target for {global_id}: {basis}")
        ifcopenshell.api.geometry.edit_object_placement(
            candidate,
            product=space,
            matrix=target.copy(),
            is_si=False,
            should_transform_children=False,
        )
        after = space_matrix(space)
        results.append(
            {
                "global_id": global_id,
                "long_name": space.LongName,
                "before_translation_mm": before[:3, 3].tolist(),
                "target_translation_mm": target[:3, 3].tolist(),
                "after_translation_mm": after[:3, 3].tolist(),
                "translation_matches_target": bool(
                    np.max(np.abs(after[:3, 3] - target[:3, 3])) <= tolerance_mm
                ),
                "basis": basis,
                "confidence": 1.0,
                "review_required": False,
            }
        )
    return results


def represented_product_mesh_digests(
    model: ifcopenshell.file,
    global_ids: list[str],
) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    """Hash exact world vertices/faces in linear time for regression safety."""

    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    records: dict[str, dict[str, Any]] = {}
    failures: dict[str, str] = {}
    for global_id in global_ids:
        product = model.by_guid(global_id)
        try:
            shape = ifcopenshell.geom.create_shape(settings, product)
            vertices = np.asarray(shape.geometry.verts, dtype=np.float64)
            faces = np.asarray(shape.geometry.faces, dtype=np.int64)
            digest = hashlib.sha256()
            digest.update(np.round(vertices, 12).astype("<f8", copy=False).tobytes())
            digest.update(faces.astype("<i8", copy=False).tobytes())
            records[global_id] = {
                "ifc_class": product.is_a(),
                "vertex_count": int(vertices.size // 3),
                "triangle_count": int(faces.size // 3),
                "world_mesh_sha256": digest.hexdigest(),
            }
        except Exception as error:  # pragma: no cover - model-specific failures
            failures[global_id] = str(error)
    return records, failures


def forward_subgraph_digests(
    model: ifcopenshell.file,
    global_ids: list[str],
) -> dict[str, str]:
    """Hash forward IFC attributes for products that have no renderable shape."""

    result = {}
    for global_id in global_ids:
        product = model.by_guid(global_id)
        digest = hashlib.sha256()
        for entity in sorted(model.traverse(product), key=lambda item: item.id()):
            digest.update(f"{entity.id()}:{entity}\n".encode("utf-8"))
        result[global_id] = digest.hexdigest()
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--baseline-ifc", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--expected-count", type=int, default=22)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.expected_count <= 0:
        raise SystemExit("--expected-count must be positive")
    if args.tolerance_mm <= 0.0:
        raise SystemExit("--tolerance-mm must be positive")
    source_path = args.input.resolve()
    baseline_path = args.baseline_ifc.resolve()
    source = ifcopenshell.open(source_path)
    baseline = ifcopenshell.open(baseline_path)
    candidate = ifcopenshell.open(source_path)
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    source_signatures = {
        space.GlobalId: space_semantic_signature(space)
        for space in source.by_type("IfcSpace")
    }
    results = apply_recovery(candidate, baseline, args.tolerance_mm)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)

    space_ids = sorted(result["global_id"] for result in results)
    candidate_vs_baseline = geometry_difference_audit(
        candidate,
        baseline,
        str(baseline_path),
        tolerance_mm=args.tolerance_mm,
        classes=[],
        global_ids=space_ids,
    )
    non_space_represented_product_ids = sorted(
        product.GlobalId
        for product in source.by_type("IfcProduct")
        if getattr(product, "GlobalId", None)
        and product.ObjectPlacement
        and product.Representation
        and not product.is_a("IfcSpace")
    )
    source_mesh_digests, source_mesh_failures = represented_product_mesh_digests(
        source, non_space_represented_product_ids
    )
    candidate_mesh_digests, candidate_mesh_failures = represented_product_mesh_digests(
        candidate, non_space_represented_product_ids
    )
    non_space_mesh_mismatch_ids = sorted(
        global_id
        for global_id in non_space_represented_product_ids
        if source_mesh_digests.get(global_id) != candidate_mesh_digests.get(global_id)
    )
    failed_mesh_ids = sorted(
        set(source_mesh_failures) | set(candidate_mesh_failures)
    )
    source_failed_subgraphs = forward_subgraph_digests(source, failed_mesh_ids)
    candidate_failed_subgraphs = forward_subgraph_digests(candidate, failed_mesh_ids)
    failed_subgraph_mismatch_ids = sorted(
        global_id
        for global_id in failed_mesh_ids
        if source_failed_subgraphs.get(global_id)
        != candidate_failed_subgraphs.get(global_id)
    )
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    candidate_signatures = {
        space.GlobalId: space_semantic_signature(space)
        for space in candidate.by_type("IfcSpace")
    }
    gates = {
        "expected_count": args.expected_count,
        "recovered_count": len(results),
        "all_translations_match_target": all(
            result["translation_matches_target"] for result in results
        ),
        "candidate_space_geometry_matches_baseline": candidate_vs_baseline[
            "over_tolerance"
        ]
        == 0,
        "non_space_represented_product_count": len(
            non_space_represented_product_ids
        ),
        "non_space_mesh_mismatch_ids": non_space_mesh_mismatch_ids,
        "source_mesh_failure_count": len(source_mesh_failures),
        "candidate_mesh_failure_count": len(candidate_mesh_failures),
        "mesh_failure_id_sets_equal": set(source_mesh_failures)
        == set(candidate_mesh_failures),
        "failed_subgraph_mismatch_ids": failed_subgraph_mismatch_ids,
        "all_non_space_world_meshes_equal": (
            not non_space_mesh_mismatch_ids
            and set(source_mesh_failures) == set(candidate_mesh_failures)
            and not failed_subgraph_mismatch_ids
        ),
        "space_semantics_equal": source_signatures == candidate_signatures,
        "schema_equal": source.schema == candidate.schema == baseline.schema,
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
        "source_entity_count": len(list(source)),
        "candidate_entity_count": len(list(candidate)),
    }
    gates["pass"] = (
        len(results) == args.expected_count
        and gates["all_translations_match_target"]
        and gates["candidate_space_geometry_matches_baseline"]
        and gates["all_non_space_world_meshes_equal"]
        and gates["space_semantics_equal"]
        and gates["schema_equal"]
        and gates["root_global_ids_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-space-placement-recovery-candidate",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "baseline": {
            "path": str(baseline_path),
            "sha256": sha256(baseline_path),
            "schema": baseline.schema,
        },
        "candidate": {
            "path": str(args.output),
            "sha256": sha256(args.output),
            "schema": candidate.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "results": results,
        "candidate_vs_baseline_space_geometry": candidate_vs_baseline,
        "non_space_world_mesh_digest_comparison": {
            "product_count": len(non_space_represented_product_ids),
            "mismatch_ids": non_space_mesh_mismatch_ids,
            "source_failures": source_mesh_failures,
            "candidate_failures": candidate_mesh_failures,
            "failed_subgraph_mismatch_ids": failed_subgraph_mismatch_ids,
        },
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
