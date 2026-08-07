#!/usr/bin/env python3
"""Compare rigid-round and geometry-preserving appliance-origin options."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.util.placement
import numpy as np

from beam_origin_reset_candidate import point_mesh_distance, reset_product_origin
from direction_noise_candidate import all_product_geometry_difference
from furniture_group_origin_candidate import (
    classify_contact_changes,
    clearance_map,
    intersection_pairs,
)
from geometry_alignment_audit import geometry_settings, sha256, world_mesh_mm
from light_fixture_anchor_candidate import apply_targets, origin_mm, renderable_context


TARGET_TYPES = {
    "1PUCikoaP5fgiYt8sJd8$6": "AC700",
    "0UOnmuAdP1MPy6p3olwiEU": "WD01",
    "3PQOXKxgj6IftqWXFXMQXG": "OV01",
    "288GLY62v8kPPydA1lAK8W": "HD01",
}
CONTEXT_CLASSES = (
    "IfcWall",
    "IfcCovering",
    "IfcSlab",
    "IfcFurniture",
    "IfcBuildingElementProxy",
    "IfcElectricAppliance",
)


def json_default(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (set, tuple)):
        return list(value)
    raise TypeError(f"unsupported JSON value {type(value).__name__}")


def current_origin_feature(
    bbox_min: np.ndarray,
    bbox_max: np.ndarray,
    tolerance_mm: float,
) -> str | None:
    """Classify local zero only when it is a bbox vertex/edge/face/centre."""
    centre = (bbox_min + bbox_max) / 2.0
    kinds = []
    for axis in range(3):
        if abs(float(centre[axis])) <= tolerance_mm:
            kinds.append("centre")
        elif min(
            abs(float(bbox_min[axis])), abs(float(bbox_max[axis]))
        ) <= tolerance_mm:
            kinds.append("boundary")
        else:
            kinds.append("other")
    boundary = kinds.count("boundary")
    centred = kinds.count("centre")
    if boundary == 3:
        return "bbox_vertex"
    if boundary == 2 and centred == 1:
        return "bbox_edge_midpoint"
    if boundary == 1 and centred == 2:
        return "bbox_face_centre"
    if centred == 3:
        return "bbox_centre"
    return None


def inventory(model: ifcopenshell.file, evidence_tolerance_mm: float) -> list[dict[str, Any]]:
    settings = geometry_settings()
    records = []
    for global_id, expected_type in TARGET_TYPES.items():
        product = model.by_guid(global_id)
        if product is None or not product.is_a("IfcElectricAppliance"):
            raise RuntimeError(f"missing target IfcElectricAppliance {global_id}")
        assigned_type = ifcopenshell.util.element.get_type(product)
        if getattr(assigned_type, "Name", None) != expected_type:
            raise RuntimeError(
                f"unexpected appliance type {global_id}: "
                f"{getattr(assigned_type, 'Name', None)}"
            )
        vertices, faces = world_mesh_mm(settings, product)
        world_points = np.array(vertices, dtype=float)
        matrix = np.array(
            ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement),
            dtype=float,
        )
        homogeneous = np.column_stack(
            (world_points, np.ones(len(world_points), dtype=float))
        )
        local_points = (
            np.linalg.inv(matrix) @ homogeneous.T
        ).T[:, :3]
        bbox_min = local_points.min(axis=0)
        bbox_max = local_points.max(axis=0)
        current = matrix[:3, 3]
        surface_distance = point_mesh_distance(current, vertices, faces)
        records.append(
            {
                "global_id": global_id,
                "type_name": expected_type,
                "type_description": getattr(assigned_type, "Description", None),
                "origin_mm": current.tolist(),
                "rounded_origin_mm": np.round(current).tolist(),
                "round_translation_mm": (np.round(current) - current).tolist(),
                "local_bbox_min_mm": bbox_min.tolist(),
                "local_bbox_max_mm": bbox_max.tolist(),
                "local_bbox_size_mm": (bbox_max - bbox_min).tolist(),
                "current_origin_feature": current_origin_feature(
                    bbox_min, bbox_max, evidence_tolerance_mm
                ),
                "current_origin_surface_distance_mm": surface_distance,
                "placement_rel_to_id": product.ObjectPlacement.PlacementRelTo.id(),
                "representation_types": [
                    representation.RepresentationType
                    for representation in product.Representation.Representations
                ],
            }
        )
    return records


def contact_state(
    targets: list[Any],
    context: list[Any],
    clearance_window_mm: float,
    tolerance_mm: float,
) -> dict[str, Any]:
    return {
        "internal_clearances": clearance_map(
            targets, targets, clearance_window_mm
        ),
        "context_clearances": clearance_map(
            targets, context, clearance_window_mm
        ),
        "internal_intersections": intersection_pairs(
            targets, targets, tolerance_mm
        ),
        "context_intersections": intersection_pairs(
            targets, context, tolerance_mm
        ),
    }


def jsonify_contact_state(state: dict[str, Any]) -> dict[str, Any]:
    return {
        "internal_clearances": [
            {"pair": list(key), "clearance_mm": value}
            for key, value in sorted(state["internal_clearances"].items())
        ],
        "context_clearances": [
            {"pair": list(key), "clearance_mm": value}
            for key, value in sorted(state["context_clearances"].items())
        ],
        "internal_intersections": [
            list(pair) for pair in sorted(state["internal_intersections"])
        ],
        "context_intersections": [
            list(pair) for pair in sorted(state["context_intersections"])
        ],
    }


def compare_option(
    source: ifcopenshell.file,
    candidate: ifcopenshell.file,
    source_state: dict[str, Any],
    target_ids: set[str],
    tolerance_mm: float,
    clearance_window_mm: float,
) -> dict[str, Any]:
    source_context, skipped_source = renderable_context(source, CONTEXT_CLASSES)
    source_context = [
        product for product in source_context if product.GlobalId not in target_ids
    ]
    candidate_context, skipped_candidate = renderable_context(
        candidate, CONTEXT_CLASSES
    )
    candidate_context = [
        product for product in candidate_context if product.GlobalId not in target_ids
    ]
    source_context_ids = {product.GlobalId for product in source_context}
    candidate_context_ids = {product.GlobalId for product in candidate_context}
    candidate_targets = [candidate.by_guid(value) for value in sorted(target_ids)]
    candidate_context_by_source = [
        candidate.by_guid(product.GlobalId) for product in source_context
    ]
    candidate_state = contact_state(
        candidate_targets,
        candidate_context_by_source,
        clearance_window_mm,
        tolerance_mm,
    )
    internal_contacts = classify_contact_changes(
        source_state["internal_clearances"],
        candidate_state["internal_clearances"],
        tolerance_mm,
        set(),
    )
    context_contacts = classify_contact_changes(
        source_state["context_clearances"],
        candidate_state["context_clearances"],
        tolerance_mm,
        set(),
    )
    return {
        "candidate_contact_state": jsonify_contact_state(candidate_state),
        "internal_contact_changes": internal_contacts,
        "context_contact_changes": context_contacts,
        "new_internal_intersections": sorted(
            candidate_state["internal_intersections"]
            - source_state["internal_intersections"]
        ),
        "lost_internal_intersections": sorted(
            source_state["internal_intersections"]
            - candidate_state["internal_intersections"]
        ),
        "new_context_intersections": sorted(
            candidate_state["context_intersections"]
            - source_state["context_intersections"]
        ),
        "lost_context_intersections": sorted(
            source_state["context_intersections"]
            - candidate_state["context_intersections"]
        ),
        "renderable_context_ids_equal": source_context_ids
        == candidate_context_ids,
        "skipped_source_context": sorted(
            record["global_id"] for record in skipped_source
        ),
        "skipped_candidate_context": sorted(
            record["global_id"] for record in skipped_candidate
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    parser.add_argument("--maximum-axis-shift-mm", type=float, default=0.5)
    parser.add_argument("--clearance-window-mm", type=float, default=2.0)
    parser.add_argument("--evidence-tolerance-mm", type=float, default=0.01)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.tolerance_mm <= args.maximum_axis_shift_mm <= 0.5:
        raise SystemExit("expected 0 < tolerance <= maximum axis shift <= 0.5 mm")
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    records = inventory(source, args.evidence_tolerance_mm)
    if len(records) != 4:
        raise RuntimeError(f"expected 4 appliance targets, got {len(records)}")
    target_ids = set(TARGET_TYPES)
    maximum_shift = max(
        max(abs(value) for value in row["round_translation_mm"])
        for row in records
    )
    if maximum_shift > args.maximum_axis_shift_mm:
        raise RuntimeError(f"appliance shift exceeds maximum: {maximum_shift}")

    source_targets = [source.by_guid(value) for value in sorted(target_ids)]
    source_context, _ = renderable_context(source, CONTEXT_CLASSES)
    source_context = [
        product for product in source_context if product.GlobalId not in target_ids
    ]
    source_state = contact_state(
        source_targets,
        source_context,
        args.clearance_window_mm,
        args.tolerance_mm,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    options = {}
    for option in ("rigid_round_legacy_origin", "pure_reset_legacy_origin"):
        candidate = ifcopenshell.open(source_path)
        rows = [
            {
                "global_id": record["global_id"],
                "target_mm": np.array(record["rounded_origin_mm"], dtype=float),
                "translation_mm": record["round_translation_mm"],
                "anchor_kind": "rounded_legacy_object_origin",
                "basis": "comparison option only; not an approved installation anchor",
                "confidence": 0.0,
            }
            for record in records
        ]
        if option == "rigid_round_legacy_origin":
            results = apply_targets(candidate, rows)
        else:
            results = [
                reset_product_origin(
                    candidate, row, args.tolerance_mm, "IfcElectricAppliance"
                )
                for row in rows
            ]
        output = args.output_dir / f"{option}.ifc"
        candidate.write(output)
        candidate = ifcopenshell.open(output)
        geometry = all_product_geometry_difference(
            candidate, source, args.tolerance_mm
        )
        geometry_by_id = {
            row["global_id"]: row for row in geometry["records"]
        }
        non_target_geometry_changes = [
            row["global_id"]
            for row in geometry["records"]
            if row["global_id"] not in target_ids and not row["within_tolerance"]
        ]
        comparison = compare_option(
            source,
            candidate,
            source_state,
            target_ids,
            args.tolerance_mm,
            args.clearance_window_mm,
        )
        options[option] = {
            "path": str(output.resolve()),
            "sha256": sha256(output),
            "results": results,
            "target_geometry": [geometry_by_id[value] for value in sorted(target_ids)],
            "non_target_geometry_changes": non_target_geometry_changes,
            **comparison,
        }
        checkpoint = {
            "source_sha256": sha256(source_path),
            "option": option,
            "candidate_sha256": sha256(output),
            "evidence": options[option],
        }
        (args.output_dir / f"{option}.json").write_text(
            json.dumps(
                checkpoint,
                ensure_ascii=False,
                indent=2,
                default=json_default,
            ),
            encoding="utf-8",
        )

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-electric-appliance-anchor-options",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "maximum_axis_shift_mm": args.maximum_axis_shift_mm,
        "clearance_window_mm": args.clearance_window_mm,
        "evidence_tolerance_mm": args.evidence_tolerance_mm,
        "inventory": records,
        "source_contact_state": jsonify_contact_state(source_state),
        "options": options,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2,
            default=json_default,
        ),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "report": str(args.report),
                "targets": len(records),
                "current_origin_features": {
                    row["type_name"]: row["current_origin_feature"]
                    for row in records
                },
                "options": {
                    name: {
                        "contact_regressions": len(
                            option["context_contact_changes"]["regressions"]
                        )
                        + len(option["internal_contact_changes"]["regressions"]),
                        "new_contacts": len(
                            option["context_contact_changes"]["new_contacts"]
                        )
                        + len(option["internal_contact_changes"]["new_contacts"]),
                        "new_intersections": len(option["new_context_intersections"])
                        + len(option["new_internal_intersections"]),
                        "lost_intersections": len(option["lost_context_intersections"])
                        + len(option["lost_internal_intersections"]),
                        "non_target_geometry_changes": len(
                            option["non_target_geometry_changes"]
                        ),
                    }
                    for name, option in options.items()
                },
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
