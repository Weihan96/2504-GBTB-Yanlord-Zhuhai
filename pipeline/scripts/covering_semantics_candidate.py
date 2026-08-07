#!/usr/bin/env python3
"""Build and verify a semantic-only IfcCovering candidate."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import numpy as np

from geometry_alignment_audit import geometry_difference_audit, sha256


def baseboard_evidence(
    name: str | None,
    predefined_type: str | None,
    dimensions_mm: list[float],
) -> tuple[bool, str]:
    """Require explicit naming plus a narrow 90-110 mm high geometry."""

    if not (name or "").startswith("Baseboard"):
        return False, "Name does not explicitly identify a baseboard"
    if predefined_type and predefined_type not in {"NOTDEFINED", "USERDEFINED"}:
        return False, f"PredefinedType is already {predefined_type}"
    x_mm, y_mm, z_mm = dimensions_mm
    if not 90.0 <= z_mm <= 110.0:
        return False, f"height {z_mm:.6f} mm is outside 90-110 mm"
    if min(x_mm, y_mm) > 25.0:
        return False, f"minimum plan thickness {min(x_mm, y_mm):.6f} mm exceeds 25 mm"
    if max(x_mm, y_mm) < 100.0:
        return False, f"plan length {max(x_mm, y_mm):.6f} mm is below 100 mm"
    return (
        True,
        "Name begins Baseboard; geometry is 90-110 mm high, at most 25 mm thick, and at least 100 mm long",
    )


def world_dimensions_mm(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> list[float]:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.array(shape.geometry.verts, dtype=float).reshape(-1, 3) * 1000.0
    return (vertices.max(axis=0) - vertices.min(axis=0)).tolist()


def apply_baseboard_semantics(model: ifcopenshell.file) -> list[dict[str, Any]]:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    results = []
    for covering in model.by_type("IfcCovering"):
        if not covering.ObjectPlacement or not covering.Representation:
            continue
        dimensions = world_dimensions_mm(settings, covering)
        accepted, basis = baseboard_evidence(
            covering.Name,
            covering.PredefinedType,
            dimensions,
        )
        if not accepted:
            continue
        before = covering.PredefinedType
        covering.PredefinedType = "SKIRTINGBOARD"
        results.append(
            {
                "global_id": covering.GlobalId,
                "name": covering.Name,
                "before": before,
                "after": covering.PredefinedType,
                "dimensions_mm": dimensions,
                "basis": basis,
                "confidence": 0.99,
                "review_required": False,
            }
        )
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--expected-count", type=int, default=19)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.expected_count <= 0:
        raise SystemExit("--expected-count must be positive")
    if args.tolerance_mm <= 0.0:
        raise SystemExit("--tolerance-mm must be positive")
    source_path = args.input.resolve()
    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    results = apply_baseboard_semantics(candidate)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(args.output)
    candidate = ifcopenshell.open(args.output)

    geometry_difference = geometry_difference_audit(
        candidate,
        source,
        str(source_path),
        tolerance_mm=args.tolerance_mm,
        classes=["IfcCovering"],
        global_ids=[],
    )
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    changed_ids = {result["global_id"] for result in results}
    written_ids = {
        covering.GlobalId
        for covering in candidate.by_type("IfcCovering")
        if covering.PredefinedType == "SKIRTINGBOARD"
    }
    gates = {
        "expected_count": args.expected_count,
        "changed_count": len(results),
        "written_count": len(written_ids),
        "changed_ids_equal_written_ids": changed_ids == written_ids,
        "geometry_over_tolerance": geometry_difference["over_tolerance"],
        "schema_equal": source.schema == candidate.schema,
        "root_global_ids_equal": source_root_ids == candidate_root_ids,
        "entity_count_equal": len(list(source)) == len(list(candidate)),
    }
    gates["pass"] = (
        len(results) == args.expected_count
        and len(written_ids) == args.expected_count
        and gates["changed_ids_equal_written_ids"]
        and gates["geometry_over_tolerance"] == 0
        and gates["schema_equal"]
        and gates["root_global_ids_equal"]
        and gates["entity_count_equal"]
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-covering-semantics-candidate",
        "source": {
            "path": str(source_path),
            "sha256": sha256(source_path),
            "schema": source.schema,
        },
        "candidate": {
            "path": str(args.output),
            "sha256": sha256(args.output),
            "schema": candidate.schema,
        },
        "tolerance_mm": args.tolerance_mm,
        "results": results,
        "geometry_difference": geometry_difference,
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
