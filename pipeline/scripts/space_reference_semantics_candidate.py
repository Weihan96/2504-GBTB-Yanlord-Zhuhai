#!/usr/bin/env python3
"""Create the read-only S003 Space Reference IFC candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api
import ifcopenshell.geom
import ifcopenshell.util.element
import numpy as np


GUID_NAMESPACE = uuid.UUID("84e3fa60-e004-4d96-9d11-92fd35e796db")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--register", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def deterministic_guid(global_id: str, role: str) -> str:
    value = uuid.uuid5(GUID_NAMESPACE, f"S003-SPACE-REFERENCE:{global_id}:{role}")
    return ifcopenshell.guid.compress(value.hex)


def read_register(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 22:
        raise RuntimeError(f"expected 22 Space Reference rows, got {len(rows)}")
    expected = [f"R{index:02d}" for index in range(1, 23)]
    if [row["candidate_reference"] for row in rows] != expected:
        raise RuntimeError("register is not the complete ordered R01-R22 sequence")
    if rows[0]["space_long_name"] != "玄关":
        raise RuntimeError("R01 is not 玄关")
    if any(row["formal_ifc_write_allowed"] != "no" for row in rows):
        raise RuntimeError("review register is not a read-only candidate")
    return rows


def current_reference(space: ifcopenshell.entity_instance) -> str:
    return str(ifcopenshell.util.element.get_psets(space).get("Pset_SpaceCommon", {}).get("Reference") or "")


def world_vertices(settings: ifcopenshell.geom.settings, product: ifcopenshell.entity_instance) -> np.ndarray:
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3)) * 1000.0
    return vertices[np.lexsort((vertices[:, 2], vertices[:, 1], vertices[:, 0]))]


def apply_references(model: ifcopenshell.file, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    applied = []
    for row in rows:
        space = model.by_guid(row["space_global_id"])
        if space is None or not space.is_a("IfcSpace"):
            raise RuntimeError(f"missing IfcSpace {row['space_global_id']}")
        if str(space.LongName or space.Name or "") != row["space_long_name"]:
            raise RuntimeError(f"Space LongName drift for {space.GlobalId}")
        if current_reference(space):
            raise RuntimeError(f"Space {space.GlobalId} already has Reference {current_reference(space)!r}")
        pset = ifcopenshell.api.run("pset.add_pset", model, product=space, name="Pset_SpaceCommon")
        pset.GlobalId = deterministic_guid(space.GlobalId, "PSET")
        relation = next(
            relationship
            for relationship in space.IsDefinedBy
            if relationship.is_a("IfcRelDefinesByProperties")
            and relationship.RelatingPropertyDefinition == pset
        )
        relation.GlobalId = deterministic_guid(space.GlobalId, "REL")
        ifcopenshell.api.run(
            "pset.edit_pset",
            model,
            pset=pset,
            properties={"Reference": row["candidate_reference"]},
        )
        applied.append(
            {
                "space_global_id": space.GlobalId,
                "space_long_name": str(space.LongName or ""),
                "reference": row["candidate_reference"],
                "pset_global_id": pset.GlobalId,
                "relation_global_id": relation.GlobalId,
            }
        )
    return applied


def main() -> None:
    args = parse_args()
    if args.tolerance_mm <= 0:
        raise RuntimeError("tolerance must be positive")
    rows = read_register(args.register)
    source_sha = sha256(args.input)
    model = ifcopenshell.open(args.input)
    if model.schema != "IFC4":
        raise RuntimeError(f"expected IFC4, got {model.schema}")

    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    source_roots = {root.GlobalId for root in model.by_type("IfcRoot")}
    source_geometry = {space.GlobalId: world_vertices(settings, space) for space in model.by_type("IfcSpace")}
    source_long_names = {space.GlobalId: str(space.LongName or "") for space in model.by_type("IfcSpace")}
    if any(current_reference(space) for space in model.by_type("IfcSpace")):
        raise RuntimeError("formal IFC already contains one or more Space References")

    applied = apply_references(model, rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    model.write(args.output)

    candidate = ifcopenshell.open(args.output)
    candidate_settings = ifcopenshell.geom.settings()
    candidate_settings.set(candidate_settings.USE_WORLD_COORDS, True)
    candidate_roots = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    expected_new_roots = {
        value
        for record in applied
        for value in (record["pset_global_id"], record["relation_global_id"])
    }
    if candidate_roots != source_roots | expected_new_roots:
        raise RuntimeError("candidate Root GlobalId boundary differs from the deterministic S003 target set")

    maximum_change = 0.0
    for row in rows:
        space = candidate.by_guid(row["space_global_id"])
        if current_reference(space) != row["candidate_reference"]:
            raise RuntimeError(f"candidate Reference mismatch for {space.GlobalId}")
        if str(space.LongName or "") != source_long_names[space.GlobalId]:
            raise RuntimeError(f"candidate LongName changed for {space.GlobalId}")
        before = source_geometry[space.GlobalId]
        after = world_vertices(candidate_settings, space)
        if before.shape != after.shape:
            raise RuntimeError(f"candidate Space geometry topology changed for {space.GlobalId}")
        maximum_change = max(maximum_change, float(np.max(np.abs(before - after), initial=0.0)))
    if maximum_change > args.tolerance_mm:
        raise RuntimeError(f"candidate Space geometry changed by {maximum_change} mm")

    report: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {"path": str(args.input), "sha256": source_sha, "schema": model.schema},
        "candidate": {"path": str(args.output), "sha256": sha256(args.output)},
        "automatic_formal_ifc_write_allowed": False,
        "tolerance_mm": args.tolerance_mm,
        "space_count": len(rows),
        "reference_count": sum(bool(current_reference(space)) for space in candidate.by_type("IfcSpace")),
        "unique_reference_count": len({row["candidate_reference"] for row in rows}),
        "r01_long_name": rows[0]["space_long_name"],
        "new_root_count": len(expected_new_roots),
        "maximum_space_world_vertex_change_mm": maximum_change,
        "long_names_preserved": True,
        "applied": applied,
        "qa": {
            "references_complete_unique": True,
            "foyer_is_r01": True,
            "root_boundary_exact": True,
            "space_geometry_within_tolerance": maximum_change <= args.tolerance_mm,
            "formal_ifc_unchanged": sha256(args.input) == source_sha,
        },
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        f"S003 candidate: 22 References, 44 deterministic new Roots, "
        f"Space geometry max change {maximum_change:.6f} mm, formal IFC unchanged"
    )


if __name__ == "__main__":
    main()
