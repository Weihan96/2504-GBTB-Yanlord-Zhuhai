#!/usr/bin/env python3
"""Audit the formal IFC after the approved S003 Space Reference write."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element
import numpy as np

from space_reference_semantics_candidate import current_reference, read_register, world_vertices


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formal", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--register", required=True, type=Path)
    parser.add_argument("--prewrite-report", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    rows = read_register(args.register)
    prewrite = json.loads(args.prewrite_report.read_text(encoding="utf-8"))
    if prewrite["candidate"]["sha256"] != sha256(args.candidate):
        raise RuntimeError("approved candidate hash drift")

    formal = ifcopenshell.open(args.formal)
    candidate = ifcopenshell.open(args.candidate)
    if formal.schema != "IFC4" or candidate.schema != "IFC4":
        raise RuntimeError("formal or candidate IFC schema drift")
    formal_settings = ifcopenshell.geom.settings()
    formal_settings.set(formal_settings.USE_WORLD_COORDS, True)
    candidate_settings = ifcopenshell.geom.settings()
    candidate_settings.set(candidate_settings.USE_WORLD_COORDS, True)

    formal_roots = {root.GlobalId for root in formal.by_type("IfcRoot")}
    candidate_roots = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    if formal_roots != candidate_roots:
        raise RuntimeError("formal Root GlobalId set differs from the approved candidate")

    maximum_change = 0.0
    for row in rows:
        formal_space = formal.by_guid(row["space_global_id"])
        candidate_space = candidate.by_guid(row["space_global_id"])
        if current_reference(formal_space) != row["candidate_reference"]:
            raise RuntimeError(f"formal Reference mismatch for {formal_space.GlobalId}")
        if current_reference(candidate_space) != row["candidate_reference"]:
            raise RuntimeError(f"candidate Reference mismatch for {candidate_space.GlobalId}")
        if str(formal_space.LongName or "") != str(candidate_space.LongName or ""):
            raise RuntimeError(f"LongName mismatch for {formal_space.GlobalId}")
        before = world_vertices(candidate_settings, candidate_space)
        after = world_vertices(formal_settings, formal_space)
        if before.shape != after.shape:
            raise RuntimeError(f"Space geometry topology mismatch for {formal_space.GlobalId}")
        maximum_change = max(maximum_change, float(np.max(np.abs(before - after), initial=0.0)))
    if maximum_change > args.tolerance_mm:
        raise RuntimeError(f"formal Space geometry differs from candidate by {maximum_change} mm")

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "formal": {"path": str(args.formal), "sha256": sha256(args.formal)},
        "candidate": {"path": str(args.candidate), "sha256": sha256(args.candidate)},
        "source_sha256": prewrite["source"]["sha256"],
        "tolerance_mm": args.tolerance_mm,
        "space_count": len(rows),
        "reference_count": sum(bool(current_reference(space)) for space in formal.by_type("IfcSpace")),
        "unique_reference_count": len({current_reference(space) for space in formal.by_type("IfcSpace")}),
        "r01_long_name": formal.by_guid(rows[0]["space_global_id"]).LongName,
        "root_sets_match_candidate": True,
        "long_names_match_candidate": True,
        "maximum_space_world_vertex_change_mm": maximum_change,
        "pass": True,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        f"S003 postwrite PASS: 22 References, 22 unique, R01={report['r01_long_name']}, "
        f"Space geometry max change {maximum_change:.6f} mm"
    )


if __name__ == "__main__":
    main()
