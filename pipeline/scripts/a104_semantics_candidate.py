#!/usr/bin/env python3
"""Build the confirmed A-104 tags and door/window groups as an IFC candidate."""

from __future__ import annotations

import argparse
import csv
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.placement
import numpy as np

from geometry_alignment_audit import geometry_difference_audit, sha256


GUID_NAMESPACE = uuid.UUID("82c77ed3-51c1-47db-ad86-7027eae91f0e")
CONFIRMED_DECISIONS = {
    "A104-DOOR-IDENTITY-001",
    "A104-RIMADESIO-SAIL-001",
    "A104-WINDOW-GROUP-DINING-001",
    "A104-WINDOW-GROUP-BED2-001",
}
NAME_UPDATES = {
    "3xKBbA2CT9mfby2$MODnzM": "次卧门",
    "1TW6$_GfnABRZusYvx0zZG": "主卧门",
    "2D5BPoo2XFSvhTdfPenCh7": "Rimadesio Sail 格栅门板",
    "0zjVS5FBbBewgUkk0fdfiv": "Rimadesio Sail 轨道",
}
GROUPS = {
    "A104-SAIL-DOOR-GROUP": {
        "name": "M05/M06 Rimadesio Sail 门组",
        "object_type": "DOOR_GROUP",
        "description": "M05 为格栅门板；M06 为轨道。保留两个现有 IfcDoor，不猜宿主。",
        "members": ("2D5BPoo2XFSvhTdfPenCh7", "0zjVS5FBbBewgUkk0fdfiv"),
    },
    "A104-DINING-BAY-WINDOW-GROUP": {
        "name": "W09/W10 餐厅飘窗窗组",
        "object_type": "WINDOW_GROUP",
        "description": "一组对称双通风窗；保留两樘 IfcWindow 及现有洞口几何。",
        "members": ("1KBBoRNWb1Svpz9__rmYKy", "1uQsUUwhT9Y8fHZJQtvJqa"),
    },
    "A104-BED2-BAY-WINDOW-GROUP": {
        "name": "W02/W03 次卧飘窗窗组",
        "object_type": "WINDOW_GROUP",
        "description": "一个窗组；保留两樘 IfcWindow 及现有洞口几何。",
        "members": ("2nisU4bAr9MfEP9b7LwE2g", "0sj5O$jhz2Huyt2vHf63WK"),
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--register", required=True, type=Path)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def deterministic_guid(key: str) -> str:
    return ifcopenshell.guid.compress(uuid.uuid5(GUID_NAMESPACE, key).hex)


def read_register(path: Path, source_sha: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 19 or any(row["source_ifc_sha256"] != source_sha for row in rows):
        raise RuntimeError("A-104 register does not match the formal IFC")
    ids = [row["candidate_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise RuntimeError("A-104 candidate identifiers are not unique")
    return rows


def require_decisions(path: Path) -> None:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        statuses = {row["decision_id"]: row["status"] for row in csv.DictReader(handle)}
    missing = sorted(decision for decision in CONFIRMED_DECISIONS if statuses.get(decision) != "confirmed")
    if missing:
        raise RuntimeError(f"A-104 decisions are not confirmed: {missing}")


def placement_matrix(product: Any) -> np.ndarray:
    return np.asarray(ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement), dtype=float)


def relation_ids(model: ifcopenshell.file, ifc_class: str) -> set[str]:
    return {relation.GlobalId for relation in model.by_type(ifc_class)}


def apply_semantics(model: ifcopenshell.file, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    records = []
    for row in rows:
        product = model.by_guid(row["global_id"])
        if product is None or product.is_a() != row["ifc_class"]:
            raise RuntimeError(f"A-104 product identity drift: {row['global_id']}")
        product.Tag = row["candidate_tag"]
        if product.GlobalId in NAME_UPDATES:
            product.Name = NAME_UPDATES[product.GlobalId]
        records.append(
            {
                "candidate_id": row["candidate_id"],
                "global_id": product.GlobalId,
                "ifc_class": product.is_a(),
                "name": str(product.Name or ""),
                "tag": str(product.Tag or ""),
            }
        )
    for key, spec in GROUPS.items():
        group = model.create_entity(
            "IfcGroup",
            GlobalId=deterministic_guid(f"{key}:GROUP"),
            Name=spec["name"],
            Description=spec["description"],
            ObjectType=spec["object_type"],
        )
        members = [model.by_guid(global_id) for global_id in spec["members"]]
        if any(member is None for member in members):
            raise RuntimeError(f"missing A-104 group member for {key}")
        model.create_entity(
            "IfcRelAssignsToGroup",
            GlobalId=deterministic_guid(f"{key}:REL"),
            RelatedObjects=members,
            RelatingGroup=group,
        )
    return records


def group_record(model: ifcopenshell.file, key: str, spec: dict[str, Any]) -> dict[str, Any]:
    group = model.by_guid(deterministic_guid(f"{key}:GROUP"))
    relations = list(group.IsGroupedBy or ()) if group else []
    members = sorted(member.GlobalId for relation in relations for member in relation.RelatedObjects)
    return {
        "global_id": group.GlobalId if group else None,
        "name": group.Name if group else None,
        "object_type": group.ObjectType if group else None,
        "members": members,
        "expected_members": sorted(spec["members"]),
        "passes": bool(group) and len(relations) == 1 and members == sorted(spec["members"]),
    }


def main() -> None:
    args = parse_args()
    if args.tolerance_mm <= 0:
        raise SystemExit("--tolerance-mm must be positive")
    source_path = args.input.resolve()
    output_path = args.output.resolve()
    if source_path == output_path:
        raise RuntimeError("candidate output must not overwrite the formal IFC")
    source_sha = sha256(source_path)
    rows = read_register(args.register, source_sha)
    require_decisions(args.decisions)
    source = ifcopenshell.open(source_path)
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    source_placements = {row["global_id"]: placement_matrix(source.by_guid(row["global_id"])) for row in rows}
    source_fills = relation_ids(source, "IfcRelFillsElement")
    source_voids = relation_ids(source, "IfcRelVoidsElement")

    candidate = ifcopenshell.open(source_path)
    written = apply_semantics(candidate, rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(output_path)
    candidate = ifcopenshell.open(output_path)

    geometry = geometry_difference_audit(
        candidate,
        source,
        str(source_path),
        tolerance_mm=args.tolerance_mm,
        classes=("IfcDoor", "IfcWindow", "IfcOpeningElement"),
        global_ids=(),
    )
    groups = [group_record(candidate, key, spec) for key, spec in GROUPS.items()]
    tag_mismatches = [
        row["global_id"]
        for row in rows
        if candidate.by_guid(row["global_id"]).Tag != row["candidate_tag"]
    ]
    name_mismatches = [
        global_id
        for global_id, name in NAME_UPDATES.items()
        if candidate.by_guid(global_id).Name != name
    ]
    placement_deltas = {
        global_id: float(np.max(np.abs(placement_matrix(candidate.by_guid(global_id)) - matrix)))
        for global_id, matrix in source_placements.items()
    }
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "schema_equal": source.schema == candidate.schema == "IFC4",
        "tag_count": len(rows) - len(tag_mismatches),
        "tag_mismatch_count": len(tag_mismatches),
        "confirmed_name_count": len(NAME_UPDATES) - len(name_mismatches),
        "group_count": len(groups),
        "group_pass_count": sum(group["passes"] for group in groups),
        "geometry_over_tolerance": geometry["over_tolerance"],
        "maximum_placement_matrix_delta": max(placement_deltas.values()),
        "fills_relationship_ids_equal": source_fills == relation_ids(candidate, "IfcRelFillsElement"),
        "voids_relationship_ids_equal": source_voids == relation_ids(candidate, "IfcRelVoidsElement"),
        "source_root_ids_preserved": source_root_ids <= candidate_root_ids,
        "new_root_count": len(candidate_root_ids - source_root_ids),
        "formal_ifc_write_allowed": False,
    }
    gates["mechanical_pass"] = (
        gates["schema_equal"]
        and gates["tag_count"] == 19
        and gates["tag_mismatch_count"] == 0
        and gates["confirmed_name_count"] == 4
        and gates["group_count"] == 3
        and gates["group_pass_count"] == 3
        and gates["geometry_over_tolerance"] == 0
        and gates["maximum_placement_matrix_delta"] <= 1e-12
        and gates["fills_relationship_ids_equal"]
        and gates["voids_relationship_ids_equal"]
        and gates["source_root_ids_preserved"]
        and gates["new_root_count"] == 6
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a104-semantics-candidate",
        "source": {"path": str(source_path), "sha256": source_sha, "schema": source.schema},
        "candidate": {"path": str(output_path), "sha256": sha256(output_path), "schema": candidate.schema},
        "tolerance_mm": args.tolerance_mm,
        "written_semantics": written,
        "groups": groups,
        "placement_deltas": placement_deltas,
        "geometry_difference": geometry,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["mechanical_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
