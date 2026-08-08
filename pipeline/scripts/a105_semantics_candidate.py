#!/usr/bin/env python3
"""Build confirmed A-105 finish, slope, and depressed-slab semantics as a candidate."""

from __future__ import annotations

import argparse
import csv
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element
import ifcopenshell.util.placement
import numpy as np

from geometry_alignment_audit import geometry_difference_audit, sha256


GUID_NAMESPACE = uuid.UUID("5b978cc9-54bd-4e93-999a-da139a13b79f")
CONFIRMED_DECISIONS = {
    "A105-FLOOR-MATERIAL-001",
    "A105-FLOOR-BUILDUP-001",
    "A105-WET-SLOPE-GRADE-001",
    "A105-WET-SLOPE-DIRECTION-001",
    "A105-DEPRESSED-SLAB-001",
}
REFERENCE_MATERIALS = {
    "1ogoq1VJP4vgBTqBCMFyCn": "地板",
    "11c$NwzxL9hASgCgwDxM$w": "银白洞石岩板",
    "3BwA5Rnkf4Bf3v3vTmeWTw": "地板",
}
DEPRESSED_SLABS = {
    "3ARl_CqPrCWQrA6_$W073W": "1Ro3zU6lXBjuo7VDM3mTa7",
    "1UixnpnQb97wuNAGrTmMJd": "2F667JcWnAm9Hchuxrf8vg",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--floor-register", required=True, type=Path)
    parser.add_argument("--floor-report", required=True, type=Path)
    parser.add_argument("--build-up-report", required=True, type=Path)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--tolerance-mm", type=float, default=0.1)
    return parser.parse_args()


def deterministic_guid(key: str) -> str:
    return ifcopenshell.guid.compress(uuid.uuid5(GUID_NAMESPACE, key).hex)


def require_reports(floor_report: Path, build_up_report: Path, source_sha: str) -> None:
    floor = json.loads(floor_report.read_text(encoding="utf-8"))
    build_up = json.loads(build_up_report.read_text(encoding="utf-8"))
    if floor.get("source", {}).get("ifc_sha256") != source_sha or not floor.get("gates", {}).get("mechanical_pass"):
        raise RuntimeError("A-105 floor report does not match the formal IFC or has not passed")
    if build_up.get("source", {}).get("ifc_sha256") != source_sha or not build_up.get("gates", {}).get("mechanical_pass"):
        raise RuntimeError("A-105 build-up report does not match the formal IFC or has not passed")


def require_decisions(path: Path) -> None:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        statuses = {row["decision_id"]: row["status"] for row in csv.DictReader(handle)}
    missing = sorted(decision for decision in CONFIRMED_DECISIONS if statuses.get(decision) != "confirmed")
    if missing:
        raise RuntimeError(f"A-105 decisions are not confirmed: {missing}")


def read_floors(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 21 or len({row["candidate_id"] for row in rows}) != 21:
        raise RuntimeError("A-105 floor register identity drift")
    for row in rows:
        row["top_plane"] = json.loads(row["top_plane"]) if row["top_plane"] else None
    return rows


def property_value(model: ifcopenshell.file, value: Any) -> Any:
    if isinstance(value, bool):
        return model.create_entity("IfcBoolean", value)
    if isinstance(value, float):
        return model.create_entity("IfcReal", value)
    return model.create_entity("IfcLabel", str(value))


def add_pset(model: ifcopenshell.file, product: Any, name: str, properties: dict[str, Any]) -> None:
    if name in ifcopenshell.util.element.get_psets(product):
        raise RuntimeError(f"candidate pset already exists on {product.GlobalId}: {name}")
    pset = model.create_entity(
        "IfcPropertySet",
        GlobalId=deterministic_guid(f"{product.GlobalId}:{name}:PSET"),
        Name=name,
        HasProperties=[
            model.create_entity("IfcPropertySingleValue", Name=key, NominalValue=property_value(model, value))
            for key, value in properties.items()
        ],
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=deterministic_guid(f"{product.GlobalId}:{name}:REL"),
        RelatedObjects=[product],
        RelatingPropertyDefinition=pset,
    )


def assign_reference_materials(model: ifcopenshell.file) -> None:
    by_name: dict[str, list[Any]] = {}
    for global_id, material_name in REFERENCE_MATERIALS.items():
        by_name.setdefault(material_name, []).append(model.by_guid(global_id))
    for material_name, products in sorted(by_name.items()):
        material = next((item for item in model.by_type("IfcMaterial") if item.Name == material_name), None)
        if material is None:
            material = model.create_entity("IfcMaterial", Name=material_name)
        model.create_entity(
            "IfcRelAssociatesMaterial",
            GlobalId=deterministic_guid(f"A105:MATERIAL:{material_name}:REL"),
            RelatedObjects=products,
            RelatingMaterial=material,
        )


def apply_semantics(model: ifcopenshell.file, rows: list[dict[str, Any]]) -> None:
    for row in rows:
        product = model.by_guid(row["global_id"])
        if product is None or not product.is_a("IfcCovering"):
            raise RuntimeError(f"missing A-105 floor product {row['global_id']}")
        product.Tag = row["candidate_id"]
        if row["kind"] == "sloped_wet_tile":
            plane = row["top_plane"]
            add_pset(
                model,
                product,
                "Pset_A105SlopeReview",
                {
                    "SlopeRatio": float(plane["slope_percent"]) / 100.0,
                    "DownhillDirection": plane["downhill_direction"],
                    "DirectionConfirmed": True,
                    "ReviewStatus": "CONFIRMED",
                },
            )
        elif product.GlobalId in REFERENCE_MATERIALS:
            add_pset(
                model,
                product,
                "Pset_A105FinishIntent",
                {
                    "GeometryRole": "FINISH_REFERENCE_PLANE",
                    "ConfirmedMaterial": REFERENCE_MATERIALS[product.GlobalId],
                    "BuildUpDepthMm": 50.0,
                    "RealLayerGeometryPending": True,
                    "ReviewStatus": "CONFIRMED",
                },
            )
        else:
            raise RuntimeError(f"unclassified A-105 floor product {product.GlobalId}")
    assign_reference_materials(model)
    for slab_id, opening_id in DEPRESSED_SLABS.items():
        slab = model.by_guid(slab_id)
        opening = model.by_guid(opening_id)
        if slab is None or opening is None or not slab.is_a("IfcSlab") or not opening.is_a("IfcOpeningElement"):
            raise RuntimeError(f"missing depressed slab evidence {slab_id}/{opening_id}")
        add_pset(
            model,
            slab,
            "Pset_A105DepressedSlabReview",
            {
                "IsDepressedSlab": True,
                "RecessDepthMm": 50.0,
                "OpeningGlobalId": opening_id,
                "ReviewStatus": "CONFIRMED",
            },
        )


def placement_matrix(product: Any) -> np.ndarray:
    return np.asarray(ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement), dtype=float)


def main() -> None:
    args = parse_args()
    if args.tolerance_mm <= 0:
        raise SystemExit("--tolerance-mm must be positive")
    source_path = args.input.resolve()
    output_path = args.output.resolve()
    if source_path == output_path:
        raise RuntimeError("candidate output must not overwrite the formal IFC")
    source_sha = sha256(source_path)
    require_reports(args.floor_report, args.build_up_report, source_sha)
    require_decisions(args.decisions)
    rows = read_floors(args.floor_register)
    source = ifcopenshell.open(source_path)
    protected_ids = [row["global_id"] for row in rows] + list(DEPRESSED_SLABS) + list(DEPRESSED_SLABS.values())
    source_placements = {global_id: placement_matrix(source.by_guid(global_id)) for global_id in protected_ids}
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    source_void_ids = {relation.GlobalId for relation in source.by_type("IfcRelVoidsElement")}

    candidate = ifcopenshell.open(source_path)
    apply_semantics(candidate, rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(output_path)
    candidate = ifcopenshell.open(output_path)

    geometry = geometry_difference_audit(
        candidate,
        source,
        str(source_path),
        tolerance_mm=args.tolerance_mm,
        classes=("IfcCovering", "IfcSlab", "IfcOpeningElement"),
        global_ids=(),
    )
    placement_deltas = {
        global_id: float(np.max(np.abs(placement_matrix(candidate.by_guid(global_id)) - matrix)))
        for global_id, matrix in source_placements.items()
    }
    tag_mismatches = [row["global_id"] for row in rows if candidate.by_guid(row["global_id"]).Tag != row["candidate_id"]]
    material_mismatches = []
    for global_id, expected in REFERENCE_MATERIALS.items():
        material = ifcopenshell.util.element.get_material(candidate.by_guid(global_id))
        if material is None or not material.is_a("IfcMaterial") or material.Name != expected:
            material_mismatches.append(global_id)
    wet_pset_count = sum(
        "Pset_A105SlopeReview" in ifcopenshell.util.element.get_psets(candidate.by_guid(row["global_id"]))
        for row in rows if row["kind"] == "sloped_wet_tile"
    )
    reference_pset_count = sum(
        "Pset_A105FinishIntent" in ifcopenshell.util.element.get_psets(candidate.by_guid(global_id))
        for global_id in REFERENCE_MATERIALS
    )
    slab_pset_count = sum(
        "Pset_A105DepressedSlabReview" in ifcopenshell.util.element.get_psets(candidate.by_guid(global_id))
        for global_id in DEPRESSED_SLABS
    )
    candidate_root_ids = {root.GlobalId for root in candidate.by_type("IfcRoot")}
    gates = {
        "schema_equal": source.schema == candidate.schema == "IFC4",
        "tag_count": 21 - len(tag_mismatches),
        "tag_mismatch_count": len(tag_mismatches),
        "wet_slope_pset_count": wet_pset_count,
        "finish_intent_pset_count": reference_pset_count,
        "depressed_slab_pset_count": slab_pset_count,
        "reference_material_count": len(REFERENCE_MATERIALS) - len(material_mismatches),
        "geometry_over_tolerance": geometry["over_tolerance"],
        "maximum_placement_matrix_delta": max(placement_deltas.values()),
        "voids_relationship_ids_equal": source_void_ids == {relation.GlobalId for relation in candidate.by_type("IfcRelVoidsElement")},
        "source_root_ids_preserved": source_root_ids <= candidate_root_ids,
        "new_root_count": len(candidate_root_ids - source_root_ids),
        "formal_ifc_write_allowed": False,
    }
    gates["mechanical_pass"] = (
        gates["schema_equal"]
        and gates["tag_count"] == 21
        and gates["tag_mismatch_count"] == 0
        and gates["wet_slope_pset_count"] == 18
        and gates["finish_intent_pset_count"] == 3
        and gates["depressed_slab_pset_count"] == 2
        and gates["reference_material_count"] == 3
        and gates["geometry_over_tolerance"] == 0
        and gates["maximum_placement_matrix_delta"] <= 1e-12
        and gates["voids_relationship_ids_equal"]
        and gates["source_root_ids_preserved"]
        and gates["new_root_count"] == 48
    )
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read-only-a105-semantics-candidate",
        "source": {"path": str(source_path), "sha256": source_sha, "schema": source.schema},
        "candidate": {"path": str(output_path), "sha256": sha256(output_path), "schema": candidate.schema},
        "tolerance_mm": args.tolerance_mm,
        "reference_materials": REFERENCE_MATERIALS,
        "depressed_slabs": DEPRESSED_SLABS,
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
