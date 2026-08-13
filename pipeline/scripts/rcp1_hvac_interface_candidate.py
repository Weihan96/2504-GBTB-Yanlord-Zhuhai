"""Validate HVAC manufacturer-interface evidence against the formal IFC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
import ifcopenshell.util.element


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def split_values(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def optional_float(value: str) -> float | None:
    return float(value) if value.strip() else None


def local_manual_path(manual_dir: Path | None, document: str) -> Path | None:
    if manual_dir is None or not document.lower().endswith("q"):
        return None
    matches = sorted(manual_dir.glob(f"**/{document}*.pdf"))
    return matches[0] if matches else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--m401-review", type=Path, required=True)
    parser.add_argument("--evidence", type=Path)
    parser.add_argument("--source-evidence", type=Path, default=Path("pipeline/decisions/source-evidence-register.csv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manual-dir", type=Path)
    args = parser.parse_args()

    formal_sha = sha256(args.input)
    m401_rows = read_csv(args.m401_review)
    m401_hashes = {row["source_ifc_sha256"] for row in m401_rows}
    if m401_hashes != {formal_sha}:
        raise RuntimeError("M-401 review does not match the formal IFC")

    model = ifcopenshell.open(args.input)
    if args.evidence:
        evidence_rows = read_csv(args.evidence)
        evidence_register = args.evidence
    else:
        canonical_sources = read_csv(args.source_evidence)
        evidence_rows = [
            json.loads(row["legacy_projection_json"])
            for row in canonical_sources
            if "rcp1-hvac-equipment-interface-evidence.csv" in split_values(row["legacy_targets"])
        ]
        evidence_register = args.source_evidence
    occurrences = []
    manual_checks = []
    type_assignments_match = True
    local_manual_hashes_match = True

    for evidence in evidence_rows:
        equipment_ids = split_values(evidence["equipment_ids"])
        global_ids = split_values(evidence["equipment_global_ids"])
        if global_ids and len(equipment_ids) != len(global_ids):
            raise RuntimeError(f"equipment/global-id count mismatch: {evidence['evidence_id']}")
        expected_type_id = evidence["ifc_type_global_id"].strip() or None
        for index, equipment_id in enumerate(equipment_ids):
            global_id = global_ids[index] if index < len(global_ids) else None
            if global_id is None:
                occurrences.append({
                    "equipment_id": equipment_id,
                    "global_id": None,
                    "formal_identity_status": "missing",
                    "type_assignment_matches": None,
                    "evidence_id": evidence["evidence_id"],
                    "interface_status": evidence["status"],
                })
                continue
            product = model.by_guid(global_id)
            if product is None:
                raise RuntimeError(f"missing formal equipment: {equipment_id} / {global_id}")
            assigned_type = ifcopenshell.util.element.get_type(product)
            actual_type_id = assigned_type.GlobalId if assigned_type is not None else None
            assignment_matches = actual_type_id == expected_type_id
            type_assignments_match = type_assignments_match and assignment_matches
            occurrences.append({
                "equipment_id": equipment_id,
                "global_id": global_id,
                "ifc_class": product.is_a(),
                "formal_identity_status": "present",
                "expected_type_global_id": expected_type_id,
                "actual_type_global_id": actual_type_id,
                "type_assignment_matches": assignment_matches,
                "ifc_type_name": evidence["ifc_type_name"] or None,
                "ifc_model_text": evidence["ifc_model_text"] or None,
                "nominal_body_mm": [
                    optional_float(evidence[key])
                    for key in ("nominal_body_width_mm", "nominal_body_depth_mm", "nominal_body_height_mm")
                ],
                "capacity_group": evidence["capacity_group"] or None,
                "manufacturer_interface": {
                    "gas_pipe_od_mm": optional_float(evidence["gas_pipe_od_mm"]),
                    "liquid_pipe_od_mm": optional_float(evidence["liquid_pipe_od_mm"]),
                    "drain_pipe_od_mm": optional_float(evidence["drain_pipe_od_mm"]),
                    "drain_slope_min": evidence["drain_slope_min"] or None,
                    "drain_slope_max": evidence["drain_slope_max"] or None,
                    "connection_side_basis": evidence["connection_side_basis"],
                    "port_coordinate_status": evidence["port_coordinate_status"],
                },
                "evidence_id": evidence["evidence_id"],
                "interface_status": evidence["status"],
                "confidence": float(evidence["confidence"]),
                "review_required": evidence["review_required"].lower() == "yes",
                "formal_ifc_write_allowed": evidence["formal_ifc_write_allowed"].lower() == "yes",
            })

        manual_path = local_manual_path(args.manual_dir, evidence["source_document"])
        expected_manual_sha = evidence["source_sha256"]
        actual_manual_sha = sha256(manual_path) if manual_path is not None else None
        hash_matches = actual_manual_sha == expected_manual_sha if actual_manual_sha is not None else None
        if hash_matches is False:
            local_manual_hashes_match = False
        manual_checks.append({
            "evidence_id": evidence["evidence_id"],
            "source_document": evidence["source_document"],
            "source_url": evidence["source_url"] or None,
            "expected_sha256": expected_manual_sha,
            "local_path": str(manual_path) if manual_path is not None else None,
            "actual_sha256": actual_manual_sha,
            "hash_matches": hash_matches,
            "pdf_page": evidence["pdf_page"] or None,
        })

    exact_interface = [row for row in occurrences if row["interface_status"].startswith("official_model_family_match")]
    geometry_only = [
        row for row in occurrences
        if row["interface_status"] == "accepted_exact_family_wildcard_geometry_only"
    ]
    missing_identity = [row for row in occurrences if row["formal_identity_status"] == "missing"]
    unresolved_model = [row for row in occurrences if row["interface_status"] == "model_text_mismatch_requires_confirmation"]
    precise_ports = [
        row for row in occurrences
        if row.get("manufacturer_interface", {}).get("port_coordinate_status") == "approved_precise_coordinates"
    ]
    formal_writes_allowed = [row for row in occurrences if row.get("formal_ifc_write_allowed")]

    output = {
        "mode": "read_only_rcp1_hvac_manufacturer_interface_candidate",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source": {
            "ifc": str(args.input),
            "ifc_sha256": formal_sha,
            "m401_review": str(args.m401_review),
            "evidence_register": str(evidence_register),
        },
        "summary": {
            "equipment_count": len(occurrences),
            "formal_identity_present_count": len(occurrences) - len(missing_identity),
            "formal_identity_missing_count": len(missing_identity),
            "official_interface_match_count": len(exact_interface),
            "official_geometry_only_match_count": len(geometry_only),
            "model_text_mismatch_count": len(unresolved_model),
            "precise_port_coordinate_ready_count": len(precise_ports),
            "formal_ifc_write_allowed_count": len(formal_writes_allowed),
        },
        "gates": {
            "source_hash_matches_m401_review": True,
            "formal_ifc_type_assignments_match": type_assignments_match,
            "available_local_manual_hashes_match": local_manual_hashes_match,
            "all_equipment_has_formal_identity": len(missing_identity) == 0,
            "all_present_model_texts_have_official_family_match": (
                len(exact_interface) + len(geometry_only) == len(occurrences) - len(missing_identity)
            ),
            "all_model_interfaces_have_official_match": len(exact_interface) == len(occurrences),
            "precise_port_coordinates_ready": len(precise_ports) == len(occurrences),
            "formal_ifc_write_allowed": len(formal_writes_allowed) == len(occurrences),
        },
        "occurrences": occurrences,
        "manual_checks": manual_checks,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output["summary"], ensure_ascii=False))
    return 0 if type_assignments_match and local_manual_hashes_match else 1


if __name__ == "__main__":
    raise SystemExit(main())
