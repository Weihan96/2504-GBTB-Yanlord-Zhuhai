#!/usr/bin/env python3
"""Classify the exact C003→INT1 proxy handoff without writing IFC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.element


HANDOFF_DECISION_ID = "COORD-HANDOFF-C003-INT1"
EXPECTED_COUNT = 16
CSV_FIELDS = [
    "global_id",
    "sheet_id",
    "ifc_class",
    "name",
    "object_type",
    "predefined_type",
    "container",
    "representation_types",
    "candidate_role",
    "bbox_min_mm",
    "bbox_max_mm",
    "dimensions_mm",
    "current_origin_mm",
    "origin_residual_mm",
    "shape_status",
    "integer_geometry_anchor",
    "basis",
    "confidence",
    "review_required",
    "review_status",
    "automatic_ifc_write_allowed",
    "source_ifc_sha256",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_report(path: Path, source_hash: str, label: str) -> dict[str, Any]:
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("source", {}).get("sha256") != source_hash:
        raise RuntimeError(f"{label} does not match the formal IFC")
    return report


def handoff_ids(rows: list[dict[str, str]]) -> list[str]:
    matches = [row for row in rows if row["decision_id"] == HANDOFF_DECISION_ID]
    if len(matches) != 1:
        raise RuntimeError("C003→INT1 handoff decision must be unique")
    row = matches[0]
    if row["status"] != "delegated":
        raise RuntimeError("C003→INT1 handoff must remain delegated until INT1 writes are approved")
    ids = [value.strip() for value in row["object_guid"].split(";") if value.strip()]
    if len(ids) != EXPECTED_COUNT or len(set(ids)) != EXPECTED_COUNT:
        raise RuntimeError(f"expected {EXPECTED_COUNT} unique C003→INT1 objects")
    return ids


def container_name(product: ifcopenshell.entity_instance) -> str:
    container = ifcopenshell.util.element.get_container(product)
    return str(getattr(container, "Name", "") or "")


def representation_types(product: ifcopenshell.entity_instance) -> list[str]:
    representation = getattr(product, "Representation", None)
    if representation is None:
        return []
    return sorted({
        "/".join(filter(None, (
            str(getattr(item, "RepresentationIdentifier", "") or ""),
            str(getattr(item, "RepresentationType", "") or ""),
        )))
        for item in representation.Representations
    })


def world_bbox(
    settings: ifcopenshell.geom.settings,
    product: ifcopenshell.entity_instance,
) -> tuple[list[float], list[float], list[float]] | None:
    try:
        shape = ifcopenshell.geom.create_shape(settings, product)
    except RuntimeError:
        return None
    values = list(shape.geometry.verts)
    if not values:
        return None
    vertices = [
        [values[index] * 1000.0, values[index + 1] * 1000.0, values[index + 2] * 1000.0]
        for index in range(0, len(values), 3)
    ]
    minimum = [min(point[axis] for point in vertices) for axis in range(3)]
    maximum = [max(point[axis] for point in vertices) for axis in range(3)]
    dimensions = [maximum[axis] - minimum[axis] for axis in range(3)]
    return tuple([round(value, 6) for value in row] for row in (minimum, maximum, dimensions))  # type: ignore[return-value]


def sheet_id(container: str) -> str:
    if container in {"KITCHEN", "VVD"}:
        return "I-501"
    if container in {"NBW", "BATHM"}:
        return "I-502"
    return "I-504"


def classify(
    product: ifcopenshell.entity_instance,
    reps: list[str],
    dimensions: list[float] | None,
) -> tuple[str, str, float, str]:
    name = str(product.Name or "")
    object_type = str(getattr(product, "ObjectType", "") or "")
    lowered = name.lower()
    if "countertop" in lowered or "worktop" in lowered:
        return (
            "named_worktop_or_countertop_candidate",
            f"existing name explicitly identifies {name}; fabrication thickness, edges and openings remain unverified",
            0.95,
            "classified_read_only",
        )
    if "pipe wall" in lowered:
        return (
            "named_bathroom_pipe_wall_candidate",
            "existing name explicitly identifies a bathroom pipe-wall component; layer build-up and access remain unverified",
            0.95,
            "classified_read_only",
        )
    if object_type == "CAD" or any("Curve3D" in item for item in reps):
        return (
            "legacy_cad_reference_needs_disposition",
            "ObjectType=CAD or Body/Curve3D identifies legacy linework/reference geometry, not approved joinery fabrication geometry",
            0.90,
            "legacy_reference_needs_disposition",
        )
    if dimensions is not None:
        x_size, y_size, z_size = dimensions
        if z_size <= 70.0 and x_size >= 300.0 and y_size >= 300.0:
            return (
                "horizontal_joinery_panel_candidate",
                f"world envelope {x_size:.3f}×{y_size:.3f}×{z_size:.3f} mm is a thin horizontal panel at the existing model position; exact countertop, shelf or other role remains unverified",
                0.80,
                "geometry_role_candidate_identity_pending",
            )
        if min(x_size, y_size) <= 70.0 and z_size >= 300.0:
            return (
                "vertical_joinery_panel_candidate",
                f"world envelope {x_size:.3f}×{y_size:.3f}×{z_size:.3f} mm is a thin vertical panel at the existing model position; exact backing, side panel or finish role remains unverified",
                0.80,
                "geometry_role_candidate_identity_pending",
            )
        if x_size >= 500.0 and y_size >= 500.0 and z_size >= 1800.0:
            return (
                "full_height_joinery_volume_candidate",
                f"world envelope {x_size:.3f}×{y_size:.3f}×{z_size:.3f} mm is a full-height joinery-sized volume; internal cabinet or enclosure role remains unverified",
                0.75,
                "geometry_role_candidate_identity_pending",
            )
    return (
        "unresolved_joinery_proxy",
        f"generic name {name or '<empty>'} and no assigned IFC type do not prove a cabinet, panel, hardware or worktop role",
        0.50,
        "identity_requires_human_review",
    )


def encode(value: Any) -> Any:
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, bool):
        return "yes" if value else "no"
    return value


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow({field: encode(record.get(field, "")) for field in CSV_FIELDS})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--p0-review", type=Path, required=True)
    parser.add_argument("--origin-review", type=Path, required=True)
    parser.add_argument("--anchor-audit", type=Path, required=True)
    parser.add_argument("--decision-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source_hash = sha256(args.input)
    model = ifcopenshell.open(args.input)
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    ids = handoff_ids(read_csv(args.p0_review))
    origin_report = read_report(args.origin_review, source_hash, "remaining-origin review")
    anchor_report = read_report(args.anchor_audit, source_hash, "integer-geometry anchor audit")
    origins = {row["global_id"]: row for row in origin_report["records"]}
    anchors = {row["global_id"]: row for row in anchor_report["records"]}
    if set(ids) - origins.keys() or set(ids) - anchors.keys():
        raise RuntimeError("C003→INT1 set is incomplete in current-hash origin reports")

    records = []
    for global_id in ids:
        product = model.by_guid(global_id)
        if product is None:
            raise RuntimeError(f"missing C003→INT1 object: {global_id}")
        origin = origins[global_id]
        anchor = anchors[global_id]
        reps = representation_types(product)
        bbox = world_bbox(settings, product) if anchor["shape_status"] == "ok" else None
        minimum, maximum, dimensions = bbox if bbox is not None else (None, None, None)
        role, basis, confidence, status = classify(product, reps, dimensions)
        container = container_name(product)
        records.append({
            "global_id": global_id,
            "sheet_id": sheet_id(container),
            "ifc_class": product.is_a(),
            "name": product.Name or "",
            "object_type": getattr(product, "ObjectType", None),
            "predefined_type": getattr(product, "PredefinedType", None),
            "container": container,
            "representation_types": reps,
            "candidate_role": role,
            "bbox_min_mm": minimum,
            "bbox_max_mm": maximum,
            "dimensions_mm": dimensions,
            "current_origin_mm": origin["current_mm"],
            "origin_residual_mm": origin["max_residual_mm"],
            "shape_status": anchor["shape_status"],
            "integer_geometry_anchor": anchor.get("anchor"),
            "basis": basis,
            "confidence": confidence,
            "review_required": True,
            "review_status": status,
            "automatic_ifc_write_allowed": False,
            "source_ifc_sha256": source_hash,
        })

    counts = Counter(record["candidate_role"] for record in records)
    sheet_counts = Counter(record["sheet_id"] for record in records)
    summary = {
        "handoff_object_count": len(records),
        "records_by_sheet": dict(sorted(sheet_counts.items())),
        "named_worktop_or_countertop_candidates": counts["named_worktop_or_countertop_candidate"],
        "named_bathroom_pipe_wall_candidates": counts["named_bathroom_pipe_wall_candidate"],
        "legacy_cad_references": counts["legacy_cad_reference_needs_disposition"],
        "horizontal_joinery_panel_candidates": counts["horizontal_joinery_panel_candidate"],
        "vertical_joinery_panel_candidates": counts["vertical_joinery_panel_candidate"],
        "full_height_joinery_volume_candidates": counts["full_height_joinery_volume_candidate"],
        "unclassified_joinery_proxies": counts["unresolved_joinery_proxy"],
        "shape_error_count": sum(record["shape_status"] != "ok" for record in records),
        "integer_geometry_anchor_count": sum(record["integer_geometry_anchor"] is not None for record in records),
    }
    gates = {
        "source_reports_match_formal_ifc": True,
        "handoff_set_matches_exactly": len(records) == EXPECTED_COUNT,
        "all_objects_classified": sum(counts.values()) == EXPECTED_COUNT,
        "all_objects_have_geometric_or_named_candidate_role": counts["unresolved_joinery_proxy"] == 0,
        "all_objects_have_fabrication_identity": False,
        "all_geometry_readable_as_solid_mesh": summary["shape_error_count"] == 0,
        "automatic_ifc_write_allowed": False,
        "fabrication_dimensions_ready": False,
    }
    report = {
        "mode": "read_only_int1_c003_handoff_candidate",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source": {
            "ifc": str(args.input),
            "ifc_sha256": source_hash,
            "p0_review": str(args.p0_review),
            "origin_review": str(args.origin_review),
            "anchor_audit": str(args.anchor_audit),
        },
        "summary": summary,
        "gates": gates,
        "records": records,
    }
    write_csv(args.decision_csv, records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": summary, "gates": gates}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
