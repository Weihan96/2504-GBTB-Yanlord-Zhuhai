#!/usr/bin/env python3
"""Aggregate final provenance for every formal IFC Drawing without regenerating geometry."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.element


OFFICIAL_ELEVATION_RE = re.compile(r"^EL-\d{2}-\d{2}-")
EXPECTED_COUNTS = {
    "plan": 5,
    "official_elevation": 36,
    "project_compiled_elevation": 8,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc")
    )
    parser.add_argument(
        "--elevation-register",
        type=Path,
        default=Path("pipeline/decisions/int1-elevation-view-register.csv"),
    )
    parser.add_argument(
        "--source-evidence",
        type=Path,
        default=Path("pipeline/decisions/source-evidence-register.csv"),
    )
    parser.add_argument(
        "--equipment-register",
        type=Path,
        default=Path("pipeline/decisions/equipment-register.csv"),
    )
    parser.add_argument(
        "--installation-requirements",
        type=Path,
        default=Path("pipeline/decisions/equipment-installation-requirements.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("build/int1/int1-drawing-integrity-report.json"),
    )
    return parser.parse_args()


def resolve(root: Path, path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def relative_path(root: Path, path: Path) -> str:
    return path.resolve().relative_to(root).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def drawing_kind(name: str, target_view: str) -> str:
    if target_view == "PLAN_VIEW":
        return "plan"
    if target_view != "ELEVATION_VIEW":
        raise RuntimeError(f"unsupported Drawing target view: {name} / {target_view}")
    if name.startswith("EL-P"):
        return "project_compiled_elevation"
    if OFFICIAL_ELEVATION_RE.match(name):
        return "official_elevation"
    raise RuntimeError(f"unclassified elevation Drawing: {name}")


def source_json_path(root: Path, kind: str, name: str) -> Path | None:
    if kind == "official_elevation":
        return root / "build/int1/native-bonsai" / f"{name}-source.json"
    if kind == "project_compiled_elevation":
        return root / "build/int1/public-native" / f"{name}-source.json"
    return None


def drawing_document_relation(model: Any, drawing: Any) -> tuple[Any, Any]:
    references = [
        relation.RelatingDocument
        for relation in model.by_type("IfcRelAssociatesDocument")
        if drawing in relation.RelatedObjects
        and relation.RelatingDocument.is_a("IfcDocumentReference")
    ]
    if len(references) != 1:
        raise RuntimeError(
            f"{drawing.Name}: expected one IfcDocumentReference, found {len(references)}"
        )
    reference = references[0]
    information = reference.ReferencedDocument
    if information is None or not information.is_a("IfcDocumentInformation"):
        raise RuntimeError(f"{drawing.Name}: document reference has no information")
    if information.Name != drawing.Name or information.Scope != "DRAWING":
        raise RuntimeError(f"{drawing.Name}: document information identity mismatch")
    return reference, information


def artifact_record(root: Path, path: Path) -> dict[str, Any]:
    return {
        "path": relative_path(root, path),
        "sha256": sha256(path),
        "passes": True,
    }


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    ifc_path = resolve(root, args.ifc)
    output_path = resolve(root, args.output)
    dependency_paths = [
        ifc_path,
        resolve(root, args.elevation_register),
        resolve(root, args.source_evidence),
        resolve(root, args.equipment_register),
        resolve(root, args.installation_requirements),
    ]
    for path in dependency_paths:
        if not path.is_file():
            raise RuntimeError(f"missing INT1 Drawing dependency: {path}")

    model = ifcopenshell.open(ifc_path)
    drawings = sorted(
        (
            drawing
            for drawing in model.by_type("IfcAnnotation")
            if drawing.ObjectType == "DRAWING"
        ),
        key=lambda drawing: drawing.Name or "",
    )
    names = [drawing.Name or "" for drawing in drawings]
    if len(names) != len(set(names)):
        raise RuntimeError("formal IFC contains duplicate Drawing names")

    records: list[dict[str, Any]] = []
    output_records: list[dict[str, Any]] = []
    for drawing in drawings:
        name = drawing.Name or ""
        pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing") or {}
        target_view = pset.get("TargetView") or ""
        kind = drawing_kind(name, target_view)
        reference, information = drawing_document_relation(model, drawing)
        if not reference.Location:
            raise RuntimeError(f"{name}: document reference has no SVG location")
        svg_path = resolve(root, Path(reference.Location))
        if not svg_path.is_file():
            raise RuntimeError(f"{name}: referenced SVG does not exist: {svg_path}")
        svg_hash = sha256(svg_path)
        output_records.append(artifact_record(root, svg_path))

        source_path = source_json_path(root, kind, name)
        source_record: dict[str, Any] | None = None
        if source_path is not None:
            if not source_path.is_file():
                raise RuntimeError(f"{name}: source JSON does not exist: {source_path}")
            source = json.loads(source_path.read_text(encoding="utf-8"))
            source_drawing = source.get("drawing") or {}
            if (
                source_drawing.get("id") != drawing.id()
                or source_drawing.get("global_id") != drawing.GlobalId
                or source_drawing.get("name") != name
                or source.get("svg") != relative_path(root, svg_path)
                or source.get("pass") is not True
            ):
                raise RuntimeError(f"{name}: source JSON identity or pass gate mismatch")
            output_records.append(artifact_record(root, source_path))
            source_record = {
                "path": relative_path(root, source_path),
                "sha256": sha256(source_path),
                "generator_svg_sha256": source.get("svg_sha256", ""),
                "generator_svg_matches_final": source.get("svg_sha256") == svg_hash,
                "generation_ifc_sha256_before_save": source.get(
                    "formal_ifc_sha256_before_save", ""
                ),
            }

        records.append(
            {
                "name": name,
                "kind": kind,
                "target_view": target_view,
                "scale": pset.get("Scale") or "",
                "drawing_ifc_id": drawing.id(),
                "drawing_global_id": drawing.GlobalId,
                "document_information_ifc_id": information.id(),
                "document_reference_ifc_id": reference.id(),
                "svg": relative_path(root, svg_path),
                "svg_sha256": svg_hash,
                "source_json": source_record,
            }
        )

    counts = Counter(record["kind"] for record in records)
    if dict(counts) != EXPECTED_COUNTS:
        raise RuntimeError(
            f"unexpected Drawing composition: expected {EXPECTED_COUNTS}, got {dict(counts)}"
        )
    if len(records) != 49 or len(output_records) != 93:
        raise RuntimeError(
            f"unexpected Drawing/artifact totals: {len(records)} / {len(output_records)}"
        )

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_int1_drawing_integrity",
        "status": True,
        "source_ifc_sha256": sha256(ifc_path),
        "source_dependencies": [
            {"path": relative_path(root, path), "sha256": sha256(path)}
            for path in dependency_paths
        ],
        "summary": {
            "drawing_count": len(records),
            "plan_count": counts["plan"],
            "official_elevation_count": counts["official_elevation"],
            "project_compiled_elevation_count": counts[
                "project_compiled_elevation"
            ],
            "drawing_document_relation_count": len(records),
            "final_svg_count": len(records),
            "source_json_count": sum(
                record["source_json"] is not None for record in records
            ),
            "generator_svg_hash_match_count": sum(
                bool(record["source_json"])
                and record["source_json"]["generator_svg_matches_final"]
                for record in records
            ),
            "final_provenance_is_aggregate_report": True,
        },
        "gates": {
            "drawing_integrity_pass": True,
            "all_document_relations_resolve": True,
            "all_final_svg_artifacts_hashed": True,
            "all_elevation_source_json_artifacts_hashed": True,
            "construction_release_ready": False,
            "fabrication_dimension_ready": False,
        },
        "output_records": output_records,
        "drawings": records,
        "blockers": [
            {
                "issue_id": "INT1-FABRICATION-INTERFACES",
                "scope": "I-501/I-502/I-503 product and project interfaces",
                "stop_condition": "close project XYZ, shop drawings, finish datums, and manufacturer/site interfaces",
            }
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "report": str(output_path),
                "summary": report["summary"],
                "gates": report["gates"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
