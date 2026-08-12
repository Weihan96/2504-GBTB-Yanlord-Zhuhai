#!/usr/bin/env python3
"""Aggregate fail-closed ELEC/RCP1 reviewed-candidate evidence without Blender."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCOPE_SHEETS = {"E-301", "E-302", "E-303", "E-304", "A-106", "M-401"}
CRITICAL_REPORTS = {
    "ELEC-PROGRAM": Path("build/elec/elec-room-program-candidate.json"),
    "ELEC-POSITIONING": Path("build/elec/elec-positioning-candidate.json"),
    "ELEC-ROUND1": Path("build/elec/elec-renovation-round1-candidate.json"),
    "ELEC-CONTROL-NETWORK": Path("build/elec/elec-control-network-candidate.json"),
    "A106": Path("build/elec/a106-ceiling-device-candidate.json"),
    "RCP1-EXISTING": Path("build/rcp1/rcp1-existing-candidate.json"),
    "M401-EXISTING": Path("build/rcp1/m401-existing-report.json"),
    "M401-SHEET": Path("build/rcp1/m401-coordination-sheet-candidate.json"),
    "RCP1-ROUTE-READINESS": Path("build/rcp1/route-readiness-candidate.json"),
}
EXCLUDED_PREVIEWS = {
    "RCP1-HVAC-INTERFACE-OLD": Path("build/rcp1/hvac-interface-candidate.json"),
    "RCP1-HVAC-PREVIEW": Path("build/rcp1/hvac-route-preview-candidate.json"),
    "RCP1-HVAC-PREVIEW-AUDIT": Path("build/rcp1/hvac-route-preview-audit.json"),
    "RCP1-PAIR-TRIAGE-OLD": Path("build/rcp1/pair-triage.json"),
}
WRITE_AUTHORITY_KEYS = {
    "automatic_ifc_write_allowed",
    "formal_ifc_write_allowed",
    "ifc_write_allowed",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc")
    )
    parser.add_argument(
        "--drawing-register",
        type=Path,
        default=Path("pipeline/decisions/drawing-register.csv"),
    )
    parser.add_argument(
        "--release-register",
        type=Path,
        default=Path("pipeline/decisions/release-report-register.csv"),
    )
    parser.add_argument(
        "--e301-e303-render",
        type=Path,
        default=Path("build/elec/E-301-E-303-renovation-round1-render.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("build/reports/critical-path-evidence-candidate.json"),
    )
    return parser.parse_args()


def resolve(root: Path, path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def inside_project(root: Path, path: Path) -> bool:
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return True


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise RuntimeError(f"missing critical-path report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def extract_ifc_hashes(value: Any) -> list[str]:
    hashes: set[str] = set()
    if isinstance(value, dict):
        for key, child in value.items():
            normalized = key.replace("-", "_").lower()
            if isinstance(child, str) and normalized in {
                "ifc_sha256",
                "source_ifc_sha256",
                "formal_ifc_sha256",
            }:
                hashes.add(child)
            hashes.update(extract_ifc_hashes(child))
    elif isinstance(value, list):
        for child in value:
            hashes.update(extract_ifc_hashes(child))
    return sorted(hashes)


def extract_write_authorities(value: Any, path: str = "") -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            child_path = f"{path}.{key}" if path else key
            if key.replace("-", "_").lower() in WRITE_AUTHORITY_KEYS:
                results.append({"path": child_path, "value": child})
            results.extend(extract_write_authorities(child, child_path))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            results.extend(extract_write_authorities(child, f"{path}[{index}]"))
    return results


def write_authority_closed(value: Any) -> bool:
    return value is False or (
        isinstance(value, str) and value.strip().lower() in {"no", "false", "0"}
    )


def report_evidence(
    root: Path, report_id: str, relative_path: Path, formal_hash: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    path = resolve(root, relative_path)
    payload = read_json(path)
    hashes = extract_ifc_hashes(payload)
    authorities = extract_write_authorities(payload)
    non_closed_authorities = [
        item for item in authorities if not write_authority_closed(item["value"])
    ]
    result = {
        "report_id": report_id,
        "path": str(path),
        "sha256": sha256(path),
        "ifc_hashes": hashes,
        "current_formal_ifc": hashes == [formal_hash],
        "write_authority_count": len(authorities),
        "non_closed_write_authorities": non_closed_authorities,
        "formal_ifc_write_closed": bool(authorities) and not non_closed_authorities,
    }
    return result, payload


def read_scope_drawings(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = {
            row["sheet_number"]: {key: (value or "").strip() for key, value in row.items()}
            for row in csv.DictReader(handle)
            if row["sheet_number"] in SCOPE_SHEETS
        }
    if set(rows) != SCOPE_SHEETS:
        raise RuntimeError(
            f"critical-path drawing register mismatch: {sorted(set(rows) ^ SCOPE_SHEETS)}"
        )
    if any(row["status"] != "candidate" or not row["publish_target"] for row in rows.values()):
        raise RuntimeError("critical-path drawings must remain explicit candidates with outputs")
    return rows


def verified_file(root: Path, declared_path: str, declared_hash: str) -> dict[str, Any]:
    path = resolve(root, Path(declared_path))
    inside = inside_project(root, path)
    exists = path.is_file()
    actual_hash = sha256(path) if exists else ""
    return {
        "path": str(path),
        "inside_project": inside,
        "exists": exists,
        "declared_sha256": declared_hash,
        "actual_sha256": actual_hash,
        "hash_match": inside and exists and actual_hash == declared_hash,
    }


def render_bundle(
    root: Path, bundle_id: str, report_path: Path, expected_pdf: Path
) -> dict[str, Any]:
    path = resolve(root, report_path)
    payload = read_json(path)
    files = {
        "source_svg": verified_file(
            root, payload["source_svg"], payload["source_svg_sha256"]
        ),
        "output_pdf": verified_file(
            root, payload["output_pdf"], payload["output_pdf_sha256"]
        ),
        "proof_png": verified_file(
            root, payload["proof_png"], payload["proof_png_sha256"]
        ),
    }
    expected = resolve(root, expected_pdf)
    pass_gate = (
        payload.get("pass") is True
        and resolve(root, Path(payload["output_pdf"])) == expected
        and all(item["hash_match"] for item in files.values())
    )
    return {
        "bundle_id": bundle_id,
        "render_report": str(path),
        "render_report_sha256": sha256(path),
        "expected_publish_target": str(expected),
        "files": files,
        "mechanical_pass": pass_gate,
    }


def m401_bundle(root: Path, payload: dict[str, Any], expected_svg: Path) -> dict[str, Any]:
    files = {
        key: verified_file(root, value["path"], value["sha256"])
        for key, value in payload["outputs"].items()
    }
    expected = resolve(root, expected_svg)
    return {
        "bundle_id": "M401",
        "expected_publish_target": str(expected),
        "files": files,
        "mechanical_pass": (
            resolve(root, Path(payload["outputs"]["svg"]["path"])) == expected
            and all(item["hash_match"] for item in files.values())
        ),
    }


def register_paths(path: Path, field: str) -> set[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return {
            (row.get(field) or "").strip()
            for row in csv.DictReader(handle)
            if (row.get(field) or "").strip()
        }


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    ifc_path = resolve(root, args.ifc)
    formal_hash = sha256(ifc_path)
    drawing_register = resolve(root, args.drawing_register)
    release_register = resolve(root, args.release_register)
    drawings = read_scope_drawings(drawing_register)

    critical_reports = []
    payloads: dict[str, dict[str, Any]] = {}
    for report_id, relative_path in CRITICAL_REPORTS.items():
        evidence, payload = report_evidence(root, report_id, relative_path, formal_hash)
        critical_reports.append(evidence)
        payloads[report_id] = payload

    publish_bundles = [
        render_bundle(
            root,
            "E301-E303",
            args.e301_e303_render,
            Path(drawings["E-301"]["publish_target"]),
        ),
        render_bundle(
            root,
            "E302-E304",
            Path("build/elec/E-302-E-304-control-network-render.json"),
            Path(drawings["E-302"]["publish_target"]),
        ),
        render_bundle(
            root,
            "A106",
            Path("build/elec/A-106-ceiling-device-render.json"),
            Path(drawings["A-106"]["publish_target"]),
        ),
        m401_bundle(
            root,
            payloads["M401-SHEET"],
            Path(drawings["M-401"]["publish_target"]),
        ),
    ]
    shared_targets_match = (
        drawings["E-301"]["publish_target"] == drawings["E-303"]["publish_target"]
        and drawings["E-302"]["publish_target"] == drawings["E-304"]["publish_target"]
    )

    registered_reports = register_paths(release_register, "report_path")
    published_targets = register_paths(drawing_register, "publish_target")
    excluded_previews = []
    for preview_id, relative_path in EXCLUDED_PREVIEWS.items():
        path = resolve(root, relative_path)
        payload = read_json(path)
        hashes = extract_ifc_hashes(payload)
        excluded_previews.append(
            {
                "preview_id": preview_id,
                "path": str(path),
                "sha256": sha256(path),
                "ifc_hashes": hashes,
                "current_formal_ifc": hashes == [formal_hash],
                "registered_for_release": str(relative_path) in registered_reports,
                "used_as_publish_target": str(relative_path) in published_targets,
                "eligible_for_release": False,
                "exclusion_reason": (
                    "preview/diagnostic evidence is not a publish artifact; "
                    "a stale IFC hash further excludes this cached copy"
                    if hashes != [formal_hash]
                    else "preview/diagnostic evidence is not a publish artifact"
                ),
            }
        )

    reports_current = all(row["current_formal_ifc"] for row in critical_reports)
    write_closed = all(row["formal_ifc_write_closed"] for row in critical_reports)
    publish_pass = all(row["mechanical_pass"] for row in publish_bundles)
    previews_excluded = all(
        not row["registered_for_release"] and not row["used_as_publish_target"]
        for row in excluded_previews
    )
    reviewed_pass = all(
        (reports_current, write_closed, publish_pass, shared_targets_match, previews_excluded)
    )
    report = {
        "schema_version": "1.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_elec_rcp1_critical_path_evidence_candidate",
        "source_ifc_sha256": formal_hash,
        "source": {
            "formal_ifc": str(ifc_path),
            "formal_ifc_sha256": formal_hash,
            "drawing_register": str(drawing_register),
            "drawing_register_sha256": sha256(drawing_register),
            "release_register": str(release_register),
            "release_register_sha256": sha256(release_register),
        },
        "summary": {
            "scope_sheet_count": len(drawings),
            "critical_report_count": len(critical_reports),
            "current_critical_report_count": sum(
                row["current_formal_ifc"] for row in critical_reports
            ),
            "publish_bundle_count": len(publish_bundles),
            "mechanical_publish_bundle_pass_count": sum(
                row["mechanical_pass"] for row in publish_bundles
            ),
            "excluded_preview_count": len(excluded_previews),
            "stale_excluded_preview_count": sum(
                not row["current_formal_ifc"] for row in excluded_previews
            ),
        },
        "critical_reports": critical_reports,
        "publish_bundles": publish_bundles,
        "excluded_previews": excluded_previews,
        "gates": {
            "all_critical_reports_match_formal_ifc": reports_current,
            "all_critical_reports_close_formal_ifc_write": write_closed,
            "all_publish_bundles_mechanically_match": publish_pass,
            "shared_e301_e303_and_e302_e304_targets_match": shared_targets_match,
            "preview_and_diagnostic_caches_excluded_from_release": previews_excluded,
            "reviewed_candidate_evidence_pass": reviewed_pass,
            "construction_release_ready": False,
            "formal_ifc_write_allowed": False,
        },
        "release_blockers": [
            "ELEC switch, socket, network and product decisions remain open in source reports",
            "RCP1 final airside, refrigerant and condensate routes remain open",
            "stale Blender preview caches must not be used as current construction evidence",
        ],
    }
    output = resolve(root, args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False, indent=2))
    if not reviewed_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
