#!/usr/bin/env python3
"""Read-only staged release gate rendered as one JSON document on stdout."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Iterable

import ifcopenshell


STAGES = ("reviewed-candidate", "construction-release-candidate")
REGISTER_REQUIRED_FIELDS = {
    "sheet_number",
    "title",
    "status",
    "publish_target",
}
IFC_HASH_KEYS = {
    "ifc_sha256",
    "ifcsha256",
    "source_ifc_sha256",
    "sourceifcsha256",
    "sourcesha256",
}
DIRECT_OUTPUT_KEYS = {
    "output",
    "output_pdf",
    "output_svg",
    "proof_png",
    "publish_target",
    "render_report",
}
CONSTRUCTION_FALSE_BLOCKER_KEYS = {
    "construction_release_pass",
    "construction_release_ready",
    "fabrication_dimension_ready",
    "final_release_pass",
    "network_positioning_complete",
    "releasable",
    "switch_positioning_complete",
    "whole_home_network_positioning_complete",
    "whole_home_socket_positioning_complete",
    "whole_home_switch_positioning_complete",
}
CONSTRUCTION_NONEMPTY_BLOCKER_KEYS = {
    "blockers",
    "blocking_input_ids",
    "open_release_items",
    "release_blocker",
    "release_blockers",
    "release_blocks",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=STAGES, required=True)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--ifc",
        type=Path,
        default=Path("2504 GBTB Yanlord Zhuhai.ifc"),
    )
    parser.add_argument(
        "--drawing-register",
        type=Path,
        default=Path("pipeline/decisions/drawing-register.csv"),
    )
    parser.add_argument(
        "--report-register",
        type=Path,
        default=Path("pipeline/decisions/release-report-register.csv"),
        help="Canonical professional report inventory for the selected stage.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        action="append",
        default=[],
        help="Optional extra professional JSON report; canonical reports cannot be omitted.",
    )
    return parser.parse_args()


def resolve_input(root: Path, path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_key(key: str) -> str:
    return key.replace("-", "_").lower()


def extract_ifc_hashes(value: Any, parent: dict[str, Any] | None = None) -> list[str]:
    hashes: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            normalized = normalized_key(key)
            if isinstance(child, str) and normalized in IFC_HASH_KEYS:
                hashes.append(child)
            elif (
                isinstance(child, str)
                and normalized == "sha256"
                and isinstance(value.get("ifc"), str)
            ):
                hashes.append(child)
            hashes.extend(extract_ifc_hashes(child, value))
    elif isinstance(value, list):
        for child in value:
            hashes.extend(extract_ifc_hashes(child, parent))
    return sorted(set(hashes))


def path_like_output(value: str) -> bool:
    return bool(value) and "://" not in value and not value.startswith("data:")


def extract_declared_outputs(
    value: Any,
    *,
    inside_outputs: bool = False,
) -> list[str]:
    outputs: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            normalized = normalized_key(key)
            child_inside_outputs = inside_outputs or normalized == "outputs"
            if isinstance(child, str):
                is_hash = "sha256" in normalized or normalized.endswith("_hash")
                is_output = normalized in DIRECT_OUTPUT_KEYS or (
                    child_inside_outputs and not is_hash
                )
                if is_output and path_like_output(child):
                    outputs.append(child)
            else:
                outputs.extend(
                    extract_declared_outputs(
                        child,
                        inside_outputs=child_inside_outputs,
                    )
                )
    elif isinstance(value, list):
        for child in value:
            if isinstance(child, str) and inside_outputs and path_like_output(child):
                outputs.append(child)
            else:
                outputs.extend(
                    extract_declared_outputs(child, inside_outputs=inside_outputs)
                )
    return sorted(set(outputs))


def construction_release_blockers(value: Any, path: str = "") -> list[str]:
    """Return explicit construction blockers without interpreting candidate-only flags."""
    blockers: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            normalized = normalized_key(key)
            child_path = f"{path}.{key}" if path else key
            if path == "" and normalized == "status" and child is False:
                blockers.append(f"{child_path}=false")
            if normalized in CONSTRUCTION_FALSE_BLOCKER_KEYS and child is False:
                blockers.append(f"{child_path}=false")
            if normalized == "status" and isinstance(child, str) and child.lower() == "block":
                blockers.append(f"{child_path}=block")
            if normalized in CONSTRUCTION_NONEMPTY_BLOCKER_KEYS and child not in (None, False, "", [], {}):
                blockers.append(f"{child_path}=nonempty")
            blockers.extend(construction_release_blockers(child, child_path))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            blockers.extend(construction_release_blockers(child, f"{path}[{index}]"))
    return sorted(set(blockers))


def make_check(
    check_id: str,
    status: str,
    message: str,
    details: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "id": check_id,
        "status": status,
        "message": message,
        "details": details or {},
    }


def read_register(path: Path) -> tuple[list[dict[str, str]], list[str], list[str]]:
    errors: list[str] = []
    if not path.is_file():
        return [], [], [f"drawing register does not exist: {path}"]
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            missing_fields = sorted(REGISTER_REQUIRED_FIELDS - set(fields))
            if missing_fields:
                errors.append(
                    "drawing register is missing fields: " + ", ".join(missing_fields)
                )
            rows = [
                {key: (value or "").strip() for key, value in row.items() if key}
                for row in reader
            ]
    except (OSError, csv.Error, UnicodeError) as exc:
        return [], [], [f"cannot read drawing register: {exc}"]

    sheet_numbers = [row.get("sheet_number", "") for row in rows]
    blank_rows = [str(index + 2) for index, sheet in enumerate(sheet_numbers) if not sheet]
    if blank_rows:
        errors.append("blank sheet_number at CSV rows: " + ", ".join(blank_rows))
    for required in ("title", "status"):
        blank = [
            str(index + 2)
            for index, row in enumerate(rows)
            if not row.get(required, "")
        ]
        if blank:
            errors.append(f"blank {required} at CSV rows: " + ", ".join(blank))
    counts = {sheet: sheet_numbers.count(sheet) for sheet in set(sheet_numbers) if sheet}
    duplicates = sorted(sheet for sheet, count in counts.items() if count > 1)
    if duplicates:
        errors.append("duplicate sheet_number values: " + ", ".join(duplicates))
    return rows, duplicates, errors


def read_report_register(path: Path, stage: str) -> tuple[list[Path], list[str]]:
    if not path.is_file():
        return [], [f"professional report register does not exist: {path}"]
    errors: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"report_id", "workstream", "report_path", "required_stage"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            return [], ["professional report register is missing fields: " + ", ".join(missing)]
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    ids = [row["report_id"] for row in rows]
    duplicates = sorted(value for value in set(ids) if ids.count(value) > 1)
    if any(not value for value in ids):
        errors.append("professional report register contains blank report_id")
    if duplicates:
        errors.append("duplicate professional report IDs: " + ", ".join(duplicates))
    allowed = set(STAGES)
    invalid = sorted({row["required_stage"] for row in rows} - allowed)
    if invalid:
        errors.append("invalid required_stage values: " + ", ".join(invalid))
    selected = [
        Path(row["report_path"])
        for row in rows
        if row["required_stage"] == "reviewed-candidate" or row["required_stage"] == stage
    ]
    if not selected:
        errors.append(f"professional report register selects no reports for {stage}")
    return selected, errors


def output_record(root: Path, declared: str, owner: str) -> dict[str, Any]:
    path = resolve_input(root, Path(declared))
    return {
        "owner": owner,
        "declared_path": declared,
        "resolved_path": str(path),
        "exists": path.is_file(),
    }


def validate_declared_artifact_records(
    root: Path,
    payload: Any,
    owner: str,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Validate a report's explicit top-level output_records contract when present."""
    if not isinstance(payload, dict) or "output_records" not in payload:
        return [], []
    raw_records = payload["output_records"]
    if not isinstance(raw_records, list):
        return [], [f"{owner}: output_records must be a list"]

    records: list[dict[str, Any]] = []
    errors: list[str] = []
    for index, raw in enumerate(raw_records):
        record_owner = f"{owner}:output_records[{index}]"
        if not isinstance(raw, dict):
            errors.append(f"{record_owner}: record must be an object")
            continue
        declared_path = raw.get("path")
        declared_hash = raw.get("sha256")
        declared_passes = raw.get("passes")
        if not isinstance(declared_path, str) or not path_like_output(declared_path):
            errors.append(f"{record_owner}: path must be a local file path")
            continue
        path = resolve_input(root, Path(declared_path))
        try:
            path.relative_to(root)
            inside_root = True
        except ValueError:
            inside_root = False
        exists = path.is_file()
        valid_hash = (
            isinstance(declared_hash, str)
            and len(declared_hash) == 64
            and all(character in "0123456789abcdef" for character in declared_hash)
        )
        actual_hash = sha256(path) if exists else ""
        hash_matches = valid_hash and exists and actual_hash == declared_hash
        record = {
            "owner": owner,
            "index": index,
            "declared_path": declared_path,
            "resolved_path": str(path),
            "inside_root": inside_root,
            "exists": exists,
            "declared_passes": declared_passes,
            "declared_sha256": declared_hash if isinstance(declared_hash, str) else "",
            "actual_sha256": actual_hash,
            "hash_matches": hash_matches,
        }
        records.append(record)
        if not inside_root:
            errors.append(f"{record_owner}: path resolves outside the project root")
        if not exists:
            errors.append(f"{record_owner}: artifact does not exist")
        if declared_passes is not True:
            errors.append(f"{record_owner}: passes must be true")
        if not valid_hash:
            errors.append(f"{record_owner}: sha256 must be 64 lowercase hexadecimal characters")
        elif exists and not hash_matches:
            errors.append(f"{record_owner}: artifact SHA-256 does not match")
    return records, errors


def evaluate(args: argparse.Namespace) -> tuple[dict[str, Any], int]:
    root = args.root.resolve()
    ifc_path = resolve_input(root, args.ifc)
    register_path = resolve_input(root, args.drawing_register)
    report_register_path = resolve_input(root, args.report_register)
    registered_reports, report_register_errors = read_report_register(report_register_path, args.stage)
    report_paths = []
    for path in [*registered_reports, *args.report]:
        resolved = resolve_input(root, path)
        if resolved not in report_paths:
            report_paths.append(resolved)
    checks: list[dict[str, Any]] = []
    errors: list[str] = []
    unresolved: list[dict[str, Any]] = []
    errors.extend(report_register_errors)
    checks.append(make_check(
        "PROFESSIONAL-REPORT-INVENTORY",
        "fail" if report_register_errors else "pass",
        "professional report inventory is complete"
        if not report_register_errors else "professional report inventory is invalid",
        {"report_count": len(report_paths), "errors": report_register_errors},
    ))

    if ifc_path.is_file():
        ifc_hash = sha256(ifc_path)
        try:
            ifc_model = ifcopenshell.open(str(ifc_path))
            entity_ids = ifc_model.wrapped_data.entity_names()
            ifc_details = {
                "schema": ifc_model.schema,
                "entity_count": len(entity_ids),
                "max_step_id": max(entity_ids, default=0),
                "root_count": len(ifc_model.by_type("IfcRoot")),
                "product_count": len(ifc_model.by_type("IfcProduct")),
            }
            checks.append(make_check("FORMAL-IFC", "pass", "formal IFC parses successfully", ifc_details))
        except Exception as exc:
            message = f"formal IFC cannot be parsed: {exc}"
            errors.append(message)
            checks.append(make_check("FORMAL-IFC", "fail", message))
    else:
        ifc_hash = ""
        message = f"formal IFC does not exist: {ifc_path}"
        errors.append(message)
        checks.append(make_check("FORMAL-IFC", "fail", message))

    rows, duplicates, register_errors = read_register(register_path)
    errors.extend(register_errors)
    checks.append(
        make_check(
            "DRAWING-NUMBERS",
            "fail" if register_errors else "pass",
            "drawing numbers are complete and unique"
            if not register_errors
            else "drawing register is incomplete or has duplicate sheet numbers",
            {"drawing_count": len(rows), "duplicates": duplicates},
        )
    )

    planned_rows = [row for row in rows if row.get("status") == "planned"]
    for row in planned_rows:
        unresolved.append(
            {
                "kind": "planned-drawing",
                "sheet_number": row.get("sheet_number", ""),
                "title": row.get("title", ""),
                "notes": row.get("notes", ""),
            }
        )
    if args.stage == "reviewed-candidate":
        checks.append(
            make_check(
                "PLANNED-DRAWINGS",
                "disclosed" if planned_rows else "pass",
                "planned drawings are disclosed and remain outside construction release"
                if planned_rows
                else "no planned drawings remain",
                {"planned_count": len(planned_rows)},
            )
        )
    else:
        if planned_rows:
            message = f"construction release has {len(planned_rows)} planned drawing(s)"
            errors.append(message)
            status = "fail"
        else:
            message = "no planned drawings remain"
            status = "pass"
        checks.append(
            make_check(
                "PLANNED-DRAWINGS",
                status,
                message,
                {"planned_count": len(planned_rows)},
            )
        )

    register_outputs: list[dict[str, Any]] = []
    for row in rows:
        target = row.get("publish_target", "")
        status = row.get("status", "")
        if target:
            register_outputs.append(
                output_record(root, target, f"drawing:{row.get('sheet_number', '')}")
            )
        elif status != "planned" or args.stage == "construction-release-candidate":
            register_outputs.append(
                {
                    "owner": f"drawing:{row.get('sheet_number', '')}",
                    "declared_path": "",
                    "resolved_path": "",
                    "exists": False,
                }
            )
    missing_register_outputs = [item for item in register_outputs if not item["exists"]]
    if missing_register_outputs:
        message = f"{len(missing_register_outputs)} required drawing output(s) are missing"
        errors.append(message)
        output_status = "fail"
    else:
        message = "all required drawing outputs exist"
        output_status = "pass"
    checks.append(
        make_check(
            "DRAWING-OUTPUTS",
            output_status,
            message,
            {
                "required_count": len(register_outputs),
                "missing": missing_register_outputs,
            },
        )
    )

    report_results: list[dict[str, Any]] = []
    report_outputs: list[dict[str, Any]] = []
    declared_artifact_records: list[dict[str, Any]] = []
    declared_artifact_errors: list[str] = []
    for report_path in report_paths:
        result: dict[str, Any] = {
            "path": str(report_path),
            "exists": report_path.is_file(),
            "ifc_hashes": [],
            "current_ifc_hash": False,
            "declared_outputs": [],
            "declared_artifact_records": [],
            "declared_artifact_errors": [],
            "construction_release_blockers": [],
        }
        if not report_path.is_file():
            result["error"] = "report does not exist"
            report_results.append(result)
            continue
        try:
            payload = json.loads(report_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            result["error"] = f"cannot read JSON report: {exc}"
            report_results.append(result)
            continue
        hashes = extract_ifc_hashes(payload)
        declared_outputs = extract_declared_outputs(payload)
        artifact_records, artifact_errors = validate_declared_artifact_records(
            root,
            payload,
            report_path.name,
        )
        result["ifc_hashes"] = hashes
        result["current_ifc_hash"] = bool(ifc_hash) and hashes == [ifc_hash]
        result["declared_outputs"] = declared_outputs
        result["declared_artifact_records"] = artifact_records
        result["declared_artifact_errors"] = artifact_errors
        result["construction_release_blockers"] = construction_release_blockers(payload)
        declared_artifact_records.extend(artifact_records)
        declared_artifact_errors.extend(artifact_errors)
        if not hashes:
            result["error"] = "report does not declare an IFC SHA-256"
        elif not result["current_ifc_hash"]:
            result["error"] = "report IFC SHA-256 does not match the formal IFC"
        for declared in declared_outputs:
            report_outputs.append(
                output_record(root, declared, f"report:{report_path.name}")
            )
        report_results.append(result)

    stale_reports = [
        result
        for result in report_results
        if not result["exists"] or not result["current_ifc_hash"]
    ]
    if stale_reports:
        message = f"{len(stale_reports)} professional report(s) are missing or stale"
        errors.append(message)
        report_status = "fail"
    else:
        message = "all professional reports match the formal IFC"
        report_status = "pass"
    checks.append(
        make_check(
            "PROFESSIONAL-REPORT-HASHES",
            report_status,
            message,
            {"reports": report_results},
        )
    )

    if declared_artifact_errors:
        message = (
            f"{len(declared_artifact_errors)} declared artifact integrity error(s) were found"
        )
        errors.append(message)
        artifact_status = "fail"
    else:
        message = "all explicitly hashed report artifacts match their declared SHA-256"
        artifact_status = "pass"
    checks.append(
        make_check(
            "REPORT-ARTIFACT-INTEGRITY",
            artifact_status,
            message,
            {
                "record_count": len(declared_artifact_records),
                "records": declared_artifact_records,
                "errors": declared_artifact_errors,
            },
        )
    )

    intrinsic_blocked_reports = [
        result
        for result in report_results
        if result["construction_release_blockers"]
    ]
    if args.stage == "construction-release-candidate":
        if intrinsic_blocked_reports:
            message = (
                f"{len(intrinsic_blocked_reports)} professional report(s) explicitly block construction release"
            )
            errors.append(message)
            readiness_status = "fail"
        else:
            message = "no professional report explicitly blocks construction release"
            readiness_status = "pass"
    else:
        message = "report-level construction blockers are disclosed but do not block reviewed candidates"
        readiness_status = "disclosed" if intrinsic_blocked_reports else "pass"
    checks.append(
        make_check(
            "PROFESSIONAL-REPORT-READINESS",
            readiness_status,
            message,
            {
                "blocked_reports": [
                    {
                        "path": item["path"],
                        "blockers": item["construction_release_blockers"],
                    }
                    for item in intrinsic_blocked_reports
                ]
            },
        )
    )

    missing_report_outputs = [item for item in report_outputs if not item["exists"]]
    if missing_report_outputs:
        message = f"{len(missing_report_outputs)} report-declared output(s) are missing"
        errors.append(message)
        report_output_status = "fail"
    else:
        message = "all report-declared outputs exist"
        report_output_status = "pass"
    checks.append(
        make_check(
            "REPORT-OUTPUTS",
            report_output_status,
            message,
            {
                "declared_count": len(report_outputs),
                "missing": missing_report_outputs,
            },
        )
    )

    passed = not errors
    construction_ready = (
        args.stage == "construction-release-candidate"
        and passed
        and not planned_rows
        and not intrinsic_blocked_reports
    )
    result = {
        "schema_version": "1.0",
        "command": "release-candidate-gate",
        "read_only": True,
        "stage": args.stage,
        "result": {
            "pass": passed,
            "reviewed_candidate_ready": passed,
            "construction_release_candidate_ready": construction_ready,
        },
        "source": {
            "root": str(root),
            "formal_ifc": str(ifc_path),
            "formal_ifc_sha256": ifc_hash,
            "drawing_register": str(register_path),
            "professional_report_register": str(report_register_path),
            "professional_reports": [str(path) for path in report_paths],
        },
        "summary": {
            "drawing_count": len(rows),
            "planned_count": len(planned_rows),
            "report_count": len(report_results),
            "fresh_report_count": sum(
                1 for item in report_results if item["current_ifc_hash"]
            ),
            "construction_blocked_report_count": len(intrinsic_blocked_reports),
            "drawing_output_count": len(register_outputs),
            "existing_drawing_output_count": sum(
                1 for item in register_outputs if item["exists"]
            ),
            "report_output_count": len(report_outputs),
            "existing_report_output_count": sum(
                1 for item in report_outputs if item["exists"]
            ),
            "declared_artifact_record_count": len(declared_artifact_records),
            "current_declared_artifact_record_count": sum(
                1 for item in declared_artifact_records if item["hash_matches"]
            ),
            "unresolved_count": len(unresolved),
            "error_count": len(errors),
        },
        "checks": checks,
        "unresolved": unresolved,
        "errors": errors,
    }
    return result, 0 if passed else 1


def main() -> int:
    args = parse_args()
    report, exit_code = evaluate(args)
    print(json.dumps(report, ensure_ascii=False, sort_keys=True))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
