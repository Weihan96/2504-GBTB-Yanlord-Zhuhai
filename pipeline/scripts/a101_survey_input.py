#!/usr/bin/env python3
"""Validate A-101 site-survey inputs and compile read-only normalized reports.

The command is stdout-only by default. It reads the formal IFC only to record
its SHA-256 and never opens or writes IFC model data. JSON/CSV reports are
written only when the caller supplies explicit output paths.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from datetime import date
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
DEFAULT_INPUT = ROOT / "pipeline/templates/a101-survey-input.csv"
SCHEMA_VERSION = "a101-survey-input/v1"

HEADERS = [
    "input_id",
    "field_type",
    "location",
    "measurement",
    "value",
    "unit",
    "status",
    "measurement_method",
    "observed_by",
    "observed_at",
    "evidence_reference",
    "confirmed_by",
    "confirmed_at",
    "notes",
]

NORMALIZED_HEADERS = [
    "input_id",
    "field_type",
    "location",
    "measurement",
    "recorded_value",
    "unit",
    "status",
    "effective_value",
    "effective_unit",
    "measurement_method",
    "observed_by",
    "observed_at",
    "evidence_reference",
    "confirmed_by",
    "confirmed_at",
    "notes",
]

STATUSES = {"pending", "observed", "confirmed"}
DIMENSION_TYPES = {
    "main_bay",
    "main_depth",
    "clear_height",
    "beam_soffit",
    "window_sill",
}
FIELD_UNITS = {field_type: "mm" for field_type in DIMENSION_TYPES} | {
    "immovable_property_condition": "n/a",
}
DIMENSION_METHODS = {"laser_distance_meter", "tape_measure", "level_instrument"}
PROPERTY_METHODS = {"visual_inspection", "property_document"}
REQUIRED_INPUTS = {
    "A101-MAIN-BAY-01": "main_bay",
    "A101-MAIN-DEPTH-01": "main_depth",
    "A101-CLEAR-HEIGHT-01": "clear_height",
    "A101-BEAM-SOFFIT-01": "beam_soffit",
    "A101-WINDOW-SILL-01": "window_sill",
    "A101-PROPERTY-CONDITION-01": "immovable_property_condition",
}
INPUT_ID_PATTERN = re.compile(r"^A101-[A-Z0-9]+(?:-[A-Z0-9]+)*$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--ifc", type=Path, default=DEFAULT_IFC)
    parser.add_argument(
        "--expected-ifc-sha256",
        help="Optional caller-frozen SHA-256; mismatch is a hard failure.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        help="Optional explicit path for the normalized JSON report.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        help="Optional explicit path for the normalized CSV report.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_input(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != HEADERS:
            raise ValueError(
                "input headers changed; expected exact schema: " + ",".join(HEADERS)
            )
        rows = [
            {header: (row.get(header) or "").strip() for header in HEADERS}
            for row in reader
            if any((value or "").strip() for value in row.values())
        ]
    return rows


def parse_iso_date(value: str, field: str, input_id: str, errors: list[str]) -> None:
    try:
        date.fromisoformat(value)
    except ValueError:
        errors.append(f"{input_id}: {field} must be YYYY-MM-DD")


def normalized_dimension(value: str, input_id: str, errors: list[str]) -> int | float | None:
    try:
        number = Decimal(value)
    except InvalidOperation:
        errors.append(f"{input_id}: dimensional value must be numeric millimetres")
        return None
    if not number.is_finite() or number <= 0:
        errors.append(f"{input_id}: dimensional value must be greater than zero")
        return None
    if number == number.to_integral_value():
        return int(number)
    return float(number)


def validate_rows(rows: list[dict[str, str]]) -> tuple[list[str], list[dict[str, Any]]]:
    errors: list[str] = []
    normalized: list[dict[str, Any]] = []
    if not rows:
        return ["input contains no survey rows"], []

    ids = [row["input_id"] for row in rows]
    duplicates = sorted({input_id for input_id in ids if ids.count(input_id) > 1})
    if duplicates:
        errors.append("duplicate input_id values: " + ", ".join(duplicates))

    by_id = {row["input_id"]: row for row in rows}
    for required_id, expected_type in REQUIRED_INPUTS.items():
        row = by_id.get(required_id)
        if row is None:
            errors.append(f"missing required survey row: {required_id}")
        elif row["field_type"] != expected_type:
            errors.append(
                f"{required_id}: field_type must remain {expected_type!r}"
            )

    for row in rows:
        input_id = row["input_id"] or "<blank>"
        field_type = row["field_type"]
        status = row["status"]

        if not INPUT_ID_PATTERN.fullmatch(row["input_id"]):
            errors.append(f"{input_id}: invalid input_id")
        if field_type not in FIELD_UNITS:
            errors.append(f"{input_id}: invalid field_type {field_type!r}")
            expected_unit = None
        else:
            expected_unit = FIELD_UNITS[field_type]
            if row["unit"] != expected_unit:
                errors.append(
                    f"{input_id}: unit must be {expected_unit!r} for {field_type}"
                )
        if not row["location"]:
            errors.append(f"{input_id}: location is required")
        if not row["measurement"]:
            errors.append(f"{input_id}: measurement is required")
        if status not in STATUSES:
            errors.append(f"{input_id}: invalid status {status!r}")

        recorded_value: str | int | float | None = None
        if row["value"]:
            if field_type in DIMENSION_TYPES:
                recorded_value = normalized_dimension(row["value"], input_id, errors)
            else:
                recorded_value = row["value"]

        observation_fields = (
            "measurement_method",
            "observed_by",
            "observed_at",
            "evidence_reference",
        )
        confirmation_fields = ("confirmed_by", "confirmed_at")

        if status == "pending":
            populated = [
                field
                for field in ("value", *observation_fields, *confirmation_fields)
                if row[field]
            ]
            if populated:
                errors.append(
                    f"{input_id}: pending row must not contain measured or confirmed fields: "
                    + ", ".join(populated)
                )
        elif status in {"observed", "confirmed"}:
            if not row["value"]:
                errors.append(f"{input_id}: {status} row requires value")
            missing_observation = [field for field in observation_fields if not row[field]]
            if missing_observation:
                errors.append(
                    f"{input_id}: {status} row missing observation fields: "
                    + ", ".join(missing_observation)
                )
            if row["observed_at"]:
                parse_iso_date(row["observed_at"], "observed_at", input_id, errors)

            allowed_methods = (
                DIMENSION_METHODS if field_type in DIMENSION_TYPES else PROPERTY_METHODS
            )
            if row["measurement_method"] not in allowed_methods:
                errors.append(
                    f"{input_id}: invalid measurement_method {row['measurement_method']!r}"
                )

            if status == "observed":
                populated_confirmation = [field for field in confirmation_fields if row[field]]
                if populated_confirmation:
                    errors.append(
                        f"{input_id}: observed row must not contain confirmation fields: "
                        + ", ".join(populated_confirmation)
                    )
            else:
                missing_confirmation = [field for field in confirmation_fields if not row[field]]
                if missing_confirmation:
                    errors.append(
                        f"{input_id}: confirmed row missing confirmation fields: "
                        + ", ".join(missing_confirmation)
                    )
                if row["confirmed_at"]:
                    parse_iso_date(row["confirmed_at"], "confirmed_at", input_id, errors)

        effective_value = recorded_value if status == "confirmed" else None
        normalized.append(
            {
                "input_id": row["input_id"],
                "field_type": field_type,
                "location": row["location"],
                "measurement": row["measurement"],
                "recorded_value": recorded_value,
                "unit": expected_unit or row["unit"],
                "status": status,
                "effective_value": effective_value,
                "effective_unit": expected_unit if effective_value is not None else None,
                "measurement_method": row["measurement_method"] or None,
                "observed_by": row["observed_by"] or None,
                "observed_at": row["observed_at"] or None,
                "evidence_reference": row["evidence_reference"] or None,
                "confirmed_by": row["confirmed_by"] or None,
                "confirmed_at": row["confirmed_at"] or None,
                "notes": row["notes"] or None,
            }
        )

    return errors, normalized


def validate_output_path(path: Path | None, suffix: str, ifc: Path) -> None:
    if path is None:
        return
    if path.suffix.lower() != suffix:
        raise ValueError(f"output path must end with {suffix}: {path}")
    if path.resolve() == ifc.resolve() or path.suffix.lower() == ".ifc":
        raise ValueError("output path must never target the formal IFC")


def build_report(
    rows: list[dict[str, Any]],
    *,
    input_path: Path,
    ifc_path: Path,
    ifc_hash: str,
    expected_hash: str | None,
) -> dict[str, Any]:
    counts = {status: sum(row["status"] == status for row in rows) for status in sorted(STATUSES)}
    blockers = [row["input_id"] for row in rows if row["status"] != "confirmed"]
    return {
        "schema_version": SCHEMA_VERSION,
        "sheet_number": "A-101",
        "report_kind": "site_survey_input",
        "source_input": str(input_path),
        "source_ifc": str(ifc_path),
        "source_ifc_sha256": ifc_hash,
        "caller_frozen_ifc_sha256": expected_hash,
        "formal_ifc_write": False,
        "stdout_only_by_default": True,
        "construction_ready": not blockers,
        "gates": {
            "source_ifc_hash_current": True,
            "pending_inputs_explicit": bool(blockers),
            "automatic_ifc_write_allowed": False,
            "construction_release_ready": not blockers,
        },
        "summary": {
            "row_count": len(rows),
            "required_row_count": len(REQUIRED_INPUTS),
            "status_counts": counts,
            "blocking_input_ids": blockers,
        },
        "inputs": rows,
    }


def csv_value(value: Any) -> str:
    return "" if value is None else str(value)


def write_normalized_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=NORMALIZED_HEADERS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({header: csv_value(row.get(header)) for header in NORMALIZED_HEADERS})


def main() -> int:
    args = parse_args()
    try:
        validate_output_path(args.output_json, ".json", args.ifc)
        validate_output_path(args.output_csv, ".csv", args.ifc)
        if not args.ifc.is_file():
            raise ValueError(f"IFC does not exist: {args.ifc}")
        ifc_hash = sha256(args.ifc)
        expected_hash = args.expected_ifc_sha256.lower() if args.expected_ifc_sha256 else None
        if expected_hash is not None:
            if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
                raise ValueError("--expected-ifc-sha256 must be 64 hexadecimal characters")
            if expected_hash != ifc_hash:
                raise ValueError(
                    f"IFC SHA-256 mismatch: expected {expected_hash}, actual {ifc_hash}"
                )

        source_rows = read_input(args.input)
        errors, normalized_rows = validate_rows(source_rows)
        if errors:
            raise ValueError("\n".join(errors))

        report = build_report(
            normalized_rows,
            input_path=args.input,
            ifc_path=args.ifc,
            ifc_hash=ifc_hash,
            expected_hash=expected_hash,
        )
        rendered = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
        if args.output_json is not None:
            args.output_json.parent.mkdir(parents=True, exist_ok=True)
            args.output_json.write_text(rendered, encoding="utf-8")
        if args.output_csv is not None:
            write_normalized_csv(args.output_csv, normalized_rows)
        sys.stdout.write(rendered)
        return 0
    except (OSError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
