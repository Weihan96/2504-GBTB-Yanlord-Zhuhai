#!/usr/bin/env python3
"""Refresh Bonsai underlays and rebuild all or selected candidate drawings."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
BUILD_REGISTER = ROOT / "pipeline/decisions/drawing-build-register.json"
UNDERLAY_SCRIPT = ROOT / "pipeline/scripts/bonsai_plan_sources.py"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_register(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1 or not isinstance(data.get("groups"), list):
        raise RuntimeError("invalid drawing build register schema")
    ownership: dict[str, int] = {}
    for index, group in enumerate(data["groups"]):
        if not group.get("sheet_numbers") or not isinstance(group.get("commands"), list):
            raise RuntimeError(f"invalid drawing build group at index {index}")
        for sheet in group["sheet_numbers"]:
            if sheet in ownership:
                raise RuntimeError(f"drawing sheet has multiple build owners: {sheet}")
            ownership[sheet] = index
        for command in group["commands"]:
            if not command or command[0] not in {"bun", "python3"}:
                raise RuntimeError(f"unsupported drawing build command: {command}")
    data["ownership"] = ownership
    return data


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sheet", action="append", default=[])
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--no-refresh-underlays", action="store_true")
    parser.add_argument("--register", type=Path, default=BUILD_REGISTER)
    args = parser.parse_args()
    register = read_register(args.register.resolve())
    requested = set(args.sheet)
    known = set(register["ownership"])
    if unknown := requested - known:
        raise RuntimeError(f"unknown drawing sheet numbers: {sorted(unknown)}")
    selected_indexes = (
        {register["ownership"][sheet] for sheet in requested}
        if requested
        else set(range(len(register["groups"])))
    )
    selected_groups = [
        group for index, group in enumerate(register["groups"]) if index in selected_indexes
    ]
    underlays = sorted(
        {name for group in selected_groups for name in group.get("underlays", [])}
    )
    commands: list[list[str]] = []
    seen: set[tuple[str, ...]] = set()
    for group in selected_groups:
        for command in group["commands"]:
            key = tuple(command)
            if key not in seen:
                seen.add(key)
                commands.append(command)

    plan = {
        "mode": "candidate_drawing_pipeline",
        "formal_ifc": str(IFC),
        "formal_ifc_sha256": sha256(IFC),
        "requested_sheets": sorted(requested) if requested else "all",
        "selected_sheets": sorted(
            {sheet for group in selected_groups for sheet in group["sheet_numbers"]}
        ),
        "underlays": underlays,
        "commands": commands,
        "refresh_underlays": not args.no_refresh_underlays,
        "dry_run": args.dry_run,
    }
    if args.dry_run:
        print(json.dumps(plan, ensure_ascii=False, indent=2))
        return 0

    formal_hash = plan["formal_ifc_sha256"]
    if not args.no_refresh_underlays and underlays:
        command = ["python3", str(UNDERLAY_SCRIPT), "refresh"]
        for name in underlays:
            command.extend(["--drawing", name])
        subprocess.run(command, cwd=ROOT, check=True)
    for command in commands:
        subprocess.run(command, cwd=ROOT, check=True)
    if sha256(IFC) != formal_hash:
        raise RuntimeError("formal IFC changed during candidate drawing build")
    print(
        json.dumps(
            {**plan, "formal_ifc_unchanged": True, "pass": True},
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
