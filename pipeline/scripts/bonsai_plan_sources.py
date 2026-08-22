#!/usr/bin/env python3
"""Check or refresh registered Bonsai plan-source SVGs without writing IFC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
DEFAULT_REGISTER = ROOT / "pipeline/decisions/bonsai-plan-source-register.csv"
BLENDER = Path("/Applications/Blender.app/Contents/MacOS/Blender")
BLENDER_SCRIPT = ROOT / "pipeline/scripts/refresh_bonsai_plan_sources.py"
XLINK_HREF = "{http://www.w3.org/1999/xlink}href"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_register(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "drawing_name", "drawing_global_id", "target_view", "scale",
        "shading_style", "svg_path", "manifest_path", "consumers",
    }
    if not rows or set(rows[0]) != required:
        raise RuntimeError("invalid Bonsai plan source register schema")
    names = [row["drawing_name"] for row in rows]
    if len(names) != len(set(names)):
        raise RuntimeError("duplicate Bonsai plan source names")
    return rows


def local_hrefs(svg_path: Path) -> list[str]:
    root = ET.parse(svg_path).getroot()
    result: set[str] = set()
    for node in root.iter():
        href = node.attrib.get(XLINK_HREF) or node.attrib.get("href")
        if href and not href.startswith(("#", "data:", "http://", "https://")):
            result.add(href)
    return sorted(result)


def inspect_row(row: dict[str, str], ifc_hash: str) -> dict[str, Any]:
    svg = ROOT / row["svg_path"]
    manifest_path = ROOT / row["manifest_path"]
    reasons: list[str] = []
    if not svg.is_file():
        reasons.append("source SVG missing")
    if not manifest_path.is_file():
        reasons.append("source manifest missing")
    manifest: dict[str, Any] = {}
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            reasons.append("source manifest invalid JSON")
    expected_manifest = {
        "formal_ifc_sha256": ifc_hash,
        "drawing_global_id": row["drawing_global_id"],
        "drawing_name": row["drawing_name"],
        "target_view": row["target_view"],
        "scale": row["scale"],
        "shading_style": row["shading_style"],
    }
    for key, expected in expected_manifest.items():
        if manifest and manifest.get(key) != expected:
            reasons.append(f"manifest {key} mismatch")
    if manifest and manifest.get("formal_ifc_write_allowed") not in {None, False}:
        reasons.append("manifest unexpectedly allows IFC write")
    if svg.is_file() and manifest and manifest.get("source_svg_sha256") != sha256(svg):
        reasons.append("source SVG hash mismatch")

    hrefs: list[str] = []
    if svg.is_file():
        try:
            hrefs = local_hrefs(svg)
        except ET.ParseError:
            reasons.append("source SVG invalid XML")
    linked = {item.get("href"): item for item in manifest.get("linked_assets", [])}
    for href in hrefs:
        if Path(href).name != href:
            reasons.append(f"nonportable linked asset: {href}")
            continue
        asset = svg.with_name(href)
        if not asset.is_file():
            reasons.append(f"linked asset missing: {href}")
            continue
        if linked:
            record = linked.get(href)
            if not record or record.get("sha256") != sha256(asset):
                reasons.append(f"linked asset hash mismatch: {href}")
        elif href == f"{row['drawing_name']}-underlay.png":
            if manifest.get("underlay_png_sha256") != sha256(asset):
                reasons.append(f"legacy underlay hash mismatch: {href}")
        else:
            reasons.append(f"linked asset absent from manifest: {href}")
    return {
        "drawing_name": row["drawing_name"],
        "svg": row["svg_path"],
        "manifest": row["manifest_path"],
        "consumers": row["consumers"].split(";") if row["consumers"] else [],
        "current": not reasons,
        "reasons": reasons,
    }


def status(ifc_path: Path, register_path: Path) -> dict[str, Any]:
    ifc_hash = sha256(ifc_path)
    rows = read_register(register_path)
    drawings = [inspect_row(row, ifc_hash) for row in rows]
    return {
        "mode": "read_only_bonsai_plan_source_status",
        "formal_ifc": str(ifc_path),
        "formal_ifc_sha256": ifc_hash,
        "drawing_count": len(drawings),
        "current_count": sum(item["current"] for item in drawings),
        "stale_count": sum(not item["current"] for item in drawings),
        "drawings": drawings,
        "pass": all(item["current"] for item in drawings),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("status", "check", "refresh"))
    parser.add_argument("--ifc", type=Path, default=DEFAULT_IFC)
    parser.add_argument("--register", type=Path, default=DEFAULT_REGISTER)
    parser.add_argument("--blender", type=Path, default=BLENDER)
    parser.add_argument("--drawing", action="append", default=[])
    args = parser.parse_args()
    ifc_path = args.ifc.resolve()
    register_path = args.register.resolve()
    before = status(ifc_path, register_path)
    selected = set(args.drawing)
    known = {item["drawing_name"] for item in before["drawings"]}
    if unknown := selected - known:
        raise RuntimeError(f"unknown registered Drawing names: {sorted(unknown)}")

    if args.command in {"status", "check"}:
        print(json.dumps(before, ensure_ascii=False, indent=2))
        return 0 if args.command == "status" or before["pass"] else 1

    stale = {
        item["drawing_name"] for item in before["drawings"] if not item["current"]
    }
    targets = sorted((selected or known) & stale)
    if not targets:
        print(json.dumps({**before, "refresh": "not_required"}, ensure_ascii=False, indent=2))
        return 0
    if not args.blender.is_file():
        raise RuntimeError(f"Blender executable missing: {args.blender}")
    run_root = ROOT / "build/bonsai-underlays"
    run_root.mkdir(parents=True, exist_ok=True)
    formal_hash = before["formal_ifc_sha256"]
    for name in targets:
        print(f"Refreshing native Bonsai Drawing: {name}", flush=True)
        candidate_dir = Path(tempfile.mkdtemp(prefix="refresh-", dir=run_root))
        command = [
            str(args.blender), "--python-exit-code", "1",
            "--python", str(BLENDER_SCRIPT), "--",
            "--ifc", str(ifc_path), "--register", str(register_path),
            "--candidate-dir", str(candidate_dir), "--drawing", name,
        ]
        try:
            process = subprocess.run(
                command,
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            markers = [
                line for line in process.stdout.splitlines()
                if line.startswith("BONSAI_PLAN_REFRESH=")
            ]
            if process.returncode or len(markers) != 1:
                diagnostic = "\n".join(
                    (process.stdout + process.stderr).splitlines()[-100:]
                )
                raise RuntimeError(
                    f"Bonsai plan refresh failed for {name} with exit "
                    f"{process.returncode}:\n{diagnostic}"
                )
        finally:
            shutil.rmtree(candidate_dir, ignore_errors=True)
        if sha256(ifc_path) != formal_hash:
            raise RuntimeError(f"formal IFC hash changed while refreshing {name}")
    after = status(ifc_path, register_path)
    if after["formal_ifc_sha256"] != formal_hash:
        raise RuntimeError("formal IFC hash changed during Bonsai plan refresh")
    selected_status = {
        item["drawing_name"]: item for item in after["drawings"]
    }
    failed = [name for name in targets if not selected_status[name]["current"]]
    if failed:
        raise RuntimeError(f"refreshed Bonsai sources remain stale: {failed}")
    print(
        json.dumps(
            {
                **after,
                "refreshed_drawings": targets,
                "formal_ifc_unchanged": True,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
