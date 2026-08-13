#!/usr/bin/env python3
"""Compile a read-only PLUM interface review matrix from canonical SSOT tables.

The output is an audit projection only. It does not infer connector coordinates,
route direction, pipe topology, or IFC relationships, and it never writes IFC.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any


FORMAL_IFC = "2504 GBTB Yanlord Zhuhai.ifc"
EQUIPMENT_REGISTER = Path("pipeline/decisions/equipment-register.csv")
REQUIREMENTS_REGISTER = Path("pipeline/decisions/equipment-installation-requirements.csv")
EVIDENCE_REGISTER = Path("pipeline/decisions/source-evidence-register.csv")
EXCLUDED_RELATIVE = Path("output/pdf/装修宝典2.0")
HASH_PATTERN = re.compile(r"^[0-9a-fA-F]{64}$")

COMPONENT_RELATIONSHIPS = (
    {
        "relationship_id": "PLUM-COMPONENT-CANDIDATE-001",
        "equipment_ids": ["SAN-001", "DRAIN-GEB-002", "DRAIN-GEB-003"],
        "relationship_kind": "candidate_component_compatibility_chain",
        "basis_requirement_ids": ["REQ-0414", "REQ-0450", "REQ-0421"],
        "source_ids": ["GEB-OFFICIAL-154446-001", "GEB-OFFICIAL-154150-001", "GEB-OFFICIAL-154298-001"],
        "review_boundary": "CleanLine50, installation set, and PVC adapter are grouped only as a review candidate; final assembly, orientation, and location remain unresolved.",
    },
    {
        "relationship_id": "PLUM-COMPONENT-CANDIDATE-002",
        "equipment_ids": ["SAN-014", "DRAIN-GEB-001"],
        "relationship_kind": "candidate_component_compatibility_pair",
        "basis_requirement_ids": ["REQ-0428", "REQ-0429", "REQ-0420", "REQ-0455"],
        "source_ids": ["GEB-OFFICIAL-224212-001", "GEB-OFFICIAL-152464-001"],
        "review_boundary": "Duofix and the PVC adapter are paired only by documented nominal interfaces; final assembly, orientation, and project occurrence assignment remain unresolved.",
    },
    {
        "relationship_id": "PLUM-COMPONENT-CANDIDATE-003",
        "equipment_ids": ["DRAIN-GEB-004"],
        "relationship_kind": "candidate_standalone_washing_machine_drain",
        "basis_requirement_ids": ["REQ-0417", "REQ-0418", "REQ-0419"],
        "source_ids": ["GEB-OFFICIAL-388013-001"],
        "review_boundary": "The product is registered only as a standalone washing-machine drain candidate; hose detail, level, slope, waterproofing, and IFC occurrence assignment remain unresolved.",
    },
)


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument("--ifc", type=Path, default=Path(FORMAL_IFC))
    parser.add_argument("--expected-ifc-sha256")
    parser.add_argument("--output-json", type=Path, default=Path("build/plum/plum-interface-review-candidate.json"))
    parser.add_argument("--output-csv", type=Path, default=Path("build/plum/plum-interface-review-candidate.csv"))
    parser.add_argument("--output-svg", type=Path, default=Path("build/plum/PLUM-interface-review-candidate.svg"))
    parser.add_argument("--output-png", type=Path, default=Path("build/plum/PLUM-interface-review-candidate.png"))
    parser.add_argument("--chrome", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def split_ids(value: str) -> list[str]:
    return [item.strip() for item in re.split(r"[;；]", value or "") if item.strip()]


def resolve(root: Path, path: Path) -> Path:
    return path if path.is_absolute() else root / path


def assert_not_excluded(root: Path, path: Path) -> None:
    excluded = (root / EXCLUDED_RELATIVE).resolve()
    candidate = path.resolve()
    if candidate == excluded or excluded in candidate.parents:
        raise RuntimeError(f"refusing to access excluded path: {candidate}")


def ensure_output(root: Path, path: Path, suffix: str, ifc_path: Path) -> None:
    if path.suffix.lower() != suffix:
        raise RuntimeError(f"output must use {suffix}: {path}")
    if path.resolve() == ifc_path.resolve():
        raise RuntimeError("derived output must not overwrite the formal IFC")
    assert_not_excluded(root, path)


def find_chrome(explicit: Path | None) -> Path:
    candidates = [
        explicit,
        Path(os.environ["CHROME_BIN"]) if os.environ.get("CHROME_BIN") else None,
        Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"),
        Path("/Applications/Chromium.app/Contents/MacOS/Chromium"),
    ]
    for name in ("google-chrome", "chromium"):
        binary = shutil.which(name)
        if binary:
            candidates.append(Path(binary))
    for candidate in candidates:
        if candidate and candidate.is_file():
            return candidate
    raise RuntimeError("Chrome/Chromium not found for mobile PNG proof")


def requirement_classification(row: dict[str, str]) -> str:
    if row["status"] == "pending" and row["blocks_release"] == "yes":
        return "pending_blocking"
    if row["status"] == "pending":
        return "pending_nonblocking"
    if row["status"] == "confirmed" and row["value_origin"] == "official_exact_model":
        return "exact_model_confirmed"
    if row["status"] == "candidate" and row["value_origin"] in {"official_model_family", "project_candidate"}:
        return "family_or_project_candidate"
    if row["status"] == "candidate":
        return "other_candidate"
    return "other_confirmed_or_observed"


def ifc_mapping_status(equipment: dict[str, str]) -> tuple[str, list[str]]:
    global_ids = split_ids(equipment["ifc_global_ids"])
    if not global_ids:
        return "no_formal_ifc_occurrence", global_ids
    if len(global_ids) == 1:
        return "linked_single_occurrence", global_ids
    return "multiple_occurrences_review_required", global_ids


def evidence_validation(
    root: Path,
    referenced_source_ids: set[str],
    evidence_by_id: dict[str, dict[str, str]],
) -> dict[str, Any]:
    missing_source_ids = sorted(referenced_source_ids - set(evidence_by_id))
    local_checks: list[dict[str, Any]] = []
    for source_id in sorted(referenced_source_ids & set(evidence_by_id)):
        source = evidence_by_id[source_id]
        local_path = source["local_path"].strip()
        registered_hash = source["sha256"].strip().lower()
        if not local_path or not HASH_PATTERN.fullmatch(registered_hash):
            continue
        path = resolve(root, Path(local_path))
        assert_not_excluded(root, path)
        exists = path.is_file()
        actual_hash = sha256(path) if exists else None
        local_checks.append({
            "source_id": source_id,
            "path": str(path),
            "exists": exists,
            "registered_sha256": registered_hash,
            "actual_sha256": actual_hash,
            "hash_matches": exists and actual_hash == registered_hash,
        })
    mismatches = [row for row in local_checks if not row["hash_matches"]]
    return {
        "referenced_source_id_count": len(referenced_source_ids),
        "all_referenced_source_ids_exist": not missing_source_ids,
        "missing_source_ids": missing_source_ids,
        "local_evidence_hash_check_count": len(local_checks),
        "local_evidence_hashes_match": not mismatches,
        "local_evidence_hash_mismatches": mismatches,
        "local_evidence_checks": local_checks,
    }


def safe_text(value: str, limit: int) -> str:
    value = " ".join((value or "").split())
    return value if len(value) <= limit else value[: max(0, limit - 1)] + "…"


def render_svg(equipment_records: list[dict[str, Any]], summary: dict[str, Any], ifc_hash: str) -> str:
    width = 1080
    row_height = 52
    header_height = 300
    height = header_height + len(equipment_records) * row_height + 110
    mapping_labels = {
        "linked_single_occurrence": "IFC 1:1",
        "no_formal_ifc_occurrence": "无 occurrence",
        "multiple_occurrences_review_required": "多实例/歧义",
    }
    rows: list[str] = []
    for index, record in enumerate(equipment_records):
        y = header_height + index * row_height
        counts = record["requirement_counts"]
        fill = "#ffffff" if index % 2 == 0 else "#f8fafc"
        rows.append(
            f'<g data-equipment-id="{html.escape(record["equipment_id"])}">'
            f'<rect x="30" y="{y}" width="1020" height="{row_height}" fill="{fill}"/>'
            f'<text class="id" x="44" y="{y + 21}">{html.escape(record["equipment_id"])}</text>'
            f'<text class="name" x="44" y="{y + 40}">{html.escape(safe_text(record["item_name"], 29))}</text>'
            f'<text class="model" x="350" y="{y + 31}">{html.escape(safe_text(record["model"] or "—", 24))}</text>'
            f'<text class="mapping {record["ifc_mapping_status"]}" x="585" y="{y + 31}">{mapping_labels[record["ifc_mapping_status"]]}</text>'
            f'<text class="count exact" x="748" y="{y + 31}">{counts.get("exact_model_confirmed", 0)}</text>'
            f'<text class="count candidate" x="830" y="{y + 31}">{counts.get("family_or_project_candidate", 0)}</text>'
            f'<text class="count pending" x="914" y="{y + 31}">{record["pending_requirement_count"]}</text>'
            f'<text class="count blocking" x="1005" y="{y + 31}">{record["blocking_requirement_count"]}</text>'
            f'<line x1="30" y1="{y + row_height}" x2="1050" y2="{y + row_height}" class="rule"/>'
            "</g>"
        )
    classification = summary["requirement_classification_counts"]
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>
.bg{{fill:#eef3f8}}.card{{fill:#fff;stroke:#cbd5e1;stroke-width:1.5}}.title{{font:700 34px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif;fill:#0f2940}}
.sub{{font:16px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif;fill:#526777}}.metric{{font:700 22px -apple-system,"PingFang SC",sans-serif;fill:#123b5d}}
.metric-label{{font:13px -apple-system,"PingFang SC",sans-serif;fill:#526777}}.warning{{font:700 14px -apple-system,"PingFang SC",sans-serif;fill:#a61e1e}}
.head{{font:700 13px -apple-system,"PingFang SC",sans-serif;fill:#344b5f}}.id{{font:700 12px ui-monospace,SFMono-Regular,monospace;fill:#23445e}}
.name{{font:13px -apple-system,"PingFang SC",sans-serif;fill:#23384a}}.model{{font:12px ui-monospace,SFMono-Regular,monospace;fill:#40586b}}
.mapping{{font:700 12px -apple-system,"PingFang SC",sans-serif}}.linked_single_occurrence{{fill:#087f5b}}.no_formal_ifc_occurrence{{fill:#b45309}}.multiple_occurrences_review_required{{fill:#7c3aed}}
.count{{font:700 14px ui-monospace,SFMono-Regular,monospace;text-anchor:middle}}.exact{{fill:#087f5b}}.candidate{{fill:#9c6a00}}.pending{{fill:#c2410c}}.blocking{{fill:#c92a2a}}
.rule{{stroke:#e2e8f0;stroke-width:1}}.footer{{font:12px ui-monospace,SFMono-Regular,monospace;fill:#64748b}}
</style>
<rect class="bg" width="{width}" height="{height}"/>
<rect class="card" x="30" y="24" width="1020" height="{height - 48}" rx="18"/>
<text class="title" x="54" y="72">PLUM 接口审核候选矩阵</text>
<text class="sub" x="54" y="103">纯派生审核视图｜不生成接口坐标、方向、管线拓扑或 IFC 写入</text>
<text class="metric" x="70" y="151">{summary["equipment_count"]}</text><text class="metric-label" x="54" y="176">设备</text>
<text class="metric" x="245" y="151">{summary["requirement_count"]}</text><text class="metric-label" x="215" y="176">PLUM 要求</text>
<text class="metric" x="440" y="151">{summary["blocking_requirement_count"]}</text><text class="metric-label" x="395" y="176">发布阻塞</text>
<text class="metric" x="630" y="151">{classification.get("exact_model_confirmed", 0)}</text><text class="metric-label" x="572" y="176">精确型号已确认</text>
<text class="metric" x="865" y="151">{classification.get("family_or_project_candidate", 0)}</text><text class="metric-label" x="785" y="176">产品族/项目候选</text>
<text class="warning" x="54" y="215">automatic_ifc_write_allowed=false ｜ construction_release_ready=false</text>
<text class="sub" x="54" y="241">三类组件关系仅作兼容性审核候选；所有方向、坐标和 occurrence 归属继续开放。</text>
<text class="head" x="44" y="282">设备 ID / 名称</text><text class="head" x="350" y="282">型号</text><text class="head" x="585" y="282">IFC 映射</text>
<text class="head" x="748" y="282" text-anchor="middle">Exact</text><text class="head" x="830" y="282" text-anchor="middle">候选</text>
<text class="head" x="914" y="282" text-anchor="middle">Pending</text><text class="head" x="1005" y="282" text-anchor="middle">阻塞</text>
{''.join(rows)}
<text class="footer" x="54" y="{height - 62}">IFC SHA-256 {ifc_hash}</text>
<text class="footer" x="54" y="{height - 38}">SSOT → deterministic review projection · no PDF · no IFC mutation</text>
</svg>
'''


def write_csv_output(path: Path, equipment_records: list[dict[str, Any]]) -> None:
    fields = [
        "equipment_id", "item_name", "manufacturer", "model", "procurement_status", "decision_status",
        "use_location", "ifc_mapping_status", "ifc_global_ids", "requirement_id", "discipline",
        "parameter_key", "value", "unit", "value_origin", "status", "classification",
        "blocks_release", "source_id", "source_locator", "source_exists", "local_evidence_hash_status",
        "automatic_ifc_write_allowed", "construction_release_ready", "notes",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for equipment in equipment_records:
            for requirement in equipment["requirements"]:
                writer.writerow({
                    "equipment_id": equipment["equipment_id"],
                    "item_name": equipment["item_name"],
                    "manufacturer": equipment["manufacturer"],
                    "model": equipment["model"],
                    "procurement_status": equipment["procurement_status"],
                    "decision_status": equipment["decision_status"],
                    "use_location": equipment["use_location"],
                    "ifc_mapping_status": equipment["ifc_mapping_status"],
                    "ifc_global_ids": ";".join(equipment["ifc_global_ids"]),
                    "requirement_id": requirement["requirement_id"],
                    "discipline": requirement["discipline"],
                    "parameter_key": requirement["parameter_key"],
                    "value": requirement["value"],
                    "unit": requirement["unit"],
                    "value_origin": requirement["value_origin"],
                    "status": requirement["status"],
                    "classification": requirement["classification"],
                    "blocks_release": "yes" if requirement["blocks_release"] else "no",
                    "source_id": requirement["source_id"],
                    "source_locator": requirement["source_locator"],
                    "source_exists": requirement["source_exists"],
                    "local_evidence_hash_status": requirement["local_evidence_hash_status"],
                    "automatic_ifc_write_allowed": "false",
                    "construction_release_ready": "false",
                    "notes": requirement["notes"],
                })


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    ifc_path = resolve(root, args.ifc)
    assert_not_excluded(root, ifc_path)
    output_json = resolve(root, args.output_json)
    output_csv = resolve(root, args.output_csv)
    output_svg = resolve(root, args.output_svg)
    output_png = resolve(root, args.output_png)
    for path, suffix in (
        (output_json, ".json"), (output_csv, ".csv"),
        (output_svg, ".svg"), (output_png, ".png"),
    ):
        ensure_output(root, path, suffix, ifc_path)

    starting_ifc_hash = sha256(ifc_path)
    if args.expected_ifc_sha256 and starting_ifc_hash != args.expected_ifc_sha256:
        raise RuntimeError(
            f"formal IFC hash changed: expected {args.expected_ifc_sha256}, got {starting_ifc_hash}"
        )

    equipment_path = root / EQUIPMENT_REGISTER
    requirements_path = root / REQUIREMENTS_REGISTER
    evidence_path = root / EVIDENCE_REGISTER
    equipment_rows = read_csv(equipment_path)
    requirement_rows = [row for row in read_csv(requirements_path) if "PLUM" in row["discipline"]]
    evidence_rows = read_csv(evidence_path)

    equipment_by_id = {row["equipment_id"]: row for row in equipment_rows}
    if len(equipment_by_id) != len(equipment_rows):
        raise RuntimeError("equipment register contains duplicate equipment IDs")
    evidence_by_id = {row["source_id"]: row for row in evidence_rows}
    if len(evidence_by_id) != len(evidence_rows):
        raise RuntimeError("source evidence register contains duplicate source IDs")
    if len({row["requirement_id"] for row in requirement_rows}) != len(requirement_rows):
        raise RuntimeError("PLUM requirement projection contains duplicate requirement IDs")

    plum_equipment_ids = sorted({row["equipment_id"] for row in requirement_rows})
    missing_equipment = sorted(set(plum_equipment_ids) - set(equipment_by_id))
    if missing_equipment:
        raise RuntimeError(f"PLUM requirements reference missing equipment: {missing_equipment}")

    referenced_source_ids = {
        row["source_id"] for row in requirement_rows if row["source_id"]
    }
    for equipment_id in plum_equipment_ids:
        referenced_source_ids.update(split_ids(equipment_by_id[equipment_id]["source_ids"]))
    for relation in COMPONENT_RELATIONSHIPS:
        referenced_source_ids.update(relation["source_ids"])
    source_validation = evidence_validation(root, referenced_source_ids, evidence_by_id)
    if not source_validation["all_referenced_source_ids_exist"]:
        raise RuntimeError(f"missing source evidence IDs: {source_validation['missing_source_ids']}")
    if not source_validation["local_evidence_hashes_match"]:
        raise RuntimeError(
            f"local source evidence hash mismatch: {source_validation['local_evidence_hash_mismatches']}"
        )
    local_hash_status = {
        row["source_id"]: "match" if row["hash_matches"] else "mismatch"
        for row in source_validation["local_evidence_checks"]
    }

    requirements_by_equipment: dict[str, list[dict[str, str]]] = {}
    for row in requirement_rows:
        requirements_by_equipment.setdefault(row["equipment_id"], []).append(row)

    equipment_records: list[dict[str, Any]] = []
    classification_counts: Counter[str] = Counter()
    for equipment_id in plum_equipment_ids:
        equipment = equipment_by_id[equipment_id]
        mapping_status, global_ids = ifc_mapping_status(equipment)
        requirements: list[dict[str, Any]] = []
        for row in requirements_by_equipment[equipment_id]:
            classification = requirement_classification(row)
            classification_counts[classification] += 1
            requirements.append({
                "requirement_id": row["requirement_id"],
                "discipline": row["discipline"],
                "parameter_key": row["parameter_key"],
                "value": row["value_number"] or row["value_text"],
                "unit": row["unit"],
                "value_origin": row["value_origin"],
                "status": row["status"],
                "classification": classification,
                "blocks_release": row["blocks_release"] == "yes",
                "source_id": row["source_id"],
                "source_locator": row["source_locator"],
                "source_exists": not row["source_id"] or row["source_id"] in evidence_by_id,
                "local_evidence_hash_status": local_hash_status.get(row["source_id"], "not_applicable"),
                "notes": row["notes"],
            })
        counts = Counter(row["classification"] for row in requirements)
        equipment_records.append({
            "equipment_id": equipment_id,
            "domain": equipment["domain"],
            "category": equipment["category"],
            "item_name": equipment["item_name"],
            "manufacturer": equipment["manufacturer"],
            "model": equipment["model"],
            "procurement_status": equipment["procurement_status"],
            "decision_status": equipment["decision_status"],
            "use_location": equipment["use_location_confirmed"] or equipment["use_location_candidate"],
            "ifc_mapping_status": mapping_status,
            "ifc_global_ids": global_ids,
            "requirement_count": len(requirements),
            "blocking_requirement_count": sum(row["blocks_release"] for row in requirements),
            "pending_requirement_count": sum(row["status"] == "pending" for row in requirements),
            "requirement_counts": dict(sorted(counts.items())),
            "requirements": requirements,
            "automatic_ifc_write_allowed": False,
        })

    requirement_ids = {row["requirement_id"] for row in requirement_rows}
    for relation in COMPONENT_RELATIONSHIPS:
        missing_ids = sorted(set(relation["equipment_ids"]) - set(equipment_by_id))
        missing_requirements = sorted(set(relation["basis_requirement_ids"]) - requirement_ids)
        if missing_ids or missing_requirements:
            raise RuntimeError(
                f"component candidate {relation['relationship_id']} is stale: "
                f"missing_equipment={missing_ids}, missing_requirements={missing_requirements}"
            )

    mapping_counts = Counter(row["ifc_mapping_status"] for row in equipment_records)
    summary = {
        "equipment_count": len(equipment_records),
        "requirement_count": len(requirement_rows),
        "blocking_requirement_count": sum(row["blocks_release"] == "yes" for row in requirement_rows),
        "ifc_mapping_counts": dict(sorted(mapping_counts.items())),
        "requirement_classification_counts": dict(sorted(classification_counts.items())),
        "component_relationship_candidate_count": len(COMPONENT_RELATIONSHIPS),
    }
    report: dict[str, Any] = {
        "mode": "read_only_plum_interface_review_candidate",
        "source_ifc_sha256": starting_ifc_hash,
        "source": {
            "ifc": {"path": str(ifc_path), "sha256": starting_ifc_hash},
            "equipment_register": {"path": str(equipment_path), "sha256": sha256(equipment_path)},
            "installation_requirements": {"path": str(requirements_path), "sha256": sha256(requirements_path)},
            "source_evidence": {"path": str(evidence_path), "sha256": sha256(evidence_path)},
        },
        "summary": summary,
        "source_validation": source_validation,
        "equipment": equipment_records,
        "component_relationship_candidates": [
            {
                **relation,
                "candidate_only": True,
                "automatic_ifc_write_allowed": False,
            }
            for relation in COMPONENT_RELATIONSHIPS
        ],
        "gates": {
            "all_plum_equipment_projected": len(equipment_records) == len(plum_equipment_ids),
            "all_plum_requirements_projected": sum(row["requirement_count"] for row in equipment_records) == len(requirement_rows),
            "all_referenced_source_ids_exist": source_validation["all_referenced_source_ids_exist"],
            "local_evidence_hashes_match": source_validation["local_evidence_hashes_match"],
            "component_relationships_are_candidate_only": True,
            "contains_connector_coordinates": False,
            "contains_route_directions": False,
            "automatic_ifc_write_allowed": False,
            "construction_release_ready": False,
        },
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_csv_output(output_csv, equipment_records)
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    svg = render_svg(equipment_records, summary, starting_ifc_hash)
    output_svg.write_text(svg, encoding="utf-8")
    subprocess.run(
        [
            str(find_chrome(args.chrome)), "--headless=new", "--disable-gpu", "--hide-scrollbars",
            "--force-device-scale-factor=1", f"--screenshot={output_png}",
            f"--window-size=1080,{300 + len(equipment_records) * 52 + 110}", output_svg.resolve().as_uri(),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if not output_png.is_file() or output_png.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError("mobile proof output is not a valid PNG")

    ending_ifc_hash = sha256(ifc_path)
    if ending_ifc_hash != starting_ifc_hash:
        raise RuntimeError("formal IFC changed during PLUM interface review generation")
    report["source"]["ifc"]["ending_sha256"] = ending_ifc_hash
    report["gates"]["ifc_unchanged_during_generation"] = True
    report["outputs"] = {
        "json": {"path": str(output_json)},
        "csv": {"path": str(output_csv), "sha256": sha256(output_csv)},
        "svg": {"path": str(output_svg), "sha256": sha256(output_svg)},
        "png": {"path": str(output_png), "sha256": sha256(output_png)},
    }
    output_json.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": summary, "gates": report["gates"], "outputs": report["outputs"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
