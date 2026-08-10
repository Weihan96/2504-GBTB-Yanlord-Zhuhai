#!/usr/bin/env python3
"""Run the project IDS with IfcTester and write a compact, reviewable report."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_blender_ifctester_path() -> str | None:
    """Expose Bonsai's pure-Python IfcTester package to the project Python."""
    candidates = sorted(
        (Path.home() / "Library/Application Support/Blender").glob(
            "*/extensions/.local/lib/python*/site-packages"
        ),
        reverse=True,
    )
    for candidate in candidates:
        if (candidate / "ifctester/__init__.py").is_file():
            candidate_text = str(candidate)
            if candidate_text not in sys.path:
                sys.path.append(candidate_text)
            return candidate_text
    return None


def compact_report(raw: dict[str, Any]) -> dict[str, Any]:
    specifications = []
    for specification in raw.get("specifications", []):
        requirements = []
        for requirement in specification.get("requirements", []):
            failed_entities = requirement.get("failed_entities", [])
            requirements.append(
                {
                    "description": requirement.get("description"),
                    "status": bool(requirement.get("status")),
                    "total_applicable": requirement.get("total_applicable", 0),
                    "total_pass": requirement.get("total_pass", 0),
                    "total_fail": requirement.get("total_fail", 0),
                    "failed_global_ids": sorted(
                        {
                            entity.get("global_id")
                            for entity in failed_entities
                            if entity.get("global_id")
                        }
                    ),
                    "failure_reasons": sorted(
                        {
                            entity.get("reason")
                            for entity in failed_entities
                            if entity.get("reason")
                        }
                    ),
                }
            )
        specifications.append(
            {
                "name": specification.get("name"),
                "status": bool(specification.get("status")),
                "total_applicable": specification.get("total_applicable", 0),
                "total_applicable_pass": specification.get("total_applicable_pass", 0),
                "total_applicable_fail": specification.get("total_applicable_fail", 0),
                "total_checks": specification.get("total_checks", 0),
                "total_checks_pass": specification.get("total_checks_pass", 0),
                "total_checks_fail": specification.get("total_checks_fail", 0),
                "requirements": requirements,
            }
        )
    return {
        "status": bool(raw.get("status")),
        "total_specifications": raw.get("total_specifications", len(specifications)),
        "total_specifications_pass": raw.get("total_specifications_pass", 0),
        "total_specifications_fail": raw.get("total_specifications_fail", 0),
        "total_checks": raw.get("total_checks", 0),
        "total_checks_pass": raw.get("total_checks_pass", 0),
        "total_checks_fail": raw.get("total_checks_fail", 0),
        "specifications": specifications,
    }


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# P0 IDS 验证报告",
        "",
        f"- 总体状态：{'通过' if report['status'] else '未通过'}",
        f"- 规格：{report['total_specifications_pass']}/{report['total_specifications']} 通过",
        f"- 检查：{report['total_checks_pass']}/{report['total_checks']} 通过",
        "",
        "| 规格 | 适用对象通过 | 检查通过 | 状态 |",
        "| --- | ---: | ---: | --- |",
    ]
    for specification in report["specifications"]:
        lines.append(
            f"| {specification['name']} | "
            f"{specification['total_applicable_pass']}/{specification['total_applicable']} | "
            f"{specification['total_checks_pass']}/{specification['total_checks']} | "
            f"{'通过' if specification['status'] else '未通过'} |"
        )
    lines.extend(["", "## 未通过要求", ""])
    failures = 0
    for specification in report["specifications"]:
        for requirement in specification["requirements"]:
            if requirement["status"]:
                continue
            failures += 1
            ids = ", ".join(f"`{value}`" for value in requirement["failed_global_ids"])
            lines.extend(
                [
                    f"### {specification['name']}",
                    "",
                    f"- 要求：{requirement['description']}",
                    f"- 失败：{requirement['total_fail']}/{requirement['total_applicable']}",
                    f"- GlobalId：{ids or '无'}",
                    "",
                ]
            )
    if not failures:
        lines.append("无。")
    return "\n".join(lines).rstrip() + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--ids", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--expected-ifc-sha256")
    parser.add_argument("--fail-on-failure", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"IFC not found: {args.input}")
    if not args.ids.is_file():
        raise SystemExit(f"IDS not found: {args.ids}")
    ifc_hash = sha256(args.input)
    if args.expected_ifc_sha256 and args.expected_ifc_sha256 != ifc_hash:
        raise SystemExit(
            "formal IFC hash differs from caller-frozen hash: "
            f"expected {args.expected_ifc_sha256}, found {ifc_hash}"
        )

    add_blender_ifctester_path()
    try:
        import ifcopenshell
        import ifctester.ids
        import ifctester.reporter
    except ModuleNotFoundError as error:
        raise SystemExit(
            "IfcTester is unavailable. Install it in the project Python or Bonsai's Blender environment."
        ) from error

    ifc = ifcopenshell.open(str(args.input))
    ids = ifctester.ids.open(str(args.ids), validate=True)
    ids.validate(ifc, should_filter_version=True)
    reporter = ifctester.reporter.Json(ids)
    report = compact_report(reporter.report())
    report.update(
        {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "ifc": str(args.input),
            "ids": str(args.ids),
            "ifc_schema": ifc.schema,
            "source_ifc_sha256": ifc_hash,
            "source": {
                "ifc": str(args.input),
                "ifc_sha256": ifc_hash,
                "ids": str(args.ids),
                "ids_sha256": sha256(args.ids),
            },
        }
    )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(markdown_report(report))

    print(
        json.dumps(
            {
                "report": str(args.report),
                "markdown": str(args.markdown) if args.markdown else None,
                "status": report["status"],
                "specifications_pass": report["total_specifications_pass"],
                "specifications_total": report["total_specifications"],
                "checks_pass": report["total_checks_pass"],
                "checks_total": report["total_checks"],
            },
            ensure_ascii=False,
        )
    )
    return 1 if args.fail_on_failure and not report["status"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
