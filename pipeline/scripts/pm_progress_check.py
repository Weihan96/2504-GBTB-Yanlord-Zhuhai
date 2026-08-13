#!/usr/bin/env python3
"""Verify that the PM summary and Gantt match the authoritative task ledger."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


DEFAULT_PM = Path("drawings/滨海湾装修施工图深化工作管理.md")
TASK_RE = re.compile(r"^- \[([ x])\] \*\*([A-Z][A-Z0-9]+)｜", re.MULTILINE)
SUMMARY_RE = re.compile(r"PM 进度（派生）.*?(\d+) 项中 (\d+) 项完成（(\d+)%）")
EQUIPMENT_SSOT_RE = re.compile(
    r"当前为设备主表 `(\d+)` 条、安装条件 `(\d+)` 条、证据 `(\d+)` 条"
)
DEFAULT_EQUIPMENT = Path("pipeline/decisions/equipment-register.csv")
DEFAULT_REQUIREMENTS = Path("pipeline/decisions/equipment-installation-requirements.csv")
DEFAULT_EVIDENCE = Path("pipeline/decisions/source-evidence-register.csv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pm", type=Path, default=DEFAULT_PM)
    parser.add_argument("--equipment-register", type=Path, default=DEFAULT_EQUIPMENT)
    parser.add_argument("--installation-requirements", type=Path, default=DEFAULT_REQUIREMENTS)
    parser.add_argument("--source-evidence", type=Path, default=DEFAULT_EVIDENCE)
    return parser.parse_args()


def section(text: str, start: str, end: str) -> str:
    start_index = text.index(start)
    end_index = text.index(end, start_index)
    return text[start_index:end_index]


def gantt_states(text: str, task_ids: set[str]) -> dict[str, bool]:
    marker = "```mermaid\ngantt"
    start_index = text.index(marker)
    end_index = text.index("```", start_index + len(marker))
    gantt = text[start_index:end_index]
    states: dict[str, bool] = {}
    for line in gantt.splitlines():
        if ":" not in line:
            continue
        tokens = [token.strip() for token in line.split(":", 1)[1].split(",")]
        matches = [token for token in tokens if token in task_ids]
        if len(matches) == 1:
            task_id = matches[0]
            states[task_id] = "done" in tokens[: tokens.index(task_id)]
    return states


def csv_record_count(path: Path) -> int:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return sum(1 for _ in csv.DictReader(stream))


def check(
    path: Path,
    equipment_register: Path,
    installation_requirements: Path,
    source_evidence: Path,
) -> dict[str, object]:
    text = path.read_text(encoding="utf-8")
    ledger = section(text, "## 5. PM 任务台账", "### 任务台账操作规则")
    tasks = [(task_id, mark == "x") for mark, task_id in TASK_RE.findall(ledger)]
    task_ids = [task_id for task_id, _ in tasks]
    total = len(tasks)
    done = sum(is_done for _, is_done in tasks)
    percent = round(done * 100 / total) if total else 0

    errors: list[str] = []
    duplicates = sorted({task_id for task_id in task_ids if task_ids.count(task_id) > 1})
    if duplicates:
        errors.append("duplicate task ids: " + ", ".join(duplicates))

    summary_match = SUMMARY_RE.search(text)
    summary = None
    if summary_match:
        summary = {
            "total": int(summary_match.group(1)),
            "done": int(summary_match.group(2)),
            "percent": int(summary_match.group(3)),
        }
        expected = {"total": total, "done": done, "percent": percent}
        if summary != expected:
            errors.append(f"PM summary {summary} does not match ledger {expected}")
    else:
        errors.append("PM progress summary row is missing")

    states = gantt_states(text, set(task_ids))
    missing = sorted(set(task_ids) - set(states))
    extra = sorted(set(states) - set(task_ids))
    if missing:
        errors.append("ledger tasks missing from Gantt: " + ", ".join(missing))
    if extra:
        errors.append("Gantt tasks missing from ledger: " + ", ".join(extra))
    mismatched = sorted(
        task_id for task_id, is_done in tasks
        if task_id in states and states[task_id] != is_done
    )
    if mismatched:
        errors.append("Gantt done state differs from ledger: " + ", ".join(mismatched))

    ssot_paths = {
        "equipment": equipment_register,
        "requirements": installation_requirements,
        "evidence": source_evidence,
    }
    missing_ssot = [name for name, source in ssot_paths.items() if not source.is_file()]
    ssot_actual = None
    ssot_declared = None
    if missing_ssot:
        errors.append("equipment SSOT files are missing: " + ", ".join(missing_ssot))
    else:
        ssot_actual = {
            name: csv_record_count(source) for name, source in ssot_paths.items()
        }
        ssot_match = EQUIPMENT_SSOT_RE.search(text)
        if not ssot_match:
            errors.append("current equipment SSOT count sentence is missing")
        else:
            ssot_declared = {
                "equipment": int(ssot_match.group(1)),
                "requirements": int(ssot_match.group(2)),
                "evidence": int(ssot_match.group(3)),
            }
            if ssot_declared != ssot_actual:
                errors.append(
                    f"PM equipment SSOT counts {ssot_declared} do not match canonical tables {ssot_actual}"
                )

    return {
        "status": not errors,
        "read_only": True,
        "pm_path": str(path),
        "summary": {"total": total, "done": done, "open": total - done, "percent": percent},
        "equipment_ssot": {"declared": ssot_declared, "actual": ssot_actual},
        "gantt_task_count": len(states),
        "errors": errors,
    }


def main() -> int:
    args = parse_args()
    result = check(
        args.pm,
        args.equipment_register,
        args.installation_requirements,
        args.source_evidence,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
