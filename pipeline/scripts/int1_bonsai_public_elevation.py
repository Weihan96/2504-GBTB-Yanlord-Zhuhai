#!/usr/bin/env python3
"""Create native Bonsai unfolded and long-section views across the public zone."""

from __future__ import annotations

import csv
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REGISTER = PROJECT_ROOT / "pipeline/decisions/int1-public-elevation-register.csv"
BASE_SCRIPT = PROJECT_ROOT / "pipeline/scripts/int1_bonsai_elevation.py"


def load_base():
    spec = importlib.util.spec_from_file_location("int1_bonsai_elevation_public_base", BASE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    return module


def run():
    base = load_base()
    rows = list(csv.DictReader(REGISTER.open(encoding="utf-8-sig")))
    reports = [base.build_view(PROJECT_ROOT, row) for row in rows]
    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "Bonsai 0.8.4 native Drawing / public true projection",
        "drawing_count": len(reports),
        "drawing_names": [report["drawing"]["name"] for report in reports],
        "complexity_exclusion_count": sum(
            report["complexity_exclusion_count"] for report in reports
        ),
        "lightweight_elevation_global_ids": sorted(
            {
                global_id
                for report in reports
                for global_id in report["lightweight_elevation_global_ids"]
            }
        ),
        "pass": len(reports) == 8
        and not any(report["complexity_exclusion_count"] for report in reports),
    }
    output = PROJECT_ROOT / base.REPORT_DIR / "native-public-elevation-summary.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))
    if not summary["pass"]:
        raise RuntimeError("public native elevation batch failed")
    return reports


if __name__ == "__main__":
    run()
