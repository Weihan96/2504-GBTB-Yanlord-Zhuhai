#!/usr/bin/env python3
"""Render an evidence-bounded P-201 water/hot-water demand endpoint candidate."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import shutil
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPECTED_SERVICE_COUNT = 24
EXPECTED_NON_SERVICE_COUNT = 3
ROOM_ORDER = ("BATHM", "BATHG", "WC", "VVD")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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
    raise RuntimeError("Chrome/Chromium not found for P-201 PNG proof")


def ensure_safe_output(path: Path, suffix: str, ifc_path: Path) -> None:
    if path.resolve() == ifc_path.resolve():
        raise RuntimeError("P-201 output must not overwrite the formal IFC")
    if path.suffix.lower() != suffix:
        raise RuntimeError(f"P-201 output must use {suffix}: {path}")


def validate_sources(
    ifc_path: Path,
    endpoints_path: Path,
    plum_report_path: Path,
    expected_ifc_sha256: str | None,
) -> tuple[str, dict[str, Any], dict[str, Any]]:
    ifc_hash = sha256(ifc_path)
    if expected_ifc_sha256 and expected_ifc_sha256 != ifc_hash:
        raise RuntimeError(
            "formal IFC hash differs from caller-frozen hash: "
            f"expected {expected_ifc_sha256}, found {ifc_hash}"
        )

    endpoints = read_json(endpoints_path)
    plum_report = read_json(plum_report_path)
    if endpoints.get("source_ifc_sha256") != ifc_hash:
        raise RuntimeError("P-201 demand endpoint evidence is stale against the formal IFC")
    if plum_report.get("source", {}).get("ifc_sha256") != ifc_hash:
        raise RuntimeError("PLUM report is stale against the formal IFC")

    records = endpoints.get("demand_endpoints", [])
    service_count = sum(bool(row.get("service_demand_candidate")) for row in records)
    non_service_count = len(records) - service_count
    if (service_count, non_service_count) != (EXPECTED_SERVICE_COUNT, EXPECTED_NON_SERVICE_COUNT):
        raise RuntimeError(
            "unexpected P-201 endpoint classification: "
            f"service={service_count}, non_service={non_service_count}"
        )
    if len({row.get("global_id") for row in records}) != len(records):
        raise RuntimeError("P-201 demand endpoint GlobalIds are not unique")
    if any(row.get("connection_requirement") != "unknown" for row in records):
        raise RuntimeError("P-201 contains a non-unknown connection requirement")
    if any(row.get("is_ifc_distribution_port") for row in records):
        raise RuntimeError("P-201 endpoint evidence unexpectedly claims an IfcDistributionPort")
    if any(row.get("candidate_is_write_authority") for row in records):
        raise RuntimeError("P-201 endpoint evidence unexpectedly grants IFC write authority")

    distribution = plum_report.get("qa", {}).get("distribution_data", {})
    expected_zero = (
        "IfcDistributionPort", "IfcSystem", "IfcDistributionSystem",
        "IfcRelConnectsPorts", "IfcRelConnectsPortToElement",
    )
    nonzero = {key: distribution.get(key) for key in expected_zero if distribution.get(key) != 0}
    if nonzero:
        raise RuntimeError(f"formal distribution topology is not absent: {nonzero}")
    return ifc_hash, endpoints, plum_report


def short_number(value: float) -> str:
    rounded = round(float(value))
    return f"{rounded:,}"


def render_endpoint_row(row: dict[str, Any], index: int, x: int, y: int, width: int) -> str:
    service = bool(row["service_demand_candidate"])
    badge_class = "service-badge" if service else "component-badge"
    status_class = "unknown" if service else "component"
    status = "冷热水需求：unknown" if service else "非服务组合构件｜不生成需求"
    origin = row["object_origin_mm"]
    type_name = row.get("type_name") or "未命名类型"
    predefined = row.get("type_predefined_type") or "NOTDEFINED"
    guid = row["global_id"]
    return (
        f'<g data-global-id="{html.escape(guid)}" data-service-demand="{str(service).lower()}">'
        f'<rect class="endpoint" x="{x}" y="{y}" width="{width}" height="54" rx="7"/>'
        f'<circle class="{badge_class}" cx="{x + 24}" cy="{y + 27}" r="15"/>'
        f'<text class="badge-text" x="{x + 24}" y="{y + 32}">{index:02d}</text>'
        f'<text class="endpoint-title" x="{x + 50}" y="{y + 20}">'
        f'{html.escape(predefined)} · {html.escape(type_name)}</text>'
        f'<text class="endpoint-meta" x="{x + 50}" y="{y + 40}">'
        f'{html.escape(guid)} · 对象原点 XYZ '
        f'{short_number(origin[0])}, {short_number(origin[1])}, {short_number(origin[2])} mm</text>'
        f'<text class="{status_class}" x="{x + width - 12}" y="{y + 31}" text-anchor="end">'
        f'{status}</text></g>'
    )


def render_room_panel(
    room: str,
    rows: list[dict[str, Any]],
    indexed: dict[str, int],
    x: int,
    y: int,
    width: int,
) -> tuple[str, int]:
    panel_height = 58 + len(rows) * 62
    service_count = sum(bool(row["service_demand_candidate"]) for row in rows)
    parts = [
        f'<rect class="panel" x="{x}" y="{y}" width="{width}" height="{panel_height}" rx="14"/>',
        f'<text class="room" x="{x + 18}" y="{y + 33}">{html.escape(room)}</text>',
        f'<text class="room-count" x="{x + width - 18}" y="{y + 32}" text-anchor="end">'
        f'{service_count} 个需求候选 · {len(rows) - service_count} 个非服务构件</text>',
    ]
    for offset, row in enumerate(rows):
        parts.append(render_endpoint_row(row, indexed[row["global_id"]], x + 12, y + 48 + offset * 62, width - 24))
    return "".join(parts), panel_height


def render_svg(records: list[dict[str, Any]], ifc_hash: str) -> str:
    room_rows = {room: [] for room in ROOM_ORDER}
    for row in records:
        room = row.get("container", {}).get("name") or "UNASSIGNED"
        room_rows.setdefault(room, []).append(row)
    unexpected_rooms = sorted(set(room_rows) - set(ROOM_ORDER))
    if unexpected_rooms:
        raise RuntimeError(f"unexpected P-201 room groups: {unexpected_rooms}")

    ordered = [row for room in ROOM_ORDER for row in room_rows[room]]
    indexed = {row["global_id"]: index for index, row in enumerate(ordered, 1)}
    left, _ = render_room_panel("BATHM", room_rows["BATHM"], indexed, 40, 180, 750)
    right_parts: list[str] = []
    right_y = 180
    for room in ("BATHG", "WC", "VVD"):
        panel, height = render_room_panel(room, room_rows[room], indexed, 810, right_y, 750)
        right_parts.append(panel)
        right_y += height + 18

    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="1220" viewBox="0 0 1600 1220">
<style>
.bg{{fill:#f1f5f9}}.panel{{fill:#fff;stroke:#cbd5e1;stroke-width:2}}.endpoint{{fill:#f8fafc;stroke:#e2e8f0}}
.title{{font:700 36px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif;fill:#0f172a}}
.meta{{font:17px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif;fill:#475569}}
.warning{{font:700 17px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif;fill:#991b1b}}
.room{{font:700 24px ui-monospace,SFMono-Regular,monospace;fill:#0f4c81}}.room-count{{font:15px -apple-system,"PingFang SC",sans-serif;fill:#475569}}
.service-badge{{fill:#dbeafe;stroke:#3b82f6}}.component-badge{{fill:#e2e8f0;stroke:#64748b}}
.badge-text{{font:700 11px ui-monospace,monospace;fill:#1e3a8a;text-anchor:middle}}
.endpoint-title{{font:700 14px -apple-system,"PingFang SC",sans-serif;fill:#1e293b}}
.endpoint-meta{{font:11px ui-monospace,SFMono-Regular,monospace;fill:#64748b}}
.unknown{{font:700 12px -apple-system,"PingFang SC",sans-serif;fill:#b45309}}.component{{font:700 12px -apple-system,"PingFang SC",sans-serif;fill:#475569}}
.footer{{font:13px ui-monospace,SFMono-Regular,monospace;fill:#64748b}}
</style>
<rect class="bg" width="1600" height="1220"/>
<text class="title" x="40" y="58">P-201 给水 / 热水需求端点候选</text>
<text class="meta" x="40" y="94">24 个服务需求候选 + 3 个非服务组合构件｜按当前 IFC 对象原点分区登记｜位置不是接口或粗装点</text>
<text class="warning" x="40" y="126">无 IfcDistributionPort / System / 正式管线拓扑；冷热水 unknown 保持 unknown；不含管径、水压、设备接口推定</text>
<text class="meta" x="40" y="154">本图是待复核需求清单，不是开发商交付参考图，也不是最终给水、热水管线施工图。</text>
{left}{''.join(right_parts)}
<text class="footer" x="40" y="1192">IFC SHA-256 {ifc_hash} · automatic_ifc_write_allowed=false · construction_release_ready=false</text>
</svg>
'''


def main() -> None:
    root_default = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root_default)
    parser.add_argument("--input-ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument("--endpoints", type=Path, default=Path("build/plum/p201-demand-endpoints.json"))
    parser.add_argument("--plum-report", type=Path, default=Path("build/plum/plum-report.json"))
    parser.add_argument("--output-svg", type=Path, default=Path("drawings/P-201-water-hotwater-demand-candidate.svg"))
    parser.add_argument("--proof-png", type=Path, default=Path("build/plum/P-201-water-hotwater-demand-candidate.png"))
    parser.add_argument("--report", type=Path, default=Path("build/plum/p201-water-hotwater-demand-candidate.json"))
    parser.add_argument("--expected-ifc-sha256")
    parser.add_argument("--chrome", type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    resolve = lambda path: path if path.is_absolute() else root / path
    ifc_path = resolve(args.input_ifc)
    endpoints_path = resolve(args.endpoints)
    plum_report_path = resolve(args.plum_report)
    svg_path = resolve(args.output_svg)
    png_path = resolve(args.proof_png)
    report_path = resolve(args.report)
    ensure_safe_output(svg_path, ".svg", ifc_path)
    ensure_safe_output(png_path, ".png", ifc_path)
    ensure_safe_output(report_path, ".json", ifc_path)

    ifc_hash, endpoints, plum_report = validate_sources(
        ifc_path, endpoints_path, plum_report_path, args.expected_ifc_sha256
    )
    records = endpoints["demand_endpoints"]
    room_counts = Counter((row.get("container") or {}).get("name") or "UNASSIGNED" for row in records)
    service_count = sum(bool(row["service_demand_candidate"]) for row in records)

    svg_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path.write_text(render_svg(records, ifc_hash), encoding="utf-8")
    png_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(find_chrome(args.chrome)), "--headless=new", "--disable-gpu", "--hide-scrollbars",
            "--force-device-scale-factor=1", f"--screenshot={png_path}", "--window-size=1600,1220",
            svg_path.resolve().as_uri(),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if not png_path.is_file() or png_path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError("P-201 proof output is not a valid PNG")

    ending_ifc_hash = sha256(ifc_path)
    if ending_ifc_hash != ifc_hash:
        raise RuntimeError("formal IFC changed during P-201 candidate generation")

    distribution = plum_report["qa"]["distribution_data"]
    report: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_p201_water_hotwater_demand_endpoint_candidate",
        "source_ifc_sha256": ifc_hash,
        "source": {
            "ifc": {"path": str(ifc_path), "sha256": ifc_hash},
            "demand_endpoints": {"path": str(endpoints_path), "sha256": sha256(endpoints_path)},
            "plum_report": {"path": str(plum_report_path), "sha256": sha256(plum_report_path)},
        },
        "summary": {
            "registered_object_count": len(records),
            "service_demand_candidate_count": service_count,
            "non_service_component_count": len(records) - service_count,
            "room_counts": dict(sorted(room_counts.items())),
            "connection_requirement_counts": dict(Counter(row["connection_requirement"] for row in records)),
            "ifc_distribution_port_count": distribution["IfcDistributionPort"],
            "ifc_system_count": distribution["IfcSystem"],
            "ifc_distribution_system_count": distribution["IfcDistributionSystem"],
            "ifc_rel_connects_ports_count": distribution["IfcRelConnectsPorts"],
            "invented_pipe_route_count": 0,
            "invented_pipe_size_count": 0,
            "invented_pressure_value_count": 0,
            "invented_equipment_interface_count": 0,
        },
        "records": records,
        "outputs": {
            "svg": {"path": str(svg_path), "sha256": sha256(svg_path)},
            "proof_png": {"path": str(png_path), "sha256": sha256(png_path)},
        },
        "gates": {
            "source_hashes_current": True,
            "caller_frozen_hash_checked": bool(args.expected_ifc_sha256),
            "all_registered_objects_drawn": True,
            "service_and_non_service_split_explicit": True,
            "unknown_water_demands_preserved": all(row["connection_requirement"] == "unknown" for row in records),
            "formal_distribution_topology_present": False,
            "ifc_unchanged_during_generation": ending_ifc_hash == ifc_hash,
            "automatic_ifc_write_allowed": False,
            "construction_release_ready": False,
        },
        "open_release_items": [
            "确认每个需求端点的冷水、热水或无需给水属性",
            "取得设备厂家接口、流量与压力要求",
            "完成正式管线路由、管径、阀件、保温和系统拓扑设计",
            "现场复核对象原点与真实接口/粗装点的偏差",
        ],
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(report_path), "summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
