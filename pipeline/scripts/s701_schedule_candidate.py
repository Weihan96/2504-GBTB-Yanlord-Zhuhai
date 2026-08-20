#!/usr/bin/env python3
"""Compile an evidence-bounded S-701 schedule candidate without writing IFC."""

from __future__ import annotations

import argparse
import csv
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

from equipment_ssot import load_canonical, requirement_map, split_ids, validate as validate_equipment_ssot


FIELDS = [
    "schedule_id", "section", "source_key", "item_name", "confirmed_scope",
    "candidate_or_observed_scope", "unresolved_for_release", "evidence_reference",
    "review_status", "automatic_ifc_write_allowed", "construction_release_ready",
    "source_ifc_sha256",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def row(
    schedule_id: str,
    section: str,
    source_key: str,
    item_name: str,
    confirmed_scope: str,
    candidate_scope: str,
    unresolved: str,
    evidence: str,
    review_status: str,
    ifc_hash: str,
) -> dict[str, str]:
    return {
        "schedule_id": schedule_id,
        "section": section,
        "source_key": source_key,
        "item_name": item_name,
        "confirmed_scope": confirmed_scope,
        "candidate_or_observed_scope": candidate_scope,
        "unresolved_for_release": unresolved,
        "evidence_reference": evidence,
        "review_status": review_status,
        "automatic_ifc_write_allowed": "false",
        "construction_release_ready": "false",
        "source_ifc_sha256": ifc_hash,
    }


def build_rows(root: Path, ifc_hash: str) -> tuple[list[dict[str, str]], dict[str, Path]]:
    sources = {
        "equipment": root / "pipeline/decisions/equipment-register.csv",
        "requirements": root / "pipeline/decisions/equipment-installation-requirements.csv",
        "evidence": root / "pipeline/decisions/source-evidence-register.csv",
        "wfin": root / "pipeline/decisions/wfin-open-issues.csv",
    }
    for name, path in sources.items():
        if not path.is_file():
            raise RuntimeError(f"missing S-701 source {name}: {path}")

    validate_equipment_ssot(root)
    _, equipment, requirements, evidence_rows = load_canonical(root)
    requirements_by_equipment = requirement_map(requirements)
    evidence_by_id = {item["source_id"]: item for item in evidence_rows}
    result: list[dict[str, str]] = []

    furniture = [item for item in equipment if item["legacy_kind"] == "furniture_product" and item["schedule_included"] == "yes"]
    for index, item in enumerate(furniture, 1):
        if item["decision_status"] != "confirmed":
            raise RuntimeError(f"furniture identity is not confirmed: {item['item_name']}")
        source_urls = [evidence_by_id[sid]["source_url"] for sid in split_ids(item["source_ids"]) if sid in evidence_by_id]
        result.append(row(
            f"S701-FUR-{index:02d}", "家具产品身份", item["selector_value"],
            " ".join(part for part in (item["manufacturer"], item["item_name"], item["variant"]) if part),
            "厂家与产品系列身份已确认", item["use_location_confirmed"] or "用途未填写",
            "本项目最终规格、数量、五金配置、安装图与现场接口", "；".join(source_urls) or item["source_ids"],
            "identity_confirmed_installation_pending", ifc_hash,
        ))

    appliances = [item for item in equipment if item["legacy_kind"] == "appliance" and item["schedule_included"] == "yes"]
    for item in appliances:
        confirmed = []
        if item["storage_location_confirmed"]:
            confirmed.append(f"存放={item['storage_location_confirmed']}")
        if item["use_location_confirmed"]:
            confirmed.append(f"使用={item['use_location_confirmed']}")
        candidate = []
        if item["storage_location_candidate"]:
            candidate.append(f"存放候选={item['storage_location_candidate']}")
        if item["use_location_candidate"]:
            candidate.append(f"使用候选={item['use_location_candidate']}")
        req = requirements_by_equipment[item["equipment_id"]]
        unresolved = [name for name, value in (
            ("型号", item["model"]), ("铭牌功率", req.get("rated_power", "")),
            ("证据", item["source_ids"]),
        ) if not value]
        unresolved.extend(
            name for name, value in (
                ("给水接口", req.get("water_required", "")), ("排水接口", req.get("drain_required", "")),
                ("燃气", req.get("gas_required", "")), ("通风", req.get("ventilation_required", "")),
            ) if value in {"", "model_dependent"}
        )
        result.append(row(
            f"S701-{item['equipment_id']}", "家电与移动厨电", item["equipment_id"], item["item_name"],
            "；".join(confirmed) or "尚无关闭项", "；".join(candidate) or "无位置候选",
            "、".join(unresolved) or "仍须厂家安装图与现场接口复核", item["source_ids"] or "source-evidence-register.csv",
            (
                "confirmed_input_installation_review_pending"
                if item["decision_status"] == "confirmed"
                else "partial_input" if item["decision_status"] == "partial" else "input_required"
            ),
            ifc_hash,
        ))

    kitchen_tool_rows = [item for item in equipment if item["legacy_kind"] == "owner_kitchen_tool" and item["schedule_included"] == "yes"]
    for item in kitchen_tool_rows:
        result.append(row(
            f"S701-{item['equipment_id']}", "家电与移动厨电", item["equipment_id"], item["item_name"],
            f"已购；存放分组={item['storage_location_confirmed']}", item["variant"],
            "到货实物 SKU／包络、抽屉内净尺寸、摆样、导轨与承重",
            item["source_ids"], "purchased_delivery_unverified_storage_group_confirmed", ifc_hash,
        ))

    door_window_rows = [item for item in equipment if item["legacy_kind"] == "a104" and item["schedule_included"] == "yes"]
    for item in door_window_rows:
        req = requirements_by_equipment[item["equipment_id"]]
        operation = req.get("operation_type", "") or "NOTDEFINED"
        unresolved = "厂家门窗表、五金型号、安装/收口与现场复核"
        if operation == "NOTDEFINED":
            unresolved = "开启方向/合页侧、" + unresolved
        result.append(row(
            f"S701-{item['legacy_id']}", "门窗与五金", item["ifc_global_ids"], item["item_name"],
            f"身份/定位；名义尺寸 {req.get('nominal_width','')}×{req.get('nominal_height','')} mm",
            f"OperationType={operation}", unresolved,
            item["source_ids"], "observed_geometry_hardware_pending", ifc_hash,
        ))

    owner_window_rows = [
        item for item in equipment
        if item["legacy_kind"] in {"owner_window_treatment", "owner_window_hardware"}
        and item["schedule_included"] == "yes"
    ]
    for item in owner_window_rows:
        blocking = [
            req["parameter_key"] for req in requirements
            if req["equipment_id"] == item["equipment_id"] and req["blocks_release"] == "yes"
        ]
        result.append(row(
            f"S701-{item['equipment_id']}", "门窗与五金", item["equipment_id"], item["item_name"],
            item["variant"], f"采购状态={item['procurement_status']}",
            "、".join(blocking) or "项目安装图与现场复核", item["source_ids"],
            "owner_direction_detail_pending", ifc_hash,
        ))

    hvac_rows = [item for item in equipment if item["legacy_kind"] == "hvac_interface" and item["ifc_type_global_id"]]
    for item in hvac_rows:
        result.append(row(
            f"S701-HVAC-{item['ifc_type_name']}", "暖通设备类型", item["ifc_type_global_id"], item["ifc_type_name"],
            f"当前 IFC 已分配类型；实例数 {item['quantity']}", item["model"],
            "最终精确型号、接口坐标、功率、风量、检修条件", item["source_ids"],
            "observed_type_not_final_selection", ifc_hash,
        ))

    hvac_insulation_rows = [item for item in equipment if item["legacy_kind"] == "hvac_insulation" and item["schedule_included"] == "yes"]
    for item in hvac_insulation_rows:
        result.append(row(
            f"S701-{item['equipment_id']}", "暖通设备类型", item["equipment_id"], item["item_name"],
            "现场观察 11 种 ID×TK 规格；华美 Class 1 产品族性能已登记",
            "现场观察值不是最终设计厚度", "M-401 逐段管径与防结露计算、日立／安装方按图复核",
            item["source_ids"], "existing_observation_final_schedule_pending", ifc_hash,
        ))

    custom_drain_rows = [item for item in equipment if item["legacy_kind"] == "owner_custom_drain" and item["schedule_included"] == "yes"]
    for item in custom_drain_rows:
        req = requirements_by_equipment[item["equipment_id"]]
        result.append(row(
            f"S701-{item['equipment_id']}", "洁具与排水", item["equipment_id"], item["item_name"],
            "定制渠道与组合方向已由业主确认", item["variant"],
            "、".join(name for name, value in (
                ("房间/数量", req.get("room_mapping", "")),
                ("准确组件型号", req.get("component_models", "")),
                ("厂家 shop drawing", req.get("shop_drawing", "")),
                ("排水/防水/完成面接口", req.get("interface_center_coordinates", "")),
            ) if not value or value == "unknown"),
            item["source_ids"], "owner_selected_custom_direction_shop_drawing_pending", ifc_hash,
        ))

    purchased_faucet_ids = {"SAN-023", "SAN-024", "SAN-025"}
    purchased_faucet_rows = [
        item for item in equipment
        if item["equipment_id"] in purchased_faucet_ids and item["schedule_included"] == "yes"
    ]
    for item in purchased_faucet_rows:
        blocking_keys = [
            requirement["parameter_key"]
            for requirement in requirements
            if requirement["equipment_id"] == item["equipment_id"]
            and requirement["blocks_release"] == "yes"
        ]
        result.append(row(
            f"S701-{item['equipment_id']}", "洁具与排水", item["equipment_id"], item["item_name"],
            f"已购买，数量 {item['quantity']}；订单款式={item['variant']}",
            f"项目位置候选={item['use_location_candidate']}；{item['model']}",
            "、".join(blocking_keys),
            item["source_ids"], "purchased_delivery_unverified_installation_pending", ifc_hash,
        ))

    owner_product_component_rows = [item for item in equipment if item["legacy_kind"] == "owner_product_component" and item["schedule_included"] == "yes"]
    for item in owner_product_component_rows:
        blocking_keys = [
            requirement["parameter_key"] for requirement in requirements
            if requirement["equipment_id"] == item["equipment_id"] and requirement["blocks_release"] == "yes"
        ]
        result.append(row(
            f"S701-{item['equipment_id']}", "洁具与排水", item["equipment_id"], item["item_name"],
            item["variant"], f"优先候选={item['model']}；未采购",
            "、".join(blocking_keys), item["source_ids"], "candidate_project_detail_external_review_pending", ifc_hash,
        ))

    for item in read_csv(sources["wfin"]):
        if item["issue_id"] not in {"WFIN-R02", "WFIN-R03", "WFIN-R05"}:
            continue
        result.append(row(
            f"S701-{item['issue_id']}", "墙面材料系统", item["issue_id"], item["scope"],
            item["current_evidence"], "材料方向候选，不代表产品确认",
            item["required_action_or_decision"], "wfin-open-issues.csv",
            item["status"], ifc_hash,
        ))
    return result, sources


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def render_svg(rows: list[dict[str, str]], ifc_hash: str) -> str:
    counts = Counter(item["section"] for item in rows)
    sections = [
        ("家具产品身份", "已确认身份；规格/安装仍待闭合"),
        ("家电与移动厨电", "存放位置与使用位置严格分列"),
        ("门窗与五金", "定位/名义尺寸可查；厂家五金未冻结"),
        ("暖通设备类型", "仅登记 IFC 既有类型，不代表最终选型"),
        ("洁具与排水", "定制渠道与组合方向已确认；准确接口仍待 shop drawing"),
        ("墙面材料系统", "Tadelakt/大白墙/中厨同材大板＋浅置物架为意图"),
    ]
    cards = []
    for index, (name, note) in enumerate(sections):
        y = 190 + index * 130
        cards.append(
            f'<rect class="card" x="55" y="{y}" width="910" height="105" rx="12"/>'
            f'<text class="count" x="95" y="{y + 64}">{counts[name]}</text>'
            f'<text class="section" x="180" y="{y + 42}">{html.escape(name)}</text>'
            f'<text class="note" x="180" y="{y + 75}">{html.escape(note)}</text>'
        )
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="1020" height="930" viewBox="0 0 1020 930">
<style>
.bg{{fill:#f4f7fb}}.card{{fill:#fff;stroke:#cbd5e1;stroke-width:2}}.title{{font:700 36px -apple-system,"PingFang SC",sans-serif;fill:#0f172a}}
.meta{{font:17px -apple-system,"PingFang SC",sans-serif;fill:#475569}}.count{{font:700 38px ui-monospace,monospace;fill:#0f4c81;text-anchor:middle}}
.section{{font:700 23px -apple-system,"PingFang SC",sans-serif;fill:#1e293b}}.note{{font:17px -apple-system,"PingFang SC",sans-serif;fill:#92400e}}
.gate{{font:700 17px ui-monospace,monospace;fill:#991b1b}}.footer{{font:14px ui-monospace,monospace;fill:#64748b}}
</style><rect class="bg" width="1020" height="930"/>
<text class="title" x="55" y="65">S-701 材料、设备及五金表候选</text>
<text class="meta" x="55" y="105">{len(rows)} 条证据化记录 · 已确认、候选、未决严格分栏 · 不是下单表</text>
<text class="gate" x="55" y="145">automatic_ifc_write_allowed=false · construction_release_ready=false</text>
{''.join(cards)}
<text class="meta" x="55" y="875">完整逐项内容见 s701-schedule-review.csv；须结合用户输入、厂家安装图及现场实测关闭。</text>
<text class="footer" x="55" y="906">IFC SHA-256 {ifc_hash}</text>
</svg>'''


def find_chrome(explicit: Path | None) -> Path:
    candidates = [explicit, Path(os.environ["CHROME_BIN"]) if os.environ.get("CHROME_BIN") else None,
                  Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")]
    for name in ("google-chrome", "chromium"):
        binary = shutil.which(name)
        if binary:
            candidates.append(Path(binary))
    for candidate in candidates:
        if candidate and candidate.is_file():
            return candidate
    raise RuntimeError("Chrome/Chromium not found")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--input-ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument("--expected-ifc-sha256")
    parser.add_argument("--review-csv", type=Path, default=Path("pipeline/decisions/s701-schedule-review.csv"))
    parser.add_argument("--output-svg", type=Path, default=Path("drawings/S-701-material-equipment-hardware-schedule-candidate.svg"))
    parser.add_argument("--proof-png", type=Path, default=Path("build/s701/S-701-material-equipment-hardware-schedule-candidate.png"))
    parser.add_argument("--report", type=Path, default=Path("build/s701/s701-report.json"))
    parser.add_argument("--chrome", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    resolve = lambda path: path if path.is_absolute() else root / path
    ifc_path = resolve(args.input_ifc)
    ifc_hash = sha256(ifc_path)
    if args.expected_ifc_sha256 and args.expected_ifc_sha256 != ifc_hash:
        raise RuntimeError(f"formal IFC hash differs from caller-frozen hash: expected {args.expected_ifc_sha256}, found {ifc_hash}")
    rows, sources = build_rows(root, ifc_hash)
    review_csv, svg, png, report_path = map(resolve, (args.review_csv, args.output_svg, args.proof_png, args.report))
    write_csv(review_csv, rows)
    svg.parent.mkdir(parents=True, exist_ok=True)
    svg.write_text(render_svg(rows, ifc_hash), encoding="utf-8")
    png.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(find_chrome(args.chrome)), "--headless=new", "--disable-gpu", "--hide-scrollbars",
                    f"--screenshot={png}", "--window-size=1020,930", svg.resolve().as_uri()], check=True, capture_output=True, text=True)
    if png.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError("S-701 proof output is not PNG")
    report: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_s701_evidence_schedule_candidate",
        "source_ifc_sha256": ifc_hash,
        "source": {"ifc_sha256": ifc_hash, "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in sources.items()}},
        "summary": {"record_count": len(rows), "section_counts": dict(Counter(item["section"] for item in rows)),
                    "automatic_ifc_write_allowed": False, "construction_release_ready": False},
        "outputs": {"review_csv": {"path": str(review_csv), "sha256": sha256(review_csv)},
                    "svg": {"path": str(svg), "sha256": sha256(svg)}, "proof_png": {"path": str(png), "sha256": sha256(png)}},
        "gates": {"source_hashes_current": True, "confirmed_candidate_unresolved_split": True,
                  "automatic_ifc_write_allowed": False, "construction_release_ready": False},
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(report_path), "summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
