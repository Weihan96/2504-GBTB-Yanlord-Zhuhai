#!/usr/bin/env python3
"""Generate a read-only D-601/D-602 variable detail index candidate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


CSV_FIELDS = [
    "node_id",
    "sheet_id",
    "category",
    "scope",
    "confirmed_evidence",
    "variable_parameters",
    "unresolved_material_or_product",
    "unresolved_manufacturer",
    "unresolved_construction",
    "review_owner",
    "source_references",
    "review_required",
    "review_status",
    "automatic_ifc_write_allowed",
    "construction_release_ready",
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


def require_text(path: Path, fragments: list[str]) -> str:
    text = path.read_text(encoding="utf-8-sig")
    missing = [fragment for fragment in fragments if fragment not in text]
    if missing:
        raise RuntimeError(f"required DET1 evidence missing from {path}: {missing}")
    return text


def find_chrome(explicit: Path | None) -> Path:
    candidates = [
        explicit,
        Path(os.environ["CHROME_BIN"]) if os.environ.get("CHROME_BIN") else None,
        Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"),
        Path("/Applications/Chromium.app/Contents/MacOS/Chromium"),
    ]
    for name in ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser"):
        binary = shutil.which(name)
        if binary:
            candidates.append(Path(binary))
    for candidate in candidates:
        if candidate and candidate.is_file():
            return candidate
    raise RuntimeError("Chrome/Chromium not found; pass --chrome or set CHROME_BIN")


def node(
    node_id: str,
    sheet_id: str,
    category: str,
    scope: str,
    confirmed_evidence: str,
    variables: str,
    material: str,
    construction: str,
    owner: str,
    references: str,
    ifc_hash: str,
    manufacturer: str = "待产品系统确认后填写；当前不指定厂家",
) -> dict[str, str]:
    return {
        "node_id": node_id,
        "sheet_id": sheet_id,
        "category": category,
        "scope": scope,
        "confirmed_evidence": confirmed_evidence,
        "variable_parameters": variables,
        "unresolved_material_or_product": material,
        "unresolved_manufacturer": manufacturer,
        "unresolved_construction": construction,
        "review_owner": owner,
        "source_references": references,
        "review_required": "yes",
        "review_status": "candidate_pending_review",
        "automatic_ifc_write_allowed": "false",
        "construction_release_ready": "false",
        "source_ifc_sha256": ifc_hash,
    }


def build_nodes(ifc_hash: str) -> list[dict[str, str]]:
    return [
        node(
            "D601-N01", "D-601", "墙面阳角", "主要通行动线与高频接触区",
            "已有需求：阳角采用圆角、倒角或护角做法；未选定其中一种",
            "位置｜转角形式｜可见面宽度｜与墙面完成面的关系",
            "护角/饰面产品与颜色样板未定",
            "基层、固定方式、转角层次及耐撞性做法未定",
            "设计/材料供应商/现场", "深化设计注意事项#9；施工图检查清单 D-601", ifc_hash,
        ),
        node(
            "D601-N02", "D-601", "地面材料交界", "F11/F21 地板与 F20 银白洞石岩板等交界",
            "F11/F21=地板、F20=银白洞石岩板；已确认完成面构造区 50 mm，但 IFC 仍为零厚度参考面",
            "交界位置｜收口型式｜缝宽｜高差｜伸缩/维修条件",
            "地板和岩板最终产品、规格、样板未定",
            "50 mm 构造区内真实分层、基层和粘结/龙骨系统未闭合",
            "设计/材料供应商/现场", "PM A-105 地坪材料与构造区；a105-floor-review.csv", ifc_hash,
        ),
        node(
            "D601-N03", "D-601", "门槛与推拉门地轨", "M01–M08 门槛及推拉门地轨联动",
            "门槛已确认齐平优先；只有机械证明溢水风险时才采用高差与 45° 倒角",
            "门号｜两侧完成面｜齐平/防水策略｜地轨型式｜门下净空",
            "门槛石/金属收口及地轨产品未定",
            "防水收口、门下净空、溢水路径及固定构造未定",
            "设计/门窗厂家/防水/现场", "a105-threshold-review.csv；PM A-105 门槛策略", ifc_hash,
        ),
        node(
            "D601-N04", "D-601", "墙地顶收口", "踢脚、墙顶交界与不同材料收口",
            "墙面候选已区分大白墙与 Tadelakt 意图；最终产品选择已延后",
            "房间/墙段｜踢脚型式｜墙顶收口｜材料分界｜阴影缝/密封策略",
            "Tadelakt 系统、大白墙涂料系统、踢脚产品和样板未定",
            "完成面总厚、基层、防水衔接与分界构造未定",
            "设计/材料供应商/现场", "wfin-open-issues.csv WFIN-R01/R02/R03；施工图检查清单 D-601", ifc_hash,
        ),
        node(
            "D601-N05", "D-601", "柜墙交界", "固定家具、台面、饰面与墙体交界",
            "INT1 现有对象只作协调包络，不作加工、开孔或粗装尺寸",
            "空间/柜组｜收口位置｜可调节缝｜密封｜拆装与检修路径",
            "柜体、台面、封板及收口产品未定",
            "墙体基层、固定件、防潮/密封和现场误差吸收方式未定",
            "室内设计/全屋定制/现场", "PM INT1；int1-drawing-report.json", ifc_hash,
        ),
        node(
            "D601-N06", "D-601", "阳台管道包覆", "阳台水管、阀门、接头与检修口",
            "检查清单要求管道包覆不得阻断阀门、接头检修",
            "管道/阀门位置｜包覆边界｜检修口｜通风/冷凝水｜拆卸顺序",
            "包覆面材、检修口和五金未定",
            "龙骨/板材、防潮、密封、管道振动隔离及可拆结构未定",
            "设计/给排水/全屋定制/现场", "施工图检查清单 D-601", ifc_hash,
        ),
        node(
            "D602-N01", "D-602", "防水范围", "卫生间、淋浴区、厨房及阳台",
            "湿区与墙面材意图已有候选；防水范围和高度尚未确认",
            "空间｜墙/地范围｜上翻高度｜终止位置｜与门窗/柜体衔接",
            "防水产品系统、底涂、增强层和密封材未定",
            "基层处理、遍数/厚度、上翻与收头构造未定",
            "防水供应商/设计/施工/监理", "施工图检查清单 D-602；wfin-open-issues.csv WFIN-R02", ifc_hash,
        ),
        node(
            "D602-N02", "D-602", "淋浴排水", "淋浴区、地漏与毛发过滤",
            "18 个湿区坡面的 1.047% 及箭头所示下坡方向已确认；毛发过滤需求已记录",
            "湿区｜地漏位置/型式｜排水路径｜毛发网取出方向｜检修与清理空间",
            "地漏、线性排水和毛发过滤产品未定",
            "地漏翼环、防水压接、找坡层和可清理构造未定",
            "给排水/防水/洁具供应商/现场", "PM 卫生间坡向；深化设计注意事项#8；a105-floor-review.csv", ifc_hash,
        ),
        node(
            "D602-N03", "D-602", "管根与穿透", "供排水、套管、墙地穿透点",
            "现有管线与洁具只有位置/包络协调证据；未得到厂家粗装图",
            "穿透点｜管径/套管｜穿墙/穿地方向｜防水收头｜防火/隔声｜检修",
            "套管、止水、密封与防火封堵系统未定",
            "预留洞、套管、防水附加层、封堵与检修顺序未定",
            "给排水/防水/洁具厂家/现场", "PM INT1 停止条件；施工图检查清单 D-602", ifc_hash,
        ),
        node(
            "D602-N04", "D-602", "湿区门槛与墙地交界", "卫生间干湿区、门口与墙地阴角",
            "门槛齐平优先已确认；防水和排水路径必须机械闭合后才能发布",
            "门口/墙段｜防水连续性｜附加层｜收头｜齐平条件｜溢水检查",
            "门槛、防水、阴角增强和密封产品未定",
            "基层阴角、附加层、门槛下连续防水及完成面衔接未定",
            "防水/设计/门窗/现场", "a105-threshold-review.csv；PM DET1 停止条件", ifc_hash,
        ),
        node(
            "D602-N05", "D-602", "厨房与阳台防水", "厨房、阳台及柜体/管道边界",
            "检查清单要求厨房与阳台纳入防水范围；当前无可发布的范围或高度",
            "空间｜用水点｜潜在溢水路径｜墙地范围｜柜体/管道收口｜排水与检修",
            "防水、密封、柜体防潮和检修口产品未定",
            "防水连续层、管根、墙地收头、柜下可视检与排水路径未定",
            "防水/厨柜/给排水/现场", "施工图检查清单 D-602；PM 阳台管道检修要求", ifc_hash,
        ),
        node(
            "D602-N06", "D-602", "试水与验收", "防水施工、保护层和覆盖前验收",
            "检查清单明确要求试水；当前无已确认的产品系统或验收参数",
            "空间｜试验阶段｜封堵边界｜水位/持续时间｜观察面｜记录与签认",
            "与最终防水系统配套的修补和保护材料未定",
            "试水条件、时长、验收标准、缺陷修补与复验流程未定",
            "防水供应商/施工/监理/业主", "施工图检查清单 D-602", ifc_hash,
        ),
    ]


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def clip(value: str, limit: int = 35) -> str:
    return value if len(value) <= limit else value[: limit - 1] + "…"


def render_svg(rows: list[dict[str, str]], ifc_hash: str) -> str:
    panels: list[str] = []
    for panel_index, sheet_id in enumerate(("D-601", "D-602")):
        x = 45 + panel_index * 775
        panel_rows = [row for row in rows if row["sheet_id"] == sheet_id]
        panels.append(f'<rect class="panel" x="{x}" y="155" width="735" height="760" rx="18"/>')
        panels.append(f'<text class="sheet" x="{x + 28}" y="205">{sheet_id}</text>')
        subtitle = "材料与收口节点" if sheet_id == "D-601" else "防水与湿区节点"
        panels.append(f'<text class="subtitle" x="{x + 150}" y="204">{subtitle}</text>')
        for index, row in enumerate(panel_rows):
            y = 235 + index * 108
            panels.extend([
                f'<rect class="row" x="{x + 20}" y="{y}" width="695" height="92" rx="10"/>',
                f'<circle class="badge" cx="{x + 55}" cy="{y + 31}" r="22"/>',
                f'<text class="badge-text" x="{x + 55}" y="{y + 36}">{html.escape(row["node_id"].split("-")[-1])}</text>',
                f'<text class="category" x="{x + 90}" y="{y + 25}">{html.escape(row["category"])}</text>',
                f'<text class="scope" x="{x + 90}" y="{y + 49}">{html.escape(clip(row["scope"], 39))}</text>',
                f'<text class="pending" x="{x + 90}" y="{y + 73}">变量：{html.escape(clip(row["variable_parameters"], 40))}</text>',
                f'<text class="status" x="{x + 610}" y="{y + 27}">待复核</text>',
                f'<text class="lock" x="{x + 610}" y="{y + 67}">不写 IFC</text>',
            ])
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="1000" viewBox="0 0 1600 1000">
  <style>
    .bg {{ fill:#f3f4f6; }} .panel {{ fill:#ffffff; stroke:#cbd5e1; stroke-width:2; }}
    .title {{ font:700 38px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif; fill:#0f172a; }}
    .meta {{ font:18px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif; fill:#475569; }}
    .sheet {{ font:700 30px ui-monospace,SFMono-Regular,monospace; fill:#0f4c81; }}
    .subtitle {{ font:600 22px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif; fill:#334155; }}
    .row {{ fill:#f8fafc; stroke:#e2e8f0; }} .badge {{ fill:#dbeafe; stroke:#60a5fa; }}
    .badge-text {{ font:700 14px ui-monospace,SFMono-Regular,monospace; fill:#1e3a8a; text-anchor:middle; }}
    .category {{ font:700 20px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif; fill:#111827; }}
    .scope {{ font:17px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif; fill:#475569; }}
    .pending {{ font:15px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif; fill:#92400e; }}
    .status {{ font:700 16px -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif; fill:#b45309; }}
    .lock {{ font:700 14px ui-monospace,SFMono-Regular,monospace; fill:#991b1b; }}
    .footer {{ font:15px ui-monospace,SFMono-Regular,monospace; fill:#64748b; }}
  </style>
  <rect class="bg" width="1600" height="1000"/>
  <text class="title" x="45" y="65">DET1 待复核节点索引候选</text>
  <text class="meta" x="45" y="103">只读索引 · 已确认意图 + 可填变量 · 不构成材料、厂家或施工尺寸确认</text>
  <text class="meta" x="45" y="132">12 个节点均为 candidate_pending_review；automatic_ifc_write_allowed=false；construction_release_ready=false</text>
  {''.join(panels)}
  <text class="footer" x="45" y="963">IFC SHA-256 {ifc_hash}</text>
</svg>
'''


def render_png(svg_path: Path, png_path: Path, chrome: Path) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            str(chrome),
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            "--force-device-scale-factor=1",
            f"--screenshot={png_path}",
            "--window-size=1600,1000",
            svg_path.resolve().as_uri(),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    if not png_path.is_file() or png_path.stat().st_size == 0:
        raise RuntimeError("DET1 proof PNG was not created")
    if png_path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError("DET1 proof output is not a PNG")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--input-ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument("--expected-ifc-sha256")
    parser.add_argument("--review-csv", type=Path, default=Path("pipeline/decisions/det1-detail-review.csv"))
    parser.add_argument("--output-svg", type=Path, default=Path("drawings/D601-D602-detail-index-candidate.svg"))
    parser.add_argument("--report", type=Path, default=Path("build/det1/detail-index-candidate.json"))
    parser.add_argument("--open-issues", type=Path, default=Path("build/det1/open-issues.md"))
    parser.add_argument("--proof-png", type=Path, default=Path("build/det1/D601-D602-detail-index-candidate.png"))
    parser.add_argument("--chrome", type=Path)
    args = parser.parse_args()

    root = args.root.resolve()
    resolve = lambda path: path if path.is_absolute() else root / path
    paths = {
        "ifc": resolve(args.input_ifc),
        "pm": root / "drawings/滨海湾装修施工图深化工作管理.md",
        "checklist": root / "drawings/滨海湾施工图纸检查清单.md",
        "attention": root / "drawings/滨海湾深化设计注意事项.md",
        "a105_floor": root / "pipeline/decisions/a105-floor-review.csv",
        "a105_threshold": root / "pipeline/decisions/a105-threshold-review.csv",
        "wfin_issues": root / "pipeline/decisions/wfin-open-issues.csv",
        "int1_report": root / "build/int1/int1-drawing-report.json",
    }
    for label, path in paths.items():
        if not path.is_file():
            raise RuntimeError(f"missing DET1 input {label}: {path}")

    require_text(paths["pm"], ["DET1｜材料收口、防水与关键节点", "构造层总厚度不闭合", "1.047%"])
    require_text(paths["checklist"], ["D-601 材料与收口节点", "D-602 防水与湿区节点", "试水要求"])
    require_text(paths["attention"], ["墙面阳角太尖", "洗澡会堵头发"])

    threshold_rows = read_csv(paths["a105_threshold"])
    if len(threshold_rows) != 8 or any(row["review_group"] != "A105-R03" for row in threshold_rows):
        raise RuntimeError("A-105 threshold evidence is incomplete")
    wfin_rows = {row["issue_id"]: row for row in read_csv(paths["wfin_issues"])}
    if not {"WFIN-R01", "WFIN-R02", "WFIN-R03"}.issubset(wfin_rows):
        raise RuntimeError("WFIN material/product open issues are incomplete")

    ifc_hash = sha256(paths["ifc"])
    if args.expected_ifc_sha256 and args.expected_ifc_sha256 != ifc_hash:
        raise RuntimeError(f"formal IFC SHA mismatch: expected {args.expected_ifc_sha256}, got {ifc_hash}")
    int1_report = json.loads(paths["int1_report"].read_text(encoding="utf-8"))
    if int1_report.get("source_ifc_sha256") != ifc_hash:
        raise RuntimeError("INT1 drawing evidence is stale against the formal IFC")

    rows = build_nodes(ifc_hash)
    review_csv = resolve(args.review_csv)
    output_svg = resolve(args.output_svg)
    report_path = resolve(args.report)
    open_issues = resolve(args.open_issues)
    proof_png = resolve(args.proof_png)
    write_csv(review_csv, rows)
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    output_svg.write_text(render_svg(rows, ifc_hash), encoding="utf-8")
    render_png(output_svg, proof_png, find_chrome(args.chrome))

    unresolved_count = sum(
        bool(row[field])
        for row in rows
        for field in ("unresolved_material_or_product", "unresolved_manufacturer", "unresolved_construction")
    )
    open_issues.parent.mkdir(parents=True, exist_ok=True)
    open_issues.write_text(
        "# DET1 待复核节点\n\n"
        f"- 当前 IFC SHA-256：`{ifc_hash}`\n"
        f"- 节点数：{len(rows)}（D-601: 6；D-602: 6）\n"
        f"- 材料/产品、厂家、构造未决字段：{unresolved_count}/{len(rows) * 3}\n"
        "- 所有尺寸仅允许来自已确认证据；索引未填造产品尺寸。\n"
        "- `automatic_ifc_write_allowed=false`\n"
        "- `construction_release_ready=false`\n\n"
        "## 复核顺序\n\n"
        "1. 先选定材料/产品系统与厂家做法。\n"
        "2. 再按厂家资料和现场实测填写构造层、收口与尺寸变量。\n"
        "3. 机械检查防水/排水路径、平立剖一致性和检修可达性后，才能进入发布门。\n",
        encoding="utf-8",
    )

    input_hashes = {label: {"path": str(path), "sha256": sha256(path)} for label, path in paths.items()}
    output_hashes = {
        "review_csv": {"path": str(review_csv), "sha256": sha256(review_csv)},
        "svg": {"path": str(output_svg), "sha256": sha256(output_svg)},
        "proof_png": {"path": str(proof_png), "sha256": sha256(proof_png)},
        "open_issues": {"path": str(open_issues), "sha256": sha256(open_issues)},
    }
    report: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_det1_variable_detail_index_candidate",
        "source": {"ifc_sha256": ifc_hash, "inputs": input_hashes},
        "summary": {
            "node_count": len(rows),
            "d601_node_count": sum(row["sheet_id"] == "D-601" for row in rows),
            "d602_node_count": sum(row["sheet_id"] == "D-602" for row in rows),
            "unresolved_field_count": unresolved_count,
            "invented_dimension_count": 0,
        },
        "records": rows,
        "outputs": output_hashes,
        "gates": {
            "input_hashes_current": True,
            "unresolved_items_explicit": unresolved_count == len(rows) * 3,
            "svg_nonempty": output_svg.stat().st_size > 0,
            "png_nonempty": proof_png.stat().st_size > 0,
            "automatic_ifc_write_allowed": False,
            "construction_release_ready": False,
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(report_path), "summary": report["summary"], "gates": report["gates"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
