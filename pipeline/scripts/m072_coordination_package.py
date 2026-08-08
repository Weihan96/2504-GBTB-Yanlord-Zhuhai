#!/usr/bin/env python3
"""Generate and mechanically audit the M072 coordination candidate package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPECTED_SHEETS = (
    "A-001", "A-101", "A-102", "A-103", "A-104", "A-105",
    "P-201", "P-202", "E-301", "E-302", "E-303", "E-304",
    "A-106", "M-401", "I-501", "I-502", "I-503", "I-504",
    "D-601", "D-602", "S-701",
)
CANDIDATE_SHEETS = ("A-001", "A-102", "A-103", "A-104", "A-105")
SOURCE_REPORTS = {
    "A-102": ("build/a102/a102-wall-plan-candidate.json", ("gates", "pass")),
    "A-103": ("build/a103/a103-report.json", ("pass",)),
    "A-104": ("build/a104/a104-report.json", ("gates", "mechanical_pass")),
    "A-105": ("build/a105/a105-report.json", ("gates", "mechanical_pass")),
}
STATUS_LABELS = {"candidate": "候选完成", "planned": "后续计划"}
STATUS_COLORS = {"candidate": "#087f5b", "planned": "#868e96"}
ISSUE_STATUS_LABELS = {
    "pending": "待确认", "delegated": "已移交", "planned": "已计划",
    "implemented": "已实施", "confirmed": "已确认", "rejected": "已驳回",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--ifc", type=Path, default=Path("2504 GBTB Yanlord Zhuhai.ifc"))
    parser.add_argument("--drawing-register", type=Path, default=Path("pipeline/decisions/drawing-register.csv"))
    parser.add_argument("--issues", type=Path, default=Path("pipeline/decisions/m072-coordination-review.csv"))
    parser.add_argument("--qa-report", type=Path, default=Path("build/reports/qa-report.json"))
    parser.add_argument("--output-svg", type=Path, default=Path("drawings/A-001-drawing-index-candidate.svg"))
    parser.add_argument("--output-pdf", type=Path, default=Path("output/pdf/A-001-drawing-index-candidate.pdf"))
    parser.add_argument("--proof-png", type=Path, default=Path("tmp/pdfs/A-001-drawing-index-candidate.png"))
    parser.add_argument("--output-markdown", type=Path, default=Path("drawings/滨海湾M072协调候选包.md"))
    parser.add_argument("--report", type=Path, default=Path("build/m072/m072-report.json"))
    parser.add_argument("--render-script", type=Path, default=Path("pipeline/scripts/render_svg_pdf.py"))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def esc(value: Any) -> str:
    return html.escape(str(value), quote=True)


def build_index_svg(rows: list[dict[str, str]], issues: list[dict[str, str]], ifc_sha: str) -> str:
    status_counts = Counter(row["status"] for row in rows)
    issue_counts = Counter(row["status"] for row in issues)
    pieces = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="500mm" height="400mm" viewBox="0 0 500 400" data-scale="NTS">',
        '<rect width="500" height="400" fill="#fff"/>',
        '<rect x="6" y="6" width="488" height="388" fill="none" stroke="#102f43" stroke-width="0.6"/>',
        '<style>@page{size:500mm 400mm;margin:0}html,body{margin:0;width:500mm;height:400mm;overflow:hidden}text{font-family:Arial,"Noto Sans CJK SC",sans-serif;fill:#102f43}.title{font-size:8px;font-weight:700}.sub{font-size:3.2px;fill:#526777}.head{font-size:3.5px;font-weight:700}.cell{font-size:3.1px}.note{font-size:3px}.small{font-size:2.6px;fill:#526777}.rule{stroke:#ced4da;stroke-width:.35}.box{fill:#f8f9fa;stroke:#adb5bd;stroke-width:.35}</style>',
        '<text class="title" x="15" y="23">A-001 图纸目录／设计说明／图例</text>',
        '<text class="sub" x="15" y="31">仁恒滨海湾｜M072 72小时协调候选｜非施工发布版</text>',
        f'<text class="sub" x="485" y="31" text-anchor="end">IFC SHA {ifc_sha[:16]}…</text>',
    ]
    columns = (rows[:11], rows[11:])
    for column_index, column in enumerate(columns):
        x = 15 + column_index * 240
        pieces.extend([
            f'<rect class="box" x="{x}" y="40" width="230" height="10"/>',
            f'<text class="head" x="{x+3}" y="47">图号</text>',
            f'<text class="head" x="{x+34}" y="47">图名</text>',
            f'<text class="head" x="{x+185}" y="47">状态</text>',
        ])
        for index, row in enumerate(column):
            y = 50 + index * 11
            color = STATUS_COLORS.get(row["status"], "#c92a2a")
            pieces.extend([
                f'<line class="rule" x1="{x}" y1="{y+11}" x2="{x+230}" y2="{y+11}"/>',
                f'<text class="cell" x="{x+3}" y="{y+7}">{esc(row["sheet_number"])}</text>',
                f'<text class="cell" x="{x+34}" y="{y+7}">{esc(row["title"])}</text>',
                f'<rect x="{x+185}" y="{y+2.5}" width="38" height="6" rx="2" fill="{color}" opacity=".14"/>',
                f'<text class="small" x="{x+204}" y="{y+7}" text-anchor="middle" fill="{color}">{esc(STATUS_LABELS.get(row["status"], row["status"]))}</text>',
            ])
    pieces.extend([
        '<rect class="box" x="15" y="185" width="230" height="112"/>',
        '<text class="head" x="21" y="196">统一设计与出图说明</text>',
        '<text class="note" x="21" y="207">1. IFC 及平面定位尺寸单位为 mm；标高以 FFL = ±0.000 为参考。</text>',
        '<text class="note" x="21" y="217">2. 装修定位默认控制完成面；结构尺寸须另行明示。</text>',
        '<text class="note" x="21" y="227">3. 带 ≈ 尺寸为候选或近似值，不得代替现场实测放线。</text>',
        '<text class="note" x="21" y="237">4. 材料封样、设备型号/尺寸/功率/接口未确认前不得下单。</text>',
        '<text class="note" x="21" y="247">5. 物业、燃气、门窗与机电限制须在封闭施工或设备下单前关闭。</text>',
        '<text class="note" x="21" y="257">6. 本包仅关闭 P0 候选协调；P1/P2 图纸已列入索引但尚未发布。</text>',
        '<text class="head" x="21" y="274">图例</text>',
        '<rect x="21" y="281" width="7" height="7" fill="#087f5b" opacity=".25"/><text class="note" x="32" y="287">候选完成</text>',
        '<rect x="83" y="281" width="7" height="7" fill="#868e96" opacity=".25"/><text class="note" x="94" y="287">后续计划</text>',
        '<rect x="145" y="281" width="7" height="7" fill="#e67700" opacity=".25"/><text class="note" x="156" y="287">现场/业主确认</text>',
        '<rect class="box" x="255" y="185" width="230" height="112"/>',
        '<text class="head" x="261" y="196">M072 协调摘要</text>',
        f'<text class="note" x="261" y="208">计划图纸 {len(rows)} 张｜候选完成 {status_counts["candidate"]} 张｜后续计划 {status_counts["planned"]} 张</text>',
        f'<text class="note" x="261" y="218">协调事项 {len(issues)} 项｜待确认 {issue_counts["pending"]} 项｜已移交 {issue_counts["delegated"]} 项</text>',
        '<text class="note" x="261" y="230">P0 候选：A-102 / A-103 / A-104 / A-105 已从同一 IFC 生成。</text>',
        '<text class="note" x="261" y="240">机械检查：图号唯一、PDF 页数/幅面、空白页、标注碰撞。</text>',
        '<text class="note" x="261" y="250">未决项：Space Reference，无宿主门，现场尺寸及设备/物业条件。</text>',
        '<text class="note" x="261" y="260">上述未决项已显式分派，不表示已关闭。</text>',
        '<text class="note" x="261" y="276">人读清单：drawings/滨海湾M072协调候选包.md</text>',
        '<text class="small" x="15" y="384">M072-CANDIDATE｜2026-08-08｜源：正式 IFC + 受版本控制决策表</text>',
        '<text class="small" x="485" y="384" text-anchor="end">A-001｜NTS｜1 / 1</text>',
        '</svg>',
    ])
    return "".join(pieces)


def nested(data: dict[str, Any], keys: tuple[str, ...]) -> Any:
    value: Any = data
    for key in keys:
        value = value[key]
    return value


def pdf_info(path: Path) -> dict[str, Any]:
    result = subprocess.run(["pdfinfo", str(path)], check=True, capture_output=True, text=True)
    pages_match = re.search(r"^Pages:\s+(\d+)$", result.stdout, re.MULTILINE)
    size_match = re.search(r"^Page size:\s+([\d.]+) x ([\d.]+) pts", result.stdout, re.MULTILINE)
    if not pages_match or not size_match:
        raise RuntimeError(f"unable to parse pdfinfo for {path}")
    width_mm = float(size_match.group(1)) * 25.4 / 72.0
    height_mm = float(size_match.group(2)) * 25.4 / 72.0
    return {"pages": int(pages_match.group(1)), "width_mm": width_mm, "height_mm": height_mm}


def parse_pgm(path: Path) -> tuple[int, int, bytes]:
    data = path.read_bytes()
    tokens: list[bytes] = []
    index = 0
    while len(tokens) < 4:
        while index < len(data) and data[index:index+1].isspace():
            index += 1
        if data[index:index+1] == b"#":
            index = data.index(b"\n", index) + 1
            continue
        end = index
        while end < len(data) and not data[end:end+1].isspace():
            end += 1
        tokens.append(data[index:end])
        index = end
    while index < len(data) and data[index:index+1].isspace():
        index += 1
    if tokens[0] != b"P5" or tokens[3] != b"255":
        raise RuntimeError(f"unsupported PGM header in {path}")
    width, height = int(tokens[1]), int(tokens[2])
    pixels = data[index:index + width * height]
    if len(pixels) != width * height:
        raise RuntimeError(f"truncated PGM pixels in {path}")
    return width, height, pixels


def raster_metrics(pdf: Path) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="m072-pdf-") as directory:
        prefix = Path(directory) / "page"
        subprocess.run(["pdftoppm", "-f", "1", "-singlefile", "-r", "20", "-gray", str(pdf), str(prefix)], check=True, capture_output=True)
        width, height, pixels = parse_pgm(prefix.with_suffix(".pgm"))
    ink = [index for index, value in enumerate(pixels) if value < 245]
    if not ink:
        return {"ink_ratio": 0.0, "ink_bbox_page_ratio": 0.0}
    xs = [index % width for index in ink]
    ys = [index // width for index in ink]
    bbox_ratio = ((max(xs) - min(xs) + 1) * (max(ys) - min(ys) + 1)) / (width * height)
    return {"ink_ratio": len(ink) / len(pixels), "ink_bbox_page_ratio": bbox_ratio}


def markdown_table(rows: list[list[str]]) -> list[str]:
    if not rows:
        return []
    rendered = ["| " + " | ".join(rows[0]) + " |"]
    rendered.append("| " + " | ".join("---" for _ in rows[0]) + " |")
    rendered.extend("| " + " | ".join(row) + " |" for row in rows[1:])
    return rendered


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    resolve = lambda path: path if path.is_absolute() else root / path
    ifc = resolve(args.ifc)
    drawing_register = resolve(args.drawing_register)
    issues_path = resolve(args.issues)
    qa_path = resolve(args.qa_report)
    output_svg = resolve(args.output_svg)
    output_pdf = resolve(args.output_pdf)
    proof_png = resolve(args.proof_png)
    output_markdown = resolve(args.output_markdown)
    report_path = resolve(args.report)
    render_script = resolve(args.render_script)

    rows = read_csv(drawing_register)
    issues = read_csv(issues_path)
    sheet_numbers = [row["sheet_number"] for row in rows]
    if tuple(sheet_numbers) != EXPECTED_SHEETS or len(set(sheet_numbers)) != len(sheet_numbers):
        raise RuntimeError("drawing register is incomplete, out of order, or contains duplicate sheet numbers")
    if any(row["status"] not in STATUS_LABELS for row in rows):
        raise RuntimeError("drawing register contains an unsupported status")
    if any(row["status"] not in ISSUE_STATUS_LABELS or not row["basis"] or not row["stop_condition"] for row in issues):
        raise RuntimeError("coordination issue rows require supported status, basis, and stop condition")

    ifc_sha = sha256(ifc)
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    output_svg.write_text(build_index_svg(rows, issues, ifc_sha), encoding="utf-8")
    subprocess.run([
        sys.executable, str(render_script), "--input-svg", str(output_svg), "--output-pdf", str(output_pdf),
        "--proof-png", str(proof_png), "--report", str(root / "build/m072/a001-pdf-render-report.json"),
        "--width-mm", "500", "--height-mm", "400",
    ], check=True, cwd=root)

    source_reports: dict[str, Any] = {}
    for sheet, (relative_path, gate_path) in SOURCE_REPORTS.items():
        data = json.loads((root / relative_path).read_text(encoding="utf-8"))
        source_reports[sheet] = {
            "path": relative_path,
            "source_sha256": data["source"]["ifc_sha256"],
            "same_formal_ifc": data["source"]["ifc_sha256"] == ifc_sha,
            "mechanical_pass": bool(nested(data, gate_path)),
        }

    pdf_records: list[dict[str, Any]] = []
    for row in rows:
        if row["status"] != "candidate":
            continue
        pdf = root / row["publish_target"]
        info = pdf_info(pdf)
        raster = raster_metrics(pdf)
        pdf_records.append({
            "sheet_number": row["sheet_number"], "path": row["publish_target"], "sha256": sha256(pdf),
            **info, **raster,
            "passes": info["pages"] == 1 and abs(info["width_mm"] - 500.0) <= 0.2
            and abs(info["height_mm"] - 400.0) <= 0.2 and raster["ink_ratio"] >= 0.002
            and raster["ink_bbox_page_ratio"] >= 0.35,
        })

    qa = json.loads(qa_path.read_text(encoding="utf-8"))
    blocks = {gate["id"] for gate in qa["gates"] if gate["status"] == "block"}
    disclosed_gates = {
        gate.strip() for issue in issues for gate in issue["qa_gate"].split(";") if gate.strip()
    }
    gates = {
        "ifc_sha256": ifc_sha,
        "drawing_register_count": len(rows),
        "drawing_register_unique": len(set(sheet_numbers)) == len(rows),
        "candidate_sheet_count": sum(row["status"] == "candidate" for row in rows),
        "planned_sheet_count": sum(row["status"] == "planned" for row in rows),
        "source_report_same_ifc_count": sum(item["same_formal_ifc"] for item in source_reports.values()),
        "source_report_mechanical_pass_count": sum(item["mechanical_pass"] for item in source_reports.values()),
        "candidate_pdf_pass_count": sum(item["passes"] for item in pdf_records),
        "candidate_pdf_count": len(pdf_records),
        "legacy_release_qa": qa["summary"],
        "legacy_block_count": len(blocks),
        "legacy_blocks_disclosed": blocks <= disclosed_gates,
        "coordination_issue_count": len(issues),
        "coordination_pending_count": sum(issue["status"] == "pending" for issue in issues),
    }
    gates["m072_mechanical_pass"] = (
        gates["drawing_register_count"] == len(EXPECTED_SHEETS)
        and gates["drawing_register_unique"]
        and gates["candidate_sheet_count"] == len(CANDIDATE_SHEETS)
        and gates["source_report_same_ifc_count"] == 4
        and gates["source_report_mechanical_pass_count"] == 4
        and gates["candidate_pdf_pass_count"] == gates["candidate_pdf_count"] == 5
        and gates["legacy_blocks_disclosed"]
    )

    candidate_rows = [["图号", "图名", "PDF", "状态"]] + [
        [row["sheet_number"], row["title"], row["publish_target"], "候选完成"]
        for row in rows if row["status"] == "candidate"
    ]
    issue_rows = [["ID", "优先级", "问题/当前状态", "责任方", "目标任务", "状态"]] + [
        [issue["issue_id"], issue["priority"], issue["current_condition"], issue["responsible_party"], issue["target_task"], ISSUE_STATUS_LABELS[issue["status"]]]
        for issue in issues if issue["status"] != "implemented"
    ]
    lines = [
        "# 滨海湾 M072 72小时协调候选包", "",
        f"- 正式 IFC SHA-256：`{ifc_sha}`",
        f"- 计划图纸：{len(rows)} 张；已完成候选：{gates['candidate_sheet_count']} 张；后续计划：{gates['planned_sheet_count']} 张",
        f"- M072 机械门：{'PASS' if gates['m072_mechanical_pass'] else 'FAIL'}",
        f"- Release QA：{qa['summary']['pass']} PASS / {qa['summary']['warn']} WARN / {qa['summary']['block']} BLOCK；本包不是施工发布版",
        "", "## 当前候选图纸", "",
        *markdown_table(candidate_rows), "", "## 机械验收", "",
        "| 检查 | 结果 |", "| --- | --- |",
        f"| A-102～A-105 生成报告与正式 IFC 同一 SHA | {gates['source_report_same_ifc_count']}/4 PASS |",
        f"| A-102～A-105 编号/尺寸闭合/标注碰撞机械门 | {gates['source_report_mechanical_pass_count']}/4 PASS |",
        f"| A-001～A-105 候选 PDF 页数、500×400 mm 幅面与非空白检查 | {gates['candidate_pdf_pass_count']}/5 PASS |",
        f"| Release QA BLOCK 是否全部显式分派 | {'PASS' if gates['legacy_blocks_disclosed'] else 'FAIL'} |",
        "", "## 协调问题与现场复核", "",
        "下表面向人阅读；机器状态源为 `pipeline/decisions/m072-coordination-review.csv`。“已移交/已计划”不等于“已解决”。", "",
        *markdown_table(issue_rows), "", "## 停止条件", "",
        "- 任一 P0 候选图不能从当前 IFC 重现。",
        "- 图号重复、PDF 空白/损坏/幅面错误，或标注碰撞重现。",
        "- Release QA 的 BLOCK 没有明确责任方、目标任务和停止条件。",
        "- 未经现场复核的候选尺寸被当作施工放线值。", "",
    ]
    output_markdown.parent.mkdir(parents=True, exist_ok=True)
    output_markdown.write_text("\n".join(lines), encoding="utf-8")
    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(), "mode": "m072-coordination-candidate",
        "source": {"ifc": str(ifc), "sha256": ifc_sha}, "source_reports": source_reports,
        "drawing_register": rows, "pdf_records": pdf_records, "issues": issues, "gates": gates,
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"svg": str(output_svg), "pdf": str(output_pdf), "markdown": str(output_markdown), "report": str(report_path), **gates}, ensure_ascii=False))
    if not gates["m072_mechanical_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
