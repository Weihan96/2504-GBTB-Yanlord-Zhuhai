#!/usr/bin/env python3
"""Mechanical guardrails for the current external confirmation forms."""

from __future__ import annotations

import csv
import hashlib
import re
import sys
from pathlib import Path
from urllib.parse import unquote, urlparse


ROOT = Path(__file__).resolve().parents[2]
FORM_DIR = ROOT / "output" / "forms" / "对外确认表"
CURRENT = [
    "01-门组做法确认表-发全屋定制.md",
    "02-日立空调接口确认表-发空调厂家.md",
    "03-给排水设备确认表-发设备与施工方.md",
    "04-燃气消防确认表-发主管单位.md",
    "05-弱电现场记录表-发现场负责人.md",
    "06-门组厂家复核表-发Rimadesio与Poliform.md",
    "07-智能面板电气接口确认表-发JINK与电气方.md",
    "08-定制家具与墙脚节点确认表-发全屋定制.md",
    "09-厨房设备与柜体深化确认表-发橱柜设备方.md",
    "10-强电箱容量现场复核表-发电气方.md",
    "11-窗改造窗饰与阳台完成面确认表-发施工方.md",
]
FORBIDDEN_FORM_JARGON = {
    "confirmed owner decision",
    "official evidence",
    "official/research conclusion",
    "research conclusion",
    "external signoff pending",
    "site measurement pending",
    "authority signoff pending",
    "project internal pending",
    "owner preference pending",
    "owner_reference_candidate_not_final",
    "unknown",
}
SOURCE_IDS = {
    "01-门组做法确认表-发全屋定制.md": "OUTBOUND-FORM-DOOR-20260817",
    "02-日立空调接口确认表-发空调厂家.md": "OUTBOUND-FORM-HVAC-20260817",
    "03-给排水设备确认表-发设备与施工方.md": "OUTBOUND-FORM-PLUM-20260817",
    "04-燃气消防确认表-发主管单位.md": "OUTBOUND-FORM-GASFIRE-20260817",
    "05-弱电现场记录表-发现场负责人.md": "OUTBOUND-FORM-SITE-20260817",
    "06-门组厂家复核表-发Rimadesio与Poliform.md": "OUTBOUND-FORM-DOOR-VENDOR-20260817",
    "07-智能面板电气接口确认表-发JINK与电气方.md": "OUTBOUND-FORM-SMART-PANEL-20260817",
    "08-定制家具与墙脚节点确认表-发全屋定制.md": "OUTBOUND-FORM-JOINERY-DETAIL-20260817",
    "09-厨房设备与柜体深化确认表-发橱柜设备方.md": "OUTBOUND-FORM-KITCHEN-20260818",
    "10-强电箱容量现场复核表-发电气方.md": "OUTBOUND-FORM-ELEC-CAPACITY-20260821",
    "11-窗改造窗饰与阳台完成面确认表-发施工方.md": "OUTBOUND-FORM-WINDOWS-20260821",
}
HISTORICAL_EXTERNAL_SOURCE_IDS = {
    "A104-SHOP-REQUEST-20260815",
    "RCP1-PLUM-REQUEST-20260815",
    "GAS-CONSULTATION-TEMPLATE-001",
    "A106-CONSULTATION-PACK-20260815",
    "E304-SITE-CHECKLIST-20260815",
}
MARKDOWN_LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
UNLINKED_ATTACHMENT = re.compile(
    r"`[^`]+\.(?:pdf|md|dwg|dxf|svg|jpe?g|png|csv|xlsx|json|webp)`",
    re.IGNORECASE,
)


def read_csv(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle)}


def cells(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def validate_tables(name: str, text: str, errors: list[str]) -> None:
    action_headers = (
        "现在请你做什么",
        "请确认的唯一问题",
        "唯一复核事项",
        "要确认的结果",
        "项目图中必须先画清的方案",
        "项目图中必须先完成",
        "项目方案和样板要求",
        "项目已经给出的方向",
        "项目方案和材料要求",
    )
    responsibility_headers = ("由谁回答", "谁回复", "谁复核")
    response_headers = ("请在这里回复", "直接填写", "外部单位只需填写")
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if not line.startswith("|"):
            continue
        header = cells(line)
        action_header = next((value for value in action_headers if value in header), None)
        if action_header is None:
            continue
        action_index = header.index(action_header)
        responsibility_header = next((value for value in responsibility_headers if value in header), None)
        response_header = next((value for value in response_headers if value in header), None)
        responsibility_index = header.index(responsibility_header) if responsibility_header else None
        response_index = header.index(response_header) if response_header else None
        if responsibility_index is None:
            errors.append(f"{name}:{index + 1}: action table is missing an explicit responsible-party column")
        if response_index is None:
            errors.append(f"{name}:{index + 1}: action table is missing a reply column")
        row_index = index + 2
        while row_index < len(lines) and lines[row_index].startswith("|"):
            row = cells(lines[row_index])
            if len(row) != len(header):
                errors.append(f"{name}:{row_index + 1}: table column count differs from header")
                row_index += 1
                continue
            if not row[action_index]:
                errors.append(f"{name}:{row_index + 1}: missing remaining action")
            if responsibility_index is not None and not row[responsibility_index]:
                errors.append(f"{name}:{row_index + 1}: missing responsible party")
            if response_index is not None and not row[response_index]:
                errors.append(f"{name}:{row_index + 1}: missing reply field")
            row_index += 1


def require(text: str, name: str, tokens: list[str], errors: list[str]) -> None:
    for token in tokens:
        if token not in text:
            errors.append(f"{name}: missing regression token {token!r}")


def validate_links(path: Path, text: str, errors: list[str]) -> None:
    """Require openable, portable links for every attachment named in a form."""
    links = MARKDOWN_LINK.findall(text)
    if not links:
        errors.append(f"{path.name}: contains no clickable attachment or evidence links")
    for raw_target in links:
        target = raw_target.strip().split(" ", 1)[0].strip("<>")
        parsed = urlparse(target)
        if parsed.scheme in {"http", "https", "mailto"} or target.startswith("#"):
            continue
        if parsed.scheme:
            errors.append(f"{path.name}: unsupported link scheme in {target!r}")
            continue
        decoded = unquote(target)
        local_target = Path(decoded)
        if local_target.is_absolute():
            errors.append(f"{path.name}: local link must be project-relative: {target!r}")
            continue
        resolved = (path.parent / local_target).resolve()
        try:
            resolved.relative_to(ROOT.resolve())
        except ValueError:
            errors.append(f"{path.name}: local link leaves the project: {target!r}")
            continue
        if not resolved.exists():
            errors.append(f"{path.name}: local link target does not exist: {target!r}")
    for match in UNLINKED_ATTACHMENT.findall(text):
        errors.append(f"{path.name}: attachment is named without a clickable link: {match}")


def main() -> int:
    errors: list[str] = []
    actual = sorted(path.name for path in FORM_DIR.glob("*.md"))
    if actual != sorted(CURRENT):
        errors.append(f"current Markdown set differs: expected {CURRENT}, found {actual}")
    if (ROOT / "drawings" / "evidence" / "EXT-外部信息最短清单-20260817.md").exists():
        errors.append("deleted external-information index has reappeared")

    texts: dict[str, str] = {}
    for name in CURRENT:
        path = FORM_DIR / name
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        texts[name] = text
        if "待填写" in text:
            errors.append(f"{name}: contains contextless 待填写")
        if "外部信息最短清单" in text:
            errors.append(f"{name}: references deleted duplicate dispatch index")
        if "当前状态" in text:
            errors.append(f"{name}: exposes the internal status column to external readers")
        for token in FORBIDDEN_FORM_JARGON:
            if token in text:
                errors.append(f"{name}: exposes internal status jargon {token!r}")
        if name != CURRENT[0] and "发送给" not in text:
            errors.append(f"{name}: does not state who receives the form")
        validate_links(path, text, errors)
        validate_tables(name, text, errors)

    require(texts.get(CURRENT[0], ""), CURRENT[0], ["业主已经确认，不用填写", "向东／图纸右侧滑开", "不要求“两道水平横档”"], errors)
    require(texts.get(CURRENT[1], ""), CURRENT[1], ["业主不用填写", "华美", "RPIZ-71FSLN5QD/P", "A06", "QD", "QDF", "冷凝水路线", "支吊点标高", "一只面板最多控制 6 台室内机", "回风口同时作为检修口"], errors)
    require(texts.get(CURRENT[2], ""), CURRENT[2], ["两个 16A 插座或两个回路", "双存水弯组合", "四分铝塑管套二分 PE 管的网络固定配方", "WTZ27510", "APP-014", "MF287", "CZ356", "CZ028", "Street", "Sorgente", "Nuna V Combo", "定制地漏"], errors)
    require(texts.get(CURRENT[3], ""), CURRENT[3], ["ER9EPA33MP", "不代表已经批准或购买", "厨房探测器的位置已经确定", "甲烷报警器＋紧急切断阀", "现有消防主机"], errors)
    require(texts.get(CURRENT[4], ""), CURRENT[4], ["现有约测净深 110 mm", "525 mm", "H+350 mm", "(4600.016, -735.369) mm", "玄关高柜右下柜格", "五孔插座", "入户临时置物位", "TL-XAP1500GE-PoE/DC", "280M-S3", "Entry A 已由业主选定", "两个可同时插入的 220 V 插座位", "PC-P1HEQ", "1～8 芯"], errors)
    require(texts.get(CURRENT[5], ""), CURRENT[5], ["可以 □　不可以 □", "本项目下单／安装图", "地面不做通长下轨"], errors)
    require(texts.get(CURRENT[6], ""), CURRENT[6], ["4 个基础功能", "4 个：客厅、书房、餐厅、照明总关", "传统有线双控", "主卧氛围 LED、主卧重点射灯"], errors)
    require(texts.get(CURRENT[7], ""), CURRENT[7], ["350～450 mm", "上下约 10 mm 双阴影缝", "全高镜柜", "材质化正投影", "谁回复"], errors)
    require(texts.get(CURRENT[8], ""), CURRENT[8], ["LS33R6VB9W/01", "895×345×873 mm", "不再让业主重复选择", "灶台后和主要操作墙使用与台面同材", "NS-01", "三个外露面各设一个", "两路独立 C16", "H70FT", "K-D01", "600 mm 石材盆", "F50 仍是基准候选"], errors)
    require(texts.get(CURRENT[9], ""), CURRENT[9], ["P230", "PX1", "63 A", "总开关", "剩余模数", "新增负荷计算"], errors)
    require(texts.get(CURRENT[10], ""), CURRENT[10], ["机械手摇", "22 mm", "K5-501", "N31-204", "LightLock", "Duolite", "阳台塑木"], errors)

    sources = read_csv(ROOT / "pipeline" / "decisions" / "source-evidence-register.csv", "source_id")
    owner_inputs = read_csv(ROOT / "pipeline" / "decisions" / "owner-input-register.csv", "input_id")
    equipment = read_csv(ROOT / "pipeline" / "decisions" / "equipment-register.csv", "equipment_id")
    requirements = read_csv(ROOT / "pipeline" / "decisions" / "equipment-installation-requirements.csv", "requirement_id")
    for name, source_id in SOURCE_IDS.items():
        if source_id not in sources:
            errors.append(f"source register missing {source_id}")
            continue
        expected = hashlib.sha256((FORM_DIR / name).read_bytes()).hexdigest()
        if sources[source_id]["sha256"] != expected:
            errors.append(f"{source_id}: stored hash differs from current Markdown")
        if "现行唯一可编辑源" not in sources[source_id]["notes"]:
            errors.append(f"{source_id}: not marked as current unique editable source")
    for source_id in HISTORICAL_EXTERNAL_SOURCE_IDS:
        row = sources.get(source_id)
        if row is None:
            errors.append(f"historical source register missing {source_id}")
            continue
        if row["status"] != "superseded_historical_evidence_current_form_linked":
            errors.append(f"{source_id}: historical request is not retired from current filling entry")
        if "现行填写入口" not in row["notes"]:
            errors.append(f"{source_id}: historical request does not point to its current form")

    ssot_checks = {
        "E304-CABINET-DIMENSIONS": ["110", "525", "H+350", "4600.016", "只缺准确净宽"],
        "RCP1-HVAC-PORTS": ["华美", "A05", "旧 IFC 紫色管线"],
        "INT1-ENTRY-PARCEL-FUNCTION": ["临时置物", "不是专用快递柜"],
        "INT1-ENTRY-PARCEL-LAYOUT": ["I-503 前", "不要求门开／门关"],
        "APP014-REPLACEABLE-SLEEVE": ["连续可抽换", "网络案例规格不作为项目规格"],
    }
    for item_id, tokens in ssot_checks.items():
        row = owner_inputs.get(item_id)
        if row is None:
            errors.append(f"owner input missing {item_id}")
            continue
        joined = " ".join(row.values())
        for token in tokens:
            if token not in joined:
                errors.append(f"{item_id}: missing SSOT token {token!r}")

    if "at least 4" not in equipment["CTRL-ENTRY-A"]["variant"] and "minimum 4" not in equipment["CTRL-ENTRY-A"]["variant"]:
        errors.append("CTRL-ENTRY-A equipment row does not preserve minimum four functions")
    if "3-Key/3-Relay（研究优选" in equipment["CTRL-ENTRY-A"]["model"]:
        errors.append("CTRL-ENTRY-A still promotes incompatible 3-key candidate")
    for req_id in [f"REQ-APP014-SLEEVE-20260818-00{index}" for index in range(1, 6)]:
        if req_id not in requirements:
            errors.append(f"installation requirements missing {req_id}")
    if "EXT-EXTERNAL-INFO-MINIMUM-20260817" in sources:
        errors.append("deleted duplicate dispatch index remains in source register")

    history = (ROOT / "drawings" / "evidence" / "E304-现场最短取证清单-20260815.md").read_text(encoding="utf-8")
    if "不是当前填写入口" not in history or "当前唯一可编辑入口" not in history:
        errors.append("historical E304 checklist is not clearly retired")

    if errors:
        print("external confirmation form validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 2
    print(f"external confirmation form validation passed: {len(CURRENT)} current forms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
