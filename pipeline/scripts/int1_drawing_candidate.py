#!/usr/bin/env python3
"""Generate readable I-501 through I-504 existing-object candidate sheets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from typing import Any


SHEETS = {
    "I-501": {
        "title": "Kitchen Coordination Candidate",
        "subtitle": "厨房平面与立面协调候选",
        "underlay": "Furniture Plan.svg",
        "pdf": "I-501-kitchen-existing-candidate.pdf",
        "svg": "I-501-kitchen-existing-candidate.svg",
    },
    "I-502": {
        "title": "Bathroom Coordination Candidate",
        "subtitle": "卫生间平面、立面与节点协调候选",
        "underlay": "Sanitary Plan.svg",
        "pdf": "I-502-bathroom-existing-candidate.pdf",
        "svg": "I-502-bathroom-existing-candidate.svg",
    },
    "I-503": {
        "title": "Entry / Laundry Candidate",
        "subtitle": "入户与家政既有对象候选",
        "underlay": "Furniture Plan.svg",
        "pdf": "I-503-entry-laundry-existing-candidate.pdf",
        "svg": "I-503-entry-laundry-existing-candidate.svg",
    },
    "I-504": {
        "title": "Fixed Furniture Candidate",
        "subtitle": "其他固定家具平面与立面候选",
        "underlay": "Furniture Plan.svg",
        "pdf": "I-504-other-fixed-furniture-existing-candidate.pdf",
        "svg": "I-504-other-fixed-furniture-existing-candidate.svg",
    },
}

ROLE_STYLE = {
    "fixed_furniture": ("#e0ad2f", "#8a6300", "Fixed furniture"),
    "loose_furniture": ("#d8d8d8", "#6c6c6c", "Loose context"),
    "fixed_equipment": ("#50b7ca", "#146777", "Fixed equipment"),
    "fixed_sanitary_fixture": ("#5b8fd8", "#234e8b", "Sanitary fixture"),
    "fixed_waste_terminal": ("#8a66bf", "#4b2c78", "Waste terminal"),
    "fixed_context": ("#9aa5ad", "#4c5961", "Wall / slab / assembly context"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_candidate(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, str]], str]:
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8-sig")))
    if not rows:
        raise RuntimeError("INT1 candidate register is empty")
    source_hashes = {row["source_ifc_sha256"] for row in rows}
    if len(source_hashes) != 1:
        raise RuntimeError(f"INT1 register has multiple IFC hashes: {source_hashes}")
    records: list[dict[str, Any]] = []
    blockers: list[dict[str, str]] = []
    for row in rows:
        if row["record_kind"] == "blocker":
            if row["review_status"] != "BLOCK":
                raise RuntimeError(f"INT1 blocker is not BLOCK: {row['object_name']}")
            blockers.append(row)
            continue
        if row["record_kind"] != "existing_object":
            raise RuntimeError(f"unknown INT1 record kind: {row['record_kind']}")
        if row["dimension_status"] != "existing_world_bbox_not_fabrication_dimension":
            raise RuntimeError(f"unsafe INT1 dimension status: {row['global_id']}")
        record: dict[str, Any] = dict(row)
        for key in ("bbox_min_mm", "bbox_max_mm", "dimensions_mm"):
            record[key] = json.loads(row[key])
            if len(record[key]) != 3:
                raise RuntimeError(f"invalid bbox field {key}: {row['global_id']}")
        records.append(record)
    return records, blockers, next(iter(source_hashes))


def read_dishwasher_interfaces(path: Path, source_hash: str) -> list[dict[str, Any]]:
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("source", {}).get("ifc_sha256") != source_hash:
        raise RuntimeError("INT1 product-interface report does not match the formal IFC")
    by_id = {
        row["equipment_id"]: row
        for row in report.get("kitchen_product_installation_requirements", [])
    }
    expected_ids = ["APP-009", "APP-010"]
    if not all(equipment_id in by_id for equipment_id in expected_ids):
        raise RuntimeError("INT1 report is missing APP-009/APP-010 interfaces")
    result = []
    for equipment_id in expected_ids:
        row = by_id[equipment_id]
        requirements = {item["parameter_key"]: item for item in row["requirements"]}
        required_keys = {
            "niche_height_min", "niche_height_max", "niche_width_min",
            "niche_width_max", "niche_depth_min", "water_connection",
            "drain_connection_od",
        }
        if not required_keys.issubset(requirements):
            raise RuntimeError(f"{equipment_id}: incomplete official interface constraints")
        for key in required_keys:
            requirement = requirements[key]
            if (
                requirement["status"] != "confirmed"
                or requirement["value_origin"] != "official_exact_model"
                or requirement["source_id"] != "APP-DW-INSTALL-001"
            ):
                raise RuntimeError(
                    f"{equipment_id}: {key} is not confirmed exact-model evidence"
                )
        if row["sheet_id"] != "I-501" or row["use_location"] != "西厨岛台":
            raise RuntimeError(f"{equipment_id}: unexpected drawing or use location")
        if row["project_interface_status"] != "unlocated":
            raise RuntimeError(f"{equipment_id}: project interface must remain unlocated")
        if row["procurement_status"] != "candidate" or row["final_product_confirmed"]:
            raise RuntimeError(f"{equipment_id}: dishwasher must remain an unpurchased candidate")
        result.append(row)
    shared_keys = (
        "niche_height_min", "niche_height_max", "niche_width_min",
        "niche_width_max", "niche_depth_min",
    )
    signatures = {
        tuple(
            (
                next(item for item in row["requirements"] if item["parameter_key"] == key)["value"],
                next(item for item in row["requirements"] if item["parameter_key"] == key)["unit"],
            )
            for key in shared_keys
        )
        for row in result
    }
    if len(signatures) != 1:
        raise RuntimeError("APP-009/APP-010 official niche constraints do not match")
    return result


def interface_requirement(interface: dict[str, Any], key: str) -> dict[str, str]:
    return next(
        item for item in interface["requirements"] if item["parameter_key"] == key
    )


def display_requirement(interface: dict[str, Any], key: str) -> str:
    requirement = interface_requirement(interface, key)
    return " ".join(
        part for part in (str(requirement["value"]), requirement["unit"]) if part
    )


def plan_rect(record: dict[str, Any]) -> tuple[float, float, float, float]:
    minimum = record["bbox_min_mm"]
    maximum = record["bbox_max_mm"]
    # Existing Bonsai plans use 1:50 around the model origin at SVG (200, 200).
    x = 200.0 + minimum[0] / 50.0
    y = 200.0 - maximum[1] / 50.0
    width = max((maximum[0] - minimum[0]) / 50.0, 0.35)
    height = max((maximum[1] - minimum[1]) / 50.0, 0.35)
    return x, y, width, height


def grouped_index(records: list[dict[str, Any]]) -> list[tuple[str, int, str]]:
    groups: Counter[tuple[str, str]] = Counter()
    for record in records:
        label = record["type_name"] or record["object_name"] or record["ifc_class"]
        groups[(record["installation_role"], label)] += 1
    return [
        (label, count, role)
        for (role, label), count in sorted(
            groups.items(), key=lambda item: (item[0][0], item[0][1], item[1])
        )
    ]


def svg_text(x: float, y: float, text: str, css_class: str = "index") -> str:
    return f'<text x="{x:.3f}" y="{y:.3f}" class="{css_class}">{escape(text)}</text>'


def make_svg(
    sheet_id: str,
    records: list[dict[str, Any]],
    blockers: list[dict[str, str]],
    source_hash: str,
    underlay_hash: str,
    dishwasher_interfaces: list[dict[str, Any]],
) -> str:
    spec = SHEETS[sheet_id]
    overlays = []
    for record in records:
        x, y, width, height = plan_rect(record)
        fill, stroke, _ = ROLE_STYLE[record["installation_role"]]
        overlays.append(
            f'<rect x="{x:.4f}" y="{y:.4f}" width="{width:.4f}" height="{height:.4f}" '
            f'fill="{fill}" fill-opacity="0.32" stroke="{stroke}" stroke-width="0.55" '
            f'data-global-id="{escape(record["global_id"])}" data-dimension-status="coordination-envelope-only"/>'
        )

    index_lines: list[str] = []
    line_y = 58.0
    for label, count, role in grouped_index(records):
        _, stroke, _ = ROLE_STYLE[role]
        index_lines.append(
            f'<rect x="383" y="{line_y - 2.4:.3f}" width="2.2" height="2.2" fill="{stroke}"/>'
        )
        display = f"{label}  x{count}"
        if len(display) > 42:
            display = display[:39] + "..."
        index_lines.append(svg_text(387.0, line_y, display))
        line_y += 4.15

    block_lines: list[str] = []
    block_y = 342.0
    for blocker in blockers:
        block_lines.append(
            svg_text(382.0, block_y, f"BLOCK {blocker['object_name']}: {blocker['type_name']}", "block")
        )
        block_y += 5.0

    interface_lines: list[str] = []
    if sheet_id == "I-501":
        interface_lines.append(svg_text(382.0, 286.0, "PRODUCT INTERFACES · NO PROJECT XYZ", "subtitle"))
        row_y = 293.0
        for interface in dishwasher_interfaces:
            model = interface["model"].replace("Siemens ", "")
            use_location = interface["use_location"]
            water = display_requirement(interface, "water_connection")
            drain = display_requirement(interface, "drain_connection_od")
            interface_lines.append(
                f'<text x="382" y="{row_y:.3f}" class="interface" '
                f'data-equipment-id="{escape(interface["equipment_id"])}" '
                f'data-project-interface-status="unlocated">'
                f'{escape(interface["equipment_id"])} · {escape(model)} · CANDIDATE / NOT PURCHASED</text>'
            )
            interface_lines.append(
                svg_text(
                    382.0,
                    row_y + 4.0,
                    f"Use: {use_location} · {water} · drain Ø{drain}",
                    "interface",
                )
            )
            row_y += 10.0
        shared = dishwasher_interfaces[0]
        niche_height_min = interface_requirement(shared, "niche_height_min")["value"]
        niche_height_max = display_requirement(shared, "niche_height_max")
        niche_width_min = interface_requirement(shared, "niche_width_min")["value"]
        niche_width_max = display_requirement(shared, "niche_width_max")
        niche_depth_min = display_requirement(shared, "niche_depth_min")
        interface_lines.append(
            svg_text(
                382.0,
                row_y + 1.0,
                f"Niche H{niche_height_min}–{niche_height_max} · "
                f"W{niche_width_min}–{niche_width_max} · D≥{niche_depth_min} "
                "(official product constraint)",
                "interface",
            )
        )
        interface_lines.append(
            svg_text(382.0, row_y + 5.0, "BLOCK: rough-in XYZ / valves / hose path / opening position", "block")
        )

    legend = []
    legend_y = 373.5
    legend_x = 12.0
    for role in sorted({record["installation_role"] for record in records}):
        fill, stroke, label = ROLE_STYLE[role]
        legend.append(
            f'<rect x="{legend_x:.2f}" y="{legend_y - 3.1:.2f}" width="3" height="3" fill="{fill}" stroke="{stroke}" stroke-width="0.35"/>'
        )
        legend.append(svg_text(legend_x + 4.3, legend_y - 0.3, label, "legend"))
        legend_x += 37.0

    return f'''<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
 width="500mm" height="400mm" viewBox="0 0 500 400" data-sheet-id="{sheet_id}"
 data-source-ifc-sha256="{source_hash}" data-underlay-sha256="{underlay_hash}">
<style>
  @page {{ size: 500mm 400mm; margin: 0; }}
  html, body {{ margin: 0; width: 500mm; height: 400mm; overflow: hidden; }}
  text {{ font-family: Arial, "PingFang SC", sans-serif; fill: #20252a; }}
  .title {{ font-size: 6px; font-weight: 700; }}
  .subtitle {{ font-size: 4px; font-weight: 600; }}
  .meta {{ font-size: 2.7px; fill: #4f5961; }}
  .index {{ font-size: 2.75px; }}
  .legend {{ font-size: 2.5px; }}
  .block {{ font-size: 2.65px; font-weight: 700; fill: #b32121; }}
  .interface {{ font-size: 2.55px; fill: #244f65; }}
</style>
<rect width="500" height="400" fill="white"/>
<rect x="7" y="7" width="366" height="366" fill="#fafafa" stroke="#30363b" stroke-width="0.5"/>
<image x="10" y="10" width="360" height="360" opacity="0.58" xlink:href="{escape(spec['underlay'])}"/>
<g transform="translate(10 10) scale(0.9)">{''.join(overlays)}</g>
<rect x="378" y="7" width="115" height="386" fill="#f7f8f9" stroke="#30363b" stroke-width="0.5"/>
{svg_text(382, 17, f"{sheet_id}  {spec['title']}", "title")}
{svg_text(382, 24, spec['subtitle'], "subtitle")}
{svg_text(382, 31, f"Status: CANDIDATE / NOT FOR CONSTRUCTION", "block")}
{svg_text(382, 37, f"Existing objects: {len(records)}", "meta")}
{svg_text(382, 42, f"IFC SHA: {source_hash[:16]}...", "meta")}
{svg_text(382, 47, "Overlay = world bbox coordination envelope", "meta")}
{svg_text(382, 51, "NOT fabrication / opening / rough-in dimensions", "block")}
{''.join(index_lines)}
{''.join(interface_lines)}
{''.join(block_lines)}
{''.join(legend)}
{svg_text(12, 384, "Exact GlobalId/object index: pipeline/decisions/int1-existing-review.csv", "meta")}
{svg_text(12, 389, "Base plan retained for context; overlay uses existing world-coordinate bbox only.", "meta")}
</svg>'''


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", required=True, type=Path)
    parser.add_argument("--source-ifc", required=True, type=Path)
    parser.add_argument("--drawings-dir", required=True, type=Path)
    parser.add_argument("--pdf-dir", required=True, type=Path)
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--existing-report", required=True, type=Path)
    parser.add_argument("--render-script", type=Path)
    parser.add_argument("--render-pdfs", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    records, blockers, source_hash = read_candidate(args.input_csv)
    actual_hash = sha256(args.source_ifc)
    if source_hash != actual_hash:
        raise RuntimeError("INT1 register does not match the formal IFC")
    if len(blockers) != 4:
        raise RuntimeError(f"expected four disclosed INT1 blockers, found {len(blockers)}")
    if Counter(record["sheet_id"] for record in records) != Counter(
        {"I-501": 70, "I-502": 33, "I-503": 1, "I-504": 25}
    ):
        raise RuntimeError("unexpected INT1 sheet object counts")
    dishwasher_interfaces = read_dishwasher_interfaces(args.existing_report, source_hash)

    args.drawings_dir.mkdir(parents=True, exist_ok=True)
    args.pdf_dir.mkdir(parents=True, exist_ok=True)
    args.build_dir.mkdir(parents=True, exist_ok=True)
    sheet_reports = []
    for sheet_id, spec in SHEETS.items():
        sheet_records = [record for record in records if record["sheet_id"] == sheet_id]
        underlay = args.drawings_dir / spec["underlay"]
        if not underlay.is_file():
            raise RuntimeError(f"missing INT1 underlay: {underlay}")
        output_svg = args.drawings_dir / spec["svg"]
        output_svg.write_text(
            make_svg(
                sheet_id,
                sheet_records,
                blockers,
                source_hash,
                sha256(underlay),
                dishwasher_interfaces,
            ) + "\n",
            encoding="utf-8",
        )
        root = ET.parse(output_svg).getroot()
        if root.attrib.get("data-sheet-id") != sheet_id:
            raise RuntimeError(f"SVG sheet ID gate failed: {sheet_id}")
        overlays = [node for node in root.iter() if node.tag.endswith("rect") and node.attrib.get("data-global-id")]
        if len(overlays) != len(sheet_records):
            raise RuntimeError(f"SVG overlay count gate failed: {sheet_id}")
        report = {
            "sheet_id": sheet_id,
            "source_ifc_sha256": source_hash,
            "source_register_sha256": sha256(args.input_csv),
            "underlay": str(underlay),
            "underlay_sha256": sha256(underlay),
            "svg": str(output_svg),
            "svg_sha256": sha256(output_svg),
            "existing_object_count": len(sheet_records),
            "overlay_count": len(overlays),
            "grouped_index_count": len(grouped_index(sheet_records)),
            "dimension_status": "existing_world_bbox_not_fabrication_dimension",
            "block_count": len(blockers),
            "installation_interface_row_count": (
                len(dishwasher_interfaces) if sheet_id == "I-501" else 0
            ),
            "mechanical_pass": True,
        }
        if args.render_pdfs:
            if not args.render_script:
                raise RuntimeError("--render-script is required with --render-pdfs")
            pdf = args.pdf_dir / spec["pdf"]
            proof = args.build_dir / f"{sheet_id}-proof.png"
            pdf_report = args.build_dir / f"{sheet_id}-pdf-report.json"
            subprocess.run(
                [
                    "python3", str(args.render_script),
                    "--input-svg", str(output_svg),
                    "--output-pdf", str(pdf),
                    "--proof-png", str(proof),
                    "--report", str(pdf_report),
                ],
                check=True,
            )
            rendered = json.loads(pdf_report.read_text(encoding="utf-8"))
            if not rendered.get("pass"):
                raise RuntimeError(f"PDF mechanical gate failed: {sheet_id}")
            report.update(
                {
                    "pdf": str(pdf),
                    "pdf_sha256": sha256(pdf),
                    "pdf_page": rendered["page"],
                    "pdf_pass": True,
                    "proof_png": str(proof),
                }
            )
        sheet_reports.append(report)

    final_report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "mode": "read_only_int1_drawing_candidate",
        "source_ifc_sha256": source_hash,
        "summary": {
            "sheet_count": 4,
            "existing_object_count": len(records),
            "block_count": len(blockers),
            "all_overlays_are_coordination_envelopes": True,
            "pdfs_rendered": args.render_pdfs,
            "installation_interface_row_count": len(dishwasher_interfaces),
            "unlocated_interface_count": sum(
                row["project_interface_status"] == "unlocated"
                for row in dishwasher_interfaces
            ),
        },
        "gates": {
            "candidate_generation_pass": all(row["mechanical_pass"] for row in sheet_reports),
            "formal_ifc_write_allowed": False,
            "fabrication_dimension_ready": False,
            "int1_completion_pass": False,
        },
        "sheets": sheet_reports,
        "blockers": [
            {
                "issue_id": row["object_name"],
                "scope": row["type_name"],
                "stop_condition": row["stop_condition"],
            }
            for row in blockers
        ],
    }
    (args.build_dir / "int1-drawing-report.json").write_text(
        json.dumps(final_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"summary": final_report["summary"], "gates": final_report["gates"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
