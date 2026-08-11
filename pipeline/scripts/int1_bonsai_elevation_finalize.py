#!/usr/bin/env python3
"""Finalize native elevation document paths and write a tracked QA manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import ifcopenshell
import ifcopenshell.util.element


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REGISTER = PROJECT_ROOT / "pipeline/decisions/int1-elevation-view-register.csv"
REPORT_DIR = PROJECT_ROOT / "build/int1/native-bonsai"
OUTPUT_DIR = PROJECT_ROOT / "drawings/elevations/native"
MANIFEST = OUTPUT_DIR / "manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


def drawing_document(model: ifcopenshell.file, drawing: ifcopenshell.entity_instance):
    relations = [
        inverse
        for inverse in model.get_inverse(drawing)
        if inverse.is_a("IfcRelAssociatesDocument")
    ]
    if len(relations) != 1:
        raise RuntimeError(f"{drawing.Name}: expected one document relation")
    return relations[0].RelatingDocument


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument(
        "--skip-ifc-write",
        action="store_true",
        help="Refresh normalized SVG hashes and the manifest without rewriting the IFC.",
    )
    arguments = parser.parse_args()
    checkpoint = arguments.checkpoint.resolve()
    model = ifcopenshell.open(checkpoint)
    rows = list(csv.DictReader(REGISTER.open(encoding="utf-8-sig")))
    expected_view_ids = {row["view_id"] for row in rows}
    drawings = [
        drawing
        for drawing in model.by_type("IfcAnnotation")
        if getattr(drawing, "ObjectType", None) == "DRAWING"
        and (drawing.Name or "").startswith("EL-")
    ]
    if len(drawings) != 36:
        raise RuntimeError(f"expected 36 native elevations, got {len(drawings)}")

    reports = []
    actual_view_ids = set()
    for drawing in sorted(drawings, key=lambda item: item.Name):
        parts = drawing.Name.split("-")
        view_id = parts[2]
        actual_view_ids.add(view_id)
        svg_path = OUTPUT_DIR / f"{drawing.Name}.svg"
        if not svg_path.is_file():
            raise FileNotFoundError(svg_path)
        normalize_svg(svg_path)
        document = drawing_document(model, drawing)
        document.Location = f"drawings/elevations/native/{drawing.Name}.svg"
        source_path = REPORT_DIR / f"{drawing.Name}-source.json"
        source = json.loads(source_path.read_text(encoding="utf-8"))
        pset = ifcopenshell.util.element.get_pset(drawing, "EPset_Drawing")
        reports.append(
            {
                "view_id": view_id,
                "sheet_id": source["sheet_id"],
                "drawing_name": drawing.Name,
                "drawing_global_id": drawing.GlobalId,
                "space_reference": source["space_reference"],
                "space_global_id": source["space_global_id"],
                "direction": source["direction"],
                "source_handle": source["source_handle"],
                "source_locator": source["source_locator"],
                "target_view": pset.get("TargetView"),
                "scale": pset.get("Scale"),
                "linework_mode": pset.get("LineworkMode"),
                "has_underlay": pset.get("HasUnderlay"),
                "has_annotation": pset.get("HasAnnotation"),
                "include_count": source["include_count"],
                "demolish_wall_count": source["demolish_wall_count"],
                "complexity_exclusions": source.get("complexity_exclusions", []),
                "integer_highlight_count": source["integer_highlight_count"],
                "major_integer_highlight_count": source[
                    "major_integer_highlight_count"
                ],
                "maximum_integer_residual_mm": source[
                    "maximum_integer_residual_mm"
                ],
                "svg": f"drawings/elevations/native/{drawing.Name}.svg",
                "svg_sha256": sha256(svg_path),
            }
        )
    if actual_view_ids != expected_view_ids:
        raise RuntimeError(
            f"view id mismatch: missing={sorted(expected_view_ids-actual_view_ids)}, "
            f"extra={sorted(actual_view_ids-expected_view_ids)}"
        )

    if arguments.skip_ifc_write:
        for drawing in drawings:
            expected = f"drawings/elevations/native/{drawing.Name}.svg"
            if drawing_document(model, drawing).Location != expected:
                raise RuntimeError(f"document path verification failed: {drawing.Name}")
    else:
        temporary = checkpoint.with_name(checkpoint.name + ".finalize")
        model.write(temporary)
        verified = ifcopenshell.open(temporary)
        for drawing in verified.by_type("IfcAnnotation"):
            if getattr(drawing, "ObjectType", None) != "DRAWING" or not (
                drawing.Name or ""
            ).startswith("EL-"):
                continue
            expected = f"drawings/elevations/native/{drawing.Name}.svg"
            if drawing_document(verified, drawing).Location != expected:
                raise RuntimeError(f"document path verification failed: {drawing.Name}")
        os.replace(temporary, checkpoint)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "Bonsai 0.8.4 native Drawing / bim.create_drawing",
        "source_ifc": "2504 GBTB Yanlord Zhuhai.ifc",
        "source_ifc_sha256_before_batch": (
            "6c2fd8da9e9ad7ddbc2b63415a27f1c979e8995b880d8fce210a2dda2ef2aab6"
        ),
        "checkpoint_sha256": sha256(checkpoint),
        "view_count": len(reports),
        "sheet_count": len({report["sheet_id"] for report in reports}),
        "linework_mode_counts": dict(
            Counter(report["linework_mode"] for report in reports)
        ),
        "demolish_wall_count": sum(
            report["demolish_wall_count"] for report in reports
        ),
        "complexity_exclusion_occurrences": sum(
            len(report["complexity_exclusions"]) for report in reports
        ),
        "complexity_exclusion_unique_global_ids": sorted(
            {
                exclusion["global_id"]
                for report in reports
                for exclusion in report["complexity_exclusions"]
            }
        ),
        "views": reports,
        "pass": True,
    }
    MANIFEST.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest | {"views": f"{len(reports)} records"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
