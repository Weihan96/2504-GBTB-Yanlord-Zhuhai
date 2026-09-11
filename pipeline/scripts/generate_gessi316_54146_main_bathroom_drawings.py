#!/usr/bin/env python3
"""Generate approved Gessi316 54146 main-bathroom Bonsai Drawings and PDF."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54146"
DRAWING_DIR = PRODUCT_DIR / "bonsai-drawings/main-bathroom"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Gessi316-54146-derived-drawing.ifc"
DERIVED_REPORT = PRODUCT_DIR / "Gessi316-54146-derived-drawing-report.json"
SESSION_BLEND = PRODUCT_DIR / "Gessi316-54146-main-bathroom-drawings.blend"
BLENDER = Path("/Applications/Blender.app/Contents/MacOS/Blender")
CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
PDFTOPPM = Path("/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/pdftoppm")
DRAWING_SCRIPT = ROOT / "pipeline/scripts/create_gessi316_54146_main_bathroom_drawings.py"
EVIDENCE = DRAWING_DIR / "GESSI316-54146-MAIN-BATH-create-drawing-evidence.json"
MANIFEST = PRODUCT_DIR / "main-bathroom-bonsai-drawing-manifest.json"
COMBINED_PDF = PRODUCT_DIR / "Gessi316-54146-main-bathroom-drawings.pdf"
VIEWS = {
    "plan": "GESSI316-54146-MAIN-BATH-PLAN",
    "front": "GESSI316-54146-MAIN-BATH-FRONT",
    "side": "GESSI316-54146-MAIN-BATH-SIDE",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def record(path: Path) -> dict:
    return {"path": relative(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def run(command):
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {command}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def svg_to_pdf(svg: Path, pdf: Path):
    with tempfile.TemporaryDirectory(prefix="gessi54146-chrome-pdf-") as profile:
        pdf.unlink(missing_ok=True)
        process = subprocess.Popen(
            [
                str(CHROME), "--headless=new", "--disable-gpu", "--disable-background-networking",
                "--disable-component-update", "--no-first-run", "--no-default-browser-check",
                "--no-pdf-header-footer", f"--user-data-dir={profile}", f"--print-to-pdf={pdf}", svg.as_uri(),
            ],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        for _ in range(300):
            if pdf.is_file() and pdf.stat().st_size:
                break
            if process.poll() is not None:
                break
            time.sleep(0.1)
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=2)
    if not pdf.is_file() or not pdf.stat().st_size:
        raise RuntimeError(f"Chrome SVG-to-PDF failed: {svg}")


def colour_gate(path: Path) -> dict:
    image = Image.open(path).convert("RGB")
    blue = grey = nonwhite = 0
    for red, green, value in image.getdata():
        if min(red, green, value) < 245:
            nonwhite += 1
        if value >= 110 and value > red * 1.20 and value > green * 1.02:
            blue += 1
        if abs(red - green) < 18 and abs(green - value) < 18 and 70 <= red <= 220:
            grey += 1
    return {
        "width": image.width,
        "height": image.height,
        "blue_pixel_count": blue,
        "grey_pixel_count": grey,
        "nonwhite_pixel_count": nonwhite,
        "blue_line_present": blue >= 100,
        "project_context_present": grey >= 100,
        "pass": blue >= 100 and grey >= 100,
    }


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-only", action="store_true")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before Gessi316 54146 Drawing package generation")
    for dependency in (DERIVED_IFC, DERIVED_REPORT, BLENDER, CHROME, PDFTOPPM, DRAWING_SCRIPT):
        if not dependency.is_file():
            raise RuntimeError(f"missing dependency: {dependency}")
    DRAWING_DIR.mkdir(parents=True, exist_ok=True)
    if not arguments.package_only:
        run([
            str(BLENDER), "--background", "--python-exit-code", "1", "--python", str(DRAWING_SCRIPT),
            "--", str(DERIVED_IFC), str(DRAWING_DIR), str(SESSION_BLEND),
        ])
    if not EVIDENCE.is_file():
        raise RuntimeError("Bonsai Create Drawing evidence missing")
    evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
    if (
        evidence.get("pass") is not True
        or evidence.get("formal_ifc_bytes_unchanged") is not True
        or evidence.get("persistence", {}).get("reload_result") != ["FINISHED"]
        or evidence.get("source_kind") != "native_dwg_review_simplification"
        or evidence.get("unaltered_official_dwg") is not False
        or evidence.get("fine_spray_nozzle_detail_path_count") != 0
    ):
        raise RuntimeError("Gessi316 54146 Create Drawing evidence gate failed")
    writer_report = json.loads(DERIVED_REPORT.read_text(encoding="utf-8"))
    writer_report["derived_ifc_sha256_before_bonsai_drawings"] = evidence["preState"]["derived_ifc_sha256"]
    writer_report["derived_ifc_sha256"] = sha256(DERIVED_IFC)
    writer_report["derived_ifc_sha256_after_bonsai_drawings"] = sha256(DERIVED_IFC)
    writer_report["bonsai_drawings_persisted"] = True
    DERIVED_REPORT.write_text(json.dumps(writer_report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    evidence_views = {item["view"]: item for item in evidence["views"]}
    view_records = []
    page_pdfs = []
    for view, drawing_name in VIEWS.items():
        svg = DRAWING_DIR / f"{drawing_name}.svg"
        page_pdf = DRAWING_DIR / f"{drawing_name}.pdf"
        if not svg.is_file():
            raise RuntimeError(f"missing Bonsai Drawing SVG: {svg}")
        svg_to_pdf(svg, page_pdf)
        page_pdfs.append(page_pdf)
        view_evidence = evidence_views[view]
        if (
            view_evidence["create_drawing"]["result"] != ["FINISHED"]
            or view_evidence["persisted_coordinate_residual_mm"] > 0.000001
            or view_evidence["svg"]["target_ifc_projection_group_count"] != 0
            or view_evidence["review_path_count"] != view_evidence["persisted_review_path_count"]
        ):
            raise RuntimeError(f"{view} persistence/no-ghost gate failed")
        view_records.append({
            "view": view,
            "drawing_name": drawing_name,
            "drawing_global_id": view_evidence["drawing"]["global_id"],
            "linework_annotation_global_id": view_evidence["linework_annotation"]["global_id"],
            "review_path_count": view_evidence["review_path_count"],
            "persisted_review_path_count": view_evidence["persisted_review_path_count"],
            "persisted_coordinate_residual_mm": view_evidence["persisted_coordinate_residual_mm"],
            "target_ifc_projection_group_count": view_evidence["svg"]["target_ifc_projection_group_count"],
            "annotation_geometry_count": view_evidence["svg"]["review_annotation_geometry_count"],
            "camera": view_evidence["camera"],
            "svg": record(svg),
            "page_pdf": record(page_pdf),
        })
    writer = PdfWriter()
    for page_pdf in page_pdfs:
        reader = PdfReader(str(page_pdf))
        if len(reader.pages) != 1:
            raise RuntimeError(f"expected one PDF page: {page_pdf}")
        writer.add_page(reader.pages[0])
    with COMBINED_PDF.open("wb") as stream:
        writer.write(stream)
    reader = PdfReader(str(COMBINED_PDF))
    if len(reader.pages) != 3:
        raise RuntimeError("combined Gessi316 54146 Drawing PDF must contain three pages")
    for page_number, view_record in enumerate(view_records, start=1):
        view = view_record["view"]
        prefix = PRODUCT_DIR / f"main-bathroom-bonsai-drawing-{view}"
        preview = prefix.with_suffix(".png")
        preview.unlink(missing_ok=True)
        run([
            str(PDFTOPPM), "-png", "-r", "180", "-f", str(page_number), "-l", str(page_number),
            "-singlefile", str(COMBINED_PDF), str(prefix),
        ])
        gate = colour_gate(preview)
        if not gate["pass"]:
            raise RuntimeError(f"{view} PDF preview colour/context gate failed: {gate}")
        view_record["pdf_rendered_preview_png"] = record(preview)
        view_record["pdf_preview_colour_gate"] = gate
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during Gessi316 54146 Drawing package generation")
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": relative(Path(__file__)),
        "task": "approved Gessi316 54146 actual main-bathroom Bonsai Plan/Front/Side Drawing package",
        "workflow": ["inspect", "plan", "execute", "persist", "reload", "verify"],
        "approval_scope": "direction-corrected Plan/Front/Side; review simplification based on official Gessi 54146 G000 native DWG outline",
        "formal_ifc": relative(FORMAL_IFC),
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(FORMAL_IFC),
        "formal_authoritative_ifc_write_allowed": False,
        "derived_ifc": record(DERIVED_IFC),
        "derived_ifc_write_report": record(DERIVED_REPORT),
        "derived_ifc_write_allowed": True,
        "source_kind": "native_dwg_review_simplification",
        "source_label_zh": "基于官方 Gessi 54146 G000 原生 DWG 轮廓的简化蓝线审核表达",
        "source_dwg_sha256": "c8ddb90f61565d5273a32574f777812bbb1d9ef33e2b98479f1e7c381e69950d",
        "unaltered_official_dwg_used_as_drawing_linework": False,
        "original_official_dwg_evidence_preserved": True,
        "fine_spray_nozzle_detail_path_count": 0,
        "target_global_id": "04DLh1Jk9Dcu9ibcaE0id8",
        "target_ifc_type_name": "Gessi316 54146",
        "room_global_id": "3a4COIs5X7lgDirMBDT4Vs",
        "room_name": "主卫湿区",
        "project_context_retained": True,
        "target_body_suppressed_from_projection": True,
        "create_drawing_operator": "bpy.ops.bim.create_drawing",
        "provider": evidence["provider"],
        "course_evidence": evidence["course_evidence"],
        "create_drawing_evidence": record(EVIDENCE),
        "session_blend": record(SESSION_BLEND),
        "combined_pdf": {**record(COMBINED_PDF), "page_count": 3, "page_order": list(VIEWS)},
        "views": view_records,
        "tests": {
            "all_create_drawing_finished": True,
            "all_annotations_persisted": all(item["review_path_count"] == item["persisted_review_path_count"] for item in view_records),
            "all_target_body_projection_groups_zero": all(item["target_ifc_projection_group_count"] == 0 for item in view_records),
            "all_pdf_previews_have_blue_line_and_project_context": all(item["pdf_preview_colour_gate"]["pass"] for item in view_records),
            "formal_ifc_hash_preserved": True,
        },
        "pass": True,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": relative(MANIFEST), "combined_pdf": relative(COMBINED_PDF), "pass": True}, indent=2))


if __name__ == "__main__":
    main()
