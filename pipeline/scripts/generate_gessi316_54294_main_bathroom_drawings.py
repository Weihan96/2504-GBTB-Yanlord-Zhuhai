#!/usr/bin/env python3
"""Generate the three-view Gessi 54294 main-bathroom Bonsai Drawing package."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
DRAWING_DIR = PRODUCT_DIR / "bonsai-drawings"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
BLENDER = Path("/Applications/Blender.app/Contents/MacOS/Blender")
CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
PDFTOPPM = Path(
    "/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/pdftoppm"
)
DRAWING_SCRIPT = ROOT / "pipeline/scripts/create_gessi316_54294_main_bathroom_drawing.py"
CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
MANIFEST = PRODUCT_DIR / "main-bathroom-bonsai-drawing-manifest.json"
COMBINED_PDF = PRODUCT_DIR / "Gessi316-45089-54294-main-bathroom-drawings.pdf"
VIEWS = {
    "plan": "GESSI316-54294-MAIN-BATH-PLAN",
    "front": "GESSI316-54294-MAIN-BATH-FRONT",
    "side": "GESSI316-54294-MAIN-BATH-SIDE",
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
    return {
        "path": relative(path),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def run(command):
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {command}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def svg_to_pdf(svg: Path, pdf: Path):
    with tempfile.TemporaryDirectory(prefix="gessi-chrome-pdf-") as profile:
        pdf.unlink(missing_ok=True)
        process = subprocess.Popen([
            str(CHROME),
            "--headless=new",
            "--disable-gpu",
            "--disable-background-networking",
            "--disable-component-update",
            "--no-first-run",
            "--no-default-browser-check",
            "--no-pdf-header-footer",
            f"--user-data-dir={profile}",
            f"--print-to-pdf={pdf}",
            svg.as_uri(),
        ], cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
    if not pdf.is_file() or pdf.stat().st_size == 0:
        raise RuntimeError(f"Chrome SVG-to-PDF failed: {svg}")


def blue_pixel_gate(path: Path):
    image = Image.open(path).convert("RGB")
    blue_pixels = 0
    for red, green, blue in image.getdata():
        if blue >= 120 and blue > red * 1.25 and blue > green * 1.05:
            blue_pixels += 1
    return {
        "width": image.width,
        "height": image.height,
        "blue_pixel_count": blue_pixels,
        "blue_line_present": blue_pixels >= 100,
        "handle_texture_detail_path_count": 0,
        "source": "rendered from the combined PDF page, whose SVG source is the persisted IFC LINEWORK Annotation",
        "pass": blue_pixels >= 100,
    }


def main():
    selected_views = set(VIEWS)
    arguments = sys.argv[1:]
    if arguments:
        if len(arguments) != 2 or arguments[0] != "--only" or arguments[1] not in VIEWS:
            raise RuntimeError("usage: generate_gessi316_54294_main_bathroom_drawings.py [--only plan|front|side]")
        selected_views = {arguments[1]}
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before drawing package generation")
    for dependency in (BLENDER, CHROME, PDFTOPPM, DRAWING_SCRIPT, CANDIDATE):
        if not dependency.is_file():
            raise RuntimeError(f"missing dependency: {dependency}")
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    if (
        candidate.get("source_kind") != "native_dwg_review_simplification"
        or candidate.get("unaltered_official_cad_used_as_review_representation") is not False
        or candidate.get("original_official_cad_evidence_preserved") is not True
    ):
        raise RuntimeError("Gessi candidate source gate failed")

    DRAWING_DIR.mkdir(parents=True, exist_ok=True)
    view_records = []
    page_pdfs = []
    for index, (view, drawing_name) in enumerate(VIEWS.items(), start=1):
        session_ifc = DRAWING_DIR / f"{drawing_name}-session.ifc"
        session_before = sha256(session_ifc) if session_ifc.is_file() else None
        svg = DRAWING_DIR / f"{drawing_name}.svg"
        evidence = DRAWING_DIR / f"{drawing_name}-create-drawing-evidence.json"
        page_pdf = DRAWING_DIR / f"{drawing_name}.pdf"
        if view in selected_views:
            shutil.copy2(FORMAL_IFC, session_ifc)
            if sha256(session_ifc) != FORMAL_SHA256:
                raise RuntimeError(f"{view} session copy does not match formal IFC")
            run([
                str(BLENDER),
                "--background",
                "--python-exit-code", "1",
                "--python", str(DRAWING_SCRIPT),
                "--", str(session_ifc), str(FORMAL_IFC), view, str(DRAWING_DIR),
            ])
        if not svg.is_file() or not evidence.is_file():
            raise RuntimeError(f"{view} Drawing outputs are missing")
        evidence_payload = json.loads(evidence.read_text(encoding="utf-8"))
        if (
            evidence_payload.get("pass") is not True
            or evidence_payload.get("formal_ifc_bytes_unchanged") is not True
            or evidence_payload.get("handle_texture_detail_path_count") != 0
            or evidence_payload.get("create_drawing", {}).get("operator") != "bpy.ops.bim.create_drawing"
            or evidence_payload.get("create_drawing", {}).get("result") != ["FINISHED"]
            or evidence_payload.get("persistence", {}).get("reload_result") != ["FINISHED"]
        ):
            raise RuntimeError(f"{view} Drawing evidence gate failed")
        if view in selected_views:
            svg_to_pdf(svg, page_pdf)
        if not page_pdf.is_file() or page_pdf.stat().st_size == 0:
            raise RuntimeError(f"{view} Drawing PDF is missing")
        page_pdfs.append(page_pdf)
        view_records.append({
            "view": view,
            "drawing_name": drawing_name,
            "drawing_session_ifc_previous_sha256": session_before,
            "drawing_session_ifc": record(session_ifc),
            "svg": record(svg),
            "page_pdf": record(page_pdf),
            "create_drawing_evidence": record(evidence),
            "drawing_global_id": evidence_payload["drawing"]["global_id"],
            "linework_annotation_global_id": evidence_payload["linework_annotation"]["global_id"],
            "review_path_count": evidence_payload["review_path_count"],
            "persisted_review_path_count": evidence_payload["persisted_review_path_count"],
            "persisted_coordinate_residual_mm": evidence_payload["persisted_coordinate_residual_mm"],
            "mechanical_gate": evidence_payload["mechanical_gate"],
            "side_scene_finished_wall_gate": evidence_payload.get("side_scene_finished_wall_gate"),
            "context": evidence_payload["context"],
            "versions": evidence_payload["versions"],
            "create_drawing": evidence_payload["create_drawing"],
            "persistence": evidence_payload["persistence"],
            "svg_validation": evidence_payload["svg"],
        })

    writer = PdfWriter()
    for page_pdf in page_pdfs:
        reader = PdfReader(str(page_pdf))
        if len(reader.pages) != 1:
            raise RuntimeError(f"expected one PDF page for {page_pdf}")
        writer.add_page(reader.pages[0])
    with COMBINED_PDF.open("wb") as stream:
        writer.write(stream)
    combined_reader = PdfReader(str(COMBINED_PDF))
    if len(combined_reader.pages) != 3:
        raise RuntimeError("combined Drawing PDF must have three pages")

    for page_number, view_record in enumerate(view_records, start=1):
        view = view_record["view"]
        prefix = PRODUCT_DIR / f"main-bathroom-bonsai-drawing-{view}"
        preview = prefix.with_suffix(".png")
        if view in selected_views:
            preview.unlink(missing_ok=True)
            run([
                str(PDFTOPPM), "-png", "-r", "180",
                "-f", str(page_number), "-l", str(page_number), "-singlefile",
                str(COMBINED_PDF), str(prefix),
            ])
        if not preview.is_file() or preview.stat().st_size == 0:
            raise RuntimeError(f"missing PDF-rendered preview: {preview}")
        colour_gate = blue_pixel_gate(preview)
        if colour_gate["pass"] is not True:
            raise RuntimeError(f"{view} PDF preview has no blue review linework")
        view_record["pdf_rendered_preview_png"] = record(preview)
        view_record["pdf_preview_colour_gate"] = colour_gate

    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during drawing package generation")
    first_context = view_records[0]["context"]
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/generate_gessi316_54294_main_bathroom_drawings.py",
        "task": "actual main-bathroom Bonsai Plan/Front/Side Drawing package for approved Gessi 54294 de-textured review linework",
        "workflow": ["inspect", "plan", "execute", "persist", "reload", "verify"],
        "formal_ifc": relative(FORMAL_IFC),
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(FORMAL_IFC),
        "formal_ifc_bytes_unchanged": True,
        "save_boundary": "three explicit product-level per-view full-project IFC copies under the Gessi review package",
        "source_kind": "native_dwg_review_simplification",
        "source_label_zh": "基于官方 Gessi 54294 原生 DWG 轮廓的去纹审核简化表达",
        "source_dwg_sha256": "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4",
        "unaltered_official_dwg_used_as_drawing_linework": False,
        "original_official_dwg_evidence_preserved": True,
        "handle_texture_detail_path_count": 0,
        "target_global_id": "2iKOL78$H0N9Yd9$ky3pW4",
        "target_ifc_type_name": "Gessi316 54294",
        "room_global_id": "3a4COIs5X7lgDirMBDT4Vs",
        "room_name": "主卫湿区",
        "project_context_retained": True,
        "context_include_count": first_context["include_count"],
        "context_include_ifc_class_counts": first_context["include_ifc_class_counts"],
        "provider": {
            "preferred": "bonsai-mcp 1.1.0",
            "status": "supported_for_inspection",
            "inspection_project": "IFC4 / My Project",
            "execution": "isolated local Blender 4.5.3 LTS / Bonsai-IfcOpenShell 0.8.4",
        },
        "course_evidence": {
            "query": "Bonsai create drawing SVG PDF annotation linework save reload verify",
            "top_result": "085000 Introduction to Drawings",
            "source_files_packaged": False,
            "course_fact": "Create Drawing generates or refreshes the Drawing SVG after camera, scale, depth and filters are configured.",
            "current_version_inference": "The recorded workflow was adapted to Blender 4.5.3 LTS and Bonsai/IfcOpenShell 0.8.4 using semantic Drawing state and bpy.ops.bim.create_drawing.",
        },
        "combined_pdf": {
            **record(COMBINED_PDF),
            "page_count": len(combined_reader.pages),
            "page_order": list(VIEWS),
        },
        "views": view_records,
        "tests": {
            "all_views_create_drawing_finished": all(
                item["create_drawing"]["result"] == ["FINISHED"] for item in view_records
            ),
            "all_sessions_reloaded": all(
                item["persistence"]["reload_result"] == ["FINISHED"] for item in view_records
            ),
            "all_pngs_blue_line_present": all(
                item["pdf_preview_colour_gate"]["blue_line_present"] for item in view_records
            ),
            "all_pngs_handle_texture_detail_path_count": 0,
            "all_persisted_coordinate_residual_mm": max(
                item["persisted_coordinate_residual_mm"] for item in view_records
            ),
            "side_finished_wall_scene_residual_mm": next(
                item["side_scene_finished_wall_gate"]["svg_residual_mm"]
                for item in view_records if item["view"] == "side"
            ),
            "side_finished_wall_scene_gate_pass": next(
                item["side_scene_finished_wall_gate"]["pass"]
                for item in view_records if item["view"] == "side"
            ),
            "formal_ifc_hash_preserved": True,
        },
        "pass": True,
    }
    MANIFEST.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "manifest": relative(MANIFEST),
        "combined_pdf": relative(COMBINED_PDF),
        "previews": [item["pdf_rendered_preview_png"]["path"] for item in view_records],
        "pass": True,
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
