#!/usr/bin/env python3
"""Package Miami Soft E09 Bonsai Drawing SVGs as PDF and previews."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/miamisoft-e09"
DRAWING_DIR = PRODUCT_DIR / "bonsai-drawings/living-area"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Baxter-Miami-Soft-E09-derived-drawing.ifc"
SESSION_BLEND = PRODUCT_DIR / "Baxter-Miami-Soft-E09-project-drawings.blend"
EVIDENCE = DRAWING_DIR / "MIAMISOFT-E09-create-drawing-evidence.json"
APPROVAL = ROOT / "pipeline/decisions/miamisoft-e09-drawing-approval.json"
MANIFEST = PRODUCT_DIR / "living-area-bonsai-drawing-manifest.json"
COMBINED_PDF = PRODUCT_DIR / "Baxter-Miami-Soft-E09-living-area-drawings.pdf"
CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
PDFTOPPM = Path(
    "/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/pdftoppm"
)
BLUE = "#1677c8"
GREY = "#a3abb3"
VIEWS = {
    "plan": "MIAMISOFT-E09-LIVING-AREA-PLAN",
    "front": "MIAMISOFT-E09-LIVING-AREA-FRONT",
    "side": "MIAMISOFT-E09-LIVING-AREA-SIDE",
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
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {command}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def svg_to_pdf(svg: Path, pdf: Path):
    with tempfile.TemporaryDirectory(prefix="miamisoft-e09-chrome-pdf-") as profile:
        pdf.unlink(missing_ok=True)
        process = subprocess.Popen(
            [
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
            ],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        for _ in range(450):
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


def image_gate(path: Path):
    image = Image.open(path).convert("RGB")
    blue_pixels = 0
    grey_pixels = 0
    ink = []
    for index, (red, green, blue) in enumerate(image.getdata()):
        x = index % image.width
        y = index // image.width
        if min(red, green, blue) < 245:
            ink.append((x, y))
        if blue >= 110 and blue > red * 1.25 and blue > green * 1.03:
            blue_pixels += 1
        if abs(red - green) <= 20 and abs(green - blue) <= 20 and 80 <= red <= 225:
            grey_pixels += 1
    if not ink:
        raise RuntimeError(f"preview has no drawing ink: {path}")
    minimum_x = min(point[0] for point in ink)
    maximum_x = max(point[0] for point in ink)
    minimum_y = min(point[1] for point in ink)
    maximum_y = max(point[1] for point in ink)
    margins = {
        "left": minimum_x,
        "right": image.width - 1 - maximum_x,
        "top": minimum_y,
        "bottom": image.height - 1 - maximum_y,
    }
    passed = blue_pixels >= 100 and grey_pixels >= 100 and min(margins.values()) >= 4
    return {
        "width": image.width,
        "height": image.height,
        "blue_pixel_count": blue_pixels,
        "grey_context_pixel_count": grey_pixels,
        "ink_bbox_px": [minimum_x, minimum_y, maximum_x, maximum_y],
        "page_margins_px": margins,
        "blue_line_present": blue_pixels >= 100,
        "grey_context_present": grey_pixels >= 100,
        "no_page_edge_clipping": min(margins.values()) >= 4,
        "pass": passed,
    }


def hidden_line_gate(svg: Path, view: str):
    root = ET.parse(svg).getroot()
    hidden = [
        element
        for element in root.iter()
        if element.attrib.get("data-line-role") == "hidden"
    ]
    dashed = [
        element
        for element in hidden
        if "stroke-dasharray:2.4,1.5" in element.attrib.get("style", "")
    ]
    expected = 259 if view == "side" else 0
    return {
        "hidden_segment_count": len(hidden),
        "dashed_hidden_segment_count": len(dashed),
        "expected_hidden_segment_count": expected,
        "pass": len(hidden) == expected and len(dashed) == expected,
    }


def main():
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before Miami E09 packaging")
    for dependency in (
        CHROME,
        PDFTOPPM,
        DERIVED_IFC,
        SESSION_BLEND,
        EVIDENCE,
        APPROVAL,
    ):
        if not dependency.is_file():
            raise RuntimeError(f"missing dependency: {dependency}")
    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
        or evidence.get("pass") is not True
        or evidence.get("formal_ifc_bytes_unchanged") is not True
        or evidence.get("tests", {}).get("all_create_drawing_finished") is not True
        or evidence.get("tests", {}).get("all_target_body_projection_counts_zero") is not True
        or evidence.get("tests", {}).get("side_hidden_path_count") != 1
        or evidence.get("tests", {}).get("side_hidden_dasharray") != "2.4,1.5"
    ):
        raise RuntimeError("Miami E09 Create Drawing evidence gate failed")

    by_view = {item["view"]: item for item in evidence["views"]}
    page_pdfs = []
    view_records = []
    for view, drawing_name in VIEWS.items():
        view_evidence = by_view[view]
        svg = DRAWING_DIR / f"{drawing_name}.svg"
        if (
            view_evidence["create_drawing"]["operator"] != "bpy.ops.bim.create_drawing"
            or view_evidence["create_drawing"]["result"] != ["FINISHED"]
            or view_evidence["create_drawing"]["linework_mode"] != "OPENCASCADE"
            or view_evidence["svg"]["target_ifc_projection_group_count"] != 0
            or view_evidence["svg"]["duplicate_target_or_annotation"] is not False
            or not svg.is_file()
            or svg.stat().st_size == 0
        ):
            raise RuntimeError(f"Miami E09 {view} SVG evidence drifted")
        hidden_gate = hidden_line_gate(svg, view)
        if hidden_gate["pass"] is not True:
            raise RuntimeError(f"Miami E09 {view} hidden-line gate failed: {hidden_gate}")
        page_pdf = DRAWING_DIR / f"{drawing_name}.pdf"
        svg_to_pdf(svg, page_pdf)
        page_pdfs.append(page_pdf)
        view_records.append(
            {
                "view": view,
                "drawing_name": drawing_name,
                "drawing_global_id": view_evidence["drawing"]["global_id"],
                "visible_annotation_global_id": view_evidence["visible_annotation_global_id"],
                "official_path_count": view_evidence["visible_path_count"],
                "official_hidden_path_count": view_evidence["hidden_path_count"],
                "persisted_official_path_count": view_evidence["persisted_visible_path_count"],
                "persisted_hidden_path_count": view_evidence["persisted_hidden_path_count"],
                "camera": view_evidence["camera"],
                "create_drawing": view_evidence["create_drawing"],
                "svg_validation": view_evidence["svg"],
                "hidden_line_gate": hidden_gate,
                "svg": record(svg),
                "page_pdf": record(page_pdf),
            }
        )

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
        raise RuntimeError("Miami E09 combined PDF must contain three pages")

    for page_number, view_record in enumerate(view_records, start=1):
        view = view_record["view"]
        prefix = PRODUCT_DIR / f"living-area-bonsai-drawing-{view}"
        preview = prefix.with_suffix(".png")
        preview.unlink(missing_ok=True)
        run(
            [
                str(PDFTOPPM),
                "-png",
                "-r",
                "180",
                "-f",
                str(page_number),
                "-l",
                str(page_number),
                "-singlefile",
                str(COMBINED_PDF),
                str(prefix),
            ]
        )
        if not preview.is_file() or preview.stat().st_size == 0:
            raise RuntimeError(f"Miami E09 {view} preview is missing")
        gate = image_gate(preview)
        if gate["pass"] is not True:
            raise RuntimeError(f"Miami E09 {view} preview gate failed: {gate}")
        view_record["pdf_rendered_preview_png"] = record(preview)
        view_record["pdf_preview_gate"] = gate

    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during Miami E09 packaging")
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/generate_miamisoft_e09_project_drawings.py",
        "task": "actual project living-area Bonsai Plan/Front/Side Drawing package for approved Baxter Miami Soft E09 native-DWG linework",
        "workflow": ["inspect", "plan", "execute", "persist", "reload", "verify"],
        "formal_ifc": relative(FORMAL_IFC),
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(FORMAL_IFC),
        "formal_ifc_bytes_unchanged": True,
        "save_boundary": "one Miami Soft E09 product-level full-project derived IFC copy",
        "derived_ifc": record(DERIVED_IFC),
        "session_blend": record(SESSION_BLEND),
        "approval_record": record(APPROVAL),
        "create_drawing_evidence": record(EVIDENCE),
        "source_kind": "native_dwg",
        "source_label_zh": "Baxter 官方 Miami Soft E09 原生二维 DWG",
        "source_dwg_sha256": "efc3315feed0a2dac8938aebbd7336856dabba6721395e4520783c39ba0cc51a",
        "target_global_id": "3osoWufdD1mhDDAM6lcix4",
        "project_highpoly_component_offsets_accepted": True,
        "component_force_fit_applied": False,
        "actual_target_body_suppressed_from_drawing_include": True,
        "official_hidden_line_source_handle": "1115B",
        "official_hidden_line_pattern": "DASHED",
        "provider": evidence["provider"],
        "course_evidence": evidence["course_evidence"],
        "combined_pdf": {
            **record(COMBINED_PDF),
            "page_count": len(combined_reader.pages),
            "page_order": list(VIEWS),
        },
        "views": view_records,
        "tests": {
            "all_create_drawing_finished": all(
                item["create_drawing"]["result"] == ["FINISHED"] for item in view_records
            ),
            "all_previews_blue_and_grey": all(
                item["pdf_preview_gate"]["blue_line_present"]
                and item["pdf_preview_gate"]["grey_context_present"]
                for item in view_records
            ),
            "all_previews_not_clipped": all(
                item["pdf_preview_gate"]["no_page_edge_clipping"] for item in view_records
            ),
            "side_hidden_line_dashed": next(
                item["hidden_line_gate"]["pass"] for item in view_records if item["view"] == "side"
            ),
            "target_projection_group_count": sum(
                item["svg_validation"]["target_ifc_projection_group_count"]
                for item in view_records
            ),
            "formal_ifc_unchanged": True,
        },
        "pass": True,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"manifest": str(MANIFEST), "combined_pdf": str(COMBINED_PDF), "pass": True}))


if __name__ == "__main__":
    main()
