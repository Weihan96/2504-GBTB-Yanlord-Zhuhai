#!/usr/bin/env python3
"""Generate the approved BED01 master-bedroom Bonsai Drawing package."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

from PIL import Image
from pypdf import PdfReader, PdfWriter


ROOT = Path(__file__).resolve().parents[2]
PRODUCT_DIR = ROOT / "output/review/highpoly-types/bed01"
DRAWING_DIR = PRODUCT_DIR / "bonsai-drawings/master-bedroom"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DERIVED_IFC = PRODUCT_DIR / "Baxter-Casablanca-BED01-derived-drawing.ifc"
BLENDER = Path("/Applications/Blender.app/Contents/MacOS/Blender")
CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
PDFTOPPM = Path(
    "/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/override/pdftoppm"
)
DRAWING_SCRIPT = ROOT / "pipeline/scripts/create_bed01_master_bedroom_drawings.py"
APPROVAL = ROOT / "pipeline/decisions/bed01-drawing-approval.json"
REFERENCE = PRODUCT_DIR / "official-dwg-review-reference.json"
REVIEW_MANIFEST = PRODUCT_DIR / "official-dwg-review-manifest.json"
EVIDENCE = DRAWING_DIR / "BED01-MASTER-BEDROOM-create-drawing-evidence.json"
MANIFEST = PRODUCT_DIR / "master-bedroom-bonsai-drawing-manifest.json"
COMBINED_PDF = PRODUCT_DIR / "Baxter-Casablanca-BED01-master-bedroom-drawings.pdf"
VIEWS = {
    "plan": "BED01-MASTER-BEDROOM-PLAN",
    "front": "BED01-MASTER-BEDROOM-FRONT",
    "side": "BED01-MASTER-BEDROOM-SIDE",
}
BLUE = "#1677c8"
GREY = "#a3abb3"


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
    with tempfile.TemporaryDirectory(prefix="bed01-chrome-pdf-") as profile:
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


def colour_and_crop_gate(path: Path):
    image = Image.open(path).convert("RGB")
    blue_pixels = 0
    grey_pixels = 0
    ink = []
    for index, (red, green, blue) in enumerate(image.getdata()):
        x = index % image.width
        y = index // image.width
        if min(red, green, blue) < 245:
            ink.append((x, y))
        if blue >= 115 and blue > red * 1.25 and blue > green * 1.04:
            blue_pixels += 1
        if abs(red - green) <= 18 and abs(green - blue) <= 18 and 90 <= red <= 220:
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
        "source": "rendered from the combined PDF page; PDF pages originate from final Bonsai SVG files",
        "pass": passed,
    }


def main():
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch before BED01 Drawing generation")
    for dependency in (
        BLENDER,
        CHROME,
        PDFTOPPM,
        DRAWING_SCRIPT,
        APPROVAL,
        REFERENCE,
        REVIEW_MANIFEST,
    ):
        if not dependency.is_file():
            raise RuntimeError(f"missing dependency: {dependency}")
    approval = json.loads(APPROVAL.read_text(encoding="utf-8"))
    if (
        approval.get("status") != "approved"
        or approval.get("derived_ifc_write_allowed") is not True
        or approval.get("formal_authoritative_ifc_write_allowed") is not False
    ):
        raise RuntimeError("BED01 product-level derived write gate is not open")

    DRAWING_DIR.mkdir(parents=True, exist_ok=True)
    derived_previous_sha = sha256(DERIVED_IFC) if DERIVED_IFC.is_file() else None
    shutil.copy2(FORMAL_IFC, DERIVED_IFC)
    if sha256(DERIVED_IFC) != FORMAL_SHA256:
        raise RuntimeError("BED01 derived IFC copy does not match formal IFC")
    blender_result = run(
        [
            str(BLENDER),
            "--background",
            "--python-exit-code",
            "1",
            "--python",
            str(DRAWING_SCRIPT),
            "--",
            str(DERIVED_IFC),
            str(FORMAL_IFC),
            str(DRAWING_DIR),
        ]
    )
    if not EVIDENCE.is_file():
        raise RuntimeError("BED01 Create Drawing evidence is missing")
    evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
    if (
        evidence.get("pass") is not True
        or evidence.get("formal_ifc_bytes_unchanged") is not True
        or evidence.get("persistence", {}).get("reload_result") != ["FINISHED"]
        or len(evidence.get("views", [])) != 3
    ):
        raise RuntimeError("BED01 Create Drawing evidence gate failed")

    view_records = []
    page_pdfs = []
    by_view = {item["view"]: item for item in evidence["views"]}
    for view, drawing_name in VIEWS.items():
        view_evidence = by_view[view]
        if (
            view_evidence["create_drawing"]["operator"] != "bpy.ops.bim.create_drawing"
            or view_evidence["create_drawing"]["result"] != ["FINISHED"]
            or view_evidence["create_drawing"]["linework_mode"] != "OPENCASCADE"
            or view_evidence["persisted_review_path_count"]
            != view_evidence["review_path_count"]
            or view_evidence["persisted_coordinate_residual_mm"] > 0.000001
            or view_evidence["svg"]["target_ifc_projection_group_count"] != 0
            or view_evidence["svg"]["blue_style_count"] == 0
            or view_evidence["svg"]["grey_style_count"] == 0
        ):
            raise RuntimeError(f"BED01 {view} Drawing evidence drifted")
        svg = DRAWING_DIR / f"{drawing_name}.svg"
        if not svg.is_file() or svg.stat().st_size == 0:
            raise RuntimeError(f"BED01 {view} SVG missing")
        page_pdf = DRAWING_DIR / f"{drawing_name}.pdf"
        svg_to_pdf(svg, page_pdf)
        page_pdfs.append(page_pdf)
        view_records.append(
            {
                "view": view,
                "drawing_name": drawing_name,
                "drawing_global_id": view_evidence["drawing"]["global_id"],
                "linework_annotation_global_id": view_evidence["linework_annotation"][
                    "global_id"
                ],
                "review_path_count": view_evidence["review_path_count"],
                "persisted_review_path_count": view_evidence[
                    "persisted_review_path_count"
                ],
                "persisted_coordinate_residual_mm": view_evidence[
                    "persisted_coordinate_residual_mm"
                ],
                "alignment": view_evidence["alignment"],
                "camera": view_evidence["camera"],
                "create_drawing": view_evidence["create_drawing"],
                "svg_validation": view_evidence["svg"],
                "svg": record(svg),
                "page_pdf": record(page_pdf),
            }
        )

    writer = PdfWriter()
    for page_pdf in page_pdfs:
        reader = PdfReader(str(page_pdf))
        if len(reader.pages) != 1:
            raise RuntimeError(f"expected one page for {page_pdf}")
        writer.add_page(reader.pages[0])
    with COMBINED_PDF.open("wb") as stream:
        writer.write(stream)
    combined_reader = PdfReader(str(COMBINED_PDF))
    if len(combined_reader.pages) != 3:
        raise RuntimeError("BED01 combined Drawing PDF must have three pages")

    for page_number, view_record in enumerate(view_records, start=1):
        view = view_record["view"]
        prefix = PRODUCT_DIR / f"master-bedroom-bonsai-drawing-{view}"
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
            raise RuntimeError(f"BED01 {view} PDF-rendered preview missing")
        gate = colour_and_crop_gate(preview)
        if gate["pass"] is not True:
            raise RuntimeError(f"BED01 {view} preview gate failed: {gate}")
        view_record["pdf_rendered_preview_png"] = record(preview)
        view_record["pdf_preview_gate"] = gate

    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during BED01 Drawing generation")
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generator": "pipeline/scripts/generate_bed01_master_bedroom_drawings.py",
        "task": "actual master-bedroom Bonsai Plan/Front/Side Drawing package for approved Baxter Casablanca native 2D DWG linework",
        "workflow": ["inspect", "plan", "execute", "persist", "reload", "verify"],
        "formal_ifc": relative(FORMAL_IFC),
        "formal_ifc_sha256": FORMAL_SHA256,
        "formal_ifc_sha256_after_generation": sha256(FORMAL_IFC),
        "formal_ifc_bytes_unchanged": True,
        "save_boundary": "one BED01 product-level full-project derived IFC copy",
        "derived_ifc_previous_sha256": derived_previous_sha,
        "derived_ifc": record(DERIVED_IFC),
        "approval_record": record(APPROVAL),
        "official_dwg_review_reference": record(REFERENCE),
        "official_dwg_review_manifest": record(REVIEW_MANIFEST),
        "source_kind": "native_dwg_review_reference",
        "source_label_zh": "Baxter Casablanca 官方独立 2D DWG 原生蓝线",
        "source_dwg_sha256": "ab697db9c27448c62a9a77537d7cc4b6286577c335328113d703184be28d9c4a",
        "source_scale": 1.0,
        "superseded_acis_candidate_excluded": True,
        "target_global_id": "3IQBEqO5vDI8Z9k1Ltge_N",
        "target_ifc_type_name": "BED01",
        "room_global_id": "3gHz6U6BfFXgV6PnRzfOf$",
        "room_name": "主卧",
        "project_context_retained": True,
        "actual_target_body_suppressed_from_drawing_include": True,
        "context_colour": GREY,
        "official_linework_colour": BLUE,
        "provider": evidence["provider"],
        "course_evidence": evidence["course_evidence"],
        "blender_stdout_tail": blender_result.stdout[-4000:],
        "create_drawing_evidence": record(EVIDENCE),
        "combined_pdf": {
            **record(COMBINED_PDF),
            "page_count": len(combined_reader.pages),
            "page_order": list(VIEWS),
        },
        "views": view_records,
        "tests": {
            "all_views_create_drawing_finished": all(
                item["create_drawing"]["result"] == ["FINISHED"]
                for item in view_records
            ),
            "all_views_opencascade": all(
                item["create_drawing"]["linework_mode"] == "OPENCASCADE"
                for item in view_records
            ),
            "all_persisted_path_counts_match": all(
                item["review_path_count"] == item["persisted_review_path_count"]
                for item in view_records
            ),
            "all_previews_blue_and_grey": all(
                item["pdf_preview_gate"]["blue_line_present"]
                and item["pdf_preview_gate"]["grey_context_present"]
                for item in view_records
            ),
            "all_previews_not_clipped": all(
                item["pdf_preview_gate"]["no_page_edge_clipping"]
                for item in view_records
            ),
            "actual_bed_body_projection_group_count": sum(
                item["svg_validation"]["target_ifc_projection_group_count"]
                for item in view_records
            ),
            "formal_ifc_hash_preserved": True,
        },
        "pass": True,
    }
    MANIFEST.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "manifest": relative(MANIFEST),
                "derived_ifc": relative(DERIVED_IFC),
                "combined_pdf": relative(COMBINED_PDF),
                "previews": [
                    item["pdf_rendered_preview_png"]["path"] for item in view_records
                ],
                "pass": True,
            },
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
