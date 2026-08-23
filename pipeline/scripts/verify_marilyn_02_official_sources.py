#!/usr/bin/env python3
"""Revalidate Baxter Marilyn pouf official web, ZIP, DWG, 3DS, SVG and PDF evidence."""

from __future__ import annotations

import argparse
import json
import shutil
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import verify_marilyn_01_official_sources as shared


ROOT = shared.ROOT
OUTPUT_DIR = ROOT / "output/review/highpoly-types/marilyn-02"
SOURCE_DIR = OUTPUT_DIR / "official-source"
OUTPUT = SOURCE_DIR / "official-source-revalidation.json"
MATTE_SVG_URL = "https://productsbook.baxter.it/models/measurements/MARIPFMN80.svg"
GLOSSY_SVG_URL = "https://productsbook.baxter.it/models/measurements/MARIPFML80.svg"
POUF_MEMBER = "Baxter_Marilyn_Armchair_2D_3D/3D/Marilyn_pouf_80x62xh45.3ds"
EXPECTED = {
    **shared.EXPECTED,
    "matte_svg": "41a56a0df8de5fa5b362f752c94923c764816525bd7b00f0499b0dadba86f509",
    "glossy_svg": "2c2a9af1f2200e83658f0b9c17bccb022d9d3ab0bec205c9b667b5d7c6f7e96c",
    "pouf_3ds": "53819a0b1321f0177cc162a4b93d62031b4c9aab77d06370e1de90c961cb97f9",
}
PRODUCT_TOKENS = [
    "MARIPFMN80",
    "MARIPFML80",
    "Baxter_Marilyn_Armchair_2D_3D.zip",
    "POUF WITH MATTE PAINTED METAL FRAME AND SWIVEL BASE",
    "POUF WITH GLOSSY LACQUERED METAL FRAME AND SWIVEL BASE",
    "<td>80</td><td>62</td><td>45</td><td>cm</td>",
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    pdftoppm = shutil.which("pdftoppm")
    if not pdftoppm:
        raise SystemExit("pdftoppm is required for Baxter technical-sheet verification")
    local = {
        "product_page": SOURCE_DIR / "baxter-marilyn-product-page.html",
        "zip": SOURCE_DIR / "Baxter_Marilyn_Armchair_2D_3D.zip",
        "technical_pdf": SOURCE_DIR / "Baxter_Marilyn_current-technical-sheet.pdf",
        "matte_svg": SOURCE_DIR / "MARIPFMN80.svg",
        "glossy_svg": SOURCE_DIR / "MARIPFML80.svg",
        "dwg": SOURCE_DIR / "Marilyn_Abaco.dwg",
        "pouf_3ds": SOURCE_DIR / "Marilyn_pouf_80x62xh45.3ds",
    }
    with tempfile.TemporaryDirectory(prefix="baxter-marilyn02-source-") as temporary:
        temporary_dir = Path(temporary)
        downloaded = {}
        for key, url, filename in (
            ("product_page", shared.PRODUCT_PAGE_URL, "product-page.html"),
            ("zip", shared.ZIP_URL, "package.zip"),
            ("technical_pdf", shared.TECHNICAL_URL, "technical.pdf"),
            ("matte_svg", MATTE_SVG_URL, "matte.svg"),
            ("glossy_svg", GLOSSY_SVG_URL, "glossy.svg"),
        ):
            path = temporary_dir / filename
            status, content_type, byte_count = shared.download(url, path)
            downloaded[key] = {
                "url": url,
                "path": path,
                "http_status": status,
                "content_type": content_type,
                "byte_count": byte_count,
                "sha256": shared.sha256_file(path),
            }

        current_html = downloaded["product_page"]["path"].read_text(encoding="utf-8")
        product_tokens = {token: token in current_html for token in PRODUCT_TOKENS}
        product_page = {
            **{key: value for key, value in downloaded["product_page"].items() if key != "path"},
            "archived_path": str(local["product_page"].relative_to(ROOT)),
            "archived_sha256": shared.sha256_file(local["product_page"]),
            "expected_archived_sha256": EXPECTED["archived_product_page"],
            "current_page_is_dynamic_and_not_required_to_be_byte_identical_to_archive": True,
            "required_identity_tokens": product_tokens,
            "all_required_identity_tokens_found": all(product_tokens.values()),
            "pass": downloaded["product_page"]["http_status"] == 200
            and shared.sha256_file(local["product_page"]) == EXPECTED["archived_product_page"]
            and all(product_tokens.values()),
        }

        zip_path = downloaded["zip"]["path"]
        with zipfile.ZipFile(zip_path) as archive:
            member_names = archive.namelist()
            dwg_bytes = archive.read(shared.DWG_MEMBER)
            pouf_bytes = archive.read(POUF_MEMBER)
        zip_bytes_match = zip_path.read_bytes() == local["zip"].read_bytes()
        native_zip = {
            **{key: value for key, value in downloaded["zip"].items() if key != "path"},
            "local_path": str(local["zip"].relative_to(ROOT)),
            "local_sha256": shared.sha256_file(local["zip"]),
            "expected_sha256": EXPECTED["zip"],
            "downloaded_bytes_match_local_archive": zip_bytes_match,
            "member_count": len(member_names),
            "dwg_member": shared.DWG_MEMBER,
            "dwg_member_sha256": shared.sha256_bytes(dwg_bytes),
            "dwg_member_matches_local_extraction": dwg_bytes == local["dwg"].read_bytes(),
            "pouf_3ds_member": POUF_MEMBER,
            "pouf_3ds_member_sha256": shared.sha256_bytes(pouf_bytes),
            "pouf_3ds_member_matches_local_extraction": pouf_bytes == local["pouf_3ds"].read_bytes(),
            "pass": downloaded["zip"]["http_status"] == 200
            and downloaded["zip"]["sha256"] == EXPECTED["zip"]
            and shared.sha256_file(local["zip"]) == EXPECTED["zip"]
            and zip_bytes_match
            and shared.sha256_bytes(dwg_bytes) == EXPECTED["dwg"]
            and dwg_bytes == local["dwg"].read_bytes()
            and shared.sha256_bytes(pouf_bytes) == EXPECTED["pouf_3ds"]
            and pouf_bytes == local["pouf_3ds"].read_bytes(),
        }

        current_render = shared.render_page(
            pdftoppm, downloaded["technical_pdf"]["path"], temporary_dir / "current-page-13", 13
        )
        archived_render = shared.render_page(
            pdftoppm, local["technical_pdf"], temporary_dir / "archived-page-13", 13
        )
        render_bytes_match = current_render.read_bytes() == archived_render.read_bytes()
        technical_pdf = {
            **{key: value for key, value in downloaded["technical_pdf"].items() if key != "path"},
            "archived_path": str(local["technical_pdf"].relative_to(ROOT)),
            "archived_sha256": shared.sha256_file(local["technical_pdf"]),
            "expected_archived_sha256": EXPECTED["archived_technical_pdf"],
            "current_pdf_is_server_regenerated_and_not_required_to_be_byte_identical_to_archive": True,
            "identity_page": 13,
            "render_dpi": 180,
            "expected_identity_page_render_sha256": EXPECTED["technical_page_13_render"],
            "current_identity_page_render_sha256": shared.sha256_file(current_render),
            "archived_identity_page_render_sha256": shared.sha256_file(archived_render),
            "current_identity_page_render_matches_archive": render_bytes_match,
            "pass": downloaded["technical_pdf"]["http_status"] == 200
            and shared.sha256_file(local["technical_pdf"]) == EXPECTED["archived_technical_pdf"]
            and shared.sha256_file(current_render) == EXPECTED["technical_page_13_render"]
            and shared.sha256_file(archived_render) == EXPECTED["technical_page_13_render"]
            and render_bytes_match,
        }

        measurement_svgs = {}
        for key in ("matte_svg", "glossy_svg"):
            bytes_match = downloaded[key]["path"].read_bytes() == local[key].read_bytes()
            measurement_svgs[key] = {
                **{field: value for field, value in downloaded[key].items() if field != "path"},
                "local_path": str(local[key].relative_to(ROOT)),
                "local_sha256": shared.sha256_file(local[key]),
                "expected_sha256": EXPECTED[key],
                "downloaded_bytes_match_local_archive": bytes_match,
                "pass": downloaded[key]["http_status"] == 200
                and downloaded[key]["sha256"] == EXPECTED[key]
                and shared.sha256_file(local[key]) == EXPECTED[key]
                and bytes_match,
            }

    report = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Baxter",
        "family": "Marilyn",
        "project_ifc_type_name": "Marilyn 02",
        "resolved_variant": "pouf with swivel base, 80 x 62 x 45 cm",
        "marilyn_01_geometry_used": False,
        "method": "fresh official web downloads; byte comparison for ZIP and exact pouf measurement SVGs; ZIP-member hash checks for DWG and pouf 3DS; page-13 raster identity for the dynamically regenerated technical PDF",
        "product_page": product_page,
        "native_zip": native_zip,
        "technical_pdf": technical_pdf,
        "measurement_svgs": measurement_svgs,
        "pass": product_page["pass"]
        and native_zip["pass"]
        and technical_pdf["pass"]
        and all(item["pass"] for item in measurement_svgs.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not report["pass"]:
        raise SystemExit("Baxter Marilyn 02 official-source revalidation failed")


if __name__ == "__main__":
    main()
