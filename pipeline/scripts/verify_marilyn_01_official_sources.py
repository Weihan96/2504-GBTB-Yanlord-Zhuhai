#!/usr/bin/env python3
"""Revalidate Baxter Marilyn bergere official web, ZIP, DWG, 3DS, SVG and PDF evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "output/review/highpoly-types/marilyn-01"
SOURCE_DIR = OUTPUT_DIR / "official-source"
OUTPUT = SOURCE_DIR / "official-source-revalidation.json"
PRODUCT_PAGE_URL = "https://www.baxter.it/en/products/marilyn-sofas-and-armchairs"
ZIP_URL = "https://dam.baxter.it/asset/9b02ac3a-7166-40df-9acc-1e8a2ad46705/Baxter_Marilyn_Armchair_2D_3D.zip"
TECHNICAL_URL = "https://productsbook.baxter.it/product-pdf/Marilyn_divani-e-poltrone_TechnicalSheet.pdf?code=MARI&lang=eng&sector=divani-e-poltrone&kind=indoor"
MATTE_SVG_URL = "https://productsbook.baxter.it/models/measurements/MARIPBMN86.svg"
GLOSSY_SVG_URL = "https://productsbook.baxter.it/models/measurements/MARIPBML86.svg"
EXPECTED = {
    "archived_product_page": "e9f630d77aff716b01d21736117df23a491984183c2e4ea0f00c87a61df85140",
    "zip": "db054a193c9c5b8dc45f0df39199b2d96dc3712a5a7e48821fe94dce1566fedb",
    "archived_technical_pdf": "8210ae5ce84b467fd6cb46a93c665fbfaed1606e806210a751ef24ea026ff581",
    "technical_page_13_render": "977e90a62a151b8a68c4abbdc41ba4beecbefd989bd41a220770debb9fb3cfc5",
    "matte_svg": "4c2b4e94f1814bab8ffd6a154481c332577669d4afed43781c6a34fa509cc7c4",
    "glossy_svg": "d5e0b04cce3465e60948c5ec6b33748666221744eac83642e71ee5e0f71e204c",
    "dwg": "cb9825eecab7c78ee28e7668e933c2fe413395a63bf75f3afab8cc4959864724",
    "bergere_3ds": "171dbd72b2298f201d834fb54a5d39bca7ca5b86a552b97a5cbc64e23dd5e385",
}
DWG_MEMBER = "Baxter_Marilyn_Armchair_2D_3D/2D/Marilyn_Abaco.dwg"
BERGERE_MEMBER = "Baxter_Marilyn_Armchair_2D_3D/3D/Marilyn_bergere_86x100xh94.3ds"
PRODUCT_TOKENS = [
    "MARIPBMN86",
    "MARIPBML86",
    "Baxter_Marilyn_Armchair_2D_3D.zip",
    "BERGÈRE ARMCHAIR WITH MATTE PAINTED METAL FRAME AND SWIVEL BASE",
    "BERGÈRE ARMCHAIR WITH GLOSSY LACQUERED METAL FRAME AND SWIVEL BASE",
    "<td>86</td><td>100</td><td>94</td><td>cm</td>",
]


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, target: Path) -> tuple[int, str | None, int]:
    request = urllib.request.Request(url, headers={"User-Agent": "Codex source verifier/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response, target.open("wb") as destination:
        shutil.copyfileobj(response, destination)
        return response.status, response.headers.get("Content-Type"), target.stat().st_size


def render_page(pdftoppm: str, source: Path, prefix: Path, page: int) -> Path:
    subprocess.run(
        [pdftoppm, "-f", str(page), "-l", str(page), "-singlefile", "-png", "-r", "180", str(source), str(prefix)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return prefix.with_suffix(".png")


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
        "technical_preview": SOURCE_DIR / "Baxter_Marilyn_technical-sheet-page-13-preview.png",
        "matte_svg": SOURCE_DIR / "MARIPBMN86.svg",
        "glossy_svg": SOURCE_DIR / "MARIPBML86.svg",
        "dwg": SOURCE_DIR / "Marilyn_Abaco.dwg",
        "bergere_3ds": SOURCE_DIR / "Marilyn_bergere_86x100xh94.3ds",
    }
    with tempfile.TemporaryDirectory(prefix="baxter-marilyn01-source-") as temporary:
        temporary_dir = Path(temporary)
        downloaded = {}
        for key, url, filename in (
            ("product_page", PRODUCT_PAGE_URL, "product-page.html"),
            ("zip", ZIP_URL, "package.zip"),
            ("technical_pdf", TECHNICAL_URL, "technical.pdf"),
            ("matte_svg", MATTE_SVG_URL, "matte.svg"),
            ("glossy_svg", GLOSSY_SVG_URL, "glossy.svg"),
        ):
            path = temporary_dir / filename
            status, content_type, byte_count = download(url, path)
            downloaded[key] = {
                "url": url,
                "path": path,
                "http_status": status,
                "content_type": content_type,
                "byte_count": byte_count,
                "sha256": sha256_file(path),
            }

        current_html = downloaded["product_page"]["path"].read_text(encoding="utf-8")
        product_tokens = {token: token in current_html for token in PRODUCT_TOKENS}
        product_page = {
            **{key: value for key, value in downloaded["product_page"].items() if key != "path"},
            "archived_path": str(local["product_page"].relative_to(ROOT)),
            "archived_sha256": sha256_file(local["product_page"]),
            "expected_archived_sha256": EXPECTED["archived_product_page"],
            "current_page_is_dynamic_and_not_required_to_be_byte_identical_to_archive": True,
            "required_identity_tokens": product_tokens,
            "all_required_identity_tokens_found": all(product_tokens.values()),
            "pass": downloaded["product_page"]["http_status"] == 200
            and sha256_file(local["product_page"]) == EXPECTED["archived_product_page"]
            and all(product_tokens.values()),
        }

        zip_path = downloaded["zip"]["path"]
        with zipfile.ZipFile(zip_path) as archive:
            member_names = archive.namelist()
            dwg_bytes = archive.read(DWG_MEMBER)
            bergere_bytes = archive.read(BERGERE_MEMBER)
        zip_bytes_match = zip_path.read_bytes() == local["zip"].read_bytes()
        native_zip = {
            **{key: value for key, value in downloaded["zip"].items() if key != "path"},
            "local_path": str(local["zip"].relative_to(ROOT)),
            "local_sha256": sha256_file(local["zip"]),
            "expected_sha256": EXPECTED["zip"],
            "downloaded_bytes_match_local_archive": zip_bytes_match,
            "member_count": len(member_names),
            "dwg_member": DWG_MEMBER,
            "dwg_member_sha256": sha256_bytes(dwg_bytes),
            "dwg_member_matches_local_extraction": dwg_bytes == local["dwg"].read_bytes(),
            "bergere_3ds_member": BERGERE_MEMBER,
            "bergere_3ds_member_sha256": sha256_bytes(bergere_bytes),
            "bergere_3ds_member_matches_local_extraction": bergere_bytes == local["bergere_3ds"].read_bytes(),
            "pass": downloaded["zip"]["http_status"] == 200
            and downloaded["zip"]["sha256"] == EXPECTED["zip"]
            and sha256_file(local["zip"]) == EXPECTED["zip"]
            and zip_bytes_match
            and sha256_bytes(dwg_bytes) == EXPECTED["dwg"]
            and dwg_bytes == local["dwg"].read_bytes()
            and sha256_bytes(bergere_bytes) == EXPECTED["bergere_3ds"]
            and bergere_bytes == local["bergere_3ds"].read_bytes(),
        }

        current_render = render_page(
            pdftoppm, downloaded["technical_pdf"]["path"], temporary_dir / "current-page-13", 13
        )
        archived_render = render_page(
            pdftoppm, local["technical_pdf"], temporary_dir / "archived-page-13", 13
        )
        render_bytes_match = current_render.read_bytes() == archived_render.read_bytes()
        technical_pdf = {
            **{key: value for key, value in downloaded["technical_pdf"].items() if key != "path"},
            "archived_path": str(local["technical_pdf"].relative_to(ROOT)),
            "archived_sha256": sha256_file(local["technical_pdf"]),
            "expected_archived_sha256": EXPECTED["archived_technical_pdf"],
            "current_pdf_is_server_regenerated_and_not_required_to_be_byte_identical_to_archive": True,
            "identity_page": 13,
            "render_dpi": 180,
            "expected_identity_page_render_sha256": EXPECTED["technical_page_13_render"],
            "current_identity_page_render_sha256": sha256_file(current_render),
            "archived_identity_page_render_sha256": sha256_file(archived_render),
            "stored_preview_path": str(local["technical_preview"].relative_to(ROOT)),
            "stored_preview_sha256": sha256_file(local["technical_preview"]),
            "current_identity_page_render_matches_archive": render_bytes_match,
            "pass": downloaded["technical_pdf"]["http_status"] == 200
            and sha256_file(local["technical_pdf"]) == EXPECTED["archived_technical_pdf"]
            and sha256_file(current_render) == EXPECTED["technical_page_13_render"]
            and sha256_file(archived_render) == EXPECTED["technical_page_13_render"]
            and sha256_file(local["technical_preview"]) == EXPECTED["technical_page_13_render"]
            and render_bytes_match,
        }

        measurement_svgs = {}
        for key, expected_key in (("matte_svg", "matte_svg"), ("glossy_svg", "glossy_svg")):
            bytes_match = downloaded[key]["path"].read_bytes() == local[key].read_bytes()
            measurement_svgs[key] = {
                **{field: value for field, value in downloaded[key].items() if field != "path"},
                "local_path": str(local[key].relative_to(ROOT)),
                "local_sha256": sha256_file(local[key]),
                "expected_sha256": EXPECTED[expected_key],
                "downloaded_bytes_match_local_archive": bytes_match,
                "pass": downloaded[key]["http_status"] == 200
                and downloaded[key]["sha256"] == EXPECTED[expected_key]
                and sha256_file(local[key]) == EXPECTED[expected_key]
                and bytes_match,
            }

    report = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Baxter",
        "family": "Marilyn",
        "project_ifc_type_name": "Marilyn 01",
        "resolved_variant": "bergere armchair with swivel base, 86 x 100 x 94 cm",
        "marilyn_02_geometry_used": False,
        "method": "fresh official web downloads; byte comparison for ZIP and exact measurement SVGs; ZIP-member hash checks for DWG and bergere 3DS; page-13 raster identity for the dynamically regenerated technical PDF",
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
        raise SystemExit("Baxter Marilyn 01 official-source revalidation failed")


if __name__ == "__main__":
    main()
