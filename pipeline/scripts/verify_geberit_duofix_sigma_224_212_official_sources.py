#!/usr/bin/env python3
"""Revalidate exact Duofix 224.212.00.2 DWGs and catalogue identity evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARTICLE = "224.212.00.2"
OUTPUT_DIR = ROOT / "output/review/highpoly-types/geberit-duofix-sigma-224-212"
SOURCE_DIR = OUTPUT_DIR / "official-source"
OUTPUT = SOURCE_DIR / "official-source-revalidation.json"
CATALOGUE_URL = "https://cdn-geberit-country-sg.prod.web.geberit.com/_assets/local-media/brochures/2025-nsea-sanitary-and-bathroom-catalogue-web.pdf"
CATALOGUE_EXTRACT = SOURCE_DIR / "geberit-224212-catalogue-page11.pdf"
CATALOGUE_EXTRACT_SHA256 = "8c3a17569e6fa53d9a769ab584fed79a865118c4a27018facd619b4a5cc5cdd3"
CATALOGUE_PAGE_PDF_INDEX_ONE_BASED = 10
EXPECTED_PAGE_RENDER_SHA256 = "757450d7fb7a85d61755fe61dce34e77d14ecc33e7522da056d8f97db3fee660"
EXPECTED_DWG = {
    "A": "f06952ac7bf91f0a0037ae85644df3df0938d549322cbfb59b7dfd9be01939d4",
    "G": "841dd3c72d594cfd2c7f921b84f10a2305ab9957b531771662fb7f03eff97207",
    "L": "e6a26c660ae9cd4fc8cb4b1effadcea38a4794b66523f0fea2c02c5e02b77f25",
    "P": "5c92a663de35357a4bf0ff7eb0838bc5e178030afe7d769ead5632f50821652f",
}


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
        [
            pdftoppm,
            "-f", str(page),
            "-l", str(page),
            "-singlefile",
            "-png",
            "-r", "180",
            str(source),
            str(prefix),
        ],
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
        raise SystemExit("pdftoppm is required for catalogue page identity verification")
    with tempfile.TemporaryDirectory(prefix="geberit-duofix-224212-source-") as temporary:
        temporary_dir = Path(temporary)
        dwg_results = {}
        for code, expected_hash in EXPECTED_DWG.items():
            filename = f"{ARTICLE}_{code}.dwg"
            url = f"https://cdn.data.geberit.com/cad/{filename}"
            downloaded_path = temporary_dir / filename
            status, content_type, downloaded_size = download(url, downloaded_path)
            local_path = SOURCE_DIR / filename
            downloaded_hash = sha256_file(downloaded_path)
            local_hash = sha256_file(local_path)
            bytes_match = downloaded_path.read_bytes() == local_path.read_bytes()
            dwg_results[code] = {
                "url": url,
                "local_path": str(local_path.relative_to(ROOT)),
                "http_status": status,
                "content_type": content_type,
                "downloaded_byte_count": downloaded_size,
                "local_byte_count": local_path.stat().st_size,
                "expected_sha256": expected_hash,
                "downloaded_sha256": downloaded_hash,
                "local_sha256": local_hash,
                "downloaded_bytes_match_local_archive": bytes_match,
                "pass": status == 200
                and downloaded_hash == expected_hash
                and local_hash == expected_hash
                and bytes_match,
            }

        downloaded_catalogue = temporary_dir / "current-official-catalogue.pdf"
        catalogue_status, catalogue_content_type, catalogue_size = download(
            CATALOGUE_URL, downloaded_catalogue
        )
        downloaded_render = render_page(
            pdftoppm,
            downloaded_catalogue,
            temporary_dir / "downloaded-article-page",
            CATALOGUE_PAGE_PDF_INDEX_ONE_BASED,
        )
        local_render = render_page(
            pdftoppm,
            CATALOGUE_EXTRACT,
            temporary_dir / "local-article-page",
            1,
        )
        downloaded_render_hash = sha256_file(downloaded_render)
        local_render_hash = sha256_file(local_render)
        catalogue_page_bytes_match = downloaded_render.read_bytes() == local_render.read_bytes()
        catalogue = {
            "url": CATALOGUE_URL,
            "http_status": catalogue_status,
            "content_type": catalogue_content_type,
            "downloaded_byte_count": catalogue_size,
            "downloaded_pdf_sha256": sha256_file(downloaded_catalogue),
            "article_page_pdf_index_one_based": CATALOGUE_PAGE_PDF_INDEX_ONE_BASED,
            "printed_page_number": 11,
            "local_extract_path": str(CATALOGUE_EXTRACT.relative_to(ROOT)),
            "local_extract_sha256": sha256_file(CATALOGUE_EXTRACT),
            "expected_local_extract_sha256": CATALOGUE_EXTRACT_SHA256,
            "render_dpi": 180,
            "expected_article_page_render_sha256": EXPECTED_PAGE_RENDER_SHA256,
            "downloaded_article_page_render_sha256": downloaded_render_hash,
            "local_extract_render_sha256": local_render_hash,
            "downloaded_article_page_render_matches_local_extract": catalogue_page_bytes_match,
            "identity_statement": "The uniquely located catalogue article row is 224.212.00.2 with B=50 cm, H=112 cm and T=12 cm.",
            "pass": catalogue_status == 200
            and catalogue_content_type == "application/pdf"
            and sha256_file(CATALOGUE_EXTRACT) == CATALOGUE_EXTRACT_SHA256
            and downloaded_render_hash == EXPECTED_PAGE_RENDER_SHA256
            and local_render_hash == EXPECTED_PAGE_RENDER_SHA256
            and catalogue_page_bytes_match,
        }

    report = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Geberit",
        "article_number": ARTICLE,
        "method": "fresh official CDN downloads with SHA-256 and byte comparison; current full official catalogue article page rendered at 180 dpi and compared byte-for-byte with the archived extract",
        "drawing_view_mapping": {"plan": "G", "front": "A", "side": "L"},
        "identity_only_code": "P",
        "dwg_results": dwg_results,
        "catalogue": catalogue,
        "all_declared_dwg_urls_accessible": all(item["http_status"] == 200 for item in dwg_results.values()),
        "all_downloaded_dwg_bytes_match_local_archive": all(
            item["downloaded_bytes_match_local_archive"] for item in dwg_results.values()
        ),
        "pass": all(item["pass"] for item in dwg_results.values()) and catalogue["pass"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not report["pass"]:
        raise SystemExit("Geberit Duofix official-source revalidation failed")


if __name__ == "__main__":
    main()
