#!/usr/bin/env python3
"""Revalidate public Geberit 154.154.00.1 identity and CAD availability."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARTICLE = "154.154.00.1"
PRODUCT_PAGE_URL = "https://catalog.geberit.co.uk/en-GB/product/PRO_199058"
ARTICLE_API_URL = (
    "https://api.prd.mbp.geberit.com/product-data/v3/articles/"
    f"{ARTICLE}?locale=en-GB&brand=GEBERIT"
)
SOURCE_DIR = ROOT / "output/review/highpoly-types/geberit-154-154-00/official-source"
ARCHIVED_PAGE = SOURCE_DIR / "geberit-154-154-00-1-product-page.html"
OUTPUT = SOURCE_DIR / "official-source-revalidation.json"
EXPECTED_EPS = {
    "perspective": (
        "https://cdn.data.geberit.com/eps/99/00/DAS_199900.eps",
        "DAS_199900-perspective.eps",
        "cee509e98b933d5857df6a5b7179eca7abbd79d859c9913132818d5f28031c84",
    ),
    "front": (
        "https://cdn.data.geberit.com/eps/99/02/DAS_199902.eps",
        "DAS_199902-front-view.eps",
        "8114001caff7c8552563c219b1b3923340db00a18c8de174ce9b061d409c10ac",
    ),
    "top": (
        "https://cdn.data.geberit.com/eps/99/04/DAS_199904.eps",
        "DAS_199904-top-view.eps",
        "3c238df45acfd43d819565106cebe0133049b6500fc012b6a4a31cbb05137291",
    ),
}


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def request_bytes(url: str) -> tuple[int, str | None, bytes]:
    request = urllib.request.Request(url, headers={"User-Agent": "Codex source verifier/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=90) as response:
            return response.status, response.headers.get("Content-Type"), response.read()
    except urllib.error.HTTPError as error:
        return error.code, error.headers.get("Content-Type"), error.read()


def page_has_exact_article_without_cad(page: bytes) -> bool:
    text = page.decode("utf-8", errors="replace")
    pattern = (
        rf'\\"id\\":\\"{re.escape(ARTICLE)}\\"'
        r'.{0,2000}?\\"cadDrawings\\":\\"\$undefined\\"'
    )
    return re.search(pattern, text, re.DOTALL) is not None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    page_status, page_content_type, current_page = request_bytes(PRODUCT_PAGE_URL)
    archived_page = ARCHIVED_PAGE.read_bytes()
    page = {
        "url": PRODUCT_PAGE_URL,
        "http_status": page_status,
        "content_type": page_content_type,
        "downloaded_byte_count": len(current_page),
        "downloaded_sha256": sha256_bytes(current_page),
        "archived_local_path": str(ARCHIVED_PAGE.relative_to(ROOT)),
        "archived_byte_count": len(archived_page),
        "archived_sha256": sha256_bytes(archived_page),
        "downloaded_exact_article_with_cadDrawings_undefined": page_has_exact_article_without_cad(current_page),
        "archived_exact_article_with_cadDrawings_undefined": page_has_exact_article_without_cad(archived_page),
    }
    page["pass"] = (
        page_status == 200
        and page_content_type is not None
        and page_content_type.startswith("text/html")
        and page["downloaded_exact_article_with_cadDrawings_undefined"]
        and page["archived_exact_article_with_cadDrawings_undefined"]
    )

    api_status, api_content_type, api_body = request_bytes(ARTICLE_API_URL)
    api = {
        "url": ARTICLE_API_URL,
        "http_status": api_status,
        "content_type": api_content_type,
        "response_byte_count": len(api_body),
        "response_sha256": sha256_bytes(api_body),
        "anonymous_access_available": api_status == 200,
        "note": "The exact article API currently requires authorization; this status is not used to infer CAD availability.",
        "pass": api_status in (200, 401, 403),
    }

    dwg_results = {}
    for code in ("A", "G", "L", "P"):
        url = f"https://cdn.data.geberit.com/cad/{ARTICLE}_{code}.dwg"
        status, content_type, body = request_bytes(url)
        dwg_results[code] = {
            "url": url,
            "http_status": status,
            "content_type": content_type,
            "response_byte_count": len(body),
            "response_sha256": sha256_bytes(body),
            "native_dwg_acquired": status == 200,
            "pass": status == 404,
        }

    eps_results = {}
    for view, (url, filename, expected_hash) in EXPECTED_EPS.items():
        status, content_type, downloaded = request_bytes(url)
        local_path = SOURCE_DIR / filename
        local = local_path.read_bytes()
        downloaded_hash = sha256_bytes(downloaded)
        local_hash = sha256_bytes(local)
        eps_results[view] = {
            "url": url,
            "local_path": str(local_path.relative_to(ROOT)),
            "http_status": status,
            "content_type": content_type,
            "expected_sha256": expected_hash,
            "downloaded_sha256": downloaded_hash,
            "local_sha256": local_hash,
            "downloaded_bytes_match_local_archive": downloaded == local,
            "role": "official vector identity and dimension evidence only; not CAD geometry",
            "pass": status == 200
            and downloaded_hash == expected_hash
            and local_hash == expected_hash
            and downloaded == local,
        }

    report = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Geberit",
        "family": "CleanLine shower channel installation set",
        "article_number": ARTICLE,
        "method": "fresh official product-page, article-API, native-DWG and EPS requests with exact identity and SHA-256 checks",
        "product_page": page,
        "article_api": api,
        "native_dwg_results": dwg_results,
        "official_eps_results": eps_results,
        "native_dwg_acquired": any(item["native_dwg_acquired"] for item in dwg_results.values()),
        "official_eps_used_as_cad_geometry": False,
        "third_party_cad_used": False,
        "drawing_geometry_source": {
            "source_kind": "geometry_derived_simplified_proxy",
            "source_label_zh": "基于原始高模几何生成的简化图纸表达",
        },
        "scope": "exact manufacturer article identity and dimensions; not a project shop drawing and not official CAD geometry",
        "pass": page["pass"]
        and api["pass"]
        and all(item["pass"] for item in dwg_results.values())
        and all(item["pass"] for item in eps_results.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not report["pass"]:
        raise SystemExit("Geberit 154.154.00.1 official-source revalidation failed")


if __name__ == "__main__":
    main()
