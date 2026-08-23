#!/usr/bin/env python3
"""Archive and revalidate exact public Gessi 54093 identity and native CAD."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
import zipfile
from datetime import datetime, timezone
from io import BytesIO
from pathlib import Path

from falper_sorgente_linework import ROOT, relative, sha256, write_json


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54093"
SOURCE_DIR = PRODUCT_DIR / "official-source"
OUTPUT = SOURCE_DIR / "official-source-revalidation.json"
API_ROOT = "https://g-ecatalogue-be-prod-we.azurewebsites.net"
DETAIL_URL = API_ROOT + "/public/product/GetProductDetails?country=it&language=en&productCode=54093"
ATTACHMENT_URL = API_ROOT + "/public/product/GetProductAttachments?country=it&language=en&productCode=PF54093"
AREA_PRO_PRODUCT = "https://areapro.gessi.com/en/product/54093"
DWG_ZIP_URL = "https://gessistorage.blob.core.windows.net/zwa/GPF5409300000G000_arc.zip"
PDF_URL = "https://gessistorage.blob.core.windows.net/zc4/GPF5409300000G000_1.pdf"
EXPECTED = {
    "zip": "615d4747639b6acfa1daca28a894124fef26e74239de89ef68455b92db31d0f9",
    "dwg": "a9308fc8498d34c8cf2f68fd28aaf90b11e62bfc59424e0a5a3ac7f0d28d5047",
    "pdf": "2eb67af1be75093a12ab8bc3a739e5bdd6e4b1c4c544948795c87cb9598515c6",
}
FILES = {
    "zip": SOURCE_DIR / "GPF5409300000G000_arc.zip",
    "dwg": SOURCE_DIR / "GPF5409300000G000_3.dwg",
    "pdf": SOURCE_DIR / "GPF5409300000G000_1.pdf",
}


def request_bytes(url: str) -> tuple[bytes, dict]:
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 Codex evidence verifier"})
    with urllib.request.urlopen(request, timeout=45) as response:
        body = response.read()
        return body, {
            "url": url,
            "http_status": response.status,
            "content_type": response.headers.get_content_type(),
            "bytes": len(body),
            "sha256": hashlib.sha256(body).hexdigest(),
        }


def select_attachment(payload: dict, original_type: str, name: str, expected_url: str) -> dict:
    matches = [
        item for item in payload["data"]["attachments"]
        if item.get("attachmentTypeOriginal") == original_type
        and item.get("attachmentName") == name
        and item.get("downloadVisible") is True
    ]
    urls = sorted({item["attachmentUrl"] for item in matches})
    if not matches or urls != [expected_url]:
        raise RuntimeError(f"exact 54093 G000 attachment drift: {original_type} / {name} / {urls}")
    return {
        "attachment_type": matches[0]["attachmentType"],
        "attachment_type_original": original_type,
        "attachment_name": name,
        "attachment_url": expected_url,
        "download_visible": True,
        "matching_configured_codes": sorted({item["configuredCode"] for item in matches}),
    }


def preserve(path: Path, body: bytes, expected_hash: str) -> None:
    actual = hashlib.sha256(body).hexdigest()
    if actual != expected_hash:
        raise RuntimeError(f"remote official source hash drift: {path.name} / {actual}")
    if path.is_file():
        if sha256(path) != expected_hash:
            raise RuntimeError(f"refusing to overwrite mismatched archived source: {path}")
        return
    path.write_bytes(body)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    detail_body, detail_request = request_bytes(DETAIL_URL)
    attachment_body, attachment_request = request_bytes(ATTACHMENT_URL)
    details, attachments = json.loads(detail_body), json.loads(attachment_body)
    product = details.get("data", {}).get("product", {})
    if (
        details.get("success") is not True
        or attachments.get("success") is not True
        or product.get("pimId") != "PF54093"
        or product.get("productId") != "54093"
        or product.get("productName") != "BASIN SPOUT"
        or product.get("productDescription") != "High spout basin"
    ):
        raise RuntimeError("official exact Gessi 54093 identity gate failed")
    configured_codes = sorted({item["configuredCode"] for item in product["productsConfigured"]})
    if not configured_codes or any(not code.startswith("GPF540930") or not code.endswith("G000") for code in configured_codes):
        raise RuntimeError("official exact Gessi 54093 configured G000 family gate failed")
    selected = {
        "native_2d_dwg_zip": select_attachment(attachments, "ZWA", "GPF5409300000G000_arc.zip", DWG_ZIP_URL),
        "technical_drawing_mm_pdf": select_attachment(attachments, "ZC4", "GPF5409300000G000_1.pdf", PDF_URL),
    }
    zip_body, zip_request = request_bytes(DWG_ZIP_URL)
    pdf_body, pdf_request = request_bytes(PDF_URL)
    preserve(FILES["zip"], zip_body, EXPECTED["zip"])
    preserve(FILES["pdf"], pdf_body, EXPECTED["pdf"])
    with zipfile.ZipFile(BytesIO(zip_body)) as archive:
        if archive.namelist() != ["GPF5409300000G000_3.dwg"]:
            raise RuntimeError(f"unexpected exact 54093 ZIP members: {archive.namelist()}")
        dwg_body = archive.read("GPF5409300000G000_3.dwg")
    preserve(FILES["dwg"], dwg_body, EXPECTED["dwg"])
    write_json(SOURCE_DIR / "gessi316-54093-product-details.json", details)
    write_json(SOURCE_DIR / "gessi316-54093-attachments.json", attachments)
    technical = {
        item["technicalName"]: item["technicalValue"]
        for item in product["technical"] if item["technicalName"] in {"Width", "Height", "Depth"}
    }
    member_hash = hashlib.sha256(dwg_body).hexdigest()
    payload = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Gessi",
        "family": "Gessi316 Meccanica",
        "resolved_article_number": "54093",
        "public_access": {
            "area_pro_product_page": AREA_PRO_PRODUCT,
            "anonymous_public_api_used": True,
            "login_or_registration_performed": False,
            "exact_54093_g000_native_dwg_publicly_downloadable": True,
        },
        "identity": {
            "pim_id": product["pimId"],
            "product_id": product["productId"],
            "product_name": product["productName"],
            "product_description": product["productDescription"],
            "technical_dimensions": technical,
            "configured_codes": configured_codes,
        },
        "api_requests": {"details": detail_request, "attachments": attachment_request},
        "selected_attachments": selected,
        "downloads": {
            "native_dwg_zip": {**zip_request, "path": relative(FILES["zip"]), "expected_sha256": EXPECTED["zip"], "pass": zip_request["sha256"] == EXPECTED["zip"]},
            "technical_pdf": {**pdf_request, "path": relative(FILES["pdf"]), "expected_sha256": EXPECTED["pdf"], "pass": pdf_request["sha256"] == EXPECTED["pdf"]},
        },
        "zip_member": {
            "zip_path": relative(FILES["zip"]), "zip_sha256": sha256(FILES["zip"]),
            "member": "GPF5409300000G000_3.dwg", "member_bytes": len(dwg_body),
            "member_sha256": member_hash, "extracted_path": relative(FILES["dwg"]),
            "extracted_sha256": sha256(FILES["dwg"]),
            "member_matches_extracted_file": member_hash == sha256(FILES["dwg"]),
            "pass": member_hash == sha256(FILES["dwg"]) == EXPECTED["dwg"],
        },
        "drawing_geometry_policy": {
            "native_dwg_used_for_three_views": "54093 G000 only",
            "g001_or_a004_variant_cad_used": False,
            "adjacent_gessi_product_cad_used": False,
            "third_party_cad_used": False,
        },
        "pass": sha256(FILES["zip"]) == EXPECTED["zip"] and sha256(FILES["dwg"]) == EXPECTED["dwg"] and sha256(FILES["pdf"]) == EXPECTED["pdf"],
    }
    write_json(args.output.resolve(), payload)
    if not payload["pass"]:
        raise RuntimeError("Gessi 54093 official source revalidation failed")
    print(json.dumps({"output": relative(args.output), "pass": True}, indent=2))


if __name__ == "__main__":
    main()
