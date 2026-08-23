#!/usr/bin/env python3
"""Revalidate exact public Gessi 54294/45089 identity and CAD attachments."""

from __future__ import annotations

import argparse
import io
import json
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from falper_sorgente_linework import ROOT, relative, sha256, write_json


PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
SOURCE_DIR = PRODUCT_DIR / "official-source"
OUTPUT = SOURCE_DIR / "official-source-revalidation.json"
API_ROOT = "https://g-ecatalogue-be-prod-we.azurewebsites.net"
DETAIL_URL = API_ROOT + "/public/product/GetProductDetails?country=it&language=en&productCode={code}"
ATTACHMENT_URL = API_ROOT + "/public/product/GetProductAttachments?country=it&language=en&productCode={pim_id}"
CATALOGUE_PAGE = "https://www.gessi.com/us/catalogs"
AREA_PRO_PRODUCT = "https://areapro.gessi.com/en/product/54294"
EXPECTED_DESCRIPTION = "External parts three-holes basin mixer with long spout, without waste."
EXPECTED = {
    "54294_zip": "fad0d98e94a83bb488b5f0703472862d628759f21700c1a87c3af343f4c1bdb9",
    "54294_dwg": "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4",
    "54294_pdf": "82058bb27752f27e6fc92029783d67c15fd443befcd393f7d572fd0147a581e3",
    "45089_zip": "b1dfbd7b379f359f969be51193d9bcc61f24a354a9c7935d72011170ff60fc30",
    "45089_dwg": "861f9a73ce4e454e3fec753c2495bcbe0e3d9f320d2b110d678790426ebd2e6c",
}
FILES = {
    "54294_zip": SOURCE_DIR / "GPF5429400000G000_arc.zip",
    "54294_dwg": SOURCE_DIR / "GPF5429400000G000_3.dwg",
    "54294_pdf": SOURCE_DIR / "GPF5429400000G000_1.pdf",
    "45089_zip": SOURCE_DIR / "GPF4508920000G000_arc.zip",
    "45089_dwg": SOURCE_DIR / "GPF4508920000G000_4.dwg",
}
URLS = {
    "54294_zip": "https://gessistorage.blob.core.windows.net/zwa/GPF5429400000G000_arc.zip",
    "54294_pdf": "https://gessistorage.blob.core.windows.net/zc4/GPF5429400000G000_1.pdf",
    "45089_zip": "https://gessistorage.blob.core.windows.net/zwa/GPF4508920000G000_arc.zip",
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
        }


def request_json(url: str) -> tuple[dict, dict]:
    body, metadata = request_bytes(url)
    payload = json.loads(body)
    metadata["sha256"] = __import__("hashlib").sha256(body).hexdigest()
    return payload, metadata


def attachment(payload: dict, *, original_type: str, name: str) -> dict:
    matches = [
        item for item in payload["data"]["attachments"]
        if item.get("attachmentTypeOriginal") == original_type
        and item.get("attachmentName") == name
        and item.get("downloadVisible") is True
    ]
    if not matches:
        raise RuntimeError(f"missing public exact attachment {original_type} / {name}")
    urls = sorted({item["attachmentUrl"] for item in matches})
    if len(urls) != 1:
        raise RuntimeError(f"attachment URL drift for {name}: {urls}")
    return {
        "attachment_type": matches[0]["attachmentType"],
        "attachment_type_original": original_type,
        "attachment_name": name,
        "attachment_url": urls[0],
        "download_visible": True,
        "matching_configured_codes": sorted({item.get("configuredCode") for item in matches if item.get("configuredCode")}),
    }


def verify_download(key: str, url: str) -> dict:
    body, metadata = request_bytes(url)
    downloaded_hash = __import__("hashlib").sha256(body).hexdigest()
    local_hash = sha256(FILES[key])
    expected_hash = EXPECTED[key]
    passed = metadata["http_status"] == 200 and downloaded_hash == local_hash == expected_hash
    return {
        **metadata,
        "path": relative(FILES[key]),
        "expected_sha256": expected_hash,
        "downloaded_sha256": downloaded_hash,
        "local_sha256": local_hash,
        "downloaded_bytes_match_local_archive": downloaded_hash == local_hash,
        "pass": passed,
    }


def verify_zip(key: str, member_key: str, member_name: str) -> dict:
    with zipfile.ZipFile(FILES[key]) as archive:
        names = archive.namelist()
        if names != [member_name]:
            raise RuntimeError(f"unexpected {key} members: {names}")
        member = archive.read(member_name)
    member_hash = __import__("hashlib").sha256(member).hexdigest()
    return {
        "zip_path": relative(FILES[key]),
        "zip_sha256": sha256(FILES[key]),
        "member": member_name,
        "member_bytes": len(member),
        "member_sha256": member_hash,
        "extracted_path": relative(FILES[member_key]),
        "extracted_sha256": sha256(FILES[member_key]),
        "member_matches_extracted_file": member_hash == sha256(FILES[member_key]),
        "pass": member_hash == sha256(FILES[member_key]) == EXPECTED[member_key],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    for key, path in FILES.items():
        if not path.is_file() or sha256(path) != EXPECTED[key]:
            raise RuntimeError(f"local Gessi evidence hash mismatch: {key}")

    details = {}
    attachments = {}
    request_records = {}
    for code, pim_id in (("54294", "PF54294"), ("45089", "PF45089")):
        detail, detail_request = request_json(DETAIL_URL.format(code=code))
        if detail.get("success") is not True or detail.get("data", {}).get("product", {}).get("pimId") != pim_id:
            raise RuntimeError(f"official exact product identity gate failed for {code}")
        if detail["data"]["product"].get("productId") != code:
            raise RuntimeError(f"official product code drift for {code}")
        attachment_payload, attachment_request = request_json(ATTACHMENT_URL.format(pim_id=pim_id))
        if attachment_payload.get("success") is not True:
            raise RuntimeError(f"official public attachment API failed for {code}")
        details[code] = detail
        attachments[code] = attachment_payload
        request_records[f"{code}_details"] = detail_request
        request_records[f"{code}_attachments"] = attachment_request
        write_json(SOURCE_DIR / f"gessi316-{code}-product-details.json", detail)
        write_json(SOURCE_DIR / f"gessi316-{code}-attachments.json", attachment_payload)

    product_54294 = details["54294"]["data"]["product"]
    product_45089 = details["45089"]["data"]["product"]
    if product_54294.get("productDescription") != EXPECTED_DESCRIPTION:
        raise RuntimeError("54294 official description drifted")
    if product_45089.get("productDescription") != "Watertight built-in part for three-hole basin mixer":
        raise RuntimeError("45089 official built-in part identity drifted")

    selected_attachments = {
        "54294_native_2d_dwg_zip": attachment(
            attachments["54294"], original_type="ZWA", name="GPF5429400000G000_arc.zip"
        ),
        "54294_technical_drawing_mm_pdf": attachment(
            attachments["54294"], original_type="ZC4", name="GPF5429400000G000_1.pdf"
        ),
        "45089_native_2d_dwg_zip": attachment(
            attachments["45089"], original_type="ZWA", name="GPF4508920000G000_arc.zip"
        ),
    }
    if selected_attachments["54294_native_2d_dwg_zip"]["attachment_url"] != URLS["54294_zip"]:
        raise RuntimeError("54294 exact official DWG URL drifted")
    if selected_attachments["54294_technical_drawing_mm_pdf"]["attachment_url"] != URLS["54294_pdf"]:
        raise RuntimeError("54294 technical PDF URL drifted")
    if selected_attachments["45089_native_2d_dwg_zip"]["attachment_url"] != URLS["45089_zip"]:
        raise RuntimeError("45089 companion DWG URL drifted")

    downloads = {key: verify_download(key, url) for key, url in URLS.items()}
    zip_members = {
        "54294": verify_zip("54294_zip", "54294_dwg", "GPF5429400000G000_3.dwg"),
        "45089": verify_zip("45089_zip", "45089_dwg", "GPF4508920000G000_4.dwg"),
    }
    payload = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Gessi",
        "family": "Gessi316 Meccanica",
        "resolved_article_number": "45089_54294",
        "external_visible_product_code": "54294",
        "companion_built_in_product_code": "45089",
        "public_access": {
            "catalogue_page": CATALOGUE_PAGE,
            "area_pro_product_page": AREA_PRO_PRODUCT,
            "anonymous_public_api_used": True,
            "login_or_registration_performed": False,
            "exact_54294_native_dwg_publicly_downloadable": True,
        },
        "identity": {
            "54294": {
                "pim_id": product_54294["pimId"],
                "product_name": product_54294["productName"],
                "product_description": product_54294["productDescription"],
                "technical_dimensions": {
                    item["technicalName"]: item["technicalValue"]
                    for item in product_54294["technical"]
                    if item["technicalName"] in {"Width", "Height", "Depth"}
                },
            },
            "45089": {
                "pim_id": product_45089["pimId"],
                "product_name": product_45089["productName"],
                "product_description": product_45089["productDescription"],
            },
        },
        "api_requests": request_records,
        "selected_attachments": selected_attachments,
        "downloads": downloads,
        "zip_members": zip_members,
        "drawing_geometry_policy": {
            "native_dwg_used_for_three_views": "54294 only",
            "companion_45089_dwg_used_as_54294_geometry": False,
            "adjacent_gessi_product_cad_used": False,
            "third_party_cad_used": False,
        },
        "all_downloaded_bytes_match_local_archive": all(item["pass"] for item in downloads.values()),
        "all_zip_members_match_extracted_files": all(item["pass"] for item in zip_members.values()),
        "pass": all(item["pass"] for item in downloads.values()) and all(item["pass"] for item in zip_members.values()),
    }
    write_json(args.output.resolve(), payload)
    if not payload["pass"]:
        raise RuntimeError("Gessi official source revalidation failed")
    print(json.dumps({"output": relative(args.output), "pass": True}, indent=2))


if __name__ == "__main__":
    main()
