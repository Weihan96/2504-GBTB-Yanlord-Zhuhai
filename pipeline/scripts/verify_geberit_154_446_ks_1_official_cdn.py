#!/usr/bin/env python3
"""Re-download exact Geberit 154.446.KS.1 DWGs and compare them with the archive."""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARTICLE = "154.446.KS.1"
OUTPUT = ROOT / "output/review/highpoly-types/geberit-154-446-ks-1/official-source/official-cdn-revalidation.json"
SOURCE_DIR = ROOT / "output/review/highpoly-types/geberit-154-446-ks-1/official-source"
EXPECTED = {
    "A": "de09a6bbf99ecbe3c832cbee7ac9d398e95b927246f274643b513bfe960586c3",
    "G": "bc41f8db7c7989de9d9089374e03c3b9deea94e2634bbb0f9418cad982b4b6e4",
    "L": "dedb47964cc69310bd9af3982379c428c8b20494b7588296aedc5d07380e3421",
    "P": "46e1d9fe9b7afda6f4a72970a5e35d7bc0935c2bdc890ccca232f0b0896a8e82",
}


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    results = {}
    with tempfile.TemporaryDirectory(prefix="geberit-154446ks1-cdn-") as temporary:
        temporary_dir = Path(temporary)
        for code, expected_hash in EXPECTED.items():
            filename = f"{ARTICLE}_{code}.dwg"
            url = f"https://cdn.data.geberit.com/cad/{filename}"
            request = urllib.request.Request(url, headers={"User-Agent": "Codex source verifier/1.0"})
            with urllib.request.urlopen(request, timeout=60) as response:
                downloaded = response.read()
                status = response.status
                content_type = response.headers.get("Content-Type")
            downloaded_path = temporary_dir / filename
            downloaded_path.write_bytes(downloaded)
            local_path = SOURCE_DIR / filename
            downloaded_hash = sha256_file(downloaded_path)
            local_hash = sha256_file(local_path)
            results[code] = {
                "url": url,
                "local_path": str(local_path.relative_to(ROOT)),
                "http_status": status,
                "content_type": content_type,
                "downloaded_byte_count": len(downloaded),
                "local_byte_count": local_path.stat().st_size,
                "expected_sha256": expected_hash,
                "downloaded_sha256": downloaded_hash,
                "local_sha256": local_hash,
                "downloaded_bytes_match_local_archive": downloaded == local_path.read_bytes(),
                "pass": status == 200
                and downloaded_hash == expected_hash
                and local_hash == expected_hash
                and downloaded == local_path.read_bytes(),
            }
    report = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Geberit",
        "article_number": ARTICLE,
        "replacement_article_number": "154.446.KS.2",
        "replacement_cad_used": False,
        "method": "fresh HTTPS download from each exact archived-article official CDN URL, SHA-256 verification and byte-for-byte local archive comparison",
        "drawing_view_mapping": {"plan": "G", "front": "A", "side": "L"},
        "identity_only_code": "P",
        "results": results,
        "all_declared_urls_accessible": all(item["http_status"] == 200 for item in results.values()),
        "all_downloaded_bytes_match_local_archive": all(
            item["downloaded_bytes_match_local_archive"] for item in results.values()
        ),
        "pass": all(item["pass"] for item in results.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not report["pass"]:
        raise SystemExit("Geberit exact archived-article official CDN revalidation failed")


if __name__ == "__main__":
    main()
