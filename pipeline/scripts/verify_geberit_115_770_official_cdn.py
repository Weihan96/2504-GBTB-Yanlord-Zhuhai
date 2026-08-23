#!/usr/bin/env python3
"""Re-download exact Geberit 115.770 DWGs and compare them with the archive."""

from __future__ import annotations

import argparse
import hashlib
import json
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARTICLE = "115.770.11.5"
OUTPUT = ROOT / "output/review/highpoly-types/geberit-115-770/official-source/official-cdn-revalidation.json"
SOURCE_DIR = ROOT / "output/review/highpoly-types/geberit-115-770/official-source"
EXPECTED = {
    "A": "88c292a019ee91e0872baebf4a74237803193662d4f2150141cb1772cef1bd8a",
    "G": "1c0e5dd30db431dbf0f890311af8fa745b84aff6f2284c5dd04addc91b07e70a",
    "L": "0f55d3204f94637ce6fe544847e045ba7a0768f1abe8cec94b6763bf2a2220ad",
    "P": "a90f58e350d4e2122e1b3c61e16e36a13d83edbba15ffbb4f1d7c2f18b86e8db",
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
    with tempfile.TemporaryDirectory(prefix="geberit-115770-cdn-") as temporary:
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
            local_bytes = local_path.read_bytes()
            results[code] = {
                "url": url,
                "local_path": str(local_path.relative_to(ROOT)),
                "http_status": status,
                "content_type": content_type,
                "downloaded_byte_count": len(downloaded),
                "local_byte_count": len(local_bytes),
                "expected_sha256": expected_hash,
                "downloaded_sha256": downloaded_hash,
                "local_sha256": local_hash,
                "downloaded_bytes_match_local_archive": downloaded == local_bytes,
                "pass": status == 200
                and downloaded_hash == expected_hash
                and local_hash == expected_hash
                and downloaded == local_bytes,
            }
    report = {
        "schema_version": 1,
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "manufacturer": "Geberit",
        "family": "Sigma01 dual-flush actuator plate",
        "article_number": ARTICLE,
        "method": "fresh HTTPS download from each declared official CDN URL, SHA-256 verification and byte-for-byte local archive comparison",
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
        raise SystemExit("Geberit 115.770 official CDN revalidation failed")


if __name__ == "__main__":
    main()
