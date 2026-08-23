#!/usr/bin/env python3
"""Generate missing review contact sheets from verified package evidence."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INVENTORY = ROOT / "pipeline/decisions/highpoly-product-inventory.json"
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
CHROME = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
WIDTH = 2400
HEIGHT = 1800


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def card(title: str, image: str, css_class: str = "") -> str:
    image_element = f'<img src="{html.escape(image)}">' if image else ""
    return (
        f'<article class="{html.escape(css_class)}"><h2>{html.escape(title)}</h2>'
        f'{image_element}</article>'
    )


def context_preview(record: dict) -> str:
    if record.get("review_preview"):
        return Path(record["review_preview"]).name
    if record.get("review_crops"):
        return Path(record["review_crops"][0]["preview"]).name
    if record.get("review_crop"):
        return Path(record["review_crop"]).name
    raise RuntimeError(f'project-context record has no review crop: {record.get("view")}')


def write_html(product: dict, folder: Path, manifest: dict, bonsai: dict, context: dict) -> Path:
    drawing_source = manifest.get("drawing_source", {})
    source_kind = (
        manifest.get("source_kind")
        or drawing_source.get("source_kind")
        or manifest.get("official_reference", {}).get("source_kind")
        or "unresolved"
    )
    official_cad = (
        manifest.get("official_cad_used") is True
        or drawing_source.get("official_cad_used") is True
        or source_kind in {"native_dwg", "native_dxf"}
        or any(record.get("blue_line_present") is True for record in manifest.get("views", []))
    )
    source_label = (
        manifest.get("source_label_zh")
        or drawing_source.get("source_label_zh")
        or ("官方原生 DWG" if source_kind == "native_dwg" else "")
        or ("官方原生 DXF" if source_kind == "native_dxf" else "")
        or source_kind
    )
    display_name = manifest.get("display_name") or product["type_name"]
    context_cards = [
        card(f'Project · {record["view"]}', context_preview(record))
        for record in context.get("views", [])[:3]
    ]
    context_cards.extend(
        card("Project context", "", "empty")
        for _ in range(3 - len(context_cards))
    )
    render_by_view = {record["view"]: Path(record["path"]).name for record in bonsai["renders"]}
    summary = f'''<article class="summary"><h2>Source / gate</h2><dl>
<dt>IFC type</dt><dd>{html.escape(str(manifest.get("ifc_type_name")))}</dd>
<dt>Line source</dt><dd>{html.escape(source_label)}</dd>
<dt>Official CAD used</dt><dd>{str(official_cad).lower()}</dd>
<dt>Formal IFC unchanged</dt><dd>{str(manifest.get("formal_ifc_bytes_unchanged") is True).lower()}</dd>
<dt>Review</dt><dd>{html.escape(str(manifest.get("review_status")))}</dd>
</dl></article>'''
    body = "".join([
        card("Three-view · Plan", "plan.svg"),
        card("Three-view · Front", "front.svg"),
        card("Three-view · Side", "side.svg"),
        summary,
        card("Bonsai · Plan camera", render_by_view["plan"], "bonsai"),
        card("Bonsai · Front camera", render_by_view["front"], "bonsai"),
        card("Bonsai · Side camera", render_by_view["side"], "bonsai"),
        card("Bonsai · Isometric camera", render_by_view["iso"], "bonsai"),
        *context_cards,
        card("Approval gate", "", "gate"),
    ])
    target = folder / "review-contact-sheet.html"
    target.write_text(
        f'''<!doctype html><html lang="en" data-generator="generate_highpoly_review_contact_sheets.py"><meta charset="utf-8"><title>{html.escape(display_name)}</title>
<style>*{{box-sizing:border-box}}html,body{{margin:0;width:{WIDTH}px;height:{HEIGHT}px;overflow:hidden;background:#e9e6df;color:#182430;font-family:Arial,sans-serif}}header{{height:90px;padding:18px 28px;background:#f8f7f3;border-bottom:2px solid #b8bec4}}h1{{margin:0 0 6px;font-size:30px}}header p{{margin:0;font-size:16px;color:#485c6f}}main{{display:grid;grid-template-columns:repeat(4,1fr);grid-template-rows:repeat(3,560px);gap:10px;padding:10px}}article{{position:relative;overflow:hidden;background:#fff;border:1px solid #c8ced4}}h2{{position:absolute;z-index:2;left:12px;top:10px;margin:0;padding:6px 9px;border-radius:4px;background:rgba(255,255,255,.92);font-size:18px}}img{{width:100%;height:100%;display:block;object-fit:contain;background:#fff}}.bonsai img{{background:#303030}}.summary{{padding:70px 24px 20px}}.summary h2{{font-size:22px}}dl{{display:grid;grid-template-columns:160px 1fr;gap:12px;font-size:18px}}dt{{font-weight:bold}}dd{{margin:0;overflow-wrap:anywhere}}.empty,.gate{{background:#f5f3ee}}.empty h2{{color:#7a8791}}.gate:after{{content:'Pending human approval · no derived IFC · no commit';position:absolute;left:30px;right:30px;top:48%;font-size:24px;text-align:center;color:#874b35}}</style>
<header><h1>{html.escape(display_name)}</h1><p>Grey = actual IFC Body · Black = simplified proxy · Blue = verified official CAD only · actual Bonsai camera renders</p></header><main>{body}</main></html>''',
        encoding="utf-8",
    )
    return target


def render(html_path: Path, png_path: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="highpoly-contact-sheet-") as profile:
        png_path.unlink(missing_ok=True)
        process = subprocess.Popen(
            [
                str(CHROME),
                "--headless=new",
                "--disable-gpu",
                "--disable-background-networking",
                "--disable-component-update",
                "--no-first-run",
                "--no-default-browser-check",
                "--hide-scrollbars",
                "--force-device-scale-factor=1",
                f"--window-size={WIDTH},{HEIGHT}",
                f"--user-data-dir={profile}",
                f"--screenshot={png_path}",
                html_path.as_uri(),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        for _ in range(300):
            if png_path.is_file() and png_path.stat().st_size:
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
    if not png_path.is_file() or png_path.stat().st_size == 0:
        raise RuntimeError(f"contact-sheet render failed for {html_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="overwrite existing pending-product sheets")
    parser.add_argument(
        "--refresh-generated",
        action="store_true",
        help="overwrite only sheets previously created by this generator's four-column layout",
    )
    parser.add_argument(
        "--slug",
        help="limit generation to one inventory folder slug",
    )
    args = parser.parse_args()
    if not CHROME.is_file():
        raise RuntimeError("Google Chrome is required")
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    inventory = load_json(INVENTORY)
    generated = []
    skipped = []
    for product in inventory["products"]:
        if product.get("status") != "review_ready_pending_approval":
            continue
        folder = ROOT / product["folder"]
        if args.slug and folder.name != args.slug:
            continue
        target = folder / "review-contact-sheet.png"
        html_path = folder / "review-contact-sheet.html"
        refresh_generated = (
            args.refresh_generated
            and html_path.is_file()
            and "grid-template-columns:repeat(4" in html_path.read_text(encoding="utf-8")
        )
        if target.is_file() and not args.all and not refresh_generated:
            skipped.append(product["type_name"])
            continue
        manifest = load_json(folder / "manifest.json")
        bonsai = load_json(folder / "bonsai-review-manifest.json")
        context = load_json(folder / "project-context-manifest.json")
        if (
            manifest.get("formal_ifc_bytes_unchanged") is not True
            or bonsai.get("mode") != "actual_bonsai_ifc_body_camera_render"
            or bonsai.get("pass") is not True
            or context.get("pass") is not True
            or {record.get("view") for record in bonsai.get("renders", [])} != {"plan", "front", "side", "iso"}
        ):
            raise RuntimeError(f"package evidence gate failed: {product['type_name']}")
        html_path = write_html(product, folder, manifest, bonsai, context)
        render(html_path, target)
        generated.append({
            "type_name": product["type_name"],
            "folder": product["folder"],
            "html": html_path.relative_to(ROOT).as_posix(),
            "png": target.relative_to(ROOT).as_posix(),
            "png_sha256": sha256(target),
        })
        print(f"generated {product['type_name']}: {generated[-1]['png_sha256']}", flush=True)
    if sha256(FORMAL_IFC) != FORMAL_SHA256:
        raise RuntimeError("formal IFC changed during contact-sheet generation")
    report = {
        "schema_version": 1,
        "formal_ifc_sha256": FORMAL_SHA256,
        "generated_count": len(generated),
        "skipped_existing_count": len(skipped),
        "generated": generated,
        "pass": True,
    }
    report_path = ROOT / "output/review/highpoly-types/contact-sheet-generation-report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(report_path.relative_to(ROOT))


if __name__ == "__main__":
    main()
