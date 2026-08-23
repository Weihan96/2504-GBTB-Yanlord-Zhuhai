#!/usr/bin/env python3
"""Overlay SIS04 semantic views on complete kitchen project drawings."""

import json
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
import render_hima01_project_context as shared


shared.PRODUCT_DIR = ROOT / "output/review/highpoly-types/sis04"
shared.CANDIDATE = shared.PRODUCT_DIR / "candidate-representations.json"
shared.PLAN_SOURCE = ROOT / "drawings/Furniture Plan.svg"
shared.FRONT_SOURCE = ROOT / "drawings/elevations/native/EL-03-07-R04-PX.svg"
shared.SIDE_SOURCE = ROOT / "drawings/elevations/native/EL-03-06-R04-PY.svg"
shared.GLOBAL_ID = "1MzM8Ms2vFo8KEm503j9w2"
shared.BBOX_TOLERANCE_SVG_UNITS = 0.08
shared.ALIGNMENT_MODE = "minimum"
shared.OVERLAY_ID_PREFIX = "sis04"
shared.GENERATOR = "pipeline/scripts/render_sis04_project_context.py"
shared.CONTEXT_VIEWS = (
    ("plan", "plan", shared.PLAN_SOURCE, "project-context-furniture-plan.svg", (10.0, 10.0), 1),
    ("front", "front", shared.FRONT_SOURCE, "project-context-r04-front-elevation.svg", (8.0, 6.0), 0),
    ("side", "side", shared.SIDE_SOURCE, "project-context-r04-side-elevation.svg", (8.0, 6.0), 0),
)


def suppress_plan_elevation_markers(path):
    """Remove only marker anchors that obscure SIS04 in the review copy."""
    content = path.read_text(encoding="utf-8")
    lines = [line for line in content.splitlines() if 'class="official-elevation-anchor"' not in line]
    if len(lines) == len(content.splitlines()):
        raise RuntimeError(f"expected official elevation anchors in {path}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def add_white_review_background(path):
    """Make transparent project paper render white in the PNG review preview."""
    content = path.read_text(encoding="utf-8")
    content, root_count = re.subn(r"<svg\s", '<svg style="background:#ffffff" ', content, count=1)
    if root_count != 1:
        raise RuntimeError(f"missing SVG root in {path}")
    match = re.search(r'viewBox="([^\"]+)"', content)
    if match is None:
        raise RuntimeError(f"missing viewBox in {path}")
    x, y, width, height = (float(value) for value in match.group(1).split())
    background = (
        f'<rect class="review-paper-background" x="{x:.6f}" y="{y:.6f}" '
        f'width="{width:.6f}" height="{height:.6f}" '
        f'style="fill:#ffffff !important;fill-opacity:1 !important;stroke:none"/>\n'
    )
    if "</defs>" not in content:
        raise RuntimeError(f"missing defs in {path}")
    path.write_text(content.replace("</defs>", "</defs>\n" + background, 1), encoding="utf-8")


def write_uncached_png_preview(source, target):
    """Render the SVG on opaque HTML paper without Quick Look's canvas tint."""
    chrome = Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
    if not chrome.is_file():
        raise RuntimeError("Google Chrome is required to render the opaque project-context preview")
    with tempfile.TemporaryDirectory(prefix="sis04-context-preview-") as temporary:
        temporary_path = Path(temporary)
        uncached_source = temporary_path / "sis04-plan-review-uncached.svg"
        wrapper = temporary_path / "sis04-plan-review.html"
        profile = temporary_path / "chrome-profile"
        shutil.copy2(source, uncached_source)
        wrapper.write_text(
            "<!doctype html><meta charset=\"utf-8\">"
            "<style>html,body{margin:0;width:100%;height:100%;overflow:hidden;background:#fff}"
            "img{display:block;width:100vw;height:100vh;object-fit:contain;background:#fff}</style>"
            f'<img src="{uncached_source.name}">',
            encoding="utf-8",
        )
        target.unlink(missing_ok=True)
        process = subprocess.Popen(
            [
                str(chrome),
                "--headless=new",
                "--disable-gpu",
                "--disable-background-networking",
                "--disable-component-update",
                "--no-first-run",
                "--no-default-browser-check",
                f"--user-data-dir={profile}",
                "--window-size=1800,1800",
                "--force-device-scale-factor=1",
                "--hide-scrollbars",
                "--default-background-color=ffffffff",
                f"--screenshot={target}",
                wrapper.as_uri(),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        for _ in range(100):
            if target.is_file() and target.stat().st_size:
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
        if not target.is_file() or not target.stat().st_size:
            raise RuntimeError("project-context preview failed")


if __name__ == "__main__":
    shared.main()
    manifest_path = shared.PRODUCT_DIR / "project-context-manifest.json"
    manifest = load_json(manifest_path)
    plan_record = next(item for item in manifest["views"] if item["view"] == "plan")
    plan_output = ROOT / plan_record["output"]
    plan_crop = ROOT / plan_record["review_crop"]
    plan_preview = ROOT / plan_record["review_preview"]
    suppress_plan_elevation_markers(plan_output)
    suppress_plan_elevation_markers(plan_crop)
    add_white_review_background(plan_crop)
    write_uncached_png_preview(plan_crop, plan_preview)
    plan_record["output_sha256"] = sha256(plan_output)
    plan_record["review_crop_sha256"] = sha256(plan_crop)
    plan_record["review_preview_sha256"] = sha256(plan_preview)
    manifest["semantic_view_mapping"] = {
        "plan": {"candidate_axes": [0, 1], "source": relative(shared.PLAN_SOURCE), "rotate_quarter_turns": 1},
        "front": {"candidate_axes": [0, 2], "source": relative(shared.FRONT_SOURCE), "project_direction": "+X"},
        "side": {"candidate_axes": [1, 2], "source": relative(shared.SIDE_SOURCE), "project_direction": "+Y"},
    }
    manifest["review_annotation_suppression"] = {
        "plan": "official-elevation-anchor groups only",
        "reason": "the EL-03 marker cluster overlaps the wall unit; walls, furniture, grids and IFC geometry remain retained",
        "geometry_removed": False,
    }
    manifest["review_preview_background"] = {
        "plan": "opaque white review-only paper background",
        "full_project_svg_remains_transparent": True,
        "drawing_geometry_changed": False,
    }
    write_json(manifest_path, manifest)
    print(json.dumps({"updated_manifest": relative(manifest_path), "pass": manifest["pass"]}, indent=2))
