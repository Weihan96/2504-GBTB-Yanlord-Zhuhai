#!/usr/bin/env python3
"""Validate the paired Bonsai raster/vector source used by review drawings."""

from __future__ import annotations

import hashlib
import json
import re
import struct
from pathlib import Path


WALL_PLAN_RASTER_UNDERLAY = re.compile(
    r'<image\b[^>]*\bxlink:href="(Wall Plan-underlay\.png)"[^>]*/>'
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def png_dimensions(path: Path) -> tuple[int, int]:
    header = path.read_bytes()[:24]
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError(f"invalid Wall Plan PNG: {path}")
    return struct.unpack(">II", header[16:24])


def validate_wall_plan_source(source: str, source_svg: Path, ifc_path: Path) -> None:
    """Require a current, camera-aligned Wall Plan SVG/PNG pair."""
    references = WALL_PLAN_RASTER_UNDERLAY.findall(source)
    if references != ["Wall Plan-underlay.png"]:
        raise RuntimeError(
            "Wall Plan source must contain exactly one paired raster underlay reference"
        )

    source_svg = source_svg.resolve()
    ifc_path = ifc_path.resolve()
    underlay = source_svg.with_name(references[0])
    manifest_path = source_svg.with_name("Wall Plan-source.json")
    if not underlay.is_file() or not manifest_path.is_file():
        raise RuntimeError("Wall Plan SVG, underlay PNG, and source manifest must travel together")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    width, height = png_dimensions(underlay)
    expected = {
        "formal_ifc_sha256": sha256(ifc_path),
        "source_svg_sha256": sha256(source_svg),
        "underlay_png_sha256": sha256(underlay),
        "underlay_width_px": width,
        "underlay_height_px": height,
    }
    mismatches = {
        key: {"manifest": manifest.get(key), "actual": value}
        for key, value in expected.items()
        if manifest.get(key) != value
    }
    if mismatches:
        raise RuntimeError(f"stale Wall Plan source manifest: {mismatches}")

    required_state = {
        "drawing_global_id": "33aMRH9An36RQi37lCJkbO",
        "drawing_name": "Wall Plan",
        "target_view": "PLAN_VIEW",
        "camera_type": "ORTHO",
        "viewport_perspectives": ["CAMERA"],
        "underlay_cache": False,
        "linework_cache": False,
        "annotation_cache": False,
        "unresolved_linked_texture_images": [],
    }
    invalid_state = {
        key: {"manifest": manifest.get(key), "required": value}
        for key, value in required_state.items()
        if manifest.get(key) != value
    }
    if invalid_state:
        raise RuntimeError(f"invalid Wall Plan camera/texture state: {invalid_state}")


def validate_electrical_coordination_source(source: str, source_svg: Path, ifc_path: Path) -> None:
    """Require a current Wall Plan derivative with visible doors and fixed joinery."""
    source_svg = source_svg.resolve()
    ifc_path = ifc_path.resolve()
    wall_plan = source_svg.with_name("Wall Plan.svg")
    if not wall_plan.is_file():
        raise RuntimeError("electrical coordination source must travel with Wall Plan.svg")
    validate_wall_plan_source(wall_plan.read_text(encoding="utf-8"), wall_plan, ifc_path)

    manifest_path = source_svg.with_name("Electrical Coordination Plan-source.json")
    if not manifest_path.is_file():
        raise RuntimeError("electrical coordination source manifest is missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {
        "formal_ifc_sha256": sha256(ifc_path),
        "wall_plan_svg_sha256": sha256(wall_plan),
        "coordination_svg_sha256": sha256(source_svg),
    }
    mismatches = {
        key: {"manifest": manifest.get(key), "actual": value}
        for key, value in expected.items()
        if manifest.get(key) != value
    }
    if mismatches:
        raise RuntimeError(f"stale electrical coordination source manifest: {mismatches}")

    counts = manifest.get("counts", {})
    required_counts = {"walls": 101, "doors": 8, "fixed_furniture": 70}
    if {key: counts.get(key) for key in required_counts} != required_counts:
        raise RuntimeError(f"electrical coordination source counts changed: {counts}")
    if source.count('data-coordination-kind="door"') != required_counts["doors"]:
        raise RuntimeError("electrical coordination source door visibility count changed")
    if source.count('data-coordination-kind="fixed_furniture"') != required_counts["fixed_furniture"]:
        raise RuntimeError("electrical coordination source fixed-furniture visibility count changed")
    if source.count('data-electrical-coordination="current-ifc"') != 1:
        raise RuntimeError("electrical coordination source identity is missing or duplicated")
