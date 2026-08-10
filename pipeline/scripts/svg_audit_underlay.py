#!/usr/bin/env python3
"""Shared safeguards for raster-free review drawings derived from Wall Plan.svg."""

from __future__ import annotations

import re


WALL_PLAN_RASTER_UNDERLAY = re.compile(
    r'\s*<image\b[^>]*\bxlink:href="Wall Plan-underlay\.png"[^>]*/>\s*'
)


def strip_wall_plan_raster_underlay(source: str, drawing_name: str) -> str:
    """Remove the single stale raster snapshot before adding review markup."""
    source, removed_underlays = WALL_PLAN_RASTER_UNDERLAY.subn("\n", source, count=1)
    if removed_underlays != 1:
        raise RuntimeError(
            f"{drawing_name} source must contain exactly one removable Wall Plan raster underlay"
        )
    return source
