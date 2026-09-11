#!/usr/bin/env python3
"""Open the approved TRAP01 derived IFC and expose the public Bonsai MCP bridge."""

from pathlib import Path
import addon_utils
import bpy


ROOT = Path(__file__).resolve().parents[2]
DERIVED_IFC = ROOT / "output/review/highpoly-types/trap01/Geberit-151.116.11.1-TRAP01-derived-drawing.ifc"

with __import__("contextlib").suppress(Exception):
    addon_utils.disable("bl_ext.user_default.project_control", default_set=False, handle_error=None)
addon_utils.enable("bonsai_bridge", default_set=False, persistent=False)
from bonsai import tool

result = bpy.ops.bim.load_project(
    filepath=str(DERIVED_IFC), should_start_fresh_session=True, use_detailed_tooltip=True
)
if result != {"FINISHED"} or not tool.Ifc.get():
    raise RuntimeError(f"failed to load TRAP01 derived IFC: {result}")
preferences = bpy.context.preferences.addons["bonsai_bridge"].preferences
preferences.allow_edits = True
bpy.ops.bonsai_mcp.start_bridge()
print(f"TRAP01_INTERNAL_DETAIL_BONSAI_READY {DERIVED_IFC}")
