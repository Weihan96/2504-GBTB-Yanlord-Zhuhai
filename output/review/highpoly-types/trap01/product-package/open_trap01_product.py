"""Task-owned read-only IFC launch and session-local public bridge on 9884."""
import sys
import bpy
import addon_utils
from bonsai import tool

args = sys.argv[sys.argv.index("--") + 1:]
assert "--no-save" in args
path = args[args.index("--ifc") + 1]
assert not bpy.data.is_saved
bpy.ops.bim.load_project(filepath=path, should_start_fresh_session=False, use_relative_path=False)
addon_utils.enable("bonsai_bridge", default_set=False, persistent=False)
import bonsai_bridge
assert bonsai_bridge.bl_info["version"] == (1, 1, 0)
prefs = bpy.context.preferences.addons["bonsai_bridge"].preferences
prefs.port = 9884
prefs.allow_edits = True
bpy.ops.bonsai_mcp.start_bridge()
assert tool.Ifc.get_path() == path
