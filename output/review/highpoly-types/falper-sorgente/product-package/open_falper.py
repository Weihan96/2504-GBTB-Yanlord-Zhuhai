"""Task-owned IFC loader; public bridge only, no preference persistence."""
import sys
import bpy
import addon_utils
from bonsai import tool

args = sys.argv[sys.argv.index('--') + 1:]
assert '--no-save' in args and not bpy.data.is_saved
path = args[args.index('--ifc') + 1]
bpy.ops.bim.load_project(filepath=path, should_start_fresh_session=False, use_relative_path=False)
addon_utils.enable('bonsai_bridge', default_set=False, persistent=False)
import bonsai_bridge
assert bonsai_bridge.bl_info['version'] == (1, 1, 0)
prefs = bpy.context.preferences.addons['bonsai_bridge'].preferences
prefs.port = 9893
prefs.allow_edits = True
bpy.ops.bonsai_mcp.start_bridge()
assert tool.Ifc.get_path() == path
target = tool.Ifc.get_object(tool.Ifc.get().by_guid('350tdaubr8QP3Cu2YMQZIN'))
for obj in bpy.context.scene.objects:
    if obj.type == 'MESH':
        obj.hide_set(obj != target)
    obj.select_set(False)
target.hide_set(False)
target.select_set(True)
bpy.context.view_layer.objects.active = target
