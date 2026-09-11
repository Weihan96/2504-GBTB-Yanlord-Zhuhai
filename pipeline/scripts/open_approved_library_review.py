"""Launcher IFC bootstrap: read-only pure IFC plus session-local library UI."""
import sys
from pathlib import Path
import importlib.util
import bpy
import addon_utils
import socket
import json
import os
import runpy

args = sys.argv[sys.argv.index("--") + 1:]
assert "--no-save" in args
path = args[args.index("--ifc") + 1]
assert not bpy.data.is_saved
assert bpy.ops.bim.load_project(filepath=path, should_start_fresh_session=False, use_relative_path=False) == {"FINISHED"}
root = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("highpoly_review_library", root / "pipeline/addons/highpoly_review_library/__init__.py")
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
module.register()
addon_utils.enable("bonsai_bridge", default_set=False, persistent=False)
import bonsai_bridge
assert bonsai_bridge.bl_info["version"] == (1, 1, 0)
prefs = bpy.context.preferences.addons["bonsai_bridge"].preferences
prefs.allow_edits = True
# This task's acceptance connection must not reuse another Blender's bridge.
# Another task may have started its own bridge since this window last closed.
# Select a free port instead of failing this IFC session or touching that task.
with socket.socket() as probe:
    try:
        probe.bind(('127.0.0.1', 9881))
    except OSError:
        probe.bind(('127.0.0.1', 0))
    prefs.port = probe.getsockname()[1]
bpy.ops.bonsai_mcp.start_bridge()
print('REVIEW_LIBRARY_BRIDGE_READY ' + json.dumps({'pid': os.getpid(), 'port': prefs.port}), flush=True)
runpy.run_path(str(root / 'pipeline/scripts/activate_approved_library.py'))
