"""Exercise the native asset grid with Blender's window event queue.

Requires a task-owned, --enable-event-simulate acceptance fixture, never a
formal project. Coordinates must come from a screenshot of that same window.
Does not save IFC or alter the global keymap. Leaves test inserts in memory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
from bonsai_mcp.blender_client import BlenderBridgeClient

parser = argparse.ArgumentParser()
parser.add_argument('--port', type=int, required=True)
parser.add_argument('--pid', type=int, required=True)
parser.add_argument('--card', type=int, nargs=2, required=True)
parser.add_argument('--drop', type=int, nargs=2, required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
out = root / 'output/review/approved-product-library'
fixture = out / 'ifc-card-ui/library-review.ifc'
formal = root / '2504 GBTB Yanlord Zhuhai.ifc'
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
report = {'status': 'running', 'pid': args.pid, 'port': args.port,
          'mechanism': 'native mouse press/move/release events; not direct insert_ifc calls',
          'fixture_sha256': digest(fixture), 'formal_sha256': digest(formal), 'events': []}
client = BlenderBridgeClient(port=args.port, timeout=60)

def execute(code):
    result = client.send('execute_code', {'code': code})
    if not result['success']:
        raise RuntimeError(result)
    return result['stdout'].strip()

def state():
    return json.loads(execute('''import bpy,json,hashlib
import bonsai.tool as tool
f=tool.Ifc.get()
print(json.dumps({'elements':[(e.GlobalId,e.Name,len(e.Representation.Representations)) for e in f.by_type('IfcElement')],
'memory_sha256':hashlib.sha256(f.to_string().encode()).hexdigest(),
'modal':[o.bl_idname for o in bpy.context.window.modal_operators]}))'''))

def event(kind, value, point):
    x, y = point
    execute(f"import bpy; bpy.context.window.event_simulate(type={kind!r},value={value!r},x={x},y={y})")
    report['events'].append([kind, value, x, y])
    time.sleep(.08)

def screenshot(name):
    path = out / name
    execute(f"import bpy; bpy.ops.wm.redraw_timer(type='DRAW_WIN_SWAP',iterations=2); bpy.ops.screen.screenshot(filepath={str(path)!r})")
    return path.name

def begin_drag():
    event('MOUSEMOVE', 'NOTHING', args.card)
    event('LEFTMOUSE', 'PRESS', args.card)
    # Cross the native threshold gradually; never invoke the drag op directly.
    for step in range(1, 16):
        xy = tuple(round(a+(b-a)*step/15) for a,b in zip(args.card,args.drop))
        event('MOUSEMOVE', 'NOTHING', xy)
    result = state()
    assert 'REVIEW_LIBRARY_OT_drag_ifc' in result['modal'], result
    return result

try:
    execute(f'''import os,bpy
import bonsai.tool as tool
assert os.getpid()=={args.pid}
assert tool.Ifc.get_path()=={str(fixture)!r}
assert bpy.app.use_event_simulate and not bpy.data.is_saved
''')
    report['before'] = state()
    event('MOUSEMOVE', 'NOTHING', args.card)
    event('LEFTMOUSE', 'PRESS', args.card)
    report['press_image'] = screenshot('drag-click-press.png')
    event('LEFTMOUSE', 'RELEASE', args.card)
    report['click_image'] = screenshot('drag-click-details.png')
    event('ESC', 'PRESS', args.card)
    event('ESC', 'RELEASE', args.card)
    assert state()['memory_sha256'] == report['before']['memory_sha256']
    begin_drag()
    event('ESC', 'PRESS', args.drop)
    event('ESC', 'RELEASE', args.drop)
    event('LEFTMOUSE', 'RELEASE', args.drop)
    report['escape'] = state()
    assert report['escape']['memory_sha256'] == report['before']['memory_sha256']
    assert 'REVIEW_LIBRARY_OT_drag_ifc' not in report['escape']['modal']
    begin_drag()
    event('MOUSEMOVE', 'NOTHING', args.card)
    event('LEFTMOUSE', 'RELEASE', args.card)
    report['sidebar_release'] = state()
    assert report['sidebar_release']['memory_sha256'] == report['before']['memory_sha256']
    begin_drag()
    report['bounds_image'] = screenshot('drag-verified-bounds.png')
    event('LEFTMOUSE', 'RELEASE', args.drop)
    report['after'] = state()
    assert len(report['after']['elements']) == len(report['before']['elements'])+1
    assert not report['after']['modal']
    report['inserted_image'] = screenshot('drag-verified-inserted.png')
    assert digest(formal) == report['formal_sha256']
    assert digest(fixture) == report['fixture_sha256']
    report['status'] = 'pass'
finally:
    (out/'drag-gesture-validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report,ensure_ascii=False))
