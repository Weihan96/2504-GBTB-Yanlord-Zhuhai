"""Read-only formal IFC inspection; live Blender UI screenshots, no saves."""
import bpy, sys, math, json, os, hashlib
from pathlib import Path
from mathutils import Vector
from bonsai import tool
args=sys.argv[sys.argv.index('--')+1:]
assert '--no-save' in args
source=Path(args[args.index('--ifc')+1])
bpy.ops.bim.load_project(filepath=str(source),should_start_fresh_session=False,use_relative_path=False)
entity=tool.Ifc.get().by_guid('3hA0vKpcn44u4Tsx4tqiUz')
obj=tool.Ifc.get_object(entity)
assert obj is not None
for other in bpy.context.view_layer.objects:
    other.hide_set(other != obj)
    other.select_set(False)
obj.select_set(True)
bpy.context.view_layer.objects.active=obj
centre=obj.matrix_world @ (sum((Vector(p) for p in obj.bound_box),Vector())/8)
out=source.parent/'output/review/highpoly-types/bst02'
state={'pid':os.getpid(),'ifc':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'global_id':entity.GlobalId,'object':obj.name,'matrix_world':[list(row) for row in obj.matrix_world],'vertices':len(obj.data.vertices),'blend_saved':bpy.data.is_saved,'screenshots':[]}
def aim(elevation):
    direction=Vector((1.8,-2.8, elevation)).normalized()
    for screen in bpy.data.screens:
        for area in screen.areas:
            if area.type=='VIEW_3D':
                space=area.spaces.active
                space.region_3d.view_rotation=direction.to_track_quat('Z','Y')
                space.region_3d.view_location=centre
                space.region_3d.view_distance=1.4
                space.region_3d.view_perspective='ORTHO'
                space.shading.type='SOLID'
                space.shading.color_type='SINGLE'
                space.shading.single_color=(0.65,0.68,0.73)
                space.overlay.show_floor=False
                area.tag_redraw()
aim(1.5)
step=0
def capture():
    global step
    target=out/('bst02-actual-ifc-bonsai-3d.png' if step==0 else 'bst02-actual-ifc-bonsai-underside.png')
    bpy.ops.screen.screenshot(filepath=str(target))
    state['screenshots'].append(str(target))
    if step==0:
        aim(-1.8);step=1;return 4.0
    (out/'bst02-actual-ifc-ui-evidence.json').write_text(json.dumps(state,indent=2))
    print('BST02_UI_SCREENSHOTS_READY '+json.dumps(state),flush=True)
    return None
bpy.app.timers.register(capture,first_interval=10.0)
