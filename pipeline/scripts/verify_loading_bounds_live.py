"""Update only the non-RNA loading method; render a display-only test overlay.

The screenshot is a synthetic 65% style check, not an actual import progress
claim. Existing IFC memory, objects, selection, camera and disk stay unchanged.
"""
import argparse
from pathlib import Path
from bonsai_mcp.blender_client import BlenderBridgeClient

parser = argparse.ArgumentParser()
parser.add_argument('--port', type=int, required=True)
parser.add_argument('--pid', type=int, required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
code = '''import ast, bpy, os, json, hashlib
from pathlib import Path
from mathutils import Vector
from bpy_extras import view3d_utils
import bonsai.tool as tool
import highpoly_review_library as lib
import highpoly_review_library.loading_feedback as feedback
from highpoly_review_library.native_assets import PlacementPreview
assert os.getpid()==EXPECTED_PID
assert feedback._active is None and not bpy.context.window.modal_operators
root=Path(ROOT)
out=root/'output/review/approved-product-library'
assert Path(tool.Ifc.get_path())==out/'ifc-card-ui/library-review.ifc'
def state():
    return {'ifc_memory':hashlib.sha256(tool.Ifc.get().to_string().encode()).hexdigest(),
      'ifc_disk':hashlib.sha256(Path(tool.Ifc.get_path()).read_bytes()).hexdigest(),
      'formal':hashlib.sha256((root/'2504 GBTB Yanlord Zhuhai.ifc').read_bytes()).hexdigest(),
      'objects':[(o.name,[v for row in o.matrix_world for v in row]) for o in bpy.data.objects],
      'selection':[o.name for o in bpy.context.selected_objects],
      'active':bpy.context.active_object.name if bpy.context.active_object else None,
      'product':bpy.context.scene.review_library.product,
      'views':[(list(a.spaces.active.region_3d.view_rotation),list(a.spaces.active.region_3d.view_location),a.spaces.active.region_3d.view_distance)
               for a in bpy.context.screen.areas if a.type=='VIEW_3D']}
report={'status':'running','pid':os.getpid(),'mode':'synthetic_style_only_no_import','before':state()}
path=out/'loading-bounds-validation.json'
path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\\n')
source=root/'pipeline/addons/highpoly_review_library/loading_feedback.py'
cls=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.ClassDef) and n.name=='LoadingFeedback')
method=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__enter__')
# No RNA unregister/re-register and no catalog/preview-cache reset.
exec(compile(ast.Module(body=[method],type_ignores=[]),str(source),'exec'),feedback.__dict__)
feedback.LoadingFeedback.__enter__=feedback.__dict__.pop('__enter__')
lib.bl_info['version']=(0,8,3)
area=next(a for a in bpy.context.screen.areas if a.type=='VIEW_3D')
region=next(r for r in area.regions if r.type=='WINDOW')
sidebar=next(r for r in area.regions if r.type=='UI')
right=sidebar.x if sidebar.width>1 else region.x+region.width
xy=((right-region.x)*.55,region.height*.50)
rv3d=area.spaces.active.region_3d
# Center the synthetic display sample in the unobstructed viewport without
# moving the user's view or relying on a floor plane behind a close-up camera.
point=view3d_utils.region_2d_to_location_3d(region,rv3d,xy,rv3d.view_location)
preview=PlacementPreview(((-.3,-.3,-.6),(.3,.3,.6)))
try:
    preview.update((area,region),point)
    with bpy.context.temp_override(area=area,region=region), feedback.LoadingFeedback(preview,'display-only-check',bpy.context) as loading:
        loading.advance(.65,'显示样式验证：蓝色边界 + 无描边绿色填充（非真实导入）')
        assert preview.outline_visible
        report['blue_outline_visible']=preview.outline_visible
        report['green_primitive']='TRIS'
        bpy.ops.screen.screenshot(filepath=str(out/'loading-bounds-blue-green.png'))
finally:
    preview.close()
    area.tag_redraw()
report['after']=state()
report['handlers_removed']=preview.handle is None and feedback._active is None and loading.handle is None
assert report['before']==report['after'] and report['handlers_removed']
report['status']='pass'
path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\\n')
print(json.dumps({'status':report['status'],'pid':os.getpid(),'version':lib.bl_info['version'],
                 'state_unchanged':report['before']==report['after'],'handlers_removed':report['handlers_removed']}))
'''.replace('EXPECTED_PID',str(args.pid)).replace('ROOT',repr(str(root)))
result=BlenderBridgeClient(port=args.port,timeout=30).send('execute_code',{'code':code})
if not result['success']:
    raise RuntimeError(result)
print(result['stdout'])
