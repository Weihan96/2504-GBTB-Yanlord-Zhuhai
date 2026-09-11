"""Duofix approved per-view sessions -> pure product and external scene recipe."""
from pathlib import Path
import json, sys, tempfile, shutil, traceback
import numpy as np
import ifcopenshell
OUT = Path(__file__).resolve().parent
PRODUCT = OUT.parent
ROOT = OUT.parents[4]
sys.path.insert(0, str(ROOT / 'pipeline/scripts'))
import review_product_package as pkg
from pure_product_package import build_pure_package, attach_recipe, spatial_scope
GUID = '3pvAlH5C14v8uVEJ1LmK8M'
FORMAL = ROOT / '2504 GBTB Yanlord Zhuhai.ifc'
FORMAL_HASH = '7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c'
SINGLE = OUT / 'DUOFIX-product.ifc'
RECIPE = OUT / 'scene-recipe.json'
REPORT = OUT / 'validation.json'
LEGACY_DIR = PRODUCT / 'bonsai-drawings/main-bathroom'
APPROVAL = ROOT / 'pipeline/decisions/geberit-duofix-sigma-224-212-drawing-approval.json'
OLD_MANIFEST = LEGACY_DIR / 'drawing-output-manifest.json'
SPECS = {
 'plan': {'drawing_guid':'0WcY_bPRH6KP5VsZwERB2F', 'annotation_guid':'0IrLefJmT6$u_j1FOarh04'},
 'front': {'drawing_guid':'0sv4QjrFH86BYBluUvo2Ma', 'annotation_guid':'0dqurPwyXAphRH_QF$$Hb1'},
 'side': {'drawing_guid':'2Dm4xfxuH18wQZamIHq1WI', 'annotation_guid':'1RLWGXL2r3Rv29jQrPNgQc'}}

def write(path, data):
 Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n')

def snapshot(model):
 t=model.by_guid(GUID)
 return {'body':pkg.body_fingerprint(t),'representations':pkg.fingerprint(t.Representation),
  'placement':pkg.placement(t).tolist(),'physical_elements':sorted(e.GlobalId for e in model.by_type('IfcElement')),
  'annotation_count':len(model.by_type('IfcAnnotation')),'unit_scale_to_m':pkg.unit_util.calculate_unit_scale(model)}

def check_pure(model):
 s=snapshot(model)
 assert s['physical_elements']==[GUID] and s['annotation_count']==0
 assert not model.by_type('IfcGroup')
 assert not [s for s in model.by_type('IfcSpatialElement') if s.Representation]
 assert not [p for p in model.by_type('IfcPropertySet') if p.Name=='EPset_Drawing']
 assert len(model.by_type('IfcRoot'))==len({e.GlobalId for e in model.by_type('IfcRoot')})
 return s

def prepare():
 assert pkg.sha256(FORMAL)==FORMAL_HASH
 manifest=json.loads(OLD_MANIFEST.read_text())
 assert json.loads(APPROVAL.read_text())['status']=='approved'
 protected=[FORMAL,APPROVAL,OLD_MANIFEST,PRODUCT/'official-source/source-access-record.json',
  ROOT/'pipeline/decisions/highpoly-subtask-completion-2026-09-05.json',
  PRODUCT/'Geberit-Duofix-Sigma-224-212-derived-drawing.ifc']
 protected.extend(p for p in (PRODUCT/'official-source').rglob('*') if p.is_file())
 combined=None; views=[]
 for view,spec in SPECS.items():
  key='left_side' if view=='side' else view
  session=LEGACY_DIR/f"duofix-main-bath-{key.replace('_','-')}-session.ifc"
  vm=manifest['views'][key]
  assert pkg.sha256(session)==vm['drawing_session_ifc_sha256']
  oldsvg=LEGACY_DIR/vm['svg']; assert pkg.sha256(oldsvg)==vm['svg_sha256']
  protected.extend([session,oldsvg])
  source=ifcopenshell.open(str(session))
  if combined is None: combined=source
  else:
   assert pkg.body_fingerprint(source.by_guid(GUID))==pkg.body_fingerprint(combined.by_guid(GUID))
   annotations=[source.by_guid(v) for v in spec.values()]
   scope=spatial_scope(source,source.by_guid(GUID))|set(annotations)
   for a in annotations: scope.update(r.RelatingGroup for r in a.HasAssignments if r.is_a('IfcRelAssignsToGroup'))
   clone=pkg.ScopedCopy(source,combined,scope,reuse_roots=True,skip_inverse_ids=[source.by_guid(GUID).id()])
   clone.memo[source.by_guid(GUID).id()]=combined.by_guid(GUID)
   for a in annotations: clone.copy(a)
   clone.inverses()
  views.append({'view':view,**spec,'old_svg':pkg.record(oldsvg),'source_session':pkg.record(session)})
 # Combining independent sessions duplicates identical IfcApplication roots
 # (non-IfcRoot entities). Share equal metadata to satisfy UR1 / UR2.
 import ifcopenshell.util.element as eu
 applications={};deduplicated_applications=[]
 for app in list(combined.by_type('IfcApplication')):
  key=pkg.fingerprint(app)
  if key in applications:
   keeper=applications[key]
   for inverse in list(combined.get_inverse(app)):eu.replace_attribute(inverse,app,keeper)
   deduplicated_applications.append({'identifier':app.ApplicationIdentifier,'version':app.Version,'fingerprint':key})
   combined.remove(app)
  else:applications[key]=app
 # A single cistern is aggregated into WC_01 / Drain Center in the project.
 # Preserve its full placement chain, but do not package physical siblings or
 # empty assembly products. Runtime matching retains the project's assemblies.
 original_assembly_relations=[]
 target=combined.by_guid(GUID)
 original_placement=pkg.placement(target).copy()
 for relation in list(target.Decomposes):
  if relation.RelatingObject.is_a('IfcElement'):
   original_assembly_relations.append({'relationship_guid':relation.GlobalId,
    'parent_guid':relation.RelatingObject.GlobalId,'parent_name':relation.RelatingObject.Name})
   remaining=[e for e in relation.RelatedObjects if e!=target]
   if remaining:relation.RelatedObjects=remaining
   else:combined.remove(relation)
 assert np.array_equal(pkg.placement(target),original_placement)
 pure,recipe,audit=build_pure_package(combined,GUID,SPECS)
 audit['external_physical_assembly_membership']=original_assembly_relations
 audit['deduplicated_identical_application_metadata']=deduplicated_applications
 audit['assembly_policy']='Pure IFC omits physical ancestor assemblies; project attachment matches the existing occurrence and retains its original assembly relationships. Full local placement dependency chain unchanged.'
 pure.write(str(SINGLE));write(RECIPE,recipe)
 state=check_pure(ifcopenshell.open(str(SINGLE)))
 dependencies=[]
 for location in sorted({s.Location for s in pure.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}):
  rel=Path(location);assert not rel.is_absolute() and '..' not in rel.parts
  source=ROOT/rel;dest=OUT/rel;dest.parent.mkdir(parents=True,exist_ok=True)
  if not dest.exists():shutil.copy2(source,dest)
  assert pkg.sha256(source)==pkg.sha256(dest)
  dependencies.append({'source':pkg.record(source),'packaged':pkg.record(dest),'location':location})
 write(REPORT,{'product':'DUOFIX','target_guid':GUID,'verdict':'pending_bonsai_scene',
  'formal_sha256_before':FORMAL_HASH,'protected_files':[pkg.record(p) for p in protected],
  'views':views,'extraction':audit,'pure_pre_bonsai':state,'single_product':pkg.record(SINGLE),
  'scene_recipe':pkg.record(RECIPE),'necessary_external_style_dependencies':dependencies,
  'courseEvidence':{'mode':'embedded-course-index','lesson':'085000','timestamps':['01:59 Create Drawing','02:13 SVG']},
  'formal_write_allowed':False,'cleanup_performed':False})
 print(json.dumps({'ifc':pkg.record(SINGLE),'extraction':audit}))

def scene():
 import bpy,bonsai_bridge as bridge
 from bonsai import tool
 import create_wd03_wardrobe_scene_drawings as context
 import create_geberit_duofix_sigma_224_212_main_bathroom_drawing as drawing
 report=json.loads(REPORT.read_text())
 assert Path(tool.Ifc.get_path()).resolve()==SINGLE and not bpy.data.is_saved
 assert pkg.sha256(FORMAL)==FORMAL_HASH
 temporary=Path(tempfile.mkdtemp(prefix='duofix-pure-package-scene-'));temp_ifc=temporary/'scene.ifc'
 report.update(temporary_project_directory=str(temporary),runtime_inputs=[str(SINGLE),str(RECIPE),str(FORMAL)],legacy_ifc_used_for_runtime=False)
 write(REPORT,report)
 try:
  report['provider']={'version':list(bridge.bl_info['version']),'port':9890,'blender':bpy.app.version_string,'ifcopenshell':ifcopenshell.version,'source':pkg.record(bridge.__file__)}
  with bpy.context.temp_override(**context.view3d_override()):
   report['product_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(SINGLE),'overwrite':True,'reload':True})
  pure=ifcopenshell.open(str(SINGLE));assert check_pure(pure)==report['pure_pre_bonsai'];report['pure_saved_reloaded']=True
  shutil.copy2(FORMAL,temp_ifc);model=ifcopenshell.open(str(temp_ifc))
  for s in model.by_type('IfcExternallyDefinedSurfaceStyle'):
   if s.Location and (ROOT/s.Location).is_file():
    dest=temporary/s.Location;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/s.Location,dest)
  bodies={e.GlobalId:pkg.body_fingerprint(e) for e in model.by_type('IfcElement') if e.Representation}
  placements={e.GlobalId:pkg.placement(e).tolist() for e in model.by_type('IfcElement')}
  report['formal_physical_elements']=len(model.by_type('IfcElement'))
  recipe=json.loads(RECIPE.read_text());report['attachment']=attach_recipe(model,pure,recipe)
  count=len(list(model));attach_recipe(model,pure,recipe);assert len(list(model))==count
  report['second_attachment_created_entities']=0
  for v in report['views']:
   entity=model.by_guid(v['drawing_guid'])
   refs=[r.RelatingDocument for r in entity.HasAssociations if r.is_a('IfcRelAssociatesDocument')];assert len(refs)==1
   refs[0].Location=str(OUT/f"DUOFIX-SCENE-{v['view'].upper()}.svg")
  model.write(str(temp_ifc))
  assert bpy.ops.bim.load_project(filepath=str(temp_ifc),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
  outputs=[]
  for v in report['views']:
   entity=tool.Ifc.get().by_guid(v['drawing_guid']);tool.Ifc.get_object(entity) or tool.Drawing.import_drawing(entity)
   with bpy.context.temp_override(**context.view3d_override()):
    assert bpy.ops.bim.activate_drawing(drawing=entity.id(),should_view_from_camera=False)=={'FINISHED'}
    props=tool.Drawing.get_document_props();props.should_use_underlay_cache=False;props.should_use_linework_cache=False;props.should_use_annotation_cache=False
    result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
   assert result=={'FINISHED'}
   svg=OUT/f"DUOFIX-SCENE-{v['view'].upper()}.svg"
   styles=drawing.add_review_highlight(svg,GUID,v['annotation_guid'],v['view'])
   layers=drawing.inspect_target_layers(svg,GUID,v['annotation_guid']);assert layers['target_ifc_projection_group_count']==0
   outputs.append({'view':v['view'],'svg':pkg.record(svg),'styles':styles,'inspection':drawing.inspect_svg(svg),'target_layers':layers,'result':sorted(result),'operator':'bpy.ops.bim.create_drawing'})
  with bpy.context.temp_override(**context.view3d_override()):
   report['temporary_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(temp_ifc),'overwrite':True,'reload':True})
  reloaded=tool.Ifc.get()
  assert len(reloaded.by_type('IfcElement'))==report['formal_physical_elements']
  assert bodies=={e.GlobalId:pkg.body_fingerprint(e) for e in reloaded.by_type('IfcElement') if e.Representation}
  assert placements=={e.GlobalId:pkg.placement(e).tolist() for e in reloaded.by_type('IfcElement')}
  report.update(scene_saved_reloaded=True,all_formal_bodies_unchanged=True,all_formal_placements_unchanged=True,target_instances=1,scene_outputs=outputs,temporary_project=pkg.record(temp_ifc),verdict='pending_independent_visual_validation')
 except Exception:
  report.update(error=traceback.format_exc(),verdict='fail');raise
 finally:
  report['formal_sha256_after']=pkg.sha256(FORMAL);report['single_product']=pkg.record(SINGLE);write(REPORT,report)
  assert report['formal_sha256_after']==FORMAL_HASH

if __name__=='__main__':prepare()
