"""Approved STREET-H sink-holder subcomponent, never the complete Street top."""
from pathlib import Path
import sys,json,shutil,tempfile,traceback
import numpy as np
import ifcopenshell,ifcopenshell.api.pset
import ifcopenshell.util.element as eu
OUT=Path(__file__).resolve().parent
PRODUCT=OUT.parent
ROOT=OUT.parents[4]
sys.path.insert(0,str(ROOT/'pipeline/scripts'))
import review_product_package as pkg
from pure_product_package import build_pure_package,attach_recipe as shared_attach_recipe
GUID='2ajpw0I9n1dBypfISg3ejX'
FORMAL=ROOT/'2504 GBTB Yanlord Zhuhai.ifc'
FORMAL_HASH='7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c'
OLD=PRODUCT/'antoniolupi-Street-H-bonsai-isolated.ifc'
SINGLE=OUT/'STREET-H-SUPPORT-product.ifc'
RECIPE=OUT/'scene-recipe.json'
REPORT=OUT/'validation.json'
APPROVAL=ROOT/'pipeline/decisions/street-h-drawing-approval.json'
SOURCE_LABEL='基于原始高模几何生成的简化图纸表达'
def write(path,data):Path(path).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
def snapshot(model):
 t=model.by_guid(GUID)
 return {'body':pkg.body_fingerprint(t),'representations':pkg.fingerprint(t.Representation),'placement':pkg.placement(t).tolist(),'physical_elements':sorted(e.GlobalId for e in model.by_type('IfcElement')),'annotation_count':len(model.by_type('IfcAnnotation')),'unit_scale_to_m':pkg.unit_util.calculate_unit_scale(model)}
def check_pure(model):
 s=snapshot(model)
 assert s['physical_elements']==[GUID] and s['annotation_count']==0
 assert not model.by_type('IfcGroup')
 assert not [e for e in model.by_type('IfcSpatialElement') if e.Representation]
 assert not [p for p in model.by_type('IfcPropertySingleValue') if p.Name in ('Include','Exclude')]
 assert not [p for p in model.by_type('IfcPropertySet') if p.Name=='EPset_Drawing']
 assert sorted(r.RepresentationIdentifier for r in model.by_guid(GUID).Representation.Representations)==['ApprovedFront','ApprovedPlan','ApprovedSide','Body','Body']
 return s
def pset(model,entity,name,values):
 p=ifcopenshell.api.pset.add_pset(model,product=entity,name=name)
 ifcopenshell.api.pset.edit_pset(model,pset=p,properties=values)
def axis(model,location,z,x):
 return model.create_entity('IfcAxis2Placement3D',Location=model.create_entity('IfcCartesianPoint',Coordinates=tuple(map(float,location))),Axis=model.create_entity('IfcDirection',DirectionRatios=tuple(map(float,z))),RefDirection=model.create_entity('IfcDirection',DirectionRatios=tuple(map(float,x))))
def align_spatial_metadata(model,formal):
 changes=[]
 for rel in model.by_type('IfcRelAggregates'):
  if not (rel.RelatingObject.is_a('IfcProject') or rel.RelatingObject.is_a('IfcSpatialElement')):continue
  children={e.GlobalId for e in rel.RelatedObjects}
  matches=[r for r in formal.by_type('IfcRelAggregates') if r.RelatingObject.GlobalId==rel.RelatingObject.GlobalId and children.issubset({e.GlobalId for e in r.RelatedObjects})]
  assert len(matches)==1,(rel,matches)
  if rel.GlobalId!=matches[0].GlobalId:
   changes.append({'old_guid':rel.GlobalId,'formal_guid':matches[0].GlobalId,'parent':rel.RelatingObject.GlobalId,'children':sorted(children)})
   rel.GlobalId=matches[0].GlobalId
 return changes
def align_target_relations(model,formal):
 changes=[]
 target=model.by_guid(GUID);original=formal.by_guid(GUID)
 for attr,owner in [('IsTypedBy','RelatingType'),('ContainedInStructure','RelatingStructure')]:
  originals=list(getattr(original,attr));relations=list(getattr(target,attr))
  assert len(originals)==1
  for rel in relations:
   assert getattr(rel,owner).GlobalId==getattr(originals[0],owner).GlobalId
   if rel.GlobalId!=originals[0].GlobalId:
    changes.append({'relationship':rel.is_a(),'old_guid':rel.GlobalId,'formal_guid':originals[0].GlobalId})
    rel.GlobalId=originals[0].GlobalId
 return changes
def reuse_added_applications(model,original_apps):
 replacements=[]
 for app in list(model.by_type('IfcApplication')):
  if app.id() in original_apps:continue
  matches=[a for a in original_apps.values() if (a.ApplicationIdentifier,a.ApplicationFullName,a.Version)==(app.ApplicationIdentifier,app.ApplicationFullName,app.Version)]
  if not matches:continue
  chosen=next((a for a in matches if pkg.fingerprint(a)==pkg.fingerprint(app)),matches[0])
  replacements.append({'imported_id':app.id(),'reused_formal_id':chosen.id(),'identifier':app.ApplicationIdentifier})
  for inverse in list(model.get_inverse(app)):eu.replace_attribute(inverse,app,chosen)
  model.remove(app)
 return replacements
def attach_recipe(model,pure,recipe):
 original_apps={a.id():a for a in model.by_type('IfcApplication')}
 result=shared_attach_recipe(model,pure,recipe)
 result['reused_formal_application_metadata']=reuse_added_applications(model,original_apps)
 return result
def prepare():
 import street_h_drawing_ifc as approved
 assert pkg.sha256(FORMAL)==FORMAL_HASH
 auth=json.loads((ROOT/'output/review/approved-product-library/migration-authorization.json').read_text())
 assert 'street-h' in auth['products'] and auth['pure_product_ifc_write_authorized']
 approval=json.loads(APPROVAL.read_text());assert approval['status']=='approved'
 assert pkg.sha256(PRODUCT/'manifest.json')==approval['candidate_manifest_sha256']
 candidate=json.loads((PRODUCT/'candidate-representations.json').read_text())
 paths=approved.candidate_paths(candidate,json.loads((PRODUCT/'official-source/source-access-record.json').read_text()))
 source=ifcopenshell.open(str(OLD));formal=ifcopenshell.open(str(FORMAL));target=source.by_guid(GUID)
 assert pkg.body_fingerprint(target)==pkg.body_fingerprint(formal.by_guid(GUID))
 assert np.array_equal(pkg.placement(target),pkg.placement(formal.by_guid(GUID)))
 protected=[pkg.record(p) for p in [FORMAL,OLD,APPROVAL,PRODUCT/'manifest.json',PRODUCT/'candidate-representations.json',*[p for p in (PRODUCT/'official-source').rglob('*') if p.is_file()],*[PRODUCT/f'{v}.svg' for v in paths]]]
 context_guids=json.loads((OUT/'context-selection.json').read_text())['guids']
 camera_specs={
  'plan':{'loc':[-6175.,-550.,2350.],'z':[0.,0.,1.],'x':[1.,0.,0.],'width':1500.,'height':1600.,'depth':3000.},
  'front':{'loc':[-6345.,-1120.,900.],'z':[0.,-1.,0.],'x':[1.,0.,0.],'width':1300.,'height':1800.,'depth':1200.},
  'side':{'loc':[-5750.,-400.,900.],'z':[1.,0.,0.],'x':[0.,1.,0.],'width':1300.,'height':1800.,'depth':1200.}}
 original=next(d for d in formal.by_type('IfcAnnotation') if d.Name=='FFL PLAN')
 specs={};views=[]
 for view,camera in camera_specs.items():
  clone=pkg.ScopedCopy(formal,source,{original},skip_inverse_ids=[original.id()])
  drawing=clone.copy(original);drawing.GlobalId=ifcopenshell.guid.new();drawing.Name=f'STREET-H-SUPPORT-SCENE-{view.upper()}'
  drawing.Description='STREET-H sink-holder support subcomponent; newly framed scene pending user acceptance'
  drawing.ObjectPlacement=source.create_entity('IfcLocalPlacement',RelativePlacement=axis(source,camera['loc'],camera['z'],camera['x']))
  block=next(e for e in source.traverse(drawing.Representation) if e.is_a('IfcBlock'))
  block.XLength=camera['width'];block.YLength=camera['height'];block.ZLength=camera['depth']
  block.Position.Location.Coordinates=(-camera['width']/2,-camera['height']/2,-camera['depth'])
  target_view='PLAN_VIEW' if view=='plan' else 'ELEVATION_VIEW'
  values={k:v for k,v in eu.get_pset(original,'EPset_Drawing').items() if k not in ('id','Include','Exclude')}
  selected_context=context_guids if view=='plan' else ['2amSzxJIb9ceddTJOslmHk','04rs0EDjn2EvxytEQSxWRB','3OVQygdDn17huGOgJJFTOY','0wbKH$_Qr2_A4hpKs0PLXw','1FgLPMw$5B4wBH2ySMkXE1']
  values.update(TargetView=target_view,Scale='1/10',HumanScale='1:10',HasAnnotation=True,GlobalReferencing=False,Include=','.join(selected_context),Exclude=','.join([GUID,'IfcSpace','IfcGrid','IfcBuildingStorey',*[a.GlobalId for a in formal.by_type('IfcAnnotation')]]))
  for key in ('Stylesheet','Markers','Symbols','Patterns','ShadingStyles'):values[key]=str(ROOT/values[key])
  pset(source,drawing,'EPset_Drawing',values)
  doc=source.create_entity('IfcDocumentReference',Location=str(OUT/f'{drawing.Name}.svg'),Identification=drawing.Name,Name=drawing.Name)
  source.create_entity('IfcRelAssociatesDocument',GlobalId=ifcopenshell.guid.new(),RelatedObjects=[drawing],RelatingDocument=doc)
  context=approved.shared.representation_context(source,'Annotation',target_view)
  rep=approved.shared.curve_representation(source,context,'Annotation',view,paths[view])
  color=source.create_entity('IfcColourRgb',Red=0.08,Green=0.08,Blue=0.08)
  style=source.create_entity('IfcCurveStyle',Name='STREET-H approved geometry-derived black',CurveWidth=source.create_entity('IfcPositiveLengthMeasure',0.25),CurveColour=color,ModelOrDraughting=True)
  source.create_entity('IfcStyledItem',Item=rep.Items[0],Styles=[style])
  ann=source.create_entity('IfcAnnotation',GlobalId=ifcopenshell.guid.new(),Name=f'STREET-H support approved {view}',ObjectType='LINEWORK',ObjectPlacement=target.ObjectPlacement,Representation=source.create_entity('IfcProductDefinitionShape',Representations=[rep]))
  pset(source,ann,'EPset_Annotation',{'Classes':'review-target-street-h geometry-derived-support','TargetGlobalId':GUID,'SourceKind':candidate['source_kind']})
  group=source.create_entity('IfcGroup',GlobalId=ifcopenshell.guid.new(),Name=drawing.Name,ObjectType='DRAWING')
  source.create_entity('IfcRelAssignsToGroup',GlobalId=ifcopenshell.guid.new(),RelatedObjects=[drawing,ann],RelatingGroup=group)
  specs[view]={'drawing_guid':drawing.GlobalId,'annotation_guid':ann.GlobalId}
  views.append({'view':view,**specs[view],'camera_matrix':pkg.placement(drawing).tolist(),'approved_path_count':len(paths[view])})
 pset(source,target,'Pset_StreetHSupportApprovedSource',{'SourceKind':candidate['source_kind'],'SourceLabelZh':SOURCE_LABEL,'Scope':'STREET-H sink-holder support subcomponent only; not complete Street basin or top','OfficialCadUsed':False,'CandidateSha256':pkg.sha256(PRODUCT/'candidate-representations.json'),'SingleProductApprovalStatus':'approved','SceneApprovalStatus':'pending','FormalIfcWriteAllowed':False})
 relation_alignment=align_spatial_metadata(source,formal)+align_target_relations(source,formal)
 pure,recipe,audit=build_pure_package(source,GUID,specs)
 audit['spatial_relation_identity_alignment']=relation_alignment
 pure.write(str(SINGLE));write(RECIPE,recipe)
 dependencies=[]
 for loc in sorted({s.Location for s in pure.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}):
  dest=OUT/loc;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/loc,dest);dependencies.append({'source':pkg.record(ROOT/loc),'packaged':pkg.record(dest)})
 write(REPORT,{'product':'street-h','display_name':'STREET-H sink-holder support subcomponent','target_guid':GUID,'verdict':'pending_bonsai_scene','formal_sha256_before':FORMAL_HASH,'protected_files':protected,'preState':snapshot(source),'pure_pre_bonsai':check_pure(ifcopenshell.open(str(SINGLE))),'extraction':audit,'views':views,'single_product':pkg.record(SINGLE),'scene_recipe':pkg.record(RECIPE),'necessary_external_style_dependencies':dependencies,'source_kind':candidate['source_kind'],'source_label_zh':SOURCE_LABEL,'single_product_approval_status':'approved','scene_approval_status':'pending','formal_write_allowed':False,'cleanup_performed':False,'courseEvidence':{'mode':'embedded-course-index','lesson':'085000','timestamp':'01:59','label':'course_fact','fact':'Create Drawing produces SVG; persisted model and visible output must both be verified.'}})
def style_svg(svg,annotation_guid,view):
 import xml.etree.ElementTree as ET
 raw=pkg.record(svg);tree=ET.parse(svg);root=tree.getroot();target=context=0
 for e in root.iter():
  if e.tag.rsplit('}',1)[-1] not in ('line','path','polyline','polygon','rect','circle','ellipse'):continue
  selected=annotation_guid in e.get('class','')
  e.set('style',f'stroke:{"#171717" if selected else "#a6adb4"};stroke-width:{0.28 if selected else 0.14};fill:none')
  if selected:e.set('data-target-global-id',GUID);target+=1
  else:context+=1
 assert target and context,(view,target,context)
 root.set('data-source-kind','geometry_derived_simplified_proxy');root.set('data-scene-approval-status','pending')
 tree.write(svg,encoding='utf-8',xml_declaration=True)
 return {'raw':raw,'black_geometry_count':target,'grey_context_geometry_count':context,'post_style_only':True,'geometry_moved_removed_or_redrawn':False}
def scene():
 import bpy,bonsai_bridge as bridge
 from bonsai import tool
 import create_wd03_wardrobe_scene_drawings as context
 report=json.loads(REPORT.read_text())
 assert Path(tool.Ifc.get_path()).resolve()==SINGLE.resolve() and not bpy.data.is_saved
 if report.get('error'):
  report.setdefault('resolved_attempts',[]).append({'temporary_project_directory':report['temporary_project_directory'],'error':report.pop('error'),'adjustment':'Retain parent basin and physical mounting walls/floor; omit secondary vanity/plumbing context that causes OpenCASCADE elevation failure'})
 temporary=Path(tempfile.mkdtemp(prefix='street-h-support-scene-'));temp_ifc=temporary/'scene.ifc'
 report.update(temporary_project_directory=str(temporary),runtime_inputs=[str(SINGLE),str(RECIPE),str(FORMAL)],legacy_ifc_used_for_runtime=False)
 write(REPORT,report)
 try:
  report['provider']={'status':'supported','version':list(bridge.bl_info['version']),'port':9896,'pid':__import__('os').getpid(),'blender':bpy.app.version_string,'ifcopenshell':ifcopenshell.version,'bridge_source':pkg.record(bridge.__file__)}
  with bpy.context.temp_override(**context.view3d_override()):report['product_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(SINGLE),'overwrite':True,'reload':True})
  pure=ifcopenshell.open(str(SINGLE));assert check_pure(pure)==report['pure_pre_bonsai'];report['pure_saved_reloaded']=True
  assert pkg.sha256(FORMAL)==FORMAL_HASH
  shutil.copy2(FORMAL,temp_ifc);model=ifcopenshell.open(str(temp_ifc))
  baseline={e.GlobalId:(pkg.body_fingerprint(e) if e.Representation else None,pkg.placement(e).tolist()) for e in model.by_type('IfcElement')};assert len(baseline)==714
  report['formal_physical_elements']=714
  for loc in sorted({s.Location for s in model.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}):
   if (ROOT/loc).is_file():dest=temporary/loc;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/loc,dest)
  recipe=json.loads(RECIPE.read_text());report['attachment']=attach_recipe(model,pure,recipe)
  count=len(list(model));attach_recipe(model,pure,recipe);assert len(list(model))==count;report['second_attachment_created_entities']=0
  model.write(str(temp_ifc))
  assert bpy.ops.bim.load_project(filepath=str(temp_ifc),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
  outputs=[]
  for v in report['views']:
   drawing=tool.Ifc.get().by_guid(v['drawing_guid']);tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
   with bpy.context.temp_override(**context.view3d_override()):
    assert bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)=={'FINISHED'}
    props=tool.Drawing.get_document_props();props.should_use_underlay_cache=props.should_use_linework_cache=props.should_use_annotation_cache=False
    result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
   assert result=={'FINISHED'}
   svg=OUT/f"STREET-H-SUPPORT-SCENE-{v['view'].upper()}.svg"
   styles=style_svg(svg,v['annotation_guid'],v['view'])
   outputs.append({'view':v['view'],'svg':pkg.record(svg),'styles':styles,'operator':'bpy.ops.bim.create_drawing','result':sorted(result)})
   report['scene_outputs']=outputs;write(REPORT,report)
  with bpy.context.temp_override(**context.view3d_override()):report['temporary_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(temp_ifc),'overwrite':True,'reload':True})
  reloaded=ifcopenshell.open(str(temp_ifc))
  assert baseline=={e.GlobalId:(pkg.body_fingerprint(e) if e.Representation else None,pkg.placement(e).tolist()) for e in reloaded.by_type('IfcElement')}
  report.update(scene_saved_reloaded=True,all_formal_bodies_unchanged=True,all_formal_placements_unchanged=True,target_instances=1,temporary_project=pkg.record(temp_ifc),verdict='pending_independent_visual_validation')
  assert bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
 except Exception:report['error']=traceback.format_exc();report['verdict']='fail';raise
 finally:
  report['formal_sha256_after']=pkg.sha256(FORMAL);report['single_product']=pkg.record(SINGLE);write(REPORT,report);assert report['formal_sha256_after']==FORMAL_HASH
def run_background():
 import bpy
 def run():
  try:scene();write(OUT/'scene-status.json',{'status':'complete'})
  except Exception:write(OUT/'scene-status.json',{'status':'failed','error':traceback.format_exc()})
  return None
 write(OUT/'scene-status.json',{'status':'running'});bpy.app.timers.register(run,first_interval=.5)
def adjust_front_camera():
 from pure_product_package import graph_from_json,graph_to_json
 recipe=json.loads(RECIPE.read_text());model=graph_from_json(recipe['graph'])
 cam=model.by_guid(recipe['views']['front']['drawing_guid'])
 cam.ObjectPlacement.RelativePlacement.Location.Coordinates=(-6345.,-1120.,900.)
 recipe['graph']=graph_to_json(model);write(RECIPE,recipe)
 report=json.loads(REPORT.read_text());report.setdefault('resolved_attempts',[]).append({'temporary_project_directory':report['temporary_project_directory'],'error':report.pop('error',None),'adjustment':'Front camera inside dry zone but ahead of basin, avoiding camera cut through parent basin high-poly Body'})
 report['scene_recipe']=pkg.record(RECIPE);write(REPORT,report)
def diagnostic_elevations(extra=()):
 import bpy
 from bonsai import tool
 import create_wd03_wardrobe_scene_drawings as context
 model=tool.Ifc.get();recipe=json.loads(RECIPE.read_text())
 keep=['2amSzxJIb9ceddTJOslmHk','04rs0EDjn2EvxytEQSxWRB','3OVQygdDn17huGOgJJFTOY','0wbKH$_Qr2_A4hpKs0PLXw',*extra]
 result={}
 for view in ['front','side']:
  drawing=model.by_guid(recipe['views'][view]['drawing_guid']);data=eu.get_pset(drawing,'EPset_Drawing')
  ifcopenshell.api.pset.edit_pset(model,pset=model.by_id(data['id']),properties={'Include':','.join(keep)})
  tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
  with bpy.context.temp_override(**context.view3d_override()):
   bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)
   props=tool.Drawing.get_document_props();props.should_use_underlay_cache=props.should_use_linework_cache=props.should_use_annotation_cache=False
   result[view]=sorted(bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False))
 write(OUT/'diagnostic-elevations.json',{'result':result,'physical_context':keep})
def apply_context_selection():
 from pure_product_package import graph_from_json,graph_to_json
 recipe=json.loads(RECIPE.read_text());model=graph_from_json(recipe['graph'])
 keep=json.loads((OUT/'diagnostic-elevations.json').read_text())['physical_context']
 for view in ['front','side']:
  drawing=model.by_guid(recipe['views'][view]['drawing_guid']);data=eu.get_pset(drawing,'EPset_Drawing')
  ifcopenshell.api.pset.edit_pset(model,pset=model.by_id(data['id']),properties={'Include':','.join(keep)})
 recipe['graph']=graph_to_json(model);write(RECIPE,recipe)
 report=json.loads(REPORT.read_text());report['scene_recipe']=pkg.record(RECIPE)
 report['elevation_context_policy']={'guids':keep,'evidence':'Both native Front and Side CreateDrawing finished with actual parent basin, mounting walls and floor; secondary vanity/plumbing context is excluded after OpenCASCADE failure, with no changes to any physical Body.','scene_approval_status':'pending'}
 write(REPORT,report)
def align_existing_package():
 from pure_product_package import graph_from_json,graph_to_json
 formal=ifcopenshell.open(str(FORMAL));pure=ifcopenshell.open(str(SINGLE));before=snapshot(pure)
 changes=align_spatial_metadata(pure,formal)+align_target_relations(pure,formal);assert snapshot(pure)==before
 pure.write(str(SINGLE));recipe=json.loads(RECIPE.read_text());graph=graph_from_json(recipe['graph'])
 recipe_changes=align_spatial_metadata(graph,formal)+align_target_relations(graph,formal);recipe['graph']=graph_to_json(graph);write(RECIPE,recipe)
 report=json.loads(REPORT.read_text());report['extraction']['spatial_relation_identity_alignment']={'pure':changes,'recipe':recipe_changes,'policy':'Use authoritative spatial relation identities to reuse existing project relations; no Body or Placement edits'}
 report['single_product']=pkg.record(SINGLE);report['scene_recipe']=pkg.record(RECIPE);write(REPORT,report)
def fix_added_scene_metadata():
 formal=ifcopenshell.open(str(FORMAL));report=json.loads(REPORT.read_text())
 path=Path(report['temporary_project']['path']);model=ifcopenshell.open(str(path));target=model.by_guid(GUID)
 before={e.GlobalId:(e.Representation.id() if e.Representation else None,e.ObjectPlacement.id() if e.ObjectPlacement else None) for e in model.by_type('IfcElement')}
 removed=[]
 for attr in ['IsTypedBy','ContainedInStructure']:
  original_ids={r.GlobalId for r in getattr(formal.by_guid(GUID),attr)}
  for rel in list(getattr(target,attr)):
   if rel.GlobalId not in original_ids:
    removed.append({'relationship':rel.is_a(),'guid':rel.GlobalId,'reason':'Duplicate isolated-source relation introduced by attachment; authoritative original retained'})
    model.remove(rel)
 original_apps={a.id():model.by_id(a.id()) for a in formal.by_type('IfcApplication')}
 assert all(pkg.fingerprint(formal.by_id(i))==pkg.fingerprint(a) for i,a in original_apps.items())
 reused=reuse_added_applications(model,original_apps)
 after={e.GlobalId:(e.Representation.id() if e.Representation else None,e.ObjectPlacement.id() if e.ObjectPlacement else None) for e in model.by_type('IfcElement')}
 assert before==after and len(after)==714
 model.write(str(path));report['temporary_project']=pkg.record(path)
 report['attachment_metadata_correction']={'removed_new_duplicate_relations':removed,'reused_formal_applications':reused,'physical_representation_and_placement_references_unchanged':True,'original_schema_defects_modified':False,'scene_svg_geometry_changed':False}
 # Rebuild once more in memory from the final three inputs to prove that the
 # product-specific adapter now prevents these metadata duplicates at source.
 pure=ifcopenshell.open(str(SINGLE));recipe=json.loads(RECIPE.read_text())
 applications_before=len(formal.by_type('IfcApplication'))
 attachment=attach_recipe(formal,pure,recipe)
 assert len(formal.by_type('IfcApplication'))==applications_before
 assert len(formal.by_guid(GUID).IsTypedBy)==len(formal.by_guid(GUID).ContainedInStructure)==1
 count=len(list(formal));attach_recipe(formal,pure,recipe);assert len(list(formal))==count
 report['final_runtime_attachment']=attachment;report['second_attachment_created_entities']=0
 report['temporary_scene_pending_metadata_reload']=True;write(REPORT,report)
def persist_corrected_metadata():
 import bpy,bonsai_bridge as bridge
 from bonsai import tool
 import create_wd03_wardrobe_scene_drawings as context
 def run():
  try:
   report=json.loads(REPORT.read_text());path=report['temporary_project']['path']
   bpy.ops.bim.load_project(filepath=path,should_start_fresh_session=False,use_relative_path=False)
   with bpy.context.temp_override(**context.view3d_override()):report['temporary_metadata_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':path,'overwrite':True,'reload':True})
   report['temporary_project']=pkg.record(path);report['temporary_scene_pending_metadata_reload']=False
   bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False)
   with bpy.context.temp_override(**context.view3d_override()):report['product_metadata_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(SINGLE),'overwrite':True,'reload':True})
   assert check_pure(tool.Ifc.get())==report['pure_pre_bonsai']
   report['single_product']=pkg.record(SINGLE);write(REPORT,report);write(OUT/'metadata-status.json',{'status':'complete'})
  except Exception:write(OUT/'metadata-status.json',{'status':'failed','error':traceback.format_exc()})
  return None
 write(OUT/'metadata-status.json',{'status':'running'});bpy.app.timers.register(run,first_interval=.5)
if __name__=='__main__':prepare()
