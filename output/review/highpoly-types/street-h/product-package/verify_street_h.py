"""Independent persisted IFC, approved points and SVG registration checks."""
import json,re,subprocess,hashlib,xml.etree.ElementTree as ET
from collections import Counter
import numpy as np
import ifcopenshell,ifcopenshell.validate
from shapely.geometry import LineString
from shapely.ops import unary_union
import migrate_street_h as task
def stable_errors(log):
 return Counter((e.get('attribute',''),re.sub(r'#\d+','#STEP',re.sub(r'\s+',' ',e['message'])),task.pkg.fingerprint(e['instance'])) for e in log.statements)
def cache_baseline():
 formal=ifcopenshell.open(str(task.FORMAL));logger=ifcopenshell.validate.json_logger();ifcopenshell.validate.validate(formal,logger)
 task.write(task.OUT/'formal-schema-baseline.json',{'formal_sha256':task.pkg.sha256(task.FORMAL),'count':len(logger.statements),'errors':[[*key,value] for key,value in stable_errors(logger).items()]})
def raster(svg,width=1300):
 png=svg.with_suffix('.png')
 subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),'--export-area-page','--export-background=white','--export-background-opacity=1',f'--export-width={width}',f'--export-filename={png}'],check=True,capture_output=True)
 return png
def verify():
 report=json.loads(task.REPORT.read_text());assert report['pure_saved_reloaded'] and report['scene_saved_reloaded']
 pure=ifcopenshell.open(str(task.SINGLE));assert task.check_pure(pure)==report['pure_pre_bonsai']
 scene=ifcopenshell.open(report['temporary_project']['path'])
 logger=ifcopenshell.validate.json_logger();ifcopenshell.validate.validate(pure,logger)
 assert not logger.statements,{'pure_schema_errors':len(logger.statements)}
 cached=json.loads((task.OUT/'formal-schema-baseline.json').read_text())
 assert cached['formal_sha256']==task.FORMAL_HASH==task.pkg.sha256(task.FORMAL)
 scene_logger=ifcopenshell.validate.json_logger();ifcopenshell.validate.validate(scene,scene_logger)
 baseline_errors=Counter({tuple(row[:3]):row[3] for row in cached['errors']});scene_errors=stable_errors(scene_logger)
 assert baseline_errors==scene_errors,{'new_errors':list((scene_errors-baseline_errors).items()),'removed_errors':list((baseline_errors-scene_errors).items())}
 report['baseline_inherited_errors']={'formal_count':cached['count'],'temporary_scene_count':len(scene_logger.statements),'new_errors':0,'comparison':'Counter of attribute, whitespace/STEP-normalized message and stable instance fingerprint','exact_baseline_match':True,'baseline_errors_modified':False,'formal_baseline_evidence':task.pkg.record(task.OUT/'formal-schema-baseline.json')}
 assert all(task.pkg.sha256(r['path'])==r['sha256'] for r in report['protected_files'])
 assert task.pkg.sha256(task.FORMAL)==task.FORMAL_HASH
 candidate=json.loads((task.PRODUCT/'candidate-representations.json').read_text())
 target=pure.by_guid(task.GUID);audits={};outputs=[];singles=[]
 for v in report['views']:
  view=v['view'];identifier=f'Approved{view.title()}'
  rep=next(r for r in target.Representation.Representations if r.RepresentationIdentifier==identifier)
  actual_points=[[list(p.Coordinates) for p in line.Points] for curves in rep.Items for line in curves.Elements]
  paths=candidate['views'][view]['proxy_paths_mm']
  expected_points=[[[float(x),float(y),0.] if view=='plan' else [float(x),0.,float(y)] if view=='front' else [0.,float(x),float(y)] for x,y in path] for path in paths]
  assert actual_points==expected_points
  audits[view]={'all_points_compared':True,'max_coordinate_error_mm':0,'paths':len(paths),'points':sum(map(len,paths)),'source':'approved candidate proxy_paths_mm','official_cad_used':False}
  svg=task.OUT/f'STREET-H-SUPPORT-SCENE-{view.upper()}.svg';root=ET.parse(svg).getroot()
  assert not [e for e in root.iter() if e.tag.endswith('image') or e.tag.endswith('use') or 'IfcSpace' in e.get('class','')]
  lines=[e for e in root.iter() if e.tag.endswith('line') and v['annotation_guid'] in e.get('class','')];assert lines
  assert not [e for e in root.iter() if e.get('{http://www.ifcopenshell.org/ns}guid')==task.GUID]
  camera=scene.by_guid(v['drawing_guid']);matrix=task.pkg.placement(camera)
  v['camera_matrix']=matrix.tolist()
  if view!='plan':assert abs(matrix[2,2])<.01
  project=np.linalg.inv(matrix)@task.pkg.placement(target)
  width=float(re.match(r'[0-9.]+',root.get('width')).group());height=float(re.match(r'[0-9.]+',root.get('height')).group())
  expected=[]
  for path in actual_points:
   points=[]
   for p in path:
    q=project@np.array([*p,1.]);points.append((q[0]/10+width/2,-q[1]/10+height/2))
   expected.append(LineString(points))
  actual=unary_union([LineString([(float(e.get('x1')),float(e.get('y1'))),(float(e.get('x2')),float(e.get('y2')))]) for e in lines])
  error=unary_union(expected).hausdorff_distance(actual)*10;assert error<.03,(view,error)
  png=raster(svg)
  outputs.append({'view':view,'svg':task.pkg.record(svg),'preview':task.pkg.record(png),'annotation_line_count':len(lines),'line_registration_error_world_mm':error,'camera_matrix_project_units':matrix.tolist(),'camera_world_positive_z':matrix[:3,2].tolist()})
  # Standalone approved drawing preview; geometry comes from saved pure IFC.
  axes=(0,1) if view=='plan' else (0,2) if view=='front' else (1,2)
  polylines=[[(p[axes[0]],p[axes[1]]) for p in line] for line in actual_points]
  flat=np.array([p for line in polylines for p in line]);lo=flat.min(0);hi=flat.max(0)
  content=[]
  for line in polylines:
   coords=' '.join(f'{p[0]-lo[0]+20},{hi[1]-p[1]+20}' for p in line)
   content.append(f'<polyline points="{coords}" fill="none" stroke="#171717" stroke-width="0.55"/>')
  single=task.OUT/f'STREET-H-SUPPORT-SINGLE-{view.upper()}.svg'
  single.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="{hi[0]-lo[0]+40}mm" height="{hi[1]-lo[1]+40}mm" viewBox="0 0 {hi[0]-lo[0]+40} {hi[1]-lo[1]+40}">'+''.join(content)+'</svg>')
  singles.append({'view':view,'svg':task.pkg.record(single),'preview':task.pkg.record(raster(single,900)),'source':'saved pure IFC approved representation'})
 report.update(independent_validation={'schema_errors':0,'pure_schema_errors':0,'temporary_scene_schema_errors':len(scene_logger.statements),'new_temporary_scene_schema_errors':0,'schema_checked_files':['pure','temporary_scene','formal_read_only_baseline'],'protected_files_unchanged':True,'approved_coordinates':audits,'approved_representation_projection':outputs},single_product_outputs=singles,scene_outputs=outputs,verdict='validated_pending_visual_check')
 report['courseEvidence'].update(source='research/bonsai-course/lessons/085000',source_sha256='85287cd1d85ea280d6c37f8ec4c6e29c32c950f6950b6f5c31ee981f6d93d04c',screenshot='research/bonsai-course/lessons/085000/screenshots/085000-01m59s-create-drawing-button.png',screenshot_sha256='a86279d39cef0f85feb9c37d5f823556e41e1f6e1ebd1cc634e564c3ebfd24b4',raw_source_available=False,screen_observation_claimed=False)
 report['provider_commit']='b0e67b1fb14ef5af93c4f747cb554b38511a844e'
 task.write(task.REPORT,report)
 # Composite is a labelled assembly of the rendered evidence, without geometry edits.
 subprocess.run(['/Users/jiaxinchen/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3',str(task.OUT/'render_preview.py'),'--with-single'],check=True)
 print(json.dumps({'schema_errors':0,'outputs':outputs,'approved_coordinates':audits},ensure_ascii=False))
def finalize():
 report=json.loads(task.REPORT.read_text());assert report['verdict']=='validated_pending_visual_check'
 log=max(task.Path('/var/folders/rz/d8p6s4y50ws53rd150nm2s3r0000gn/T/codex-task-blender-501').glob('1f7b5768c7c9b46d71df-*.log'),key=lambda p:p.stat().st_mtime)
 ready=json.loads(next(line.split('TASK_BLENDER_IFC_READY ',1)[1] for line in log.read_text().splitlines() if line.startswith('TASK_BLENDER_IFC_READY ')))
 assert ready['ifc']==str(task.SINGLE) and ready['sha256']==task.pkg.sha256(task.SINGLE)
 assert not ready['blend_saved'] and not ready['has_blend_warning']
 fresh={'status':'ready','launcher':'bonsai-launcher','task_owner':'street-h-pure-product','owner_marker':'--codex-task-owner=1f7b5768c7c9b46d71df','bridge_port':9896,'log':str(log),**ready}
 task.write(task.OUT/'fresh-session.json',fresh);report['fresh_session']=fresh
 report['workspace_safety']={'protected_staged_files':937,'protected_git_index_entries_sha256':hashlib.sha256(subprocess.check_output(['git','ls-files','--stage','-z'],cwd=task.ROOT)).hexdigest(),'new_files_only':True}
 assert report['workspace_safety']['protected_git_index_entries_sha256']=='06230b8b18fbb1fb5a1e4bd3e7749bc2086dda09ab4dadc3d78abee58712e67f'
 report['verdict']='pass';report['visual_verification']={'reviewed':True,'file':task.pkg.record(task.OUT/'STREET-H-SUPPORT-contact-sheet.png'),'historical_symbols_or_spatial_fills':False,'black_support_linework_and_grey_physical_context':True,'scope':'STREET-H sink-holder support only; complete Street top is grey scene context','scene_approval_status':'pending'}
 task.write(task.REPORT,report)
 manifest={'schema_version':1,'profile_key':'street-h','display_name':'STREET-H sink-holder support subcomponent','status':'complete','validation_verdict':'pass','approval_status':'approved','single_product_approval_status':'approved','scene_approval_status':'pending','label_zh':'单品已通过、场景待验收','target_global_id':task.GUID,'persistent_ifc':task.pkg.record(task.SINGLE),'source_kind':'geometry_derived_simplified_proxy','source_label_zh':task.SOURCE_LABEL,'official_cad_used':False,'official_product_url':'https://www.antoniolupi.it/en/products/sinks/street','approval_record':task.pkg.record(task.APPROVAL),'candidate':task.pkg.record(task.PRODUCT/'candidate-representations.json'),'geometry':{'physical_products':1,'original_body_representation_count':2,'approved_views':['plan','front','side'],'annotations':0,'drawing_cameras':0,'spatial_geometry':0,'placement_matrix_project_units':report['pure_pre_bonsai']['placement'],'units_to_m':.001},'scene_recipe':task.pkg.record(task.RECIPE),'scene_outputs':report['scene_outputs'],'single_product_outputs':report['single_product_outputs'],'contact_sheet':task.pkg.record(task.OUT/'STREET-H-SUPPORT-contact-sheet.png'),'external_style_dependencies':report['necessary_external_style_dependencies'],'formal_write_allowed':False,'formal_sha256':task.FORMAL_HASH,'legacy_source_retained':True,'validation':'validation.json','scope':'STREET-H support subcomponent, not the complete STREET basin/top; official parent-family DXF is reference only.'}
 task.write(task.OUT/'manifest.json',manifest)
 task.write(task.OUT/'operator-result.json',{
  'task':{'target':task.GUID,'scope':'STREET-H sink-holder support subcomponent only','artifact':str(task.SINGLE),'versions':report['provider']},
  'courseEvidence':report['courseEvidence'],
  'plan':['Inspect approved black-line candidate and source IFC','Preserve Body and coordinate dependencies in pure IFC','Externalize camera and scene filters to JSON recipe','Attach only approved expressions to identity-matched formal copy','Create actual Bonsai drawings, save and reload','Compare persisted geometry, schema baseline, SVG registration and rendered evidence'],
  'preState':{'unchanged_source_ifc':task.pkg.record(task.OLD),'source_state':task.snapshot(ifcopenshell.open(str(task.OLD))),'formal_sha256':task.FORMAL_HASH,'save_boundary':'Only new product-package files and disposable temporary full-project copies'},
  'execution':{'native_drawings':report['scene_outputs'],'attachment':report['final_runtime_attachment'],'metadata_correction':report['attachment_metadata_correction']},
  'persistence':{'pure':report['product_metadata_bonsai_save_reload'],'temporary_scene':report['temporary_metadata_bonsai_save_reload'],'fresh_pure_session':fresh},
  'postState':task.check_pure(ifcopenshell.open(str(task.SINGLE))),
  'outputs':{'pure':task.pkg.record(task.SINGLE),'recipe':task.pkg.record(task.RECIPE),'manifest':task.pkg.record(task.OUT/'manifest.json'),'validation':task.pkg.record(task.REPORT)},
  'visual':report['visual_verification'],'verdict':'pass','scene_approval_status':'pending','remaining_uncertainty':'New scene cameras and installation view remain awaiting user acceptance; formal schema baseline defects are inherited, not repaired.'})
 task.write(task.OUT/'cleanup-proposal.json',{'status':'proposal_only_no_deletion','requires_user_confirmation':True,'legacy_files':[task.pkg.record(task.OLD)],'temporary_directories':[str(p) for p in task.Path(report['temporary_project_directory']).parent.glob('street-h-support-scene-*') if p.is_dir()],'keep':[str(task.SINGLE),str(task.RECIPE),'scene SVG/PNG','single-product SVG/PNG','source and approval evidence','validation files']})
 task.write(task.OUT/'handoff.json',{'product':'street-h','status':'complete','validation_verdict':'pass','scene_approval_status':'pending','label_zh':'单品已通过、场景待验收','package':'manifest.json','ifc':'STREET-H-SUPPORT-product.ifc','scene_recipe':'scene-recipe.json','validation':'validation.json','fresh_session':'fresh-session.json','pid':ready['pid'],'new_files_only':True,'staged_baseline_modified':False,'legacy_deleted':False,'formal_modified':False,'bridge_port':9896,'scope':'STREET-H sink-holder support subcomponent only'})
if __name__=='__main__':verify()
