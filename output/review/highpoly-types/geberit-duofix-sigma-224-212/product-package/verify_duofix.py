"""Independent IFC/schema/approved SVG geometry verification and preview export."""
import json, copy, subprocess, xml.etree.ElementTree as ET
import ifcopenshell, ifcopenshell.validate
from shapely.geometry import LineString
from shapely.ops import unary_union
import migrate_duofix as task
pkg=task.pkg

def lines(path,guid):
 root=ET.parse(path).getroot()
 result=[e for e in root.iter() if e.tag.endswith('line') and guid in e.get('class','')]
 assert result and not [e for e in root.iter() if e.tag.endswith('image')]
 return result,unary_union([LineString([(float(e.get('x1')),float(e.get('y1'))),(float(e.get('x2')),float(e.get('y2')))]) for e in result])

def render(svg,width=1400):
 png=svg.with_suffix('.png')
 subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),'--export-area-page','--export-background=white','--export-background-opacity=1',f'--export-width={width}',f'--export-filename={png}'],check=True,capture_output=True)
 return pkg.record(png)

def main():
 report=json.loads(task.REPORT.read_text());assert report['pure_saved_reloaded'] and report['scene_saved_reloaded']
 model=ifcopenshell.open(str(task.SINGLE));state=task.check_pure(model);assert state==report['pure_pre_bonsai']
 logger=ifcopenshell.validate.json_logger();ifcopenshell.validate.validate(model,logger)
 assert not logger.statements,logger.statements
 assert all(pkg.sha256(x['path'])==x['sha256'] for x in report['protected_files'])
 from pure_product_package import graph_from_json
 recipe=json.loads(task.RECIPE.read_text());graph=graph_from_json(recipe['graph'])
 assert all(not e.Representation for e in graph.by_type('IfcElement'))
 assert all(not graph.by_guid(s['annotation_guid']).Representation for s in recipe['views'].values())
 formal=ifcopenshell.open(str(task.FORMAL));temporary=ifcopenshell.open(report['temporary_project']['path'])
 def assembly_memberships(m):
  return sorted((r.GlobalId,r.RelatingObject.GlobalId) for r in m.by_guid(task.GUID).Decomposes)
 assert assembly_memberships(formal)==assembly_memberships(temporary)
 comparisons=[]
 for v in report['views']:
  svg=task.OUT/f"DUOFIX-SCENE-{v['view'].upper()}.svg"
  oldlines,old=lines(v['old_svg']['path'],v['annotation_guid']);newlines,new=lines(svg,v['annotation_guid'])
  oldroot=ET.parse(v['old_svg']['path']).getroot();scale=oldroot.get('data-scale')
  # Bonsai stores the paper:model ratio, e.g. 1:25; camera EPset is independently recorded.
  import ifcopenshell.util.element as eu
  camera_scale=eu.get_psets(graph.by_guid(v['drawing_guid']))['EPset_Drawing']['Scale']
  numerator,denominator=map(float,camera_scale.split('/'));scale_factor=denominator/numerator
  error=old.segmentize(.1).hausdorff_distance(new.segmentize(.1))*scale_factor
  assert error<.01,{'view':v['view'],'world_mm':error}
  preview=render(svg)
  thumb=ET.Element('{http://www.w3.org/2000/svg}svg',width='1000',height='760')
  minx,miny,maxx,maxy=new.bounds;margin=max(maxx-minx,maxy-miny)*.08
  thumb.set('viewBox',f'{minx-margin} {miny-margin} {maxx-minx+margin*2} {maxy-miny+margin*2}')
  for e in newlines:thumb.append(copy.deepcopy(e))
  single=task.OUT/f"DUOFIX-SINGLE-{v['view'].upper()}.svg"
  ET.ElementTree(thumb).write(single,encoding='utf-8',xml_declaration=True)
  comparisons.append({'view':v['view'],'approved_edges':len(oldlines),'new_edges':len(newlines),
   'max_linework_hausdorff_error_world_mm':error,'tolerance_mm':.01,'scale_metadata':scale,
   'svg':pkg.record(svg),'preview':preview,'single_product_svg':pkg.record(single),'single_product_thumbnail':render(single,1000)})
 report['independent_validation']={'pure_state':state,'schema_errors':0,'approved_scene_comparison':comparisons,
  'protected_files_unchanged':True,'formal_sha256':pkg.sha256(task.FORMAL),
  'recipe_has_no_body_or_linework_geometry':True,'project_assembly_memberships_unchanged':assembly_memberships(temporary)}
 report['verdict']='validated_pending_visual_check';task.write(task.REPORT,report)
 source=task.PRODUCT/'official-source/source-access-record.json';source_data=json.loads(source.read_text())
 task.write(task.OUT/'manifest.json',{'schema_version':1,'profile_key':'geberit-duofix-sigma-224-212',
  'display_name':'Geberit Duofix Sigma 224.212.00.2','status':'validated_pending_visual_check',
  'target_global_id':task.GUID,'persistent_ifc':pkg.record(task.SINGLE),'source_kind':'native_dwg',
  'source_label_zh':'Geberit 精确型号 224.212.00.2 官方 G/A/L 原生 DWG 蓝线',
  'official_native_dwg':source_data['official_native_dwg'],'source_access_record':pkg.record(source),
  'approval_record':pkg.record(task.APPROVAL),'scene_approval_record':pkg.record(task.ROOT/'pipeline/decisions/highpoly-subtask-completion-2026-09-05.json'),
  'geometry':{'physical_products':1,'body_representations':3,'original_body_preserved':True,
   'approved_views':['plan','front','side'],'annotations':0,'drawing_cameras':0,'spatial_geometry':0,
   'placement_matrix_project_units':state['placement'],'units_to_m':state['unit_scale_to_m']},
  'scene_recipe':pkg.record(task.RECIPE),'scene_outputs':comparisons,
  'external_style_dependencies':report['necessary_external_style_dependencies'],'validation':'validation.json',
  'legacy_source_retained':True,'formal_write_allowed':False,'formal_sha256':task.FORMAL_HASH,
  'known_limitations':['Exact manufacturer article reference, not a project shop drawing.','Pure package targets the accepted main-bathroom occurrence, not the earlier representative instance.']})
 task.write(task.OUT/'cleanup-proposal.json',{'status':'proposal_only_no_deletion','requires_user_confirmation':True,
  'legacy_files':[v['source_session'] for v in report['views']],
  'temporary_directories':[report['temporary_project_directory']],
  'keep':['DUOFIX-product.ifc','scene-recipe.json','SVG/PNG','official DWG sources and approvals']})
 task.write(task.OUT/'handoff.json',{'product':'geberit-duofix-sigma-224-212','status':'validated_pending_visual_check',
  'package':'manifest.json','ifc':'DUOFIX-product.ifc','validation':'validation.json','scene_recipe':'scene-recipe.json',
  'new_files_only':True,'staged_baseline_modified':False,'legacy_deleted':False,'formal_modified':False,'bridge_port':9890})
 print(json.dumps(comparisons))

if __name__=='__main__':main()
