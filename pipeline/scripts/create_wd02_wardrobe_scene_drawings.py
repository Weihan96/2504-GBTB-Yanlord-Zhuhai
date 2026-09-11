#!/usr/bin/env python3
"""Persist approved WD02 semantic views and generate native scene drawings."""
import contextlib
import json
import sys
import traceback
from pathlib import Path
from datetime import datetime, timezone

import bpy
import addon_utils
import numpy as np
import ifcopenshell
import ifcopenshell.api.geometry
import ifcopenshell.util.placement
from mathutils import Matrix, Vector
from bonsai import tool
import bonsai_bridge as bridge

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'pipeline/scripts'))
import create_wd03_wardrobe_scene_drawings as base

PRODUCT = ROOT / 'output/review/highpoly-types/wd02'
OUT = PRODUCT / 'bonsai-drawings/wardrobe'
IFC = PRODUCT / 'Poliform-Senzafine-WD02-derived-drawing.ifc'
BLEND = PRODUCT / 'Poliform-Senzafine-WD02-project-drawings.blend'
APPROVAL = ROOT / 'pipeline/decisions/wd02-drawing-approval.json'
EVIDENCE = OUT / 'WD02-WARDROBE-create-drawing-evidence.json'
GUID = '3yuXF4PtnBHgIXDlmXFJ$7'
COUNTS = {'plan':15, 'front':41, 'side':7}
QUOTE = '单品SVG可以验收，只差给我一个场景SVG'

def record(path):
    return {'path': str(path), 'sha256':base.sha256(path), 'bytes':path.stat().st_size}

def write(path, data):
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,default=str)+'\n')

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assert base.sha256(base.FORMAL_IFC) == base.FORMAL_SHA256
    assert bridge.bl_info['version'] == (1,1,0)
    candidate = json.loads((PRODUCT/'candidate-representations.json').read_text())
    assert candidate['representative_global_id']==GUID and not candidate['official_cad_used']
    assert {v:len(x['proxy_paths_mm']) for v,x in candidate['views'].items()}==COUNTS
    frozen = {str(PRODUCT/f):base.sha256(PRODUCT/f) for f in ['candidate-representations.json','wd02-semantic-segmentation.json','plan.svg','front.svg','side.svg','wd02-semantic-diagnostic-plan.svg','wd02-semantic-diagnostic-front.svg','wd02-semantic-diagnostic-side.svg']}
    approval = json.loads(APPROVAL.read_text())
    approval.update(status='approved',reviewer='user',review_date='2026-09-06',approved_views=list(COUNTS),derived_ifc_write_allowed=True,formal_authoritative_ifc_write_allowed=False,approval_evidence=QUOTE,approved_artifact_sha256=frozen,scene_review_status='pending',note='Approved semantic single-product views; product-derived IFC and native scene drawings requested. Formal authority remains closed.')
    write(APPROVAL, approval)
    with contextlib.suppress(Exception):
        addon_utils.disable('bl_ext.user_default.project_control',default_set=False,handle_error=None)
    assert bpy.ops.bim.load_project(filepath=str(base.FORMAL_IFC),should_start_fresh_session=True)=={'FINISHED'}
    with bpy.context.temp_override(**base.view3d_override()):
        before = bridge._h_get_scene_info({})
        copy_save = bridge._h_save_ifc_file({'output_path':str(IFC),'overwrite':True,'reload':True})
    model = tool.Ifc.get()
    target = model.by_guid(GUID)
    target_obj = tool.Ifc.get_object(target)
    placement = np.array(ifcopenshell.util.placement.get_local_placement(target.ObjectPlacement))
    target_matrix = target_obj.matrix_world.copy()
    original_body = [r.id() for r in target.Representation.Representations if r.RepresentationIdentifier=='Body']
    base.TARGET_GLOBAL_ID=GUID
    base.EXPECTED_PATH_COUNTS=COUNTS
    base.VIEW_DEFINITIONS = {v:{'drawing_name':f'WD02-WARDROBE-{v.upper()}','target_view':'PLAN_VIEW' if v=='plan' else 'ELEVATION_VIEW','location_hint':'EAST' if v=='side' else 'SOUTH'} for v in COUNTS}
    target_bbox = base.shared.world_bbox(target_obj)
    centre = [(target_bbox[0][a]+target_bbox[1][a])/2 for a in range(3)]
    scene_bbox=((centre[0]-1.6,centre[1]-1.6,-.2),(centre[0]+1.6,centre[1]+1.6,3))
    context = [e for e,*_ in base.shared.room_elements(scene_bbox) if e.GlobalId!=GUID]
    override=base.view3d_override()
    views=[]
    # Local camera basis columns are right, up and outward (camera looks -Z).
    orientations = {'plan':Matrix.Identity(4), 'front':Matrix(((1,0,0,0),(0,0,-1,0),(0,1,0,0),(0,0,0,1))), 'side':Matrix(((0,0,-1,0),(-1,0,0,0),(0,1,0,0),(0,0,0,1)))}
    camera_centres = {'plan':(.274,.28,3.1),'front':(.274,-1.15,1.2),'side':(-1.15,.28,1.2)}
    for view,definition in base.VIEW_DEFINITIONS.items():
        svg=OUT/f'{definition["drawing_name"]}.svg'
        drawing,camera,*_=base.shared.add_drawing(model,definition,scene_bbox,context,svg)
        local=orientations[view].copy();local.translation=camera_centres[view]
        camera.matrix_world=target_matrix@local
        camera.data.clip_end=4.5 if view=='plan' else 3.4
        props=tool.Drawing.get_camera_props(camera)
        props.update_props=False;props.width=3.0 if view=='plan' else 2.4;props.height=3.0 if view=='plan' else 3.2
        props.update_camera_resolution();props.update_props=True
        ifcopenshell.api.geometry.edit_object_placement(model,product=drawing,matrix=np.array(camera.matrix_world),is_si=True,should_transform_children=False)
        drawing.Description=f'WD02 {view}: actual project placement; approved semantic linework; Front from local -Y glass-door side.'
        with bpy.context.temp_override(**override):
            assert bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)=={'FINISHED'}
        annotation=base.add_semantic_annotation(model,drawing,target,target_obj,view,candidate['views'][view]['proxy_paths_mm'])
        annotation.Name=annotation.Name.replace('WD03','WD02')
        p=ifcopenshell.util.element.get_pset(annotation,'EPset_Annotation')
        ifcopenshell.api.pset.edit_pset(model,pset=model.by_id(p['id']),properties={'Classes':'review-target-wd02 geometry-derived semantic-linework approved'})
        # Persist reusable product representations in their correct local planes.
        rep=base.curve_representation(model,annotation.Representation.Representations[0].ContextOfItems,view,candidate['views'][view]['proxy_paths_mm'])
        rep.RepresentationIdentifier=f'Wd02{view.title()}'
        target.Representation.Representations=tuple(target.Representation.Representations)+(rep,)
        props.has_annotation=True
        p=ifcopenshell.util.element.get_pset(drawing,'EPset_Drawing')
        ifcopenshell.api.pset.edit_pset(model,pset=model.by_id(p['id']),properties={'HasAnnotation':True})
        dp=tool.Drawing.get_document_props();dp.should_use_underlay_cache=False;dp.should_use_linework_cache=False;dp.should_use_annotation_cache=False
        with bpy.context.temp_override(**override):
            result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
        assert result=={'FINISHED'} and svg.is_file()
        validation=base.style_and_inspect(svg,annotation.GlobalId)
        views.append({'view':view,'drawing_global_id':drawing.GlobalId,'annotation_global_id':annotation.GlobalId,'path_count':COUNTS[view],'camera_world_matrix_m':[list(row) for row in camera.matrix_world],'view_direction_world':list(-(camera.matrix_world.to_3x3()@Vector((0,0,1)))),'create_drawing_result':sorted(result),'svg':{**record(svg),**validation}})
    p=ifcopenshell.api.pset.add_pset(model,product=target,name='Pset_Wd02DrawingSource')
    ifcopenshell.api.pset.edit_pset(model,pset=p,properties={'SourceKind':candidate['source_kind'],'SourceLabelZh':candidate['source_label_zh'],'OfficialCadUsed':False,'ApprovalEvidence':QUOTE,'CandidateSha256':frozen[str(PRODUCT/'candidate-representations.json')],'ReviewStatus':'APPROVED_SINGLE_PRODUCT_SCENE_PENDING'})
    with bpy.context.temp_override(**base.view3d_override()):
        persistence=bridge._h_save_ifc_file({'output_path':str(IFC),'overwrite':True,'reload':True})
    model=tool.Ifc.get();target=model.by_guid(GUID)
    assert np.allclose(placement,ifcopenshell.util.placement.get_local_placement(target.ObjectPlacement),atol=1e-9)
    assert original_body==[r.id() for r in target.Representation.Representations if r.RepresentationIdentifier=='Body']
    for item in views:
        a=model.by_guid(item['annotation_global_id']);d=model.by_guid(item['drawing_global_id'])
        assert base.persisted_path_count(a)==item['path_count']
        annotation_matrix=ifcopenshell.util.placement.get_local_placement(a.ObjectPlacement)
        installation_error=float(np.max(np.abs(placement[:3,3]-annotation_matrix[:3,3])))
        assert installation_error<0.001 and np.allclose(placement[:3,:3],annotation_matrix[:3,:3],atol=1e-6)
        persisted=ifcopenshell.util.placement.get_local_placement(d.ObjectPlacement);persisted[:3,3]/=1000
        assert np.allclose(persisted,item['camera_world_matrix_m'],atol=1e-6)
        paths=a.Representation.Representations[0].Items[0].Elements
        error=max(abs(x-y) for poly,path in zip(paths,candidate['views'][item['view']]['proxy_paths_mm']) for point,pair in zip(poly.Points,path) for x,y in zip(point.Coordinates,base.coordinates_mm(item['view'],*pair)))
        assert error<0.001
        item['reload_verification']={'path_count':len(paths),'local_coordinate_error_mm':error,'installation_placement_error_mm':installation_error,'camera_persisted':True,'shared_product_origin_world_mm':placement[:3,3].tolist(),'persistence_tolerance_mm':0.001}
    bpy.ops.wm.save_as_mainfile(filepath=str(BLEND),check_existing=False)
    assert frozen=={p:base.sha256(Path(p)) for p in frozen}
    assert base.sha256(base.FORMAL_IFC)==base.FORMAL_SHA256
    write(EVIDENCE,{'task':'Approved WD02 semantic single-product linework to native scene SVG','courseEvidence':{'mode':'embedded-course-index','lesson':'085000','timestamps':['01:59 Create Drawing','02:13 SVG','02:41 Element Filters']},'plan':['copy formal to product-derived IFC','persist approved local semantic curves','orient cameras from product door semantics','native Create Drawing','public handler save and reload','verify coordinates and render SVG'],'preState':{'formal_sha256':base.FORMAL_SHA256,'approved_hashes':frozen,'target_placement_mm':placement.tolist(),'provider_scene':before},'execution':{'provider':'installed public bonsai_bridge handlers, isolated Blender process','provider_version':bridge.bl_info['version'],'blender':bpy.app.version_string,'ifcopenshell':ifcopenshell.version,'generator':'bpy.ops.bim.create_drawing','linework_mode':'OPENCASCADE','front_semantics':'local -Y glass door looking +Y, transformed to world +X'},'persistence':{'copy':copy_save,'final':persistence},'postState':{'target_placement_mm':placement.tolist(),'body_representation_ids_unchanged':original_body,'candidate_bytes_unchanged':True,'formal_sha256':base.FORMAL_SHA256},'outputs':{'views':views,'derived_ifc':record(IFC),'blend':record(BLEND)},'visual':{'status':'awaiting_render_inspection'},'review_status':'scene_pending_user_review','verdict':'awaiting_visual_verification'})
    print('WD02_NATIVE_SCENES_GENERATED',str(EVIDENCE))

if __name__=='__main__':
    try:main()
    except Exception:
        (PRODUCT/'wd02-scene-error.log').write_text(traceback.format_exc());raise
    finally:
        bpy.app.timers.register(lambda:bpy.ops.wm.quit_blender() and None,first_interval=2)
