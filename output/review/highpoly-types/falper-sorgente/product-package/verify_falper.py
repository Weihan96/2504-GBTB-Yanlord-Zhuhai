"""Independent schema, source, registration, SVG and resource verification."""
import json, subprocess, re, shutil
import xml.etree.ElementTree as ET
import numpy as np
import ifcopenshell, ifcopenshell.validate
from shapely.geometry import LineString
from shapely.ops import unary_union
import migrate_falper as task

def main():
    report=json.loads(task.REPORT.read_text())
    assert report['pure_saved_reloaded'] and report['scene_saved_reloaded']
    model=ifcopenshell.open(str(task.SINGLE))
    assert task.check_pure(model)==report['pure_pre_bonsai']
    logger=ifcopenshell.validate.json_logger()
    ifcopenshell.validate.validate(model,logger)
    assert not logger.statements,logger.statements
    assert all(task.pkg.sha256(r['path'])==r['sha256'] for r in report['protected_files'])
    assert task.pkg.sha256(task.FORMAL)==task.FORMAL_HASH
    archived=[]
    for source in [task.ROOT/json.loads(task.REGISTER.read_text())['variants']['WFB']['source_dwg'],task.APPROVAL,task.REGISTER,task.CORRECTED,task.PRODUCT/'INVALIDATED-proxy-contaminated.json']:
        dest=task.OUT/'source-evidence'/source.name
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,dest)
        archived.append({'source':task.pkg.record(source),'packaged':task.pkg.record(dest)})
    report['packaged_source_evidence']=archived
    scene=ifcopenshell.open(report['temporary_project']['path'])
    target=model.by_guid(task.GUID)
    outputs=[]
    for v in report['views']:
        svg=task.OUT/f"FALPER-WFB-SCENE-{v['view'].upper()}.svg"
        root=ET.parse(svg).getroot()
        assert not [e for e in root.iter() if e.tag.endswith('image')]
        assert not [e for e in root.iter() if e.tag.endswith('use') or 'IfcSpace' in e.get('class','')], 'Unrelated reference or room geometry remains'
        lines=[e for e in root.iter() if e.tag.endswith('line') and v['annotation_guid'] in e.get('class','')]
        assert lines
        assert not [e for e in root.iter() if e.get('{http://www.ifcopenshell.org/ns}guid')==task.GUID]
        camera=scene.by_guid(v['drawing_guid'])
        matrix=task.pkg.placement(camera)
        if v['view']!='plan': assert abs(matrix[2,2])<.01
        projector=np.linalg.inv(matrix)@task.pkg.placement(target)
        rep=next(r for r in target.Representation.Representations if r.RepresentationIdentifier==v['representation_identifier'])
        width=float(re.match(r'[0-9.]+',root.get('width')).group())
        height=float(re.match(r'[0-9.]+',root.get('height')).group())
        expected=[]
        for curve_set in rep.Items:
            for curve in curve_set.Elements:
                points=[]
                for p in curve.Points:
                    q=projector@np.array([*p.Coordinates,1.])
                    points.append((q[0]/25+width/2,-q[1]/25+height/2))
                expected.append(LineString(points))
        actual=unary_union([LineString([(float(e.get('x1')),float(e.get('y1'))),(float(e.get('x2')),float(e.get('y2')))]) for e in lines])
        error=unary_union(expected).hausdorff_distance(actual)*25
        assert error<.02,{'view':v['view'],'registration_error_mm':error,'camera':matrix.tolist(),'width':width,'height':height}
        png=svg.with_suffix('.png')
        subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),'--export-area-page','--export-background=white','--export-background-opacity=1','--export-width=1400',f'--export-filename={png}'],check=True,capture_output=True)
        outputs.append({'view':v['view'],'svg':task.pkg.record(svg),'preview':task.pkg.record(png),
                        'annotation_line_count':len(lines),'line_registration_error_world_mm':error,
                        'camera_matrix_project_units':matrix.tolist(),'camera_world_positive_z':matrix[:3,2].tolist()})
    report.update({'independent_validation':{'schema_errors':0,'protected_files_unchanged':True,'approved_representation_projection':outputs},
                   'verdict':'validated_pending_visual_check','provider_commit':'b0e67b1fb14ef5af93c4f747cb554b38511a844e',
                   'provider_bridge_sha256':'5902c69c146aff5f53d627b6151fdc53725433a044b96c199f7720a3dbccbd75'})
    report['scene_recipe']=task.pkg.record(task.RECIPE)
    report['courseEvidence'].update({'raw_source_available':False,'source':'research/bonsai-course/lessons/085000',
        'source_sha256':'85287cd1d85ea280d6c37f8ec4c6e29c32c950f6950b6f5c31ee981f6d93d04c',
        'screenshots':[{'timestamp':'01:59','path':'research/bonsai-course/lessons/085000/screenshots/085000-01m59s-create-drawing-button.png','sha256':'a86279d39cef0f85feb9c37d5f823556e41e1f6e1ebd1cc634e564c3ebfd24b4'},
                       {'timestamp':'02:13','path':'research/bonsai-course/lessons/085000/screenshots/085000-02m13s-svg-in-browser.png','sha256':'d559033d0259bd3cf3192b3ccbcb9b8bb8607d206c5e8c5b9e8a402d62771d62'}],
        'screen_observation_claimed':False})
    task.write(task.REPORT,report)
    task.write(task.OUT/'manifest.json',{'schema_version':1,'profile_key':'falper-sorgente','display_name':'Falper Sorgente WFB',
        'status':'validated_pending_visual_check','validation_verdict':'pending_visual_check','approval_status':'approved','scene_approval_status':'pending',
        'label_zh':'单品已通过、场景待验收','target_global_id':task.GUID,'persistent_ifc':task.pkg.record(task.SINGLE),
        'source_kind':'official_native_dwg_paths_mm','source_dwg':report['source_dwg'],'official_download_url':report['official_download_url'],
        'official_product_url':report['product_url'],'approval_record':report['approval_record'],'corrected_manifest':report['corrected_manifest'],
        'invalidated_proxy_source_used':False,'geometry':{'physical_products':1,'original_body_representation_count':2,'approved_views':['plan','front','side'],'annotations':0,'drawing_cameras':0,'spatial_geometry':0,'placement_matrix_project_units':report['pure_pre_bonsai']['placement'],'units_to_m':report['pure_pre_bonsai']['unit_scale_to_m']},
        'scene_recipe':task.pkg.record(task.RECIPE),'scene_outputs':outputs,'external_style_dependencies':report['necessary_external_style_dependencies'],
        'formal_write_allowed':False,'formal_sha256':task.FORMAL_HASH,'legacy_source_retained':True,'validation':'validation.json',
        'scope':'Approved official family reference; newly generated scene remains pending user acceptance.'})
    task.write(task.OUT/'cleanup-proposal.json',{'status':'proposal_only_no_deletion','requires_user_confirmation':True,
        'legacy_files':[task.pkg.record(task.OLD)],'missing_legacy_corrected_full_ifc':str(task.PRODUCT/'Falper-Sorgente-WFB-derived-drawing-corrected.ifc'),
        'temporary_directories':[report['temporary_project_directory'],*report.get('previous_temporary_attempts',[])],'keep':[str(task.SINGLE),str(task.RECIPE),'scene SVG/PNG','source CAD and approval evidence','validation files']})
    task.write(task.OUT/'handoff.json',{'product':'falper-sorgente','status':'validated_pending_visual_check','validation_verdict':'pending_visual_check',
        'scene_approval_status':'pending','label_zh':'单品已通过、场景待验收','package':'manifest.json','ifc':'FALPER-WFB-product.ifc',
        'scene_recipe':'scene-recipe.json','validation':'validation.json','new_files_only':True,'staged_baseline_modified':False,'legacy_deleted':False,'formal_modified':False,'bridge_port':9893})
    print(json.dumps({'schema_errors':0,'outputs':outputs},ensure_ascii=False))

if __name__=='__main__': main()
