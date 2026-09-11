"""Independent native-curve, schema, scope and rendered-output verification."""
import json, subprocess, shutil
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np
import ifcopenshell, ifcopenshell.validate
import migrate_bed02 as task
pkg=task.pkg

def render(svg):
    png=svg.with_suffix('.png')
    subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),'--export-area-page',
        '--export-background=white','--export-background-opacity=1','--export-width=1400',f'--export-filename={png}'],
        check=True,capture_output=True)
    return png

def paths(rep):
    return [[list(p.Coordinates) for p in line.Points] for item in rep.Items for line in item.Elements]

def single_views(model):
    reps=model.by_guid(task.GUID).Representation.Representations
    outputs={}
    for view,axes in [('plan',(0,1)),('front',(0,2)),('side',(1,2))]:
        rep=next(r for r in reps if r.RepresentationIdentifier=='Approved'+view.title())
        curves=[[(p[axes[0]],-p[axes[1]]) for p in line] for line in paths(rep)]
        xs=[p[0] for l in curves for p in l];ys=[p[1] for l in curves for p in l]
        width=max(xs)-min(xs);height=max(ys)-min(ys);size=max(width,height)
        cx=(max(xs)+min(xs))/2;cy=(max(ys)+min(ys))/2
        w=max(width*1.15,height*1.15*1.4);h=w/1.4
        root=ET.Element('svg',{'xmlns':'http://www.w3.org/2000/svg','width':'1400','height':'1000',
            'viewBox':f'{cx-w/2} {cy-h/2} {w} {h}','data-source-kind':'native_dwg_review_reference',
            'data-global-id':task.GUID,'data-view':view,'data-geometry-source':'BED02-product.ifc',
            'data-path-count':str(len(curves)),'data-source-scale':'1.0','data-exact-project-configuration':'false'})
        ET.SubElement(root,'title').text=f'Baxter Viktor 160x200 family / BED02 / {view}; project width differs'
        for line in curves:
            ET.SubElement(root,'polyline',{'points':' '.join(f'{x:.9f},{y:.9f}' for x,y in line),
                'fill':'none','stroke':'#1677c8','stroke-width':str(size/800),'stroke-linejoin':'round'})
        svg=task.OUT/f'BED02-SINGLE-{view.upper()}.svg';ET.ElementTree(root).write(svg,encoding='utf-8',xml_declaration=True)
        outputs[view]={'svg':pkg.record(svg),'png':pkg.record(render(svg)),
            'source_representation':rep.RepresentationIdentifier,'path_count':len(curves),
            'source':'Exact approved native DWG curves read from pure IFC at scale 1.0'}
    iso=task.PRODUCT/'bonsai-camera-iso.png'
    if iso.is_file():
        target=task.OUT/'BED02-BODY-ISO.png';shutil.copy2(iso,target)
        outputs['iso']={'png':pkg.record(target),'source':pkg.record(iso),
            'method':'Existing actual IFC Body review image; original Body fingerprint preserved'}
    return outputs

def main():
    report=json.loads(task.REPORT.read_text())
    assert report['pure_saved_reloaded'] and report['scene_saved_reloaded']
    model=ifcopenshell.open(str(task.SINGLE));state=task.check_pure(model)
    assert state==report['pure_pre_bonsai']
    logger=ifcopenshell.validate.json_logger();ifcopenshell.validate.validate(model,logger)
    assert not logger.statements,logger.statements
    assert all(pkg.sha256(p['path'])==p['sha256'] for p in report['protected_files'])
    from pure_product_package import graph_from_json
    recipe=json.loads(task.RECIPE.read_text());graph=graph_from_json(recipe['graph'])
    saved_scene=ifcopenshell.open(report['temporary_project']['path'])
    report['scene_recipe']=pkg.record(task.RECIPE)
    report['bedroom_context_boundary']=recipe['bedroom_context_boundary']
    assert not recipe['body_geometry_included'] and not recipe['linework_geometry_included']
    assert all(not e.Representation for e in graph.by_type('IfcElement'))
    assert all(not graph.by_guid(v['annotation_guid']).Representation for v in report['views'])
    reference=json.loads(task.REFERENCE.read_text());residuals={};cameras={};scene_outputs=[]
    reps=model.by_guid(task.GUID).Representation.Representations
    for v in report['views']:
        view=v['view'];a=report['approved_alignment'][view];source=reference['views'][view]['paths_mm']
        lo=min(p[0] for l in source for p in l);hi=max(p[0] for l in source for p in l)
        expected=[]
        for line in source:
            curve=[]
            for x,y in line:
                if a['view_direction_reflection_x']: x=lo+hi-x
                x+=a['translation_mm'][0];y+=a['translation_mm'][1]
                curve.append([x,y,0.] if view=='plan' else ([x,0.,y] if view=='front' else [0.,x,y]))
            expected.append(curve)
        rep=next(r for r in reps if r.RepresentationIdentifier=='Approved'+view.title())
        actual=paths(rep);assert len(actual)==len(expected)==task.COUNTS[view]
        error=max(float(np.max(np.abs(np.array(x)-np.array(y)))) for x,y in zip(actual,expected))
        assert error<1e-8
        residuals[view]={'native_path_count':len(actual),'max_point_residual_mm':error,'uniform_scale':1.0}
        camera=graph.by_guid(v['drawing_guid']);matrix=pkg.placement(camera)
        camera_roundtrip_error=float(np.abs(matrix-pkg.placement(saved_scene.by_guid(v['drawing_guid']))).max())
        assert camera_roundtrip_error<1e-8
        import ifcopenshell.util.element as eu
        assert eu.get_pset(camera,'EPset_Drawing','Exclude')==eu.get_pset(saved_scene.by_guid(v['drawing_guid']),'EPset_Drawing','Exclude')
        cameras[view]={'placement_matrix_project_units':matrix.tolist(),'world_positive_z':matrix[:3,2].tolist(),
            'recipe_to_saved_camera_max_numeric_error':camera_roundtrip_error}
        if view=='front':
            # BED02 front is local -Y, transformed to world +X.
            expected_front=-pkg.placement(model.by_guid(task.GUID))[:3,1]
            assert np.dot(matrix[:3,2],expected_front)>0.999999
            assert abs(matrix[2,2])<.01
        svg=task.OUT/f'BED02-SCENE-{view.upper()}.svg';root=ET.parse(svg).getroot()
        geometry=[e for e in root.iter() if e.tag.rsplit('}',1)[-1] in ('line','path','polyline','polygon','circle','ellipse','rect')]
        blue=[e for e in geometry if '#1677c8' in e.get('style','')]
        grey=[e for e in geometry if '#a3abb3' in e.get('style','')]
        assert blue and grey and not [e for e in root.iter() if e.tag.endswith('image')]
        assert not [e for e in root.iter() if e.tag.rsplit('}',1)[-1]=='use']
        assert not [e for e in root.iter() if 'IfcSpace' in e.get('class','')]
        assert not [e for e in root.iter() if recipe['bedroom_context_boundary']['excluded_guid'] in e.get('class','')]
        annotation_guids={part[9:] for e in root.iter() for part in e.get('class','').split() if part.startswith('GlobalId-') and 'IfcAnnotation' in e.get('class','')}
        assert annotation_guids=={v['annotation_guid']}, annotation_guids
        assert all(e.tag.rsplit('}',1)[-1]=='line' for e in blue)
        assert len(blue)==sum(len(line)-1 for line in actual)
        paper=list(map(float,root.get('viewBox').split()))
        transform=np.linalg.inv(matrix)@pkg.placement(model.by_guid(task.GUID))
        projected=[]
        for line in actual:
            for p in line:
                point=transform@np.array([*p,1.])
                projected.append([point[0]/25+paper[2]/2,-point[1]/25+paper[3]/2])
        observed=np.array([[float(e.get(x)),float(e.get(y))] for e in blue for x,y in [('x1','y1'),('x2','y2')]])
        expected_points=np.array(projected)
        distances=np.sqrt(((observed[:,None,:]-expected_points[None,:,:])**2).sum(axis=2))
        registration_error=float(max(distances.min(axis=1).max(),distances.min(axis=0).max())*25)
        assert registration_error<.02, {'view':view,'registration_error_mm':registration_error}
        assert not [e for e in root.iter() if task.GUID in e.get('class','') and 'projection' in e.get('class','')]
        scene_outputs.append({'view':view,'svg':pkg.record(svg),'png':pkg.record(render(svg)),
            'blue_geometry_count':len(blue),'grey_context_geometry_count':len(grey),
            'pure_ifc_to_bonsai_svg_registration_error_mm':registration_error,'scene_approval_status':'pending'})
    previews=single_views(model)
    index_after=task.index_snapshot();assert index_after==report['index_before']
    report.update(independent_validation={'schema_errors':0,'pure_state':state,'native_geometry_comparison':residuals,
        'recipe_has_no_body_or_line_geometry':True,'protected_files_unchanged':True,'camera_views':cameras},
        library_previews=previews,scene_previews=scene_outputs,index_after=index_after,verdict='validated_pending_visual_check')
    task.write(task.REPORT,report)
    task.write(task.OUT/'manifest.json',{'schema_version':1,'profile_key':'bed02','display_name':'Baxter Viktor 160×200 family / BED02',
        'status':'validated_pending_visual_check','target_global_id':task.GUID,'persistent_ifc':pkg.record(task.SINGLE),
        'single_product_approval_status':'approved','scene_approval_status':'pending','approval_label':'单品已通过、场景待验收',
        'source_kind':'native_dwg_review_reference','source_label_zh':'Baxter Viktor 官方160×200家族原生DWG（非精确项目配置）',
        'official_cad_used':True,'source_scale':1.0,'exact_project_configuration':False,
        'approval_record':pkg.record(task.APPROVAL),'migration_authorization':report['derived_write_authority'],
        'source_records':[report['source_dwg'],pkg.record(task.REFERENCE),pkg.record(task.MANIFEST)],
        'geometry':{'physical_products':1,'body_representations':3,'approved_views':list(task.COUNTS),
            'annotations':0,'groups':0,'drawing_cameras':0,'spatial_geometry':0,
            'placement_matrix_project_units':state['placement'],'units_to_m':state['unit_scale_to_m'],
            'native_path_counts':task.COUNTS},'scene_recipe':pkg.record(task.RECIPE),'scene_outputs':scene_outputs,
        'library_previews':previews,'validation':'validation.json','formal_write_allowed':False,'formal_sha256':task.BASELINE})
    task.write(task.OUT/'cleanup-proposal.json',{'status':'proposal_only_no_deletion','requires_user_confirmation':True,
        'legacy_files':[pkg.record(p) for p in [task.OLD,task.PRODUCT/'Baxter-Viktor-BED02-bonsai-review.blend'] if p.is_file()],
        'temporary_directories':[str(Path(report['migration_temporary_ifc']['path']).parent),report['temporary_project_directory'],*report.get('previous_temporary_project_directories',[]),
            '/var/folders/rz/d8p6s4y50ws53rd150nm2s3r0000gn/T/bed02-migration-input-tjjly75v'],
        'keep':['BED02-product.ifc','scene-recipe.json','scene/single SVG and PNG','source/approval records','official-source/**']})
    print(json.dumps({'state':state,'native_geometry_comparison':residuals,'scene_outputs':scene_outputs}))

if __name__=='__main__':main()
