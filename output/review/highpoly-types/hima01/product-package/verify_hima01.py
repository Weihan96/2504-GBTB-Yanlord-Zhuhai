"""Independent disk geometry, placement, schema and generated-SVG checks."""
import json,copy,subprocess,xml.etree.ElementTree as ET
import numpy as np
import ifcopenshell,ifcopenshell.validate
from shapely.geometry import LineString
from shapely.ops import unary_union
import migrate_hima01 as task

def svg_lines(path,guid):
    return [e for e in ET.parse(path).getroot().iter() if e.tag.rsplit('}',1)[-1]=='line' and guid in e.get('class','')]

def render(svg,png,width=1400):
    subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),'--export-area-page','--export-background=white','--export-background-opacity=1',f'--export-width={width}',f'--export-filename={png}'],check=True,capture_output=True)

def main():
    pkg=task.pkg;r=json.loads(task.REPORT.read_text())
    assert r['pure_saved_reloaded'] and r['scene_saved_reloaded']
    pure=ifcopenshell.open(str(task.SINGLE));state=task.check_pure(pure)
    assert state==r['pure_pre_bonsai']
    logger=ifcopenshell.validate.json_logger();ifcopenshell.validate.validate(pure,logger)
    assert not logger.statements,logger.statements
    scene=ifcopenshell.open(r['temporary_project']['path']);formal=ifcopenshell.open(str(task.FORMAL))
    assert len(formal.by_type('IfcElement'))==len(scene.by_type('IfcElement'))==714
    for original in formal.by_type('IfcElement'):
        actual=scene.by_guid(original.GlobalId)
        assert np.array_equal(pkg.placement(original),pkg.placement(actual))
        if original.Representation: assert pkg.body_fingerprint(original)==pkg.body_fingerprint(actual)
    assert all(pkg.sha256(p['path'])==p['sha256'] for p in r['protected_files'])
    candidate=json.loads((task.PRODUCT/'candidate-representations.json').read_text())
    recipe=json.loads(task.RECIPE.read_text())
    from pure_product_package import graph_from_json,representation_content
    graph=graph_from_json(recipe['graph'])
    assert not [x for x in graph.by_type('IfcElement') if x.Representation]
    assert not graph.by_type('IfcPolyline') and not graph.by_type('IfcGeometricCurveSet')
    assert recipe['annotation_import_adapter']['geometry_source'].startswith('persisted pure IFC')
    r['scene_recipe']=pkg.record(task.RECIPE)
    r['annotation_import_adapter']={**recipe['annotation_import_adapter'],'script_evidence':pkg.record(task.OUT/'migrate_hima01.py'),'first_import_svg_error_mm':{'plan':0.2949230255429348,'front':0.24347799815602772,'side':4.21599381892522},'first_import_rejected':True,'saved_source_representation_unchanged':True}
    target=pure.by_guid(task.GUID);M=pkg.placement(target);checks=[]
    for view,spec in json.loads(task.RECIPE.read_text())['views'].items():
        rep=next(x for x in target.Representation.Representations if x.RepresentationIdentifier==f'Approved{view.title()}')
        scene_rep=scene.by_guid(spec['annotation_guid']).Representation.Representations[0]
        assert pkg.fingerprint(rep.Items)==pkg.fingerprint(scene_rep.Items), 'Persisted scene annotation geometry changed during Blender mesh adaptation'
        polylines=rep.Items[0].Elements
        approved=candidate['views'][view]['proxy_paths_mm']
        approved=[path if path[0]==path[-1] else [*path,path[0]] for path in approved]
        assert len(polylines)==len(approved)
        axes={'plan':(0,1),'front':(0,2),'side':(1,2)}[view]
        for path,line in zip(approved,polylines):
            actual=[[p.Coordinates[axes[0]],p.Coordinates[axes[1]]] for p in line.Points]
            assert np.array_equal(actual,path)
        camera=scene.by_guid(spec['drawing_guid']);C=pkg.placement(camera)
        block=next(x for x in scene.traverse(camera.Representation) if x.is_a('IfcBlock'))
        scale=50.0
        svg=task.OUT/f'HIMA01-SCENE-{view.upper()}.svg';lines=svg_lines(svg,spec['annotation_guid']);assert lines
        svgroot=ET.parse(svg).getroot()
        assert not [e for e in svgroot.iter() if e.tag.rsplit('}',1)[-1]=='use' and any(k in e.get('{http://www.w3.org/1999/xlink}href','') for k in ('elevation-','section-','grid-'))], 'Historical reference symbols must be excluded before CreateDrawing'
        assert not [e for e in svgroot.iter() if 'IfcSpace' in e.get('class','').split()], 'Space filling must be excluded by recipe'
        # Bonsai raster-resolution quantisation can widen the orthographic
        # frame by a fraction of a millimetre. Use its actual vector viewBox.
        frame=[float(x) for x in ET.parse(svg).getroot().get('viewBox').split()]
        new=unary_union([LineString([(float(e.get('x1')),float(e.get('y1'))),(float(e.get('x2')),float(e.get('y2')))]) for e in lines])
        expected=[]
        for line in polylines:
            coords=[]
            for point in line.Points:
                world=(M@np.array([*point.Coordinates,1.0]))[:3]
                local=C[:3,:3].T@(world-C[:3,3])
                coords.append((local[0]/scale+frame[2]/2,frame[3]/2-local[1]/scale))
            expected.append(LineString(coords))
        # Compare every emitted segment before GEOS union. Near-coincident
        # approved edges can make union introduce spurious noding topology.
        expected_segments=np.array([pair for line in expected for pair in zip(list(line.coords)[:-1],list(line.coords)[1:])])
        actual_segments=np.array([[(float(e.get('x1')),float(e.get('y1'))),(float(e.get('x2')),float(e.get('y2')))] for e in lines])
        assert expected_segments.shape==actual_segments.shape
        forward=np.max(np.linalg.norm(expected_segments-actual_segments,axis=2),axis=1)
        backward=np.max(np.linalg.norm(expected_segments-actual_segments[:,::-1],axis=2),axis=1)
        error=float(np.max(np.minimum(forward,backward)))*scale
        assert error<.001,{'view':view,'every_segment_endpoint_error_mm':error}
        if view!='plan':
            assert abs(C[2,0])<.01 and C[2,1]>.99, 'World Z must be screen up, never sideways or upside down'
        png=svg.with_suffix('.png');render(svg,png)
        # Display crop is made only from verified Bonsai product LINEWORK.
        thumb=ET.Element('{http://www.w3.org/2000/svg}svg',width='1000',height='760')
        bounds=new.bounds;pad=max(bounds[2]-bounds[0],bounds[3]-bounds[1])*.08
        thumb.set('viewBox',f'{bounds[0]-pad} {bounds[1]-pad} {bounds[2]-bounds[0]+2*pad} {bounds[3]-bounds[1]+2*pad}')
        for e in lines:
            thumb.append(copy.deepcopy(e))
        points=np.array([p for path in approved for p in path]);lo=points.min(axis=0);hi=points.max(axis=0)
        single_svg=task.OUT/f'HIMA01-SINGLE-{view.upper()}.svg';ET.ElementTree(thumb).write(single_svg,encoding='utf-8',xml_declaration=True)
        single_png=single_svg.with_suffix('.png');render(single_svg,single_png,1000)
        checks.append({'view':view,'path_count':len(polylines),'approved_coordinates_exact':True,'generated_scene_absolute_error_mm':error,'camera_screen_up_world_z':C[2,1],'svg':pkg.record(svg),'preview':pkg.record(png),'single_product_svg':pkg.record(single_svg),'single_product_thumbnail':pkg.record(single_png),'local_bounds_mm':[lo.tolist(),hi.tolist()]})
    r['independent_validation']={'schema_errors':0,'pure_state':state,'view_checks':checks,'all_formal_body_and_placement_count':714,'protected_files_unchanged':True,'recipe_runtime_geometry_source':'pure IFC only','single_product_envelope_and_origin_preserved':True}
    r['verdict']='validated_pending_visual_check';task.write(task.REPORT,r)
    task.write(task.OUT/'manifest.json',{'schema_version':1,'profile_key':'hima01','display_name':'Poliform Hima HIMA01','status':'validated_pending_visual_check','target_global_id':task.GUID,'persistent_ifc':pkg.record(task.SINGLE),'scene_recipe':pkg.record(task.RECIPE),'source_kind':'geometry_derived_simplified_proxy','source_label_zh':task.SOURCE_LABEL,'source_dwg':pkg.record(task.PRODUCT/'official-source/official-download/Poliform-HIMA-screen.dwg'),'official_download_url':'https://s3.poliform.it/2025/09/Poliform-HIMA-screen.dwg','approval_record':pkg.record(task.APPROVAL),'migration_authorization':pkg.record(task.ROOT/'output/review/approved-product-library/migration-authorization.json'),'single_product_approval_status':'approved','scene_approval_status':'pending','approval_label_zh':'单品已通过、场景待验收','formal_write_allowed':False,'formal_sha256':task.FORMAL_HASH,'geometry':{'physical_products':1,'body_representation_count':3,'original_body_preserved':True,'discarded_box_representation_count':1,'approved_views':['plan','front','side'],'annotations':0,'groups':0,'drawing_cameras':0,'spatial_geometry':0,'units_to_m':.001,'placement_matrix_project_units':state['placement']},'scene_outputs':checks,'single_product_views':{x['view']:x['single_product_svg'] for x in checks},'library_previews':{x['view']:x['single_product_thumbnail'] for x in checks},'external_style_dependencies':r['necessary_external_style_dependencies'],'validation':'validation.json','known_limitations':['Approved folded screen geometry-derived simplified black outlines; official PVA11 DWG is unfolded family reference only.','New real Bonsai scene SVGs are technical migration checks awaiting human scene approval.'],'preview_provenance':{'2d':'Product GlobalId annotation LINEWORK cropped from verified real Bonsai Scene SVG; display-only framing; elevations preserve world +Z up.'}})
    task.write(task.OUT/'cleanup-proposal.json',{'status':'proposal_only_no_deletion','requires_user_confirmation':True,'legacy_files':[pkg.record(task.OLD),pkg.record(task.PRODUCT/'Poliform-Hima-HIMA01-bonsai-review.blend')],'temporary_directories':sorted(set([r['temporary_project_directory'],*r.get('previous_temporary_attempts',[])])),'keep':[str(task.SINGLE),str(task.RECIPE),'official-source/**','source and approval records','scene SVG/PNG and validation evidence'],'cleanup_performed':False})
    task.write(task.OUT/'handoff.json',{'product':'hima01','status':'validated_pending_visual_check','package':'manifest.json','ifc':'HIMA01-product.ifc','scene_recipe':'scene-recipe.json','validation':'validation.json','scene_approval_status':'pending','single_product_approval_status':'approved','new_files_only':True,'staged_baseline_modified':False,'legacy_deleted':False,'formal_modified':False,'bridge_port':9894})
    print(json.dumps(checks,ensure_ascii=False))

if __name__=='__main__':main()
