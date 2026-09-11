"""BED02 pure product migration; runtime geometry comes only from the package."""
from pathlib import Path
import json, sys, shutil, tempfile, traceback, subprocess, hashlib
import numpy as np
import ifcopenshell
OUT = Path(__file__).resolve().parent
PRODUCT = OUT.parent
ROOT = OUT.parents[4]
sys.path.insert(0, str(ROOT / 'pipeline/scripts'))
import review_product_package as pkg
FORMAL = ROOT / '2504 GBTB Yanlord Zhuhai.ifc'
BASELINE = '7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c'
GUID = '1i_pqgLv9A7uuV7MjaArBW'
OLD = PRODUCT / 'Baxter-Viktor-BED02-bonsai-isolated.ifc'
APPROVAL = ROOT / 'pipeline/decisions/bed02-drawing-approval.json'
REFERENCE = PRODUCT / 'official-dwg-review-reference.json'
MANIFEST = PRODUCT / 'official-dwg-review-manifest.json'
SINGLE = OUT / 'BED02-product.ifc'
RECIPE = OUT / 'scene-recipe.json'
REPORT = OUT / 'validation.json'
COUNTS = {'plan':51, 'front':99, 'side':84}

def write(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')

def index_snapshot():
    entries = subprocess.check_output(['git','ls-files','--stage','-z'], cwd=ROOT)
    return {'entries_sha256':hashlib.sha256(entries).hexdigest(),
            'staged_files':len(subprocess.check_output(['git','diff','--cached','--name-only'], cwd=ROOT).splitlines()),
            'tracked_unstaged_files':subprocess.check_output(['git','diff','--name-only'], cwd=ROOT).decode().splitlines()}

def snapshot(model):
    target = model.by_guid(GUID)
    return {'body':pkg.body_fingerprint(target), 'representations':pkg.fingerprint(target.Representation),
        'placement':pkg.placement(target).tolist(), 'physical_elements':sorted(e.GlobalId for e in model.by_type('IfcElement')),
        'annotation_count':len(model.by_type('IfcAnnotation')), 'group_count':len(model.by_type('IfcGroup')),
        'unit_scale_to_m':pkg.unit_util.calculate_unit_scale(model)}

def check_pure(model):
    result=snapshot(model)
    assert result['physical_elements']==[GUID]
    assert result['annotation_count']==result['group_count']==0
    assert not [e for e in model.by_type('IfcSpatialElement') if e.Representation]
    assert not [e for e in model.by_type('IfcPropertySet') if e.Name=='EPset_Drawing']
    assert not [e for e in model.by_type('IfcPropertySingleValue') if e.Name in ('Include','Exclude')]
    reps=model.by_guid(GUID).Representation.Representations
    assert len([r for r in reps if r.RepresentationIdentifier=='Body'])==3
    assert len(reps)==6
    for view,count in COUNTS.items():
        rep=next(r for r in reps if r.RepresentationIdentifier=='Approved'+view.title())
        assert sum(len(item.Elements) for item in rep.Items)==count
    return result

def copy_resources(model, dest):
    records=[]
    for location in sorted({s.Location for s in model.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}):
        rel=Path(location)
        assert not rel.is_absolute() and '..' not in rel.parts
        src=ROOT/rel
        if not src.is_file():
            records.append({'source_path':str(src),'status':'already_missing_from_formal_source'})
            continue
        target=dest/rel
        target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists(): shutil.copy2(src,target)
        assert pkg.sha256(src)==pkg.sha256(target)
        records.append({'source':pkg.record(src),'packaged':pkg.record(target)})
    return records

def preflight():
    assert not REPORT.exists() and not SINGLE.exists()
    assert pkg.sha256(FORMAL)==BASELINE
    approval=json.loads(APPROVAL.read_text());manifest=json.loads(MANIFEST.read_text());ref=json.loads(REFERENCE.read_text())
    assert approval['status']=='approved' and approval['approved_views']==list(COUNTS)
    # The user approved the explicitly listed official SVGs. The historical
    # candidate hash predates manifest metadata updates; keep that record intact.
    assert approval['approval_evidence']['approved_artifacts']==[str((PRODUCT/f'official-dwg-{v}.svg').relative_to(ROOT)) for v in COUNTS]
    assert pkg.sha256(REFERENCE)==manifest['official_dwg_review_reference_sha256']
    assert pkg.sha256(ROOT/ref['source_dwg'])==ref['source_dwg_sha256']
    protected=[FORMAL,OLD,APPROVAL,REFERENCE,MANIFEST,PRODUCT/'manifest.json',ROOT/ref['source_dwg']]
    for view in manifest['views']:
        assert pkg.sha256(ROOT/view['svg'])==view['svg_sha256']
        assert view['alignment']['uniform_scale']==1.0
        assert not view['alignment']['source_geometry_deformed']
        assert len(ref['views'][view['view']]['paths_mm'])==COUNTS[view['view']]
        protected.append(ROOT/view['svg'])
    source=ifcopenshell.open(str(OLD))
    temporary=Path(tempfile.mkdtemp(prefix='bed02-migration-input-'))
    temp_ifc=temporary/'scene.ifc';shutil.copy2(FORMAL,temp_ifc)
    formal=ifcopenshell.open(str(temp_ifc))
    assert pkg.body_fingerprint(source.by_guid(GUID))==pkg.body_fingerprint(formal.by_guid(GUID))
    assert np.array_equal(pkg.placement(source.by_guid(GUID)),pkg.placement(formal.by_guid(GUID)))
    resources=copy_resources(formal,temporary)
    report={'product':'BED02','target_guid':GUID,'verdict':'pending_bonsai_prepare',
        'formal_sha256_before':BASELINE,'index_before':index_snapshot(),
        'protected_files':[pkg.record(p) for p in protected], 'source_pre_state':snapshot(source),
        'migration_temporary_ifc':pkg.record(temp_ifc), 'migration_resources':resources,
        'single_product_approval_status':'approved','scene_approval_status':'pending',
        'derived_write_authority':pkg.record(ROOT/'output/review/approved-product-library/migration-authorization.json'),
        'historical_candidate_manifest_hash':{'recorded':approval['candidate_manifest_sha256'],'current_official_manifest':pkg.sha256(MANIFEST),'match':False,'resolution':'Explicit approved artifact list checked against current official manifest SVG hashes; original approval retained unchanged.'},
        'formal_write_allowed':False, 'source_kind':'native_dwg_review_reference',
        'family_configuration':ref['configuration'],'exact_project_configuration':False,
        'source_dwg':pkg.record(ROOT/ref['source_dwg']),'official_download_url':ref['official_download_url'],
        'approved_alignment':{v['view']:v['alignment'] for v in manifest['views']},
        'source_dimensions':ref['official_dimension_cross_check_mm'],'project_body_bounds_mm':manifest['bounds_mm'],
        'courseEvidence':{'mode':'embedded-course-index','lesson':'085000','timestamps':['01:59 Create Drawing','02:13 SVG']}}
    write(REPORT,report)
    print(str(temp_ifc))

def view3d_override():
    import bpy
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type=='VIEW_3D':
                return {'window':window,'screen':window.screen,'area':area,
                        'region':next(r for r in area.regions if r.type=='WINDOW')}
    raise RuntimeError('No VIEW_3D available')

def prepare():
    import bpy, bonsai_bridge as bridge
    from bonsai import tool
    import ifcopenshell.api
    import create_bed01_master_bedroom_drawings as draw
    from pure_product_package import build_pure_package
    report=json.loads(REPORT.read_text())
    assert pkg.sha256(FORMAL)==BASELINE and not SINGLE.exists()
    assert Path(tool.Ifc.get_path()).resolve()==Path(report['migration_temporary_ifc']['path']).resolve()
    try:
        model=tool.Ifc.get();target=model.by_guid(GUID);obj=tool.Ifc.get_object(target)
        assert pkg.body_fingerprint(target)==report['source_pre_state']['body']
        bbox=draw.shared.world_bbox(obj); centre=[(bbox[0][i]+bbox[1][i])/2 for i in range(3)]
        rooms=[]
        for room in model.by_type('IfcSpace'):
            room_obj=tool.Ifc.get_object(room)
            if room_obj is None: continue
            rb=draw.shared.world_bbox(room_obj)
            if all(rb[0][i]<=centre[i]<=rb[1][i] for i in (0,1)):
                rooms.append((room,rb))
        assert len(rooms)==1, [(r.Name,r.LongName,b) for r,b in rooms]
        room,room_bbox=rooms[0]
        context=[x[0] for x in draw.shared.room_elements(room_bbox) if x[0]!=target]
        report['room']={'guid':room.GlobalId,'name':room.LongName,'bounds_m':[list(v) for v in room_bbox]}
        report['target_world_bbox_m']=[list(v) for v in bbox]
        report['context_elements']=[{'guid':e.GlobalId,'class':e.is_a(),'name':e.Name} for e in context]
        draw.EXPECTED_PATH_COUNTS=COUNTS;draw.TARGET_GLOBAL_ID=GUID
        draw.SOURCE_LABEL_ZH='Baxter Viktor 官方160×200家族原生DWG（非精确项目配置）'
        draw.SOURCE_DWG_SHA256=report['source_dwg']['sha256']
        draw.shared.ROOM_NAME=room.LongName;draw.shared.ARTICLE='Viktor 160x200 family'
        draw.shared.TARGET_GLOBAL_ID=GUID
        reference=json.loads(REFERENCE.read_text())
        views={};records=[]
        with bpy.context.temp_override(**view3d_override()):
            for view in COUNTS:
                definition={'drawing_name':f'BED02-SCENE-{view.upper()}',
                    'target_view':'PLAN_VIEW' if view=='plan' else 'ELEVATION_VIEW',
                    'location_hint':'PLAN' if view=='plan' else ('EAST' if view=='front' else 'SOUTH')}
                drawing,camera,*_=draw.shared.add_drawing(model,definition,room_bbox,context,OUT/f'BED02-SCENE-{view.upper()}.svg')
                assert bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)=={'FINISHED'}
                alignment=report['approved_alignment'][view]
                paths=draw.align_paths(view,reference['views'][view]['paths_mm'],alignment)
                annotation,rep,path_count,edge_count=draw.add_official_annotation(model,drawing,target,obj,view,paths,alignment)
                annotation.Name=f'Baxter Viktor official 160x200 family / {view}'
                annotation.Description='Approved family native DWG, scale 1.0; project Body configuration differs.'
                pset=ifcopenshell.util.element.get_pset(annotation,'EPset_Annotation')
                ifcopenshell.api.pset.edit_pset(model,pset=model.by_id(pset['id']),properties={'Classes':'review-target-bed02 official-native-dwg','ExactProjectConfiguration':False})
                for item in rep.Items:
                    for styled in item.StyledByItem:
                        for style in styled.Styles: style.Name='Baxter Viktor official native 2D DWG linework'
                camera_props=tool.Drawing.get_camera_props(camera);camera_props.has_annotation=True
                pset=ifcopenshell.util.element.get_pset(drawing,'EPset_Drawing')
                ifcopenshell.api.pset.edit_pset(model,pset=model.by_id(pset['id']),properties={'HasAnnotation':True})
                # Persist camera geometry and placement into the external recipe.
                tool.Geometry.record_object_position(camera)
                ifcopenshell.api.geometry.edit_object_placement(model,product=drawing,matrix=np.array(camera.matrix_world),is_si=True)
                views[view]={'annotation_guid':annotation.GlobalId,'drawing_guid':drawing.GlobalId}
                records.append({'view':view,**views[view],'path_count':path_count,'edge_count':edge_count})
        pure,recipe,audit=build_pure_package(model,GUID,views)
        for view,spec in recipe['views'].items():
            spec['representation_identifier']=recipe['pure_view_representations'][view]['representation_identifier']
        pure.write(str(SINGLE));write(RECIPE,recipe)
        report.update(extraction=audit,views=records,pure_pre_bonsai=check_pure(ifcopenshell.open(str(SINGLE))),
            necessary_external_style_dependencies=copy_resources(pure,OUT),single_product=pkg.record(SINGLE),
            scene_recipe=pkg.record(RECIPE),verdict='pending_bonsai_scene',
            provider={'version':list(bridge.bl_info['version']),'port':9891,'pid':__import__('os').getpid(),
                      'blender':bpy.app.version_string,'ifcopenshell':ifcopenshell.version,'source':pkg.record(bridge.__file__)})
        assert bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
    except Exception:
        report.update(error=traceback.format_exc(),verdict='fail');raise
    finally:
        report['formal_sha256_after']=pkg.sha256(FORMAL);write(REPORT,report)
        assert report['formal_sha256_after']==BASELINE

def style_svg(path,annotation_guid,view):
    import xml.etree.ElementTree as ET
    original=pkg.sha256(path);tree=ET.parse(path);root=tree.getroot()
    def local(s): return s.rsplit('}',1)[-1]
    geometry={'path','polyline','polygon','line','circle','ellipse','rect'}
    groups=[e for e in root.iter() if any(local(k)=='guid' and v==annotation_guid for k,v in e.attrib.items())
            or f'GlobalId-{annotation_guid}' in e.get('class','').split()]
    assert groups, 'Missing official annotation in Bonsai SVG'
    blue={id(e) for g in groups for e in g.iter() if local(e.tag) in geometry}
    before=[(local(e.tag),{k:v for k,v in e.attrib.items() if k not in ('style','class')}) for e in root.iter() if local(e.tag) in geometry]
    grey=0
    for e in root.iter():
        if local(e.tag) not in geometry: continue
        if id(e) in blue: style='stroke:#1677c8;stroke-width:0.35;fill:none'
        else: style='stroke:#a3abb3;stroke-width:0.22;fill:none;stroke-opacity:0.72';grey+=1
        e.set('style',e.get('style','').rstrip(';')+';'+style)
    after=[(local(e.tag),{k:v for k,v in e.attrib.items() if k not in ('style','class')}) for e in root.iter() if local(e.tag) in geometry]
    assert before==after and blue and grey
    for g in groups:
        g.set('data-source-kind','native_dwg_review_reference');g.set('data-source-scale','1.0');g.set('data-exact-project-configuration','false')
    tree.write(path,encoding='utf-8',xml_declaration=True)
    return {'raw_bonsai_sha256':original,'blue_geometry_count':len(blue),'grey_context_geometry_count':grey,'geometry_attributes_unchanged':True}

def resume_scene_after_style_adapter_fix():
    """Resume inspected existing runtime; its PLAN already exists unmodified."""
    import bpy, bonsai_bridge as bridge
    from bonsai import tool
    report=json.loads(REPORT.read_text());temp_ifc=Path(report['temporary_project_directory'])/'scene.ifc'
    assert Path(tool.Ifc.get_path()).resolve()==temp_ifc.resolve()
    assert report['pure_saved_reloaded'] and report['second_attachment_created_entities']==0
    report.setdefault('resolved_attempt_errors',[]).append(report.pop('error'))
    report['verdict']='bonsai_scene_resumed_after_svg_class_adapter_fix';write(REPORT,report)
    try:
        formal=ifcopenshell.open(str(FORMAL))
        bodies={e.GlobalId:pkg.body_fingerprint(e) for e in formal.by_type('IfcElement') if e.Representation}
        coords={e.GlobalId:pkg.placement(e).tolist() for e in formal.by_type('IfcElement')}
        outputs=[]
        for v in report['views']:
            svg=OUT/f"BED02-SCENE-{v['view'].upper()}.svg"
            if v['view']!='plan':
                drawing=tool.Ifc.get().by_guid(v['drawing_guid']);tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
                with bpy.context.temp_override(**view3d_override()):
                    assert bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)=={'FINISHED'}
                    props=tool.Drawing.get_document_props();props.should_use_underlay_cache=False;props.should_use_linework_cache=False;props.should_use_annotation_cache=False
                    assert bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)=={'FINISHED'}
            style=style_svg(svg,v['annotation_guid'],v['view'])
            outputs.append({'view':v['view'],'svg':pkg.record(svg),'style':style,'operator':'bpy.ops.bim.create_drawing','result':['FINISHED']})
            report['scene_outputs']=outputs;write(REPORT,report)
        with bpy.context.temp_override(**view3d_override()):
            report['temporary_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(temp_ifc),'overwrite':True,'reload':True})
        reloaded=tool.Ifc.get()
        assert len(reloaded.by_type('IfcElement'))==714
        assert bodies=={e.GlobalId:pkg.body_fingerprint(e) for e in reloaded.by_type('IfcElement') if e.Representation}
        assert coords=={e.GlobalId:pkg.placement(e).tolist() for e in reloaded.by_type('IfcElement')}
        report.update(scene_saved_reloaded=True,all_formal_bodies_unchanged=True,all_formal_placements_unchanged=True,
            target_instances=1,temporary_project=pkg.record(temp_ifc),verdict='pending_independent_visual_validation')
        assert bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
    except Exception:
        report.update(error=traceback.format_exc(),verdict='fail');raise
    finally:
        report['formal_sha256_after']=pkg.sha256(FORMAL);report['single_product']=pkg.record(SINGLE);write(REPORT,report)
        assert report['formal_sha256_after']==BASELINE

def clean_scene_context():
    """Exclude historical references in JSON, then reconstruct from formal."""
    import bpy
    from bonsai import tool
    import ifcopenshell.api
    from pure_product_package import graph_from_json, graph_to_json
    report=json.loads(REPORT.read_text())
    assert report['scene_saved_reloaded'] and Path(tool.Ifc.get_path()).resolve()==SINGLE.resolve()
    assert pkg.sha256(FORMAL)==BASELINE
    report['superseded_scene_outputs_with_global_reference_markers']=report['scene_outputs']
    report['previous_temporary_project_directories']=[report['temporary_project_directory']]
    recipe=json.loads(RECIPE.read_text());graph=graph_from_json(recipe['graph'])
    formal=ifcopenshell.open(str(FORMAL))
    historical_guids=sorted(a.GlobalId for a in formal.by_type('IfcAnnotation'))
    exclude=','.join(['IfcSpace','IfcGrid','IfcBuildingStorey',*historical_guids])
    for spec in recipe['views'].values():
        drawing=graph.by_guid(spec['drawing_guid']);pset=ifcopenshell.util.element.get_pset(drawing,'EPset_Drawing')
        ifcopenshell.api.pset.edit_pset(graph,pset=graph.by_id(pset['id']),properties={'GlobalReferencing':False,'Exclude':exclude})
    recipe['graph']=graph_to_json(graph)
    recipe['scene_filter_policy']={'excluded_reference_classes':['IfcSpace','IfcGrid','IfcBuildingStorey'],
        'excluded_historical_annotation_guids':historical_guids,'deleted_ifc_entities':False,
        'new_product_linework_annotations_excluded':False}
    write(RECIPE,recipe)
    report.update(scene_recipe=pkg.record(RECIPE),verdict='reconstructing_scene_with_external_reference_filters',
        scene_contains_only_requested_annotation_graphics=True)
    write(REPORT,report)
    scene()

def limit_bedroom_context():
    """Remove a confirmed guest-bath fixture from the bedroom view recipe."""
    import ifcopenshell.api, ifcopenshell.geom
    import ifcopenshell.util.element as eu
    from pure_product_package import graph_from_json, graph_to_json
    report=json.loads(REPORT.read_text());recipe=json.loads(RECIPE.read_text())
    graph=graph_from_json(recipe['graph']);formal=ifcopenshell.open(str(FORMAL))
    assert pkg.sha256(FORMAL)==BASELINE
    offender='245NU$zZL0d9tYTVwwBdk$';element=formal.by_guid(offender)
    settings=ifcopenshell.geom.settings();settings.set(settings.USE_WORLD_COORDS,True)
    shape=ifcopenshell.geom.create_shape(settings,element);points=np.array(shape.geometry.verts).reshape(-1,3)
    lo=points.min(0);hi=points.max(0);centre=(lo+hi)/2
    bedroom_max_x=report['room']['bounds_m'][1][0]
    assert centre[0]>bedroom_max_x and centre[0]<bedroom_max_x+.05
    for spec in recipe['views'].values():
        drawing=graph.by_guid(spec['drawing_guid']);pset=eu.get_pset(drawing,'EPset_Drawing')
        included=pset['Include'].split(',');assert offender in included;included.remove(offender)
        excluded=pset['Exclude'].split(',');excluded.append(offender)
        ifcopenshell.api.run('pset.edit_pset',graph,pset=graph.by_id(pset['id']),properties={'Include':','.join(included),'Exclude':','.join(excluded)})
    front=graph.by_guid(recipe['views']['front']['drawing_guid']);matrix=pkg.placement(front);before=matrix.tolist()
    matrix[0,3]=-1900.0
    ifcopenshell.api.run('geometry.edit_object_placement',graph,product=front,matrix=matrix,is_si=False)
    block=front.Representation.Representations[0].Items[0].TreeRootExpression
    assert block.is_a('IfcBlock');old_depth=block.ZLength
    depth=2900.0;block.ZLength=depth
    point=block.Position.Location;point.Coordinates=(*point.Coordinates[:2],-depth)
    recipe['graph']=graph_to_json(graph)
    audit={'excluded_guid':offender,'ifc_class':element.is_a(),'type_name':eu.get_type(element).Name,
        'container':eu.get_container(element).Name,'owner_space_guid':'3bF9ub5u5FFgmXQtSgfWnW','owner_space_name':'客卫',
        'bounds_m':[lo.tolist(),hi.tolist()],'centre_m':centre.tolist(),'bedroom_east_boundary_x_m':bedroom_max_x,
        'cause':'50 mm centre tolerance admitted a fixture from the adjoining guest bathroom; original camera was outside the bedroom.',
        'fix':'Exclude guest-bath fixture in all views; Front camera is now inside bedroom and stops at its west inner wall.',
        'front_camera_before_project_units':before,'front_camera_after_project_units':pkg.placement(front).tolist(),
        'front_depth_before_mm':old_depth,'front_depth_after_mm':depth,
        'front_x_visibility_interval_m':[-4.8,-1.9],'product_geometry_or_position_changed':False}
    recipe['bedroom_context_boundary']=audit
    write(RECIPE,recipe)
    report.setdefault('previous_temporary_project_directories',[]).append(report['temporary_project_directory'])
    report.update(bedroom_context_boundary=audit,scene_recipe=pkg.record(RECIPE),verdict='pending_room_boundary_regeneration')
    write(REPORT,report)

def scene():
    import bpy, bonsai_bridge as bridge
    from bonsai import tool
    from pure_product_package import attach_recipe
    report=json.loads(REPORT.read_text())
    assert Path(tool.Ifc.get_path()).resolve()==SINGLE.resolve()
    assert pkg.sha256(FORMAL)==BASELINE and not bpy.data.is_saved
    temporary=Path(tempfile.mkdtemp(prefix='bed02-pure-package-scene-'));temp_ifc=temporary/'scene.ifc'
    report.update(temporary_project_directory=str(temporary),runtime_inputs=[str(SINGLE),str(RECIPE),str(FORMAL)],
                  legacy_ifc_used_for_runtime=False,verdict='bonsai_scene_running')
    write(REPORT,report)
    try:
        with bpy.context.temp_override(**view3d_override()):
            report['product_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(SINGLE),'overwrite':True,'reload':True})
        pure=ifcopenshell.open(str(SINGLE))
        assert check_pure(pure)==report['pure_pre_bonsai']
        report['pure_saved_reloaded']=True
        shutil.copy2(FORMAL,temp_ifc);model=ifcopenshell.open(str(temp_ifc));copy_resources(model,temporary)
        bodies={e.GlobalId:pkg.body_fingerprint(e) for e in model.by_type('IfcElement') if e.Representation}
        coordinates={e.GlobalId:pkg.placement(e).tolist() for e in model.by_type('IfcElement')}
        report['formal_physical_elements']=len(model.by_type('IfcElement'))
        recipe=json.loads(RECIPE.read_text());report['attachment']=attach_recipe(model,pure,recipe)
        count=len(list(model));attach_recipe(model,pure,recipe);assert len(list(model))==count
        report['second_attachment_created_entities']=0
        for v in report['views']:
            drawing=model.by_guid(v['drawing_guid']);refs=[r.RelatingDocument for r in drawing.HasAssociations if r.is_a('IfcRelAssociatesDocument')]
            assert len(refs)==1;refs[0].Location=str(OUT/f"BED02-SCENE-{v['view'].upper()}.svg")
        model.write(str(temp_ifc))
        assert bpy.ops.bim.load_project(filepath=str(temp_ifc),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
        outputs=[]
        for v in report['views']:
            drawing=tool.Ifc.get().by_guid(v['drawing_guid']);tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
            with bpy.context.temp_override(**view3d_override()):
                assert bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)=={'FINISHED'}
                props=tool.Drawing.get_document_props();props.should_use_underlay_cache=False;props.should_use_linework_cache=False;props.should_use_annotation_cache=False
                result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
            assert result=={'FINISHED'}
            svg=OUT/f"BED02-SCENE-{v['view'].upper()}.svg"
            style=style_svg(svg,v['annotation_guid'],v['view'])
            outputs.append({'view':v['view'],'svg':pkg.record(svg),'style':style,'operator':'bpy.ops.bim.create_drawing','result':sorted(result)})
            report['scene_outputs']=outputs;write(REPORT,report)
        with bpy.context.temp_override(**view3d_override()):
            report['temporary_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(temp_ifc),'overwrite':True,'reload':True})
        reloaded=tool.Ifc.get()
        assert len(reloaded.by_type('IfcElement'))==report['formal_physical_elements']==714
        assert bodies=={e.GlobalId:pkg.body_fingerprint(e) for e in reloaded.by_type('IfcElement') if e.Representation}
        assert coordinates=={e.GlobalId:pkg.placement(e).tolist() for e in reloaded.by_type('IfcElement')}
        assert len([e for e in reloaded.by_type('IfcElement') if e.GlobalId==GUID])==1
        report.update(scene_saved_reloaded=True,all_formal_bodies_unchanged=True,all_formal_placements_unchanged=True,
            target_instances=1,temporary_project=pkg.record(temp_ifc),verdict='pending_independent_visual_validation')
        assert bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
    except Exception:
        report.update(error=traceback.format_exc(),verdict='fail');raise
    finally:
        report['formal_sha256_after']=pkg.sha256(FORMAL);report['single_product']=pkg.record(SINGLE);write(REPORT,report)
        assert report['formal_sha256_after']==BASELINE

if __name__=='__main__': preflight()
