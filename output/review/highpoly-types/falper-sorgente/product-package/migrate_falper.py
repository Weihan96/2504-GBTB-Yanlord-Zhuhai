"""Falper official-source migration, with no formal or legacy writes."""
from pathlib import Path
import sys, json, shutil, tempfile, traceback
import numpy as np
import ifcopenshell
import ifcopenshell.util.element as eu

OUT = Path(__file__).resolve().parent
PRODUCT = OUT.parent
ROOT = OUT.parents[4]
sys.path.insert(0, str(ROOT / 'pipeline/scripts'))
import review_product_package as pkg
from pure_product_package import build_pure_package, attach_recipe
GUID = '350tdaubr8QP3Cu2YMQZIN'
FORMAL = ROOT / '2504 GBTB Yanlord Zhuhai.ifc'
FORMAL_HASH = '7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c'
OLD = PRODUCT / 'Falper-Sorgente-WFB-bonsai-isolated.ifc'
APPROVAL = ROOT / 'pipeline/decisions/falper-sorgente-drawing-approval.json'
REGISTER = ROOT / 'pipeline/decisions/falper-sorgente-official-dwg-linework.json'
CORRECTED = PRODUCT / 'derived-ifc-corrected-manifest.json'
SINGLE = OUT / 'FALPER-WFB-product.ifc'
RECIPE = OUT / 'scene-recipe.json'
REPORT = OUT / 'validation.json'
IDENTIFIERS = {'plan': 'FalperWFBPlan', 'front': 'FalperWFBFront', 'side': 'FalperWFBSide'}

def write(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')

def snapshot(model):
    target = model.by_guid(GUID)
    return {'body': pkg.body_fingerprint(target), 'representations': pkg.fingerprint(target.Representation),
            'placement': pkg.placement(target).tolist(), 'physical_elements': sorted(e.GlobalId for e in model.by_type('IfcElement')),
            'annotation_count': len(model.by_type('IfcAnnotation')), 'group_count': len(model.by_type('IfcGroup')),
            'drawing_pset_count': len([p for p in model.by_type('IfcPropertySet') if p.Name == 'EPset_Drawing']),
            'unit_scale_to_m': pkg.unit_util.calculate_unit_scale(model)}

def check_pure(model):
    state = snapshot(model)
    assert state['physical_elements'] == [GUID]
    assert state['annotation_count'] == state['group_count'] == state['drawing_pset_count'] == 0
    assert not [e for e in model.by_type('IfcSpatialElement') if e.Representation]
    assert not [p for p in model.by_type('IfcPropertySingleValue') if p.Name in ('Include', 'Exclude')]
    reps = model.by_guid(GUID).Representation.Representations
    assert sorted(r.RepresentationIdentifier for r in reps) == ['Body','Body','FalperWFBFront','FalperWFBPlan','FalperWFBSide']
    return state

def source_audit():
    import falper_sorgente_drawing_ifc as official
    corrected = json.loads(CORRECTED.read_text())
    assert pkg.sha256(FORMAL) == FORMAL_HASH
    for path, key in [(PRODUCT/'manifest.json','candidate_manifest_sha256'),(APPROVAL,'approval_record_sha256'),(REGISTER,'official_linework_register_sha256')]:
        assert pkg.sha256(path) == corrected[key]
    approval = json.loads(APPROVAL.read_text())
    assert approval['status'] == 'approved' and set(approval['approved_views']) == set(IDENTIFIERS)
    register = json.loads(REGISTER.read_text())
    dwg = ROOT / register['variants']['WFB']['source_dwg']
    assert pkg.sha256(dwg) == corrected['source_dwg_sha256']
    model = ifcopenshell.open(str(OLD))
    target = model.by_guid(GUID)
    expected = official.official_native_dwg_view_paths(json.loads((PRODUCT/'manifest.json').read_text()), register)
    audits = {}
    for view, identifier in IDENTIFIERS.items():
        rep = next(r for r in target.Representation.Representations if r.RepresentationIdentifier == identifier)
        actual = [[list(p.Coordinates) for p in line.Points] for curves in rep.Items for line in curves.Elements]
        paths = [[list((float(x),float(y),0.) if view=='plan' else ((float(x),0.,float(y)) if view=='front' else (0.,float(x),float(y)))) for x,y in path] for path in expected[view]]
        assert len(actual) == len(paths) == corrected['representation_path_counts'][view]
        assert [len(p) for p in actual] == [len(p) for p in paths]
        error = max(float(np.max(np.abs(np.array(a)-np.array(b)))) for a,b in zip(actual,paths))
        assert error < 1e-9
        audits[view] = {'paths':len(actual),'points':sum(map(len,actual)), 'max_coordinate_error_mm':error,
                        'tolerance_mm':1e-9, 'approved_representation':identifier, 'all_points_compared':True}
    formal = ifcopenshell.open(str(FORMAL))
    assert pkg.body_fingerprint(target) == pkg.body_fingerprint(formal.by_guid(GUID))
    protected = [pkg.record(p) for p in [FORMAL,OLD,APPROVAL,REGISTER,CORRECTED,PRODUCT/'manifest.json',PRODUCT/'INVALIDATED-proxy-contaminated.json',dwg]]
    return {'official_coordinate_comparison':audits,'protected_files':protected,'formal_sha256_before':FORMAL_HASH,
            'source_dwg':pkg.record(dwg),'approval_record':pkg.record(APPROVAL),'corrected_manifest':pkg.record(CORRECTED),
            'official_download_url':register['autocad_2d_zip_url'],'product_url':register['product_url'],
            'invalidated_source_used':False,'legacy_corrected_full_ifc_missing':not (ROOT/corrected['derived_ifc']).exists(),
            'approval_status':'approved','scene_approval_status':'pending','label_zh':'单品已通过、场景待验收',
            'formal_write_allowed':False,'cleanup_performed':False}

def override():
    import create_wd03_wardrobe_scene_drawings as context
    return context.view3d_override()

def prepare():
    import bpy, ifcopenshell.api, inspect, importlib
    from bonsai import tool
    from bonsai.core import drawing as core
    import create_gessi316_54294_main_bathroom_drawing as shared
    importlib.reload(shared)
    # Isolated source retains BATHG/FUR storeys, not the formal FFL storey.
    # Adapt only the in-memory function's initial placement hint; preserve all
    # spatial entity names and final explicit camera coordinates.
    code = inspect.getsource(shared.add_drawing).replace('if item.Name == "FFL"', 'if item.Name == "BATHG"')
    exec(code, shared.__dict__)
    report = source_audit()
    assert not SINGLE.exists()
    assert Path(tool.Ifc.get_path()).resolve() == OLD.resolve() and not bpy.data.is_saved
    # Previous failed preparation has no persisted output; reconstruct only
    # this task-owned in-memory candidate from its unchanged source.
    if tool.Ifc.get().by_type('IfcAnnotation'):
        bpy.ops.bim.load_project(filepath=str(OLD),should_start_fresh_session=False,use_relative_path=False)
    source = tool.Ifc.get()
    target = source.by_guid(GUID)
    formal = ifcopenshell.open(str(FORMAL))
    context_guids = [e.GlobalId for e in formal.by_type('IfcElement') if e.GlobalId != GUID and e.Representation and
                     (e.is_a('IfcWall') or e.is_a('IfcSlab') or np.linalg.norm(pkg.placement(e)[:2,3]-pkg.placement(target)[:2,3]) < 2000)]
    # The original approved product coordinates are retained. Cameras are new,
    # unapproved scene review framing; front faces local -Y and side local +X.
    bbox = ([-1.65,-.75,0.], [.55,1.45,2.8])
    specs = {}
    for view, identifier in IDENTIFIERS.items():
        definition = {'drawing_name':f'FALPER-WFB-SCENE-{view.upper()}',
                      'target_view':'PLAN_VIEW' if view=='plan' else 'ELEVATION_VIEW',
                      'location_hint':'PLAN' if view=='plan' else ('SOUTH' if view=='front' else 'EAST')}
        with bpy.context.temp_override(**override()):
            drawing, camera, *_ = shared.add_drawing(source,definition,bbox,[],OUT/f'FALPER-WFB-SCENE-{view.upper()}.svg')
            bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)
            context = tool.Drawing.get_annotation_context(definition['target_view']) or tool.Drawing.create_annotation_context(definition['target_view'])
            obj = core.add_annotation(tool.Ifc,tool.Collector,tool.Drawing,drawing=drawing,object_type='LINEWORK',relating_type=None,enable_editing=False)
            ann = tool.Ifc.get_entity(obj)
            ann.Name = f'Falper Sorgente WFB official native DWG / {view}'
            ann.ObjectPlacement = target.ObjectPlacement
            obj.matrix_world = tool.Ifc.get_object(target).matrix_world
            original = next(r for r in target.Representation.Representations if r.RepresentationIdentifier == identifier)
            rep = source.create_entity('IfcShapeRepresentation',ContextOfItems=context,RepresentationIdentifier='Annotation',RepresentationType=original.RepresentationType,Items=original.Items)
            ann.Representation = source.create_entity('IfcProductDefinitionShape',Representations=[rep])
            cprops = tool.Drawing.get_camera_props(camera)
            cprops.has_annotation = True
            data = eu.get_pset(drawing,'EPset_Drawing')
            ifcopenshell.api.pset.edit_pset(source,pset=source.by_id(data['id']),properties={'HasAnnotation':True,'Include':','.join([*context_guids,ann.GlobalId]),'Exclude':GUID})
            drawing.Description = 'Falper Sorgente WFB scene review candidate; product linework approved, scene pending user acceptance'
            # Write camera placement and orthographic extent into the recipe.
            import bonsai.core.geometry as core_geometry
            core_geometry.edit_object_placement(tool.Ifc,tool.Geometry,tool.Surveyor,obj=camera)
            cprops.update_representation(camera.matrix_world)
            bpy.ops.bim.update_representation(obj=camera.name,ifc_representation_class='')
        specs[view] = {'drawing_guid':drawing.GlobalId,'annotation_guid':ann.GlobalId,'representation_identifier':identifier}
    pure,recipe,audit = build_pure_package(source,GUID,specs)
    pure.write(str(SINGLE))
    write(RECIPE,recipe)
    resources = []
    for location in {s.Location for s in pure.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}:
        src,dest=ROOT/location,OUT/location
        assert src.is_file()
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dest)
        resources.append({'source':pkg.record(src),'packaged':pkg.record(dest)})
    report.update({'product':'falper-sorgente','target_guid':GUID,'verdict':'pending_bonsai_scene','extraction':audit,
                   'views':[{'view':k,**v} for k,v in specs.items()],'pure_pre_bonsai':check_pure(ifcopenshell.open(str(SINGLE))),
                   'necessary_external_style_dependencies':resources,'single_product':pkg.record(SINGLE),'scene_recipe':pkg.record(RECIPE),
                   'scene_context_guids':context_guids,'new_camera_scope':'Product surroundings within 2.2m square, original project coordinates',
                   'courseEvidence':{'mode':'embedded-course-index','lesson':'085000','timestamps':['01:59 Create Drawing','02:13 SVG'],
                                     'course_fact':'Drawing camera, scale, depth and filters govern regenerated SVG; validate persisted output.'}})
    write(REPORT,report)
    assert pkg.sha256(OLD) == next(r['sha256'] for r in report['protected_files'] if r['path']==str(OLD))
    bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False)

def style_svg(path,ann_guid,view):
    import xml.etree.ElementTree as ET
    tree=ET.parse(path); root=tree.getroot(); blue=grey=0
    geometry={'path','polyline','polygon','line','circle','ellipse','rect'}
    for e in root.iter():
        if e.tag.rsplit('}',1)[-1] not in geometry: continue
        target=ann_guid in e.get('class','')
        colour='#1677c8' if target else '#a3abb3'
        e.set('style',f'stroke:{colour};stroke-width:{0.35 if target else 0.22};fill:none')
        blue+=target; grey+=not target
    assert blue and grey,(view,blue,grey)
    root.set('data-source-kind','official_native_dwg_paths_mm');root.set('data-scene-approval','pending')
    tree.write(path,encoding='utf-8',xml_declaration=True)
    return {'blue_geometry':blue,'grey_context_geometry':grey,'style_only':True}

def scene():
    import bpy, bonsai_bridge as bridge
    from bonsai import tool
    report=json.loads(REPORT.read_text())
    if report.get('error'):
        report.setdefault('resolved_attempt_errors',[]).append(report.pop('error'))
        report.setdefault('previous_temporary_attempts',[]).append(report['temporary_project_directory'])
    assert Path(tool.Ifc.get_path()).resolve()==SINGLE.resolve() and not bpy.data.is_saved
    assert pkg.sha256(FORMAL)==FORMAL_HASH
    temporary=Path(tempfile.mkdtemp(prefix='falper-pure-package-scene-'))
    temp_ifc=temporary/'scene.ifc'
    report['temporary_project_directory']=str(temporary)
    report['runtime_inputs']=[str(SINGLE),str(RECIPE),str(FORMAL)]
    report['legacy_ifc_used_for_runtime']=False
    write(REPORT,report)
    try:
        report['provider']={'version':list(bridge.bl_info['version']),'port':9893,'blender':bpy.app.version_string,'ifcopenshell':ifcopenshell.version,'bridge_source':pkg.record(bridge.__file__)}
        with bpy.context.temp_override(**override()):
            report['product_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(SINGLE),'overwrite':True,'reload':True})
        pure=ifcopenshell.open(str(SINGLE))
        assert check_pure(pure)==report['pure_pre_bonsai']
        report['pure_saved_reloaded']=True
        shutil.copy2(FORMAL,temp_ifc)
        model=ifcopenshell.open(str(temp_ifc))
        for location in {s.Location for s in model.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}:
            src,dest=ROOT/location,temporary/location
            if src.is_file():
                dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
        baseline={e.GlobalId:{'body':pkg.body_fingerprint(e) if e.Representation else None,'placement':pkg.placement(e).tolist()} for e in model.by_type('IfcElement')}
        assert len(baseline)==714
        recipe=json.loads(RECIPE.read_text())
        report['attachment']=attach_recipe(model,pure,recipe)
        count=len(list(model));attach_recipe(model,pure,recipe);assert len(list(model))==count
        report['second_attachment_created_entities']=0
        for v in report['views']:
            drawing=model.by_guid(v['drawing_guid'])
            refs=[r.RelatingDocument for r in drawing.HasAssociations if r.is_a('IfcRelAssociatesDocument')]
            assert len(refs)==1
            refs[0].Location=str(OUT/f"FALPER-WFB-SCENE-{v['view'].upper()}.svg")
        model.write(str(temp_ifc))
        bpy.ops.bim.load_project(filepath=str(temp_ifc),should_start_fresh_session=False,use_relative_path=False)
        outputs=[]
        for v in report['views']:
            drawing=tool.Ifc.get().by_guid(v['drawing_guid'])
            tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
            with bpy.context.temp_override(**override()):
                assert bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)=={'FINISHED'}
                props=tool.Drawing.get_document_props()
                props.should_use_underlay_cache=props.should_use_linework_cache=props.should_use_annotation_cache=False
                result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
            assert result=={'FINISHED'}
            svg=OUT/f"FALPER-WFB-SCENE-{v['view'].upper()}.svg"
            raw=pkg.record(svg); styling=style_svg(svg,v['annotation_guid'],v['view'])
            outputs.append({'view':v['view'],'svg':pkg.record(svg),'raw_bonsai_svg':raw,'styles':styling,'operator':'bpy.ops.bim.create_drawing','result':sorted(result)})
            report['scene_outputs']=outputs;write(REPORT,report)
        with bpy.context.temp_override(**override()):
            report['temporary_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(temp_ifc),'overwrite':True,'reload':True})
        reloaded=tool.Ifc.get()
        after={e.GlobalId:{'body':pkg.body_fingerprint(e) if e.Representation else None,'placement':pkg.placement(e).tolist()} for e in reloaded.by_type('IfcElement')}
        assert baseline==after
        report.update({'scene_saved_reloaded':True,'all_formal_bodies_unchanged':True,'all_formal_placements_unchanged':True,'formal_physical_elements':714,'target_instances':1,'temporary_project':pkg.record(temp_ifc),'verdict':'pending_independent_visual_validation'})
        bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False)
    except Exception:
        report['error']=traceback.format_exc();report['verdict']='fail';raise
    finally:
        report['formal_sha256_after']=pkg.sha256(FORMAL);report['single_product']=pkg.record(SINGLE);write(REPORT,report)
        assert report['formal_sha256_after']==FORMAL_HASH

def run_background(operation):
    import bpy
    def run():
        try:
            globals()[operation]()
            write(OUT/f'{operation}-status.json',{'status':'complete'})
        except Exception:
            write(OUT/f'{operation}-status.json',{'status':'failed','error':traceback.format_exc()})
        return None
    write(OUT/f'{operation}-status.json',{'status':'running'})
    bpy.app.timers.register(run,first_interval=0.5)

def exclude_unrelated_references():
    import ifcopenshell.api.pset
    from pure_product_package import graph_from_json, graph_to_json
    report=json.loads(REPORT.read_text())
    formal=ifcopenshell.open(str(FORMAL))
    recipe=json.loads(RECIPE.read_text())
    graph=graph_from_json(recipe['graph'])
    excluded={e.GlobalId for cls in ('IfcAnnotation','IfcGrid','IfcBuildingStorey') for e in formal.by_type(cls)}
    excluded.update(v['drawing_guid'] for v in recipe['views'].values())
    excluded.add(GUID)
    excluded.add('IfcSpace')
    for v in recipe['views'].values():
        drawing=graph.by_guid(v['drawing_guid'])
        pset=eu.get_pset(drawing,'EPset_Drawing')
        ifcopenshell.api.pset.edit_pset(graph,pset=graph.by_id(pset['id']),properties={'Exclude':','.join(sorted(excluded)),'GlobalReferencing':False})
    recipe['graph']=graph_to_json(graph)
    recipe['reference_filter_policy']='Explicitly exclude historical Drawing, Annotation, Grid and Storey references; retain target LINEWORK and physical context.'
    write(RECIPE,recipe)
    report['scene_recipe']=pkg.record(RECIPE)
    report['reference_filter_policy']=recipe['reference_filter_policy']
    report.setdefault('previous_temporary_attempts',[]).append(report['temporary_project_directory'])
    write(REPORT,report)
