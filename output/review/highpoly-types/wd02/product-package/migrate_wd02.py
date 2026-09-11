"""WD02 pure geometry package; runtime never reads a legacy project IFC."""
from pathlib import Path
import json, sys, shutil, tempfile, traceback
import ifcopenshell
import ifcopenshell.validate
import numpy as np
OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[4]
PRODUCT = OUT.parent
sys.path.insert(0, str(ROOT/'pipeline/scripts'))
import review_product_package as pkg
from pure_product_package import build_pure_package, attach_recipe
GUID = '3yuXF4PtnBHgIXDlmXFJ$7'
FORMAL = ROOT/'2504 GBTB Yanlord Zhuhai.ifc'
BASELINE = '7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c'
OLD = PRODUCT/'Poliform-Senzafine-WD02-derived-drawing.ifc'
SINGLE = OUT/'WD02-product.ifc'
RECIPE = OUT/'scene-recipe.json'
REPORT = OUT/'validation.json'
APPROVAL = ROOT/'pipeline/decisions/wd02-drawing-approval.json'
EVIDENCE = PRODUCT/'bonsai-drawings/wardrobe/WD02-WARDROBE-create-drawing-evidence.json'

def write(path, value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')

def state(model):
    p = model.by_guid(GUID)
    return dict(body=pkg.body_fingerprint(p),representations=pkg.fingerprint(p.Representation),
        placement=pkg.placement(p).tolist(),unit_scale_to_m=pkg.unit_util.calculate_unit_scale(model))

def check(model):
    assert [e.GlobalId for e in model.by_type('IfcElement')] == [GUID]
    assert not model.by_type('IfcAnnotation') and not model.by_type('IfcGroup')
    assert not [e for e in model.by_type('IfcSpatialElement') if e.Representation]
    assert not [p for p in model.by_type('IfcPropertySet') if p.Name=='EPset_Drawing' or any(x.Name in ('Include','Exclude') for x in p.HasProperties)]
    reps=model.by_guid(GUID).Representation.Representations
    assert len([r for r in reps if r.RepresentationIdentifier=='Body'])==2
    assert {r.RepresentationIdentifier for r in reps}=={'Body','Wd02Plan','Wd02Front','Wd02Side'}
    logger=ifcopenshell.validate.json_logger()
    ifcopenshell.validate.validate(model,logger)
    assert not logger.statements,logger.statements
    return state(model)

def prepare():
    assert pkg.sha256(FORMAL)==BASELINE
    approval=json.loads(APPROVAL.read_text())
    assert approval['scene_svg_approved'] and approval['derived_ifc_write_allowed']
    assert not SINGLE.exists()
    for path,sha in approval['approved_artifact_sha256'].items():
        assert pkg.sha256(path)==sha
    views=json.loads(EVIDENCE.read_text())['outputs']['views']
    source=ifcopenshell.open(str(OLD))
    specs={v['view']:dict(drawing_guid=v['drawing_global_id'],annotation_guid=v['annotation_global_id'],
        representation_identifier='Wd02'+v['view'].title()) for v in views}
    pure,recipe,audit=build_pure_package(source,GUID,specs)
    pure.write(str(SINGLE))
    write(RECIPE,recipe)
    assert check(ifcopenshell.open(str(SINGLE)))==state(source)
    protected=[FORMAL,OLD,APPROVAL,EVIDENCE,*[Path(v['svg']['path']) for v in views],
        *[p for p in (PRODUCT/'official-source').rglob('*') if p.is_file()],
        *[PRODUCT/f'{v}.svg' for v in ('plan','front','side')]]
    dependencies=[]
    for loc in sorted({s.Location for s in pure.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}):
        rel=Path(loc)
        assert not rel.is_absolute() and '..' not in rel.parts
        if (ROOT/rel).is_file():
            (OUT/rel).parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(ROOT/rel,OUT/rel)
            dependencies.append(dict(source=pkg.record(ROOT/rel),packaged=pkg.record(OUT/rel)))
        else:
            dependencies.append(dict(location=loc,status='already_missing_at_formal_source'))
    write(REPORT,dict(profile_key='wd02',status='extracted',pure_pre_bonsai=state(pure),extraction=audit,
        protected_files=[pkg.record(p) for p in protected],approved_views=views,source_full_project=pkg.record(OLD),
        single_product=pkg.record(SINGLE),scene_recipe=pkg.record(RECIPE),style_dependencies=dependencies,
        formal_sha256_before=BASELINE,formal_ifc_write_allowed=False,cleanup_performed=False,
        courseEvidence=dict(mode='embedded-course-index',lesson='085000',timestamps=['01:59 Create Drawing','02:13 SVG'])) )
    print(json.dumps(pkg.record(SINGLE)))

def scene():
    import bpy, bonsai_bridge
    from bonsai import tool
    import create_wd03_wardrobe_scene_drawings as draw
    report=json.loads(REPORT.read_text())
    assert Path(tool.Ifc.get_path()).resolve()==SINGLE
    assert not bpy.data.is_saved
    assert bonsai_bridge.bl_info['version']==(1,1,0)
    assert pkg.sha256(FORMAL)==BASELINE
    temp=Path(tempfile.mkdtemp(prefix='wd02-pure-scene-'))/'scene.ifc'
    report.update(status='scene_running',temporary_project_ifc=str(temp),runtime_inputs=[str(SINGLE),str(RECIPE),str(FORMAL)],legacy_ifc_used_for_runtime=False)
    write(REPORT,report)
    try:
        report['provider']=dict(version=list(bonsai_bridge.bl_info['version']),port=9887,pid=__import__('os').getpid(),
            blender=bpy.app.version_string,ifcopenshell=ifcopenshell.version,source=pkg.record(bonsai_bridge.__file__))
        with bpy.context.temp_override(**draw.view3d_override()):
            report['product_bonsai_save_reload']=bonsai_bridge._h_save_ifc_file(dict(output_path=str(SINGLE),overwrite=True,reload=True))
        pure=ifcopenshell.open(str(SINGLE))
        assert check(pure)==report['pure_pre_bonsai']
        report['pure_saved_reloaded']=True
        shutil.copy2(FORMAL,temp)
        model=ifcopenshell.open(str(temp))
        bodies={e.GlobalId:pkg.body_fingerprint(e) for e in model.by_type('IfcElement') if e.Representation}
        coords={e.GlobalId:pkg.placement(e).tolist() for e in model.by_type('IfcElement')}
        report['formal_physical_elements']=len(model.by_type('IfcElement'))
        recipe=json.loads(RECIPE.read_text())
        assert not recipe['body_geometry_included'] and not recipe['linework_geometry_included']
        report['attachment']=attach_recipe(model,pure,recipe)
        count=len(list(model));attach_recipe(model,pure,recipe)
        assert len(list(model))==count
        report['second_attachment_created_entities']=0
        for loc in sorted({s.Location for s in model.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}):
            rel=Path(loc)
            assert not rel.is_absolute() and '..' not in rel.parts
            if (ROOT/rel).is_file():
                (temp.parent/rel).parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(ROOT/rel,temp.parent/rel)
        for v in report['approved_views']:
            d=model.by_guid(v['drawing_global_id'])
            refs=[r.RelatingDocument for r in d.HasAssociations if r.is_a('IfcRelAssociatesDocument')]
            assert len(refs)==1
            refs[0].Location=str(OUT/f"WD02-SCENE-{v['view'].upper()}.svg")
        model.write(str(temp))
        assert bpy.ops.bim.load_project(filepath=str(temp),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
        draw.TARGET_GLOBAL_ID=GUID
        draw.EXPECTED_PATH_COUNTS={'plan':15,'front':41,'side':7}
        outputs=[]
        for v in report['approved_views']:
            d=tool.Ifc.get().by_guid(v['drawing_global_id'])
            tool.Ifc.get_object(d) or tool.Drawing.import_drawing(d)
            with bpy.context.temp_override(**draw.view3d_override()):
                assert bpy.ops.bim.activate_drawing(drawing=d.id(),should_view_from_camera=False)=={'FINISHED'}
                props=tool.Drawing.get_document_props()
                props.should_use_underlay_cache=props.should_use_linework_cache=props.should_use_annotation_cache=False
                result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
            assert result=={'FINISHED'}
            svg=OUT/f"WD02-SCENE-{v['view'].upper()}.svg"
            inspection=draw.style_and_inspect(svg,v['annotation_global_id'])
            outputs.append(dict(view=v['view'],svg=pkg.record(svg),inspection=inspection,generator='bpy.ops.bim.create_drawing'))
        with bpy.context.temp_override(**draw.view3d_override()):
            report['temporary_save_reload']=bonsai_bridge._h_save_ifc_file(dict(output_path=str(temp),overwrite=True,reload=True))
        reload=tool.Ifc.get()
        assert bodies=={e.GlobalId:pkg.body_fingerprint(e) for e in reload.by_type('IfcElement') if e.Representation}
        assert coords=={e.GlobalId:pkg.placement(e).tolist() for e in reload.by_type('IfcElement')}
        assert len(reload.by_type('IfcElement'))==report['formal_physical_elements']==714
        report.update(scene_reload_verified=True,all_formal_bodies_unchanged=True,all_formal_placements_unchanged=True,
            target_instances=1,scene_outputs=outputs,status='scene_generated_pending_independent_validation')
    except Exception:
        report.update(status='failed',error=traceback.format_exc())
        raise
    finally:
        report['formal_sha256_after']=pkg.sha256(FORMAL)
        report['single_product']=pkg.record(SINGLE)
        write(REPORT,report)
        assert report['formal_sha256_after']==BASELINE

def finalize_visual():
    """Called after independent checks and inspecting all three rendered PNGs."""
    import bpy
    from bonsai import tool
    import create_wd03_wardrobe_scene_drawings as draw
    assert Path(tool.Ifc.get_path()).resolve()==SINGLE
    assert len(tool.Ifc.get().by_type('IfcElement'))==1
    assert not tool.Ifc.get().by_type('IfcAnnotation')
    assert not [o for o in bpy.data.objects if o.type=='CAMERA']
    assert not bpy.data.is_saved and not tool.Blender.get_bim_props().has_blend_warning
    meshes=[o for o in bpy.data.objects if o.type=='MESH']
    assert all(tool.Ifc.get_entity(o) and (tool.Ifc.get_entity(o).GlobalId==GUID or tool.Ifc.get_entity(o).is_a('IfcTypeProduct')) for o in meshes)
    obj=tool.Ifc.get_object(tool.Ifc.get().by_guid(GUID))
    bpy.ops.object.select_all(action='DESELECT')
    obj.hide_set(False);obj.select_set(True);bpy.context.view_layer.objects.active=obj
    with bpy.context.temp_override(**draw.view3d_override()):
        bpy.ops.view3d.view_selected(use_all_regions=False)
    r=json.loads(REPORT.read_text())
    assert r['status']=='validated_pending_visual_check'
    r['final_live_state']=dict(path=tool.Ifc.get_path(),physical_products=1,annotations=0,cameras=0,
        mesh_objects=[o.name for o in meshes],no_scene_orphan_meshes=True,blend_saved=False,has_blend_warning=False,
        fresh_session_reload=True,pid=__import__('os').getpid(),port=9887)
    r['visual_review']=dict(inspected=['WD02-SCENE-PLAN.png','WD02-SCENE-FRONT.png','WD02-SCENE-SIDE.png'],
        observations=['Black semantic wardrobe linework present in all three views',
            'Plan closed cabinet perimeter and Front hanging rail/drawers preserved',
            'Grey scene context preserved; no duplicate product Body overlay',
            'Side and Front floor contact and wall registration match approved source'],verdict='pass')
    r['status']='pass';r['verdict']='pass'
    assert pkg.sha256(FORMAL)==BASELINE
    write(REPORT,r)
    for name in ('manifest.json','handoff.json'):
        p=OUT/name;data=json.loads(p.read_text());data['status']='pass';data['verdict']='pass';write(p,data)
    print(json.dumps(r['final_live_state']))

if __name__=='__main__':
    prepare()
