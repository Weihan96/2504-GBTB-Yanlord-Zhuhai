"""SXB010 migration; legacy session IFCs are extraction inputs only."""
from pathlib import Path
import json, sys, shutil, tempfile, traceback
import numpy as np
import ifcopenshell
OUT = Path(__file__).resolve().parent
PRODUCT = OUT.parent
ROOT = OUT.parents[4]
sys.path.insert(0, str(ROOT / 'pipeline/scripts'))
import review_product_package as pkg
FORMAL = ROOT / '2504 GBTB Yanlord Zhuhai.ifc'
BASELINE = '7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c'
GUID = '1O9JRXCI56VRUbpuLJy86Z'
APPROVAL = ROOT / 'pipeline/decisions/sxb010-drawing-approval.json'
OLD_MANIFEST = PRODUCT / 'bonsai-drawings/dining-bay/drawing-output-manifest.json'
SINGLE = OUT / 'SXB010-product.ifc'
RECIPE = OUT / 'scene-recipe.json'
REPORT = OUT / 'validation.json'

def write(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')

def snapshot(model):
    target = model.by_guid(GUID)
    return {'body':pkg.body_fingerprint(target), 'representations':pkg.fingerprint(target.Representation),
        'placement':pkg.placement(target).tolist(), 'physical_elements':sorted(e.GlobalId for e in model.by_type('IfcElement')),
        'annotation_count':len(model.by_type('IfcAnnotation')), 'group_count':len(model.by_type('IfcGroup')),
        'unit_scale_to_m':pkg.unit_util.calculate_unit_scale(model)}

def check_pure(model):
    result = snapshot(model)
    assert result['physical_elements'] == [GUID]
    assert result['annotation_count'] == result['group_count'] == 0
    assert not [e for e in model.by_type('IfcSpatialElement') if e.Representation]
    assert not [e for e in model.by_type('IfcPropertySet') if e.Name == 'EPset_Drawing']
    assert not [e for e in model.by_type('IfcPropertySingleValue') if e.Name in ('Include','Exclude')]
    reps = model.by_guid(GUID).Representation.Representations
    assert len([r for r in reps if r.RepresentationIdentifier == 'Body']) == 2
    for view, count in [('Plan',5),('Front',51),('Side',55)]:
        rep = next(r for r in reps if r.RepresentationIdentifier == 'Approved'+view)
        assert sum(len(item.Elements) for item in rep.Items) == count
    return result

def prepare(retry=False):
    from pure_product_package import build_pure_package
    assert retry or not SINGLE.exists()
    approval = json.loads(APPROVAL.read_text())
    assert approval['scene_svg_approved'] and not approval['formal_authoritative_ifc_write_allowed']
    assert pkg.sha256(FORMAL) == BASELINE
    assert pkg.sha256(OLD_MANIFEST) == approval['scene_svg_approval_evidence']['drawing_output_manifest_sha256']
    old = json.loads(OLD_MANIFEST.read_text())
    combined = None
    specs, views, protected = {}, [], [FORMAL, APPROVAL, OLD_MANIFEST]
    for v in old['views']:
        session = ROOT/v['session_ifc']
        assert pkg.sha256(session) == v['session_ifc_sha256']
        assert pkg.sha256(ROOT/v['svg']) == v['svg_sha256']
        source = ifcopenshell.open(str(session))
        selected = [v['drawing_global_id'], v['linework_annotation_global_id']]
        scoped = pkg.extract(source, GUID, selected)
        if combined is None:
            combined = scoped
        else:
            pkg.attach(combined, scoped, GUID, selected)
        specs[v['view']] = {'drawing_guid':selected[0], 'annotation_guid':selected[1]}
        views.append({'view':v['view'],'drawing_guid':selected[0], 'annotation_guid':selected[1],
            'old_svg':pkg.record(ROOT/v['svg']), 'source_session':pkg.record(session)})
        protected += [session, ROOT/v['svg']]
    # Independently saved sessions share the same application identity. IFC
    # UR1/UR2 requires one record, not one copy per session; no author is changed.
    from ifcopenshell.util.element import replace_attribute, remove_deep2
    applications = {}
    merged_applications = 0
    for app in list(combined.by_type('IfcApplication')):
        key = pkg.fingerprint(app)
        if key not in applications:
            applications[key] = app
            continue
        for inverse in combined.get_inverse(app):
            replace_attribute(inverse, app, applications[key])
        remove_deep2(combined,app)
        merged_applications += 1
    pure, recipe, audit = build_pure_package(combined, GUID, specs)
    audit['identical_cross_session_application_records_merged'] = merged_applications
    pure.write(str(SINGLE))
    write(RECIPE, recipe)
    protected += [p for p in (PRODUCT/'official-source').rglob('*') if p.is_file()]
    protected += [PRODUCT/'sxb010-derived-drawing.ifc', PRODUCT/'line-simplification-audit.json']
    report = {'product':'SXB010', 'target_guid':GUID,'verdict':'pending_bonsai_scene',
        'views':views, 'protected_files':[pkg.record(p) for p in protected], 'extraction':audit,
        'pure_pre_bonsai':check_pure(ifcopenshell.open(str(SINGLE))), 'single_product':pkg.record(SINGLE),
        'scene_recipe':pkg.record(RECIPE), 'formal_sha256_before':BASELINE, 'formal_write_allowed':False,
        'source_kind':'geometry_derived_simplified_proxy', 'source_label_zh':'基于原始高模几何生成的简化图纸表达',
        'colour_policy':'Preserve approved scene blue review highlighting; geometry-derived simplified proxy, NOT official CAD',
        'legacy_merge':'Three final approved per-view sessions scoped and combined in memory only',
        'courseEvidence':{'mode':'embedded-course-index','lesson':'085000','timestamps':['01:59','02:13','02:41']}}
    dependencies = []
    for location in sorted({s.Location for s in pure.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}):
        rel = Path(location)
        assert not rel.is_absolute() and '..' not in rel.parts
        source = ROOT/rel
        if source.is_file():
            dest = OUT/rel; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(source,dest)
            dependencies.append({'source':pkg.record(source),'packaged':pkg.record(dest)})
        else:
            dependencies.append({'source_path':str(source),'status':'already_missing_from_legacy_source'})
    report['necessary_external_style_dependencies'] = dependencies
    write(REPORT, report)
    print(json.dumps({'ifc':report['single_product'],'state':report['pure_pre_bonsai']}))

def scene():
    import bpy, bonsai_bridge as bridge
    from bonsai import tool
    from pure_product_package import attach_recipe
    import create_sxb010_dining_bay_drawing as draw
    report = json.loads(REPORT.read_text())
    assert Path(tool.Ifc.get_path()).resolve() == SINGLE.resolve()
    assert pkg.sha256(FORMAL) == BASELINE and not bpy.data.is_saved
    temporary = Path(tempfile.mkdtemp(prefix='sxb010-pure-package-scene-'))
    temp_ifc = temporary/'scene.ifc'
    report['temporary_project_directory'] = str(temporary)
    report['runtime_inputs'] = [str(SINGLE),str(RECIPE),str(FORMAL)]
    report['legacy_ifc_used_for_runtime'] = False
    report['verdict'] = 'bonsai_scene_running'
    write(REPORT,report)
    try:
        report['provider'] = {'version':list(bridge.bl_info['version']),'port':9888,
            'pid':__import__('os').getpid(),'blender':bpy.app.version_string,'ifcopenshell':ifcopenshell.version,
            'source':pkg.record(bridge.__file__)}
        with bpy.context.temp_override(**draw.view3d_override()):
            report['product_bonsai_save_reload'] = bridge._h_save_ifc_file({'output_path':str(SINGLE),'overwrite':True,'reload':True})
        pure = ifcopenshell.open(str(SINGLE))
        assert check_pure(pure) == report['pure_pre_bonsai']
        report['pure_saved_reloaded'] = True
        shutil.copy2(FORMAL,temp_ifc)
        model = ifcopenshell.open(str(temp_ifc))
        for location in sorted({s.Location for s in model.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}):
            rel = Path(location)
            assert not rel.is_absolute() and '..' not in rel.parts
            if (ROOT/rel).is_file():
                dest=temporary/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,dest)
        bodies = {e.GlobalId:pkg.body_fingerprint(e) for e in model.by_type('IfcElement') if e.Representation}
        coordinates = {e.GlobalId:pkg.placement(e).tolist() for e in model.by_type('IfcElement')}
        report['formal_physical_elements'] = len(model.by_type('IfcElement'))
        recipe = json.loads(RECIPE.read_text())
        report['attachment'] = attach_recipe(model,pure,recipe)
        count = len(list(model));attach_recipe(model,pure,recipe);assert len(list(model)) == count
        report['second_attachment_created_entities'] = 0
        for v in report['views']:
            drawing = model.by_guid(v['drawing_guid'])
            refs=[r.RelatingDocument for r in drawing.HasAssociations if r.is_a('IfcRelAssociatesDocument')]
            assert len(refs)==1
            refs[0].Location=str(OUT/f"SXB010-SCENE-{v['view'].upper()}.svg")
        model.write(str(temp_ifc))
        assert bpy.ops.bim.load_project(filepath=str(temp_ifc),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
        outputs=[]
        for v in report['views']:
            drawing=tool.Ifc.get().by_guid(v['drawing_guid'])
            tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
            with bpy.context.temp_override(**draw.view3d_override()):
                assert bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)=={'FINISHED'}
                props=tool.Drawing.get_document_props()
                props.should_use_underlay_cache=False;props.should_use_linework_cache=False;props.should_use_annotation_cache=False
                result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
            assert result=={'FINISHED'}
            svg=OUT/f"SXB010-SCENE-{v['view'].upper()}.svg"
            style=draw.style_review_svg(svg,v['annotation_guid'],v['view'])
            outputs.append({'view':v['view'],'svg':pkg.record(svg),'style':style,'inspection':draw.inspect_svg(svg),
                'operator':'bpy.ops.bim.create_drawing','result':sorted(result)})
            report['scene_outputs']=outputs;write(REPORT,report)
        with bpy.context.temp_override(**draw.view3d_override()):
            report['temporary_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(temp_ifc),'overwrite':True,'reload':True})
        reloaded=tool.Ifc.get()
        assert len(reloaded.by_type('IfcElement'))==report['formal_physical_elements']
        assert bodies=={e.GlobalId:pkg.body_fingerprint(e) for e in reloaded.by_type('IfcElement') if e.Representation}
        assert coordinates=={e.GlobalId:pkg.placement(e).tolist() for e in reloaded.by_type('IfcElement')}
        assert len([e for e in reloaded.by_type('IfcElement') if e.GlobalId==GUID])==1
        report.update(scene_saved_reloaded=True,all_formal_bodies_unchanged=True,all_formal_placements_unchanged=True,
            target_instances=1,temporary_project=pkg.record(temp_ifc),verdict='pending_independent_visual_validation')
        assert bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
    except Exception:
        report.update(error=traceback.format_exc(),verdict='fail')
        raise
    finally:
        report['formal_sha256_after']=pkg.sha256(FORMAL)
        report['single_product']=pkg.record(SINGLE)
        write(REPORT,report)
        assert report['formal_sha256_after']==BASELINE

if __name__=='__main__':
    prepare()
