"""WD03-only migration adapter. Legacy inputs are used during extraction only."""
from pathlib import Path
import sys, json
ROOT = Path(__file__).resolve().parents[5]
sys.path.insert(0, str(ROOT / 'pipeline/scripts'))
import ifcopenshell
import ifcopenshell.validate
import review_product_package as pkg

OUT = Path(__file__).resolve().parent
PRODUCT = OUT.parent
GUID = '3cmikd9MTB$egM5KQaNgUf'
FORMAL = ROOT / '2504 GBTB Yanlord Zhuhai.ifc'
BASELINE = '7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c'
OLD = PRODUCT / 'Poliform-Senzafine-WD03-derived-drawing.ifc'
SINGLE = OUT / 'WD03-product.ifc'
APPROVAL = ROOT / 'pipeline/decisions/wd03-drawing-approval.json'
EVIDENCE = PRODUCT / 'bonsai-drawings/wardrobe/WD03-WARDROBE-create-drawing-evidence.json'

def write(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')

def prepare():
    assert pkg.sha256(FORMAL) == BASELINE
    approval = json.loads(APPROVAL.read_text())
    assert approval['scene_svg_approved'] and approval['derived_ifc_write_allowed']
    assert not approval['formal_authoritative_ifc_write_allowed']
    assert pkg.sha256(OLD) == approval['write_execution']['derived_ifc_sha256']
    views = json.loads(EVIDENCE.read_text())['outputs']['views']
    for view in views:
        assert pkg.sha256(view['svg']['path']) == approval['scene_svg_sha256'][view['view']]
    source = ifcopenshell.open(str(OLD))
    pure = pkg.extract(source, GUID, [])
    repairs = pkg.repair_missing_metadata(pure)
    assert not pure.by_type('IfcAnnotation')
    assert not pure.by_type('IfcGroup')
    assert not [p for p in pure.by_type('IfcProduct') if p.Representation and p.GlobalId != GUID]
    assert not SINGLE.exists(), 'Existing output must be examined before replacement'
    pure.write(str(SINGLE))
    reloaded = ifcopenshell.open(str(SINGLE))
    validated = pkg.validate(source, reloaded, GUID, [])
    logger = ifcopenshell.validate.json_logger()
    ifcopenshell.validate.validate(reloaded, logger)
    assert not logger.statements, logger.statements
    protected = [FORMAL, OLD, APPROVAL, EVIDENCE, *[Path(v['svg']['path']) for v in views],
                 *[p for p in (PRODUCT/'official-source').rglob('*') if p.is_file()]]
    write(OUT/'handoff.json', {'profile_key':'wd03', 'status':'pure_extracted_scene_pending',
        'single_product':pkg.record(SINGLE), 'source_full_project':pkg.record(OLD),
        'validation':validated, 'schema_errors':0, 'inherited_metadata_repairs':repairs,
        'drawing_annotations':0, 'scene_camera_objects':0, 'represented_products':1,
        'protected_files':[pkg.record(p) for p in protected], 'approved_views':views,
        'formal_sha256_before':BASELINE, 'formal_sha256_after':pkg.sha256(FORMAL),
        'legacy_files_retained':True, 'formal_ifc_write_allowed':False,
        'courseEvidence':{'mode':'embedded-course-index','lesson':'085000','timestamps':['01:59','02:13','02:41']}})

def roundtrip():
    import bpy, bonsai_bridge
    from bonsai import tool
    import create_wd03_wardrobe_scene_drawings as draw
    report = json.loads((OUT/'handoff.json').read_text())
    assert Path(tool.Ifc.get_path()).resolve() == SINGLE
    assert len(tool.Ifc.get().by_type('IfcElement')) == 1
    assert not tool.Ifc.get().by_type('IfcAnnotation')
    with bpy.context.temp_override(**draw.view3d_override()):
        result = bonsai_bridge._h_save_ifc_file({'output_path':str(SINGLE),'overwrite':True,'reload':True})
    package = ifcopenshell.open(str(SINGLE))
    report['bonsai_save_reload'] = result
    report['validation'] = pkg.validate(ifcopenshell.open(str(OLD)),package,GUID,[])
    assert tool.Ifc.get_object(tool.Ifc.get().by_guid(GUID)) is not None
    report['provider'] = {'version':list(bonsai_bridge.bl_info['version']), 'path':bonsai_bridge.__file__,
        'sha256':pkg.sha256(bonsai_bridge.__file__), 'blender':bpy.app.version_string,
        'ifcopenshell':ifcopenshell.version, 'pid':__import__('os').getpid(),'port':9883}
    report['single_product'] = pkg.record(SINGLE)
    report['status'] = 'pure_bonsai_roundtrip_verified_scene_pending'
    assert pkg.sha256(FORMAL) == BASELINE
    write(OUT/'handoff.json', report)

def verify_pure():
    report = json.loads((OUT/'handoff.json').read_text())
    package = ifcopenshell.open(str(SINGLE))
    assert not package.by_type('IfcAnnotation')
    assert not package.by_type('IfcGroup')
    for pset in package.by_type('IfcPropertySet'):
        assert pset.Name != 'EPset_Drawing'
        assert not any(p.Name in ('Include','Exclude') for p in pset.HasProperties)
    represented = [p for p in package.by_type('IfcProduct') if p.Representation]
    assert len(represented) == 1 and represented[0].GlobalId == GUID
    expected = {'Wd03Plan':5, 'Wd03Front':9, 'Wd03Side':3}
    observed = {}
    for rep in package.by_guid(GUID).Representation.Representations:
        if rep.RepresentationIdentifier in expected:
            observed[rep.RepresentationIdentifier] = sum(len(i.Elements) for i in rep.Items)
    assert observed == expected, observed
    logger = ifcopenshell.validate.json_logger()
    ifcopenshell.validate.validate(package, logger)
    assert not logger.statements, logger.statements
    assert all(pkg.sha256(p['path']) == p['sha256'] for p in report['protected_files'])
    report['independent_pure_validation'] = {'one_represented_product':True,
        'drawing_or_camera_annotations':0,'scene_filter_properties':0,
        'three_view_path_counts':observed,'schema_errors':0,
        'formal_and_approved_legacy_sources_unchanged':True}
    report['single_product'] = pkg.record(SINGLE)
    write(OUT/'handoff.json', report)
    legacy = [OLD, PRODUCT/'Poliform-Senzafine-WD03-project-drawings.blend']
    write(OUT/'cleanup-proposal.json', {'status':'proposal_only_no_deletion',
        'files':[dict(pkg.record(p),requires_user_confirmation=True) for p in legacy if p.exists()],
        'new_package_scene_validation_required_before_cleanup':True,
        'keep':['WD03-product.ifc','scene-recipe.json','scene SVG/PNG','validation records',
                'official-source/**',str(APPROVAL)]})

def export_recipe():
    import pure_product_package as helper
    report = json.loads((OUT/'handoff.json').read_text())
    views = {v['view']:{'annotation_guid':v['annotation_global_id'],
        'drawing_guid':v['drawing_global_id'], 'representation_identifier':'Wd03'+v['view'].title()}
        for v in report['approved_views']}
    source = ifcopenshell.open(str(OLD))
    pure, recipe, audit = helper.build_pure_package(source,GUID,views)
    # Existing pure IFC has already passed Bonsai roundtrip; assert the new
    # shared extractor returns the same product graph before reusing it.
    pkg.validate(source,pure,GUID,[])
    pkg.validate(source,ifcopenshell.open(str(SINGLE)),GUID,[])
    write(OUT/'scene-recipe.json',recipe)
    report['pure_package_audit'] = audit
    report['scene_recipe'] = pkg.record(OUT/'scene-recipe.json')
    write(OUT/'handoff.json',report)

def scene(helper_module):
    """Runtime dependencies: pure IFC, external recipe, formal baseline only."""
    import importlib, shutil, tempfile, traceback
    import bpy, bonsai_bridge
    import numpy as np
    from bonsai import tool
    import create_wd03_wardrobe_scene_drawings as draw
    helper = importlib.import_module(helper_module)
    report = json.loads((OUT/'handoff.json').read_text())
    recipe = json.loads((OUT/'scene-recipe.json').read_text())
    assert recipe['body_geometry_included'] is False
    assert recipe['linework_geometry_included'] is False
    assert pkg.sha256(FORMAL) == BASELINE
    package = ifcopenshell.open(str(SINGLE))
    assert not package.by_type('IfcAnnotation')
    temporary = Path(tempfile.mkdtemp(prefix='wd03-pure-scene-'))
    temp_ifc = temporary/'scene.ifc'
    report.setdefault('temporary_attempts',[]).append(str(temp_ifc))
    report['status'] = 'scene_running'
    report['scene_adapter'] = {'file':pkg.record(__file__),'helper':pkg.record(helper.__file__)}
    report['scene_reload_verified'] = False
    report.pop('scene_outputs',None)
    report.pop('error',None)
    write(OUT/'handoff.json',report)
    try:
        shutil.copy2(FORMAL,temp_ifc)
        model = ifcopenshell.open(str(temp_ifc))
        formal_body = {e.GlobalId:pkg.body_fingerprint(e) for e in model.by_type('IfcElement') if e.Representation}
        report['attachment'] = helper.attach_recipe(model,package,recipe)
        count = len(list(model))
        helper.attach_recipe(model,package,recipe)
        assert len(list(model)) == count
        report['second_attachment_created_entities'] = 0
        assert formal_body == {e.GlobalId:pkg.body_fingerprint(e) for e in model.by_type('IfcElement') if e.Representation}
        for view in report['approved_views']:
            drawing = model.by_guid(view['drawing_global_id'])
            refs = [r.RelatingDocument for r in drawing.HasAssociations if r.is_a('IfcRelAssociatesDocument')]
            assert len(refs) == 1
            refs[0].Location = str(OUT/f"WD03-SCENE-{view['view'].upper()}.svg")
        resources = []
        for loc in sorted({s.Location for s in model.by_type('IfcExternallyDefinedSurfaceStyle') if s.Location}):
            relative = Path(loc)
            assert not relative.is_absolute() and '..' not in relative.parts
            resource = ROOT/relative
            if resource.is_file():
                dest = temporary/relative
                dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(resource,dest)
                resources.append({'source':pkg.record(resource),'temporary':str(dest)})
            else:
                resources.append({'location':loc,'status':'already_missing_at_formal_source'})
        report['temporary_resources'] = resources
        model.write(str(temp_ifc))
        assert bpy.ops.bim.load_project(filepath=str(temp_ifc),should_start_fresh_session=False,use_relative_path=False) == {'FINISHED'}
        outputs = []
        scene_anns = {}
        for view in report['approved_views']:
            drawing = tool.Ifc.get().by_guid(view['drawing_global_id'])
            annotation = tool.Ifc.get().by_guid(view['annotation_global_id'])
            scene_anns[annotation.GlobalId] = {'representation':pkg.fingerprint(annotation.Representation),
                'placement':pkg.placement(annotation).tolist()}
            tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
            with bpy.context.temp_override(**draw.view3d_override()):
                assert bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False) == {'FINISHED'}
                props = tool.Drawing.get_document_props()
                props.should_use_underlay_cache = False
                props.should_use_linework_cache = False
                props.should_use_annotation_cache = False
                result = bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
            assert result == {'FINISHED'}
            svg = OUT/f"WD03-SCENE-{view['view'].upper()}.svg"
            inspection = draw.style_and_inspect(svg,annotation.GlobalId)
            outputs.append({'view':view['view'],'svg':pkg.record(svg),'inspection':inspection,
                'generator':'bpy.ops.bim.create_drawing','annotation_global_id':annotation.GlobalId,
                'drawing_global_id':drawing.GlobalId})
        with bpy.context.temp_override(**draw.view3d_override()):
            report['temporary_save_reload'] = bonsai_bridge._h_save_ifc_file({'output_path':str(temp_ifc),'overwrite':True,'reload':True})
        reload = tool.Ifc.get()
        assert formal_body == {e.GlobalId:pkg.body_fingerprint(e) for e in reload.by_type('IfcElement') if e.Representation}
        assert len([e for e in reload.by_type('IfcElement') if e.GlobalId == GUID]) == 1
        assert np.array_equal(pkg.placement(reload.by_guid(GUID)),pkg.placement(package.by_guid(GUID)))
        for guid,expected in scene_anns.items():
            assert pkg.fingerprint(reload.by_guid(guid).Representation) == expected['representation']
            assert pkg.placement(reload.by_guid(guid)).tolist() == expected['placement']
        report['scene_outputs'] = outputs
        report['scene_reload_verified'] = True
        report['all_formal_bodies_unchanged'] = True
        report['runtime_input_dependencies'] = [str(SINGLE),str(OUT/'scene-recipe.json'),str(FORMAL)]
        report['status'] = 'scene_generated_pending_independent_visual_check'
        assert bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False) == {'FINISHED'}
    except Exception:
        report['status'] = 'scene_failed'
        report['error'] = traceback.format_exc()
        raise
    finally:
        report['formal_sha256_after_scene'] = pkg.sha256(FORMAL)
        write(OUT/'handoff.json',report)
        assert report['formal_sha256_after_scene'] == BASELINE

def verify_scene():
    import subprocess
    from verify_wd02_product_storage import line_geometry
    report = json.loads((OUT/'handoff.json').read_text())
    assert report['status'] == 'scene_generated_pending_independent_visual_check'
    assert report['scene_reload_verified']
    assert report['attachment']['linework_geometry_source'] == 'pure_product_ifc_representations'
    recipe = json.loads((OUT/'scene-recipe.json').read_text())
    assert not recipe['body_geometry_included'] and not recipe['linework_geometry_included']
    import pure_product_package as helper
    reconstructed_recipe = helper.graph_from_json(recipe['graph'])
    assert all(not reconstructed_recipe.by_guid(v['annotation_guid']).Representation for v in recipe['views'].values())
    compare = []
    for view in report['approved_views']:
        name = view['view']
        svg = OUT/f'WD03-SCENE-{name.upper()}.svg'
        old = line_geometry(Path(view['svg']['path']),view['annotation_global_id'])
        new = line_geometry(svg,view['annotation_global_id'])
        error = old.segmentize(.1).hausdorff_distance(new.segmentize(.1))*25
        assert error < .01, (name,error)
        png = svg.with_suffix('.png')
        subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),
            '--export-area-page','--export-background=white','--export-background-opacity=1',
            '--export-width=1400',f'--export-filename={png}'],check=True,capture_output=True)
        compare.append({'view':name,'max_linework_hausdorff_error_world_mm':error,
            'svg':pkg.record(svg),'png':pkg.record(png)})
    assert all(pkg.sha256(p['path']) == p['sha256'] for p in report['protected_files'])
    report['approved_scene_comparison'] = compare
    report['formal_sha256_after_independent_validation'] = pkg.sha256(FORMAL)
    assert report['formal_sha256_after_independent_validation'] == BASELINE
    report['status'] = 'independently_validated_pending_visual_check'
    write(OUT/'handoff.json',report)
    cleanup = json.loads((OUT/'cleanup-proposal.json').read_text())
    cleanup['temporary_scene_directories'] = [{'directory':str(Path(p).parent),
        'requires_user_confirmation':True} for p in report['temporary_attempts']]
    write(OUT/'cleanup-proposal.json',cleanup)
    (OUT/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>WD03 pure IFC validation</title>'
        '<style>body{font-family:system-ui;margin:32px}img{width:100%;max-width:1000px}</style>'
        '<h1>WD03 · 单品 IFC 迁移</h1><p>3D + 已批准 Plan / Front / Side；场景配方外置。正式 IFC 未改，旧副本未删除。</p>'
        '<p><a href="WD03-product.ifc">单品 IFC</a> · <a href="scene-recipe.json">场景配方</a> · '
        '<a href="handoff.json">验证记录</a> · <a href="cleanup-proposal.json">待清理清单</a></p>'
        + ''.join(f'<h2>{x["view"].title()} 场景 SVG</h2><a href="WD03-SCENE-{x["view"].upper()}.svg">'
          f'<img src="WD03-SCENE-{x["view"].upper()}.png"></a>' for x in compare))

if __name__ == '__main__':
    prepare()
