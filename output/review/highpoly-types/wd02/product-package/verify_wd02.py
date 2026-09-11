"""Independent persisted IFC and approved-scene geometry verification."""
from pathlib import Path
import json, subprocess
import ifcopenshell
import migrate_wd02 as task
from verify_wd02_product_storage import line_geometry
from pure_product_package import graph_from_json
pkg=task.pkg

def main():
    r=json.loads(task.REPORT.read_text())
    assert r['scene_reload_verified'] and r['pure_saved_reloaded']
    assert task.check(ifcopenshell.open(str(task.SINGLE)))==r['pure_pre_bonsai']
    recipe=json.loads(task.RECIPE.read_text())
    assert not recipe['body_geometry_included'] and not recipe['linework_geometry_included']
    graph=graph_from_json(recipe['graph'])
    assert not graph.by_guid(task.GUID).Representation
    assert all(not graph.by_guid(v['annotation_guid']).Representation for v in recipe['views'].values())
    compare=[]
    for v in r['approved_views']:
        svg=task.OUT/f"WD02-SCENE-{v['view'].upper()}.svg"
        old=line_geometry(Path(v['svg']['path']),v['annotation_global_id'])
        new=line_geometry(svg,v['annotation_global_id'])
        error=old.segmentize(.1).hausdorff_distance(new.segmentize(.1))*25
        assert error<.01,(v['view'],error)
        png=svg.with_suffix('.png')
        subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),
            '--export-area-page','--export-background=white','--export-background-opacity=1',
            '--export-width=1400',f'--export-filename={png}'],check=True,capture_output=True)
        compare.append(dict(view=v['view'],svg=pkg.record(svg),preview=pkg.record(png),
            max_linework_hausdorff_error_world_mm=error,tolerance_mm=.01))
    assert all(pkg.sha256(x['path'])==x['sha256'] for x in r['protected_files'])
    assert pkg.sha256(task.FORMAL)==task.BASELINE
    r.update(independent_validation=dict(schema_errors=0,recipe_contains_geometry=False,
        all_protected_files_unchanged=True,approved_scene_comparison=compare),status='validated_pending_visual_check')
    task.write(task.REPORT,r)
    previews={v:pkg.record(task.PRODUCT/f'{v}-preview.png') for v in ('plan','front','side')}
    previews['iso']=pkg.record(task.PRODUCT/'bonsai-camera-iso.png')
    task.write(task.OUT/'manifest.json',dict(schema_version=1,profile_key='wd02',
        display_name='Poliform Senzafine glass wardrobe / WD02',status=r['status'],
        persistent_ifc=pkg.record(task.SINGLE),target_global_id=task.GUID,
        source_kind='geometry_derived_simplified_proxy',source_label_zh='基于原始高模几何生成的简化图纸表达',
        approval_record=pkg.record(task.APPROVAL),source_access_record=pkg.record(task.PRODUCT/'official-source/source-access-record.json'),
        geometry=dict(physical_products=1,body_representations=2,annotations=0,drawing_cameras=0,
            spatial_geometry=0,approved_views=['plan','front','side'],placement_matrix_project_units=r['pure_pre_bonsai']['placement']),
        scene_recipe=pkg.record(task.RECIPE),scene_outputs=compare,style_dependencies=r['style_dependencies'],
        single_product_svgs={v:pkg.record(task.PRODUCT/f'{v}.svg') for v in ('plan','front','side')},
        library_previews=previews,preview_provenance='Existing approved semantic SVG previews and actual Bonsai Body iso camera render; source Body fingerprint unchanged',
        formal_write_allowed=False,formal_sha256=task.BASELINE,legacy_source_retained=True,
        validation='validation.json'))
    task.write(task.OUT/'cleanup-proposal.json',dict(status='proposal_only_no_deletion',requires_user_confirmation=True,
        legacy_files=[pkg.record(p) for p in [task.OLD,task.PRODUCT/'Poliform-Senzafine-WD02-project-drawings.blend',
            task.PRODUCT/'product-storage-pilot/WD02-product.ifc'] if p.is_file()],
        temporary_directories=[str(Path(r['temporary_project_ifc']).parent)],
        keep=['WD02-product.ifc','scene-recipe.json','SVG/PNG','Materials.blend','approval and official sources']))
    task.write(task.OUT/'handoff.json',dict(product='wd02',status=r['status'],package='manifest.json',
        ifc='WD02-product.ifc',validation='validation.json',scene_recipe='scene-recipe.json',
        new_files_only=True,legacy_deleted=False,formal_modified=False,staged_baseline_modified=False,bridge_port=9887))
    print(json.dumps(compare))

if __name__=='__main__':
    main()
