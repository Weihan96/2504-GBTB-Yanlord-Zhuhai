"""Independent IFC schema, approved SVG registration and preview validation."""
import json, subprocess, shutil
from pathlib import Path
from xml.etree import ElementTree as ET
import ifcopenshell, ifcopenshell.validate
from shapely.geometry import LineString
from shapely.ops import unary_union
import migrate_sxb010 as task
pkg=task.pkg

def line_geometry(path,guid):
    root=ET.parse(path).getroot()
    lines=[e for e in root.iter() if guid in e.get('class','') and e.tag.endswith('line')]
    assert lines
    assert not [e for e in root.iter() if e.tag.endswith('image')]
    return unary_union([LineString([(float(e.get('x1')),float(e.get('y1'))),(float(e.get('x2')),float(e.get('y2')))]) for e in lines])

def render(svg):
    png=svg.with_suffix('.png')
    subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),'--export-area-page',
        '--export-background=white','--export-background-opacity=1','--export-width=1400',f'--export-filename={png}'],
        check=True,capture_output=True)
    return png

def single_views(model):
    """Project exact persisted pure-IFC approved paths; do not rebuild geometry."""
    reps=model.by_guid(task.GUID).Representation.Representations
    outputs={}
    for view,axes in [('plan',(0,2)),('front',(0,1)),('side',(2,1))]:
        rep=next(r for r in reps if r.RepresentationIdentifier=='Approved'+view.title())
        paths=[[(float(p.Coordinates[axes[0]]),-float(p.Coordinates[axes[1]])) for p in line.Points]
            for item in rep.Items for line in item.Elements]
        xs=[p[0] for line in paths for p in line];ys=[p[1] for line in paths for p in line]
        width=max(xs)-min(xs);height=max(ys)-min(ys);size=max(width,height)
        cx=(max(xs)+min(xs))/2;cy=(max(ys)+min(ys))/2
        w=max(width*1.15,height*1.15*1.4);h=w/1.4
        root=ET.Element('svg',{'xmlns':'http://www.w3.org/2000/svg','width':'1400','height':'1000',
            'viewBox':f'{cx-w/2} {cy-h/2} {w} {h}','data-source-kind':'geometry_derived_simplified_proxy',
            'data-source-label-zh':'基于原始高模几何生成的简化图纸表达','data-global-id':task.GUID,
            'data-view':view,'data-geometry-source':'SXB010-product.ifc','data-path-count':str(len(paths))})
        ET.SubElement(root,'title').text=f'SXB010 / {view.title()} 单品 SVG / 批准的百叶中心线'
        for line in paths:
            ET.SubElement(root,'polyline',{'points':' '.join(f'{x:.9f},{y:.9f}' for x,y in line),
                'fill':'none','stroke':'#202832','stroke-width':str(size/800),'stroke-linejoin':'round'})
        svg=task.OUT/f'SXB010-SINGLE-{view.upper()}.svg'
        ET.ElementTree(root).write(svg,encoding='utf-8',xml_declaration=True)
        outputs[view]={'svg':pkg.record(svg),'png':pkg.record(render(svg)),'source_representation':rep.RepresentationIdentifier,
            'path_count':len(paths),'source':'Exact approved curves read from pure IFC, orthographic semantic axes'}
    source=task.PRODUCT/'bonsai-camera-iso.png'
    camera=json.loads((task.PRODUCT/'bonsai-review-manifest.json').read_text())
    expected=next(v for v in camera['renders'] if v['view']=='iso')
    assert pkg.sha256(source)==expected['sha256']
    iso=task.OUT/'SXB010-BODY-ISO.png';shutil.copy2(source,iso)
    outputs['iso']={'png':pkg.record(iso),'source':pkg.record(source),
        'source_manifest':pkg.record(task.PRODUCT/'bonsai-review-manifest.json'),
        'method':'Reused approved actual IFC Body camera render; current Body fingerprint unchanged'}
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
    assert not recipe['body_geometry_included'] and not recipe['linework_geometry_included']
    assert all(not e.Representation for e in graph.by_type('IfcElement'))
    assert all(not graph.by_guid(v['annotation_guid']).Representation for v in report['views'])
    comparisons=[]
    for v in report['views']:
        svg=task.OUT/f"SXB010-SCENE-{v['view'].upper()}.svg"
        old=line_geometry(Path(v['old_svg']['path']),v['annotation_guid']);new=line_geometry(svg,v['annotation_guid'])
        root=ET.parse(svg).getroot();scale=root.get('data-scale','1:25')
        # Bonsai data-scale is a numeric ratio such as 1/25; source uses 1:25.
        scale_world=25.0
        if scale and '/' in scale: scale_world=float(scale.split('/')[1])/float(scale.split('/')[0])
        elif scale and ':' in scale: scale_world=float(scale.split(':')[1])/float(scale.split(':')[0])
        error=old.segmentize(.1).hausdorff_distance(new.segmentize(.1))*scale_world
        assert error<.01,{'view':v['view'],'error_mm':error}
        comparisons.append({'view':v['view'],'max_linework_hausdorff_error_world_mm':error,'tolerance_mm':.01,
            'svg':pkg.record(svg),'preview':pkg.record(render(svg)),'target_bounds_approved':list(old.bounds),'target_bounds_new':list(new.bounds)})
    previews=single_views(model)
    assert pkg.sha256(task.FORMAL)==task.BASELINE
    report['independent_validation']={'schema_errors':0,'pure_state':state,'approved_scene_comparison':comparisons,
        'recipe_has_no_body_or_line_geometry':True,'protected_files_unchanged':True,'formal_sha256':task.BASELINE}
    report['library_previews']=previews
    report['verdict']='validated_pending_visual_check'
    task.write(task.REPORT,report)
    task.write(task.OUT/'manifest.json',{'schema_version':1,'profile_key':'sxb010','display_name':'Hunter Douglas 25 mm 百叶窗帘 / SXB010',
        'status':'validated_pending_visual_check','target_global_id':task.GUID,'persistent_ifc':pkg.record(task.SINGLE),
        'source_kind':'geometry_derived_simplified_proxy','source_label_zh':'基于原始高模几何生成的简化图纸表达',
        'official_cad_used':False,'colour_policy':'Single views black; scene retained approved blue review highlight, not official CAD',
        'approval_record':pkg.record(task.APPROVAL),'source_records':[r for r in report['protected_files'] if '/official-source/' in r['path']],
        'geometry':{'physical_products':1,'body_representations':2,'approved_views':['plan','front','side'],'annotations':0,'groups':0,
            'drawing_cameras':0,'spatial_geometry':0,'placement_matrix_project_units':state['placement'],'units_to_m':state['unit_scale_to_m'],
            'semantic_path_counts':{'plan':5,'front':51,'side':55},'slat_count':45},
        'scene_recipe':pkg.record(task.RECIPE),'scene_outputs':comparisons,'library_previews':previews,
        'validation':'validation.json','legacy_source_retained':True,'formal_write_allowed':False,'formal_sha256':task.BASELINE})
    legacy=[task.PRODUCT/'sxb010-derived-drawing.ifc',task.PRODUCT/'sxb010-bonsai-review.blend']
    legacy += [Path(v['source_session']['path']) for v in report['views']]
    task.write(task.OUT/'cleanup-proposal.json',{'status':'proposal_only_no_deletion','requires_user_confirmation':True,
        'legacy_files':[pkg.record(p) for p in legacy if p.is_file()],
        'temporary_directories':[report['temporary_project_directory']],
        'keep':['SXB010-product.ifc','scene-recipe.json','scene/single SVG and PNG','source/approval records','official-source/**']})
    task.write(task.OUT/'handoff.json',{'product':'sxb010','status':'validated_pending_visual_check','package':'manifest.json',
        'ifc':'SXB010-product.ifc','scene_recipe':'scene-recipe.json','validation':'validation.json',
        'new_files_only':True,'staged_baseline_modified':False,'legacy_deleted':False,'formal_modified':False,'bridge_port':9888})
    print(json.dumps({'state':state,'comparison':comparisons}))

if __name__=='__main__':main()
