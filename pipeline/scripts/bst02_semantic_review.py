#!/usr/bin/env python3
"""BST02: semantic visible linework from actual connected Body components."""
import base64
import json
import math
import struct
from pathlib import Path
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone

import ifcopenshell
import ifcopenshell.geom
from shapely.geometry import Polygon, LineString
from shapely.ops import unary_union, linemerge, polygonize_full

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from int1_highpoly_type_review import projected_raw_edges, display_edge_sample
from render_sis04_project_context import write_uncached_png_preview as shared_png_preview
from wd03_semantic_review import path_d, render_contact_svg

DIR = ROOT / 'output/review/highpoly-types/bst02'
IFC = DIR / 'Baxter-Beside-BST02-bonsai-isolated.ifc'
FORMAL = ROOT / '2504 GBTB Yanlord Zhuhai.ifc'
BASELINE = '7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c'
GUID = '3hA0vKpcn44u4Tsx4tqiUz'
SOURCE = 'geometry_derived_simplified_proxy'
LABEL = '基于原始高模几何生成的简化图纸表达'
AXES = {'plan': (0, 1), 'front': (0, 2), 'side': (1, 2)}


def write_uncached_png_preview(source, target):
    try:
        shared_png_preview(source, target)
    except subprocess.TimeoutExpired:
        # Shared renderer may time out only while reaping Chrome after its PNG
        # has already been produced. Accept only a complete newly rendered PNG.
        from PIL import Image
        with Image.open(target) as preview:
            preview.verify()


def body():
    model = ifcopenshell.open(IFC)
    shape = ifcopenshell.geom.create_shape(ifcopenshell.geom.settings(), model.by_guid(GUID))
    g = shape.geometry
    vertices = [tuple(c * 1000 for c in g.verts[i:i+3]) for i in range(0, len(g.verts), 3)]
    faces = [tuple(g.faces[i:i+3]) for i in range(0, len(g.faces), 3)]
    adj = defaultdict(set)
    for f in faces:
        for i in f:
            adj[i].update(f)
    unseen = set(adj)
    records, groups = [], {}
    while unseen:
        ids, stack = set(), [min(unseen)]
        while stack:
            i = stack.pop()
            if i in ids:
                continue
            ids.add(i)
            stack.extend(adj[i] - ids)
        unseen -= ids
        lo = [min(vertices[i][a] for i in ids) for a in range(3)]
        hi = [max(vertices[i][a] for i in ids) for a in range(3)]
        name = 'vertical_handle' if len(ids) == 12 else 'cylindrical_cabinet_body' if len(ids) == 256 else 'stepped_top_cap'
        assert len(ids) in (12, 256, 512)
        selected = [n for n, f in enumerate(faces) if f[0] in ids]
        item_ids = sorted({g.item_ids[n] for n in selected})
        groups[name] = [faces[n] for n in selected]
        records.append(dict(semantic=name, item_ids=item_ids, vertex_count=len(ids), face_count=len(selected),
                            minimum_mm=lo, maximum_mm=hi, size_mm=[hi[a]-lo[a] for a in range(3)],
                            segmentation='edge-connected triangle component inside actual IfcPolygonalFaceSet'))
    assert len(records) == 3 and len(vertices) == 780 and len(faces) == 1418
    return vertices, faces, groups, records


def polygon(vertices, faces, axes):
    triangles = [Polygon([(round(vertices[i][axes[0]],6),round(vertices[i][axes[1]],6)) for i in f]) for f in faces]
    p = unary_union([t for t in triangles if t.area > 1e-8])
    assert p.geom_type == 'Polygon'
    assert sum(Polygon(r).area for r in p.interiors) < 1e-6
    p = Polygon(p.exterior)
    s = p.simplify(0.05, preserve_topology=True)
    assert p.boundary.hausdorff_distance(s.boundary) <= 0.050001
    return s


def lines(g):
    if g.is_empty:
        return []
    if g.geom_type == 'LineString':
        return [list(g.coords)]
    return [p for child in g.geoms for p in lines(child)]


def views(vertices, groups):
    result, proof = {}, {}
    for view, axes in AXES.items():
        projected = {name: polygon(vertices, fs, axes) for name, fs in groups.items()}
        # Product-local camera positions: +Z (Plan), -Y (Front), -X (Side).
        if view == 'plan':
            order = ['stepped_top_cap', 'vertical_handle', 'cylindrical_cabinet_body']
        elif view == 'front':
            order = ['vertical_handle', 'stepped_top_cap', 'cylindrical_cabinet_body']
        else:
            order = ['stepped_top_cap', 'cylindrical_cabinet_body', 'vertical_handle']
        occluders, drawn, entries = [], [], []
        for name in order:
            p = projected[name]
            visible = p.boundary
            if occluders:
                # Remove edges behind opaque components; contact edges are emitted once.
                visible = visible.difference(unary_union(occluders))
            if drawn:
                visible = visible.difference(unary_union(drawn))
            if visible.geom_type == 'MultiLineString':
                visible = linemerge(visible)
            for n, path in enumerate(lines(visible)):
                if LineString(path).length < 0.01:
                    continue
                entries.append(dict(id=f'{view}.{name}.{n}', component=name,
                    kind='visible_component_boundary', path_mm=[[round(x,6),round(y,6)] for x,y in path],
                    closed=path[0]==path[-1], projection_boundary_error_mm=0.05))
            occluders.append(p)
            drawn.append(visible)
        if view != 'plan':
            # Projection union preserves the outside step, but the recessed collar
            # is also a real visible depth change across the middle of the cap.
            top_points = {i for f in groups['stepped_top_cap'] for i in f}
            collar = [vertices[i] for i in top_points
                      if abs(vertices[i][2]-230)<1e-5
                      and (vertices[i][0]**2+vertices[i][1]**2)**0.5 < 270]
            a = axes[0]
            extent = [round(min(p[a] for p in collar),6),round(max(p[a] for p in collar),6)]
            entries.append(dict(id=f'{view}.stepped_top_cap.depth_seam',component='stepped_top_cap',
                kind='visible_depth_discontinuity',path_mm=[[extent[0],230],[extent[1],230]],closed=False,
                evidence='Actual radius-266 collar meets radius-280 cap at z=230; cap front surface is 14 mm nearer.',
                projection_boundary_error_mm=0))
        result[view] = entries
        proof[view] = dict(front_to_back_components=order,
            fully_hidden_components=[n for n in order if not any(e['component']==n for e in entries)],
            rule='Subtract nearer opaque component projections from farther component boundaries; deduplicate shared contact edges.',
            total_visible_length_mm=sum(LineString(e['path_mm']).length for e in entries),
            individual_component_projection_areas_mm2={n:p.area for n,p in projected.items()})
    return result, proof


def render(view, entries, raw=None):
    points = [p for e in entries for p in e['path_mm']]
    lo = [min(p[a] for p in points) for a in (0,1)]
    hi = [max(p[a] for p in points) for a in (0,1)]
    scale = min(1060/(hi[0]-lo[0]),650/(hi[1]-lo[1]))
    tr = lambda p: (700+(p[0]-(lo[0]+hi[0])/2)*scale,555-(p[1]-(lo[1]+hi[1])/2)*scale)
    grey = '' if raw is None else f'<path class="original-highpoly" d="{path_d(raw,tr)}" stroke="#888" stroke-opacity="0.22" stroke-width="0.6" fill="none"/>'
    markup = ''.join(f'<path data-semantic-id="{e["id"]}" d="{path_d([e["path_mm"]],tr)}"/>' for e in entries)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="980" viewBox="0 0 1400 980"><rect width="1400" height="980" fill="white"/>
<g font-family="Arial" fill="#20252b"><text x="50" y="55" font-size="30">Baxter Beside / BST02 — {view.upper()} semantic single product</text><text x="50" y="95" font-size="19">{LABEL}</text><text x="50" y="132" font-size="18">{'Grey: actual Body / ' if raw else ''}Black: visible semantic boundaries · pending approval</text></g>{grey}
<g class="simplified-proxy-silhouette geometry-derived" data-source-kind="{SOURCE}" fill="none" stroke="#000" stroke-width="2" stroke-linejoin="round">{markup}</g>
<g font-family="Arial" font-size="17"><text x="50" y="925">Actual Body: 560 × 577.889 × 250 mm. Official family height: 350 mm; 100 mm difference preserved.</text><text x="50" y="955">Cylinder / stepped top cap / vertical handle · {GUID}</text></g></svg>'''


def validate():
    candidate=load_json(DIR/'candidate-representations.json')
    manifest=load_json(DIR/'manifest.json')
    records={}
    for view in AXES:
        linework=unary_union([LineString(p) for p in candidate['views'][view]['proxy_paths_mm']])
        regions,cuts,dangles,invalid=polygonize_full(linework)
        assert cuts.length < 1e-6 and dangles.length < 1e-6 and invalid.length < 1e-6
        records[view]=dict(closed_semantic_region_count=len(regions.geoms),
            dangling_line_length_mm=dangles.length,cut_line_length_mm=cuts.length,
            invalid_ring_length_mm=invalid.length,svg_sha256=sha256(DIR/f'{view}.svg'))
    assert sha256(DIR/'candidate-representations.json')==manifest['candidate_representations_sha256']
    assert sha256(FORMAL)==BASELINE
    write_json(DIR/'bst02-semantic-validation.json',dict(
        generated_at=datetime.now(timezone.utc).isoformat(),views=records,**{'pass': True},
        candidate_sha256=sha256(DIR/'candidate-representations.json'),
        manifest_sha256=sha256(DIR/'manifest.json'),formal_ifc_sha256=BASELINE,
        no_ifc_write=True,writer_followup='Existing bst02_drawing_ifc.py still expects old 6/1/1 paths; update only after new semantic candidate approval.'))


def main():
    assert sha256(FORMAL)==BASELINE
    vertices, faces, groups, records = body()
    project_semantics, _ = views(vertices, groups)
    handle = [vertices[i] for i in {i for f in groups['vertical_handle'] for i in f}]
    centre = [sum(p[a] for p in handle)/len(handle) for a in (0,1)]
    angle = -math.atan2(centre[0], -centre[1])
    c, s = math.cos(angle), math.sin(angle)
    vertices = [(c*x-s*y,s*x+c*y,z) for x,y,z in vertices]
    aligned = [c*centre[0]-s*centre[1],s*centre[0]+c*centre[1]]
    assert abs(aligned[0]) < 1e-8
    semantics, visibility = views(vertices,groups)
    audit_path = DIR/'bst02-semantic-segmentation.json'
    audit = dict(schema_version=1,generated_at=datetime.now(timezone.utc).isoformat(),representative_global_id=GUID,
        source_kind=SOURCE,source_label_zh=LABEL,isolated_ifc=relative(IFC),isolated_ifc_sha256=sha256(IFC),
        representation_item_count=2,connected_component_count=3,mesh_vertex_count=780,mesh_face_count=1418,
        components=records,views=semantics,visibility=visibility,
        coordinate_frame='Handle-aligned review XYZ in mm; actual IFC placement remains unchanged. Plan +Z, Front -Y, Side -X.',
        orientation_audit=dict(review_rotation_z_degrees=math.degrees(angle),handle_centre_original_mm=centre,
            handle_centre_aligned_mm=aligned,method='Mean of the 12 handle vertices relative to the cylindrical axis at local (0,0).',
            project_context_uses_original_local_projection=True,ifc_placement_modified=False),
        base_audit=dict(official_measurement_svg='official-source/BESICOBS55.svg',
            official_measurement_svg_sha256=sha256(DIR/'official-source/BESICOBS55.svg'),
            observation='Official front and side show an inset bottom plinth. Actual IFC contains only three connected components and no separate plinth below the cylinder at z=0.',
            conclusion='Base is absent from this actual Body. Its exact dimensions and the allocation of the 100 mm total-height difference are not established.',
            invented_base_added=False,requires_upstream_geometry_correction=True),
        semantic_notes=['220 mm cylindrical body; no separately modelled doors, drawers, base or hinges.',
            'Top cap contains a 10 mm recessed collar at z=220..230, radius 266 mm, and a 20 mm top at z=230..250, radius 280 mm.',
            'The circular opaque top occludes all internal triangulation and most of the handle; only its external protrusion is shown in Plan.',
            'Side: body occludes the part of the handle within the body projection; Front: handle is closer than body.',
            'Open component edge paths end at genuine occlusion or shared contact; combined outside contour is closed.'],
        unsupported_invented_features=[],meaningless_coplanar_internal_lines=0,
        actual_height_mm=250,official_family_height_mm=350,height_discrepancy_mm=100,
        review_status='visual_review_pending',derived_ifc_write_allowed=False,formal_ifc_write_allowed=False)
    asset=Path('/Users/jiaxinchen/Documents/Projects/仁恒-滨海湾/Blender文件/assets/Baxter Models/Imported/Beside/3D/Beside_55x58xh35.glb')
    with asset.open('rb') as stream:
        stream.read(12)
        size,_=struct.unpack('<II',stream.read(8))
        glb=json.loads(stream.read(size))
    audit['base_audit'].update(local_asset_path=str(asset),local_asset_sha256=sha256(asset),
        local_asset_provenance='User local imported Baxter model; official-download origin not independently verified.',
        local_asset_components=[dict(name=m.get('name'),bounds=[dict(minimum=glb['accessors'][p['attributes']['POSITION']]['min'],maximum=glb['accessors'][p['attributes']['POSITION']]['max']) for p in m['primitives']]) for m in glb['meshes']],
        base_in_local_asset_mm=dict(diameter=475,height=30,minimum_z=0,maximum_z=30),
        units_evidence='Local asset overall height 35 matches official 350 mm, hence listed mesh coordinates interpreted as centimetres.',
        conclusion='Official dimension SVG and user local GLB both show a base. GLB Base is diameter475 x height30 mm; project Body has no base. Therefore the 100 mm height difference is not solely a missing base. Candidate remains actual Body only pending choice about upstream model repair.')
    write_json(audit_path,audit)
    candidate_path=DIR/'candidate-representations.json'
    candidate=load_json(candidate_path)
    candidate['project_local_views'] = {v:dict(proxy_paths_mm=[e['path_mm'] for e in es],official_cad_paths_mm=[]) for v,es in project_semantics.items()}
    candidate['review_rotation_z_degrees'] = math.degrees(angle)
    entries=[]
    for view,axes in AXES.items():
        raw=display_edge_sample(projected_raw_edges(vertices,faces,axes),maximum=1800)
        for suffix, edges in [('',None),('-body-comparison',raw)]:
            svg=DIR/f'{view}{suffix}.svg'
            svg.write_text(render(view,semantics[view],edges),encoding='utf8')
            write_uncached_png_preview(svg,DIR/f'{view}{suffix}-preview.png')
        candidate['views'][view].update(proxy_paths_mm=[e['path_mm'] for e in semantics[view]],
            semantic_paths=semantics[view],official_cad_paths_mm=[])
        entries.append(dict(view=view,svg=relative(DIR/f'{view}.svg'),svg_sha256=sha256(DIR/f'{view}.svg'),
            preview=relative(DIR/f'{view}-preview.png'),body_comparison=relative(DIR/f'{view}-body-comparison.svg'),
            projection_axes=list(axes),silhouette_path_count=len(semantics[view]),semantic_path_count=len(semantics[view]),
            drawing_line_source_kind=SOURCE,drawing_line_source_label_zh=LABEL,official_cad_path_count=0,blue_line_present=False))
    candidate.update(semantic_segmentation=relative(audit_path),semantic_segmentation_sha256=sha256(audit_path),
        derived_ifc_write_allowed=False,formal_ifc_write_allowed=False,review_status='visual_review_pending')
    write_json(candidate_path,candidate)
    render_contact_svg(DIR/'semantic-contact-sheet.svg','BST02 semantic single-product review',[DIR/f'{v}-preview.png' for v in AXES])
    write_uncached_png_preview(DIR/'semantic-contact-sheet.svg',DIR/'semantic-contact-sheet.png')
    subprocess.run([sys.executable,str(ROOT/'pipeline/scripts/render_bst02_project_context.py')],cwd=ROOT,check=True)
    panels=[]
    for n,view in enumerate(AXES):
        encoded=base64.b64encode((DIR/f'{view}.svg').read_bytes()).decode('ascii')
        panels.append(f'<image x="{n*600}" y="90" width="600" height="430" href="data:image/svg+xml;base64,{encoded}"/>')
    for n,camera in enumerate(['plan','front-elevation','side-elevation','iso']):
        encoded=base64.b64encode((DIR/f'bonsai-camera-{camera}.png').read_bytes()).decode('ascii')
        panels.append(f'<text x="{n*450+20}" y="570" font-size="22">Body {camera}</text><image x="{n*450}" y="590" width="450" height="340" href="data:image/png;base64,{encoded}"/>')
    for n,name in enumerate(['furniture-plan','side-elevation']):
        encoded=base64.b64encode((DIR/f'project-context-{name}-review-preview.png').read_bytes()).decode('ascii')
        panels.append(f'<text x="{n*900+30}" y="980" font-size="22">Project {name} review overlay</text><image x="{n*900}" y="1000" width="900" height="710" href="data:image/png;base64,{encoded}"/>')
    (DIR/'review-contact-sheet.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg" width="1800" height="1740" viewBox="0 0 1800 1740"><rect width="1800" height="1740" fill="white"/><g font-family="Arial"><text x="30" y="50" font-size="32">BST02 semantic candidate — pending review</text>'+''.join(panels)+'</g></svg>')
    write_uncached_png_preview(DIR/'review-contact-sheet.svg',DIR/'review-contact-sheet.png')
    manifest_path=DIR/'manifest.json'
    manifest=load_json(manifest_path)
    manifest.update(generator='pipeline/scripts/bst02_semantic_review.py',generated_at=audit['generated_at'],views=entries,
        semantic_segmentation=dict(path=relative(audit_path),sha256=sha256(audit_path),classified_component_count=3),
        candidate_representations_sha256=sha256(candidate_path),review_status='visual_review_pending',
        approved_for_drawing_ifc=False,derived_ifc_write_allowed=False,formal_ifc_write='not performed',
        semantic_contact_sheet=relative(DIR/'semantic-contact-sheet.png'),
        review_contact_sheet=relative(DIR/'review-contact-sheet.png'),review_contact_sheet_sha256=sha256(DIR/'review-contact-sheet.png'))
    manifest['project_context']['manifest_sha256']=sha256(DIR/'project-context-manifest.json')
    write_json(manifest_path,manifest)
    approval_path=ROOT/'pipeline/decisions/bst02-drawing-approval.json'
    approval=load_json(approval_path)
    approval.update(candidate_manifest_sha256=sha256(manifest_path),status='pending',approved_views=[],
        derived_ifc_write_allowed=False,formal_authoritative_ifc_write_allowed=False,
        note='User requested semantic reconstruction. Three actual connected components, projection and occlusion audited. Current 250 mm height preserved. New candidate awaits approval; no IFC write.')
    write_json(approval_path,approval)
    import bst02_review
    bst02_review.write_index(manifest)
    index=DIR/'index.html'
    extra='<h2>Semantic single-product review</h2><p>基于原始高模几何生成的简化图纸表达。审核坐标摆正 -7.031244°，把手居中；项目placement未改变。官方尺寸图和本地GLB均有内缩底座，本地模型底座Ø475×30mm；实际Body未建底座。100mm总高差并非仅缺底座。保留实际现高250mm，未添加补件。</p><a href="bst02-semantic-segmentation.json">Semantic visibility audit</a><a href="semantic-contact-sheet.png">Latest semantic contact sheet</a>'
    extra+=''.join(f'<a href="{v}-body-comparison.svg">{v.title()} grey Body comparison</a>' for v in AXES)
    index.write_text(index.read_text().replace('Black line = drawing representation from the exact Gessi 54294 official native DWG.', 'Black line = semantic boundaries derived from actual BST02 Body.').replace('<h2>Three-view review</h2>',extra+'<h2>Three-view review</h2>'))
    assert sha256(FORMAL)==BASELINE
    validate()
    print(json.dumps({v:len(es) for v,es in semantics.items()}))


if __name__=='__main__':
    main()
