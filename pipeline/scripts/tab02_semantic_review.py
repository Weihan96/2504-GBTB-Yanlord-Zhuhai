#!/usr/bin/env python3
"""TAB02 semantic single-product candidate from the three actual Body items."""
import html
import base64
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone

import ifcopenshell
import ifcopenshell.geom
from shapely.geometry import Polygon
from shapely.ops import unary_union

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from int1_highpoly_type_review import projected_raw_edges, display_edge_sample
from render_sis04_project_context import write_uncached_png_preview
from wd03_semantic_review import path_d, render_contact_svg

DIR = ROOT / 'output/review/highpoly-types/tab02'
IFC = DIR / 'RODA-Bernardo-367-TAB02-bonsai-isolated.ifc'
FORMAL = ROOT / '2504 GBTB Yanlord Zhuhai.ifc'
BASELINE = '7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c'
GUID = '2eX84IyLr8_e34nuOoRPLQ'
SOURCE = 'geometry_derived_simplified_proxy'
LABEL = '基于原始高模几何生成的简化图纸表达'
AXES = {'plan': (0, 1), 'front': (0, 2), 'side': (1, 2)}


def body():
    model = ifcopenshell.open(IFC)
    shape = ifcopenshell.geom.create_shape(ifcopenshell.geom.settings(), model.by_guid(GUID))
    g = shape.geometry
    vertices = [tuple(c * 1000 for c in g.verts[i:i+3]) for i in range(0, len(g.verts), 3)]
    faces = [g.faces[i:i+3] for i in range(0, len(g.faces), 3)]
    groups = defaultdict(list)
    for face, item in zip(faces, g.item_ids):
        groups[item].append(face)
    records = []
    for item, item_faces in groups.items():
        indices = {i for f in item_faces for i in f}
        lo = [min(vertices[i][a] for i in indices) for a in range(3)]
        hi = [max(vertices[i][a] for i in indices) for a in range(3)]
        size = [hi[a] - lo[a] for a in range(3)]
        if abs(size[0]-300)<0.01 and abs(size[2]-40)<0.01:
            semantic = 'stone_base'
        elif abs(size[0]-500)<0.01 and abs(size[2]-3)<0.01:
            semantic = 'rounded_tabletop'
        elif abs(size[0]-50)<0.01 and abs(size[2]-627)<0.01:
            semantic = 'central_column'
        else:
            raise RuntimeError(f'Unclassified Body item {item}: {size}')
        records.append(dict(item_id=item, ifc_entity=model.by_id(item).is_a(), semantic=semantic,
                            face_count=len(item_faces), minimum_mm=lo, maximum_mm=hi, size_mm=size))
    assert len(records)==3 and len(faces)==1860 and len(vertices)==936
    return vertices, faces, groups, records


def projected_component(vertices, faces, axes):
    triangles = [Polygon([(round(vertices[i][axes[0]],6),round(vertices[i][axes[1]],6)) for i in face]) for face in faces]
    union = unary_union([p for p in triangles if p.area > 1e-8])
    assert union.geom_type == 'Polygon'
    # Floating-point triangle intersections can leave zero-area pinholes.
    assert sum(Polygon(ring).area for ring in union.interiors) < 1e-8
    union = Polygon(union.exterior)
    simplified = union.simplify(0.05, preserve_topology=True)
    assert union.boundary.hausdorff_distance(simplified.boundary) <= 0.050001
    return list(simplified.exterior.coords), union.boundary.hausdorff_distance(simplified.boundary)


def make_views(vertices, groups, records):
    result = {}
    for view, axes in AXES.items():
        entries = []
        for record in records:
            semantic = record['semantic']
            path, error = projected_component(vertices, groups[record['item_id']], axes)
            if view == 'plan' and semantic != 'rounded_tabletop':
                continue  # Opaque tabletop fully occludes the column and smaller base.
            paths = [path]
            if view != 'plan' and semantic == 'central_column':
                # Contact edges coincide with the tabletop underside and base top.
                # Keep column sides only so the shared contact lines are drawn once.
                lo, hi = record['minimum_mm'], record['maximum_mm']
                a = axes[0]
                paths = [[[lo[a],lo[2]],[lo[a],hi[2]]], [[hi[a],lo[2]],[hi[a],hi[2]]]]
            for n, p in enumerate(paths):
                entries.append(dict(id=f'{view}.{semantic}.{n}', component=semantic,
                    item_id=record['item_id'], kind='visible_component_boundary',
                    path_mm=[[round(x,6),round(y,6)] for x,y in p],
                    projection_boundary_error_mm=error,
                    sides=[semantic, 'exterior_void_or_adjacent_contact_component']))
        result[view] = entries
    return result


def render(view, entries, raw=None):
    pts = [p for e in entries for p in e['path_mm']]
    x0,x1 = min(p[0] for p in pts),max(p[0] for p in pts)
    y0,y1 = min(p[1] for p in pts),max(p[1] for p in pts)
    scale = min(860/(x1-x0),660/(y1-y0))
    transform = lambda p: (700+(p[0]-(x0+x1)/2)*scale, 560-(p[1]-(y0+y1)/2)*scale)
    raw_markup = '' if raw is None else f'<path class="original-highpoly" d="{path_d(raw, transform)}" fill="none" stroke="#999" stroke-opacity="0.25" stroke-width="0.6"/>'
    paths = ''.join(f'<path data-semantic-id="{e["id"]}" data-ifc-item-id="{e["item_id"]}" d="{path_d([e["path_mm"]],transform)}"/>' for e in entries)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="980" viewBox="0 0 1400 980">
<rect width="1400" height="980" fill="white"/><text x="55" y="55" font-family="Arial" font-size="30">RODA Bernardo 367 / TAB02 — {view.upper()} semantic single product</text>
<text x="55" y="96" font-family="Arial" font-size="19">{LABEL}</text>
<text x="55" y="131" font-family="Arial" font-size="17">{'Grey: actual Body comparison / ' if raw is not None else ''}Black: visible semantic boundaries · pending approval</text>
{raw_markup}<g class="simplified-proxy-silhouette geometry-derived" data-source-kind="{SOURCE}" fill="none" stroke="#000" stroke-width="2" stroke-linejoin="round">{paths}</g>
<text x="55" y="940" font-family="Arial" font-size="16">3 mm rounded tabletop / diameter 50 mm column / diameter 300 x 40 mm stone base · {GUID}</text></svg>'''


def main():
    assert sha256(FORMAL)==BASELINE
    vertices, faces, groups, records = body()
    semantics = make_views(vertices,groups,records)
    audit_path = DIR/'tab02-semantic-segmentation.json'
    audit = dict(schema_version=1, generated_at=datetime.now(timezone.utc).isoformat(),
        representative_global_id=GUID, source_kind=SOURCE, source_label_zh=LABEL,
        isolated_ifc=relative(IFC), isolated_ifc_sha256=sha256(IFC),
        mesh_vertex_count=len(vertices), mesh_face_count=len(faces), representation_item_count=3,
        components=records, views=semantics, meaningless_coplanar_internal_lines=0,
        coordinate_frame='IFC product local XYZ in mm; no translation, reflection, fitting or scaling',
        visibility=dict(plan='Opaque tabletop occludes column and base completely.',
                        front='Camera from -Y; all three components visible.', side='Camera from -X; all three components visible.'),
        omissions=dict(separate_mounting_hardware='No separate mounting hardware exists in the actual Body; none invented.',
                       column_contact_edges='Coincide with tabletop underside/base top; drawn once by those components.'),
        derived_ifc_write_allowed=False, formal_ifc_write_allowed=False, review_status='visual_review_pending',pass_=True)
    audit['pass'] = audit.pop('pass_')
    write_json(audit_path,audit)
    candidate_path=DIR/'candidate-representations.json'
    candidate=load_json(candidate_path)
    views=[]
    previews=[]
    for view,axes in AXES.items():
        entries=semantics[view]
        raw=display_edge_sample(projected_raw_edges(vertices,faces,axes),maximum=975)
        for suffix, edges in [('',None),('-body-comparison',raw)]:
            svg=DIR/f'{view}{suffix}.svg'
            svg.write_text(render(view,entries,edges),encoding='utf8')
            write_uncached_png_preview(svg,DIR/f'{view}{suffix}-preview.png')
        previews.append(DIR/f'{view}-preview.png')
        shutil.copyfile(DIR/f'{view}-preview.png',DIR/f'{view}-review-preview.png')
        candidate['views'][view].update(proxy_paths_mm=[e['path_mm'] for e in entries],semantic_paths=entries,official_cad_paths_mm=[])
        views.append(dict(view=view,svg=relative(DIR/f'{view}.svg'),svg_sha256=sha256(DIR/f'{view}.svg'),
            preview=relative(DIR/f'{view}-preview.png'),body_comparison=relative(DIR/f'{view}-body-comparison.svg'),
            projection_axes=list(axes),silhouette_path_count=1,semantic_path_count=len(entries),
            drawing_line_source_kind=SOURCE,drawing_line_source_label_zh=LABEL,official_cad_path_count=0,blue_line_present=False))
    candidate.update(semantic_segmentation=relative(audit_path),semantic_segmentation_sha256=sha256(audit_path),
                     derived_ifc_write_allowed=False,formal_ifc_write_allowed=False,review_status='visual_review_pending')
    write_json(candidate_path,candidate)
    render_contact_svg(DIR/'semantic-contact-sheet.svg','TAB02 semantic single-product review',previews)
    write_uncached_png_preview(DIR/'semantic-contact-sheet.svg',DIR/'semantic-contact-sheet.png')
    subprocess.run([sys.executable,str(ROOT/'pipeline/scripts/render_tab02_project_context.py')],cwd=ROOT,check=True)
    # Refresh the existing review-sheet entry point as well as the focused sheet.
    panels=[]
    for n,view in enumerate(AXES):
        encoded=base64.b64encode((DIR/f'{view}.svg').read_bytes()).decode('ascii')
        panels.append(f'<image x="{n*600}" y="90" width="600" height="430" href="data:image/svg+xml;base64,{encoded}"/>')
    for n,camera in enumerate(['plan','front-elevation','side-elevation','iso']):
        encoded=base64.b64encode((DIR/f'bonsai-camera-{camera}.png').read_bytes()).decode('ascii')
        panels.append(f'<text x="{n*450+20}" y="570" font-size="22">Body {camera}</text><image x="{n*450}" y="590" width="450" height="340" href="data:image/png;base64,{encoded}"/>')
    encoded=base64.b64encode((DIR/'project-context-furniture-plan-review-preview.png').read_bytes()).decode('ascii')
    panels.append(f'<text x="30" y="980" font-size="24">Project Plan review overlay (not a new Bonsai scene drawing)</text><image x="300" y="1000" width="1200" height="710" href="data:image/png;base64,{encoded}"/>')
    (DIR/'review-contact-sheet.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg" width="1800" height="1740" viewBox="0 0 1800 1740"><rect width="1800" height="1740" fill="white"/><g font-family="Arial"><text x="30" y="50" font-size="32">TAB02 semantic single product — latest pending candidate</text>'+''.join(panels)+'</g></svg>',encoding='utf8')
    write_uncached_png_preview(DIR/'review-contact-sheet.svg',DIR/'review-contact-sheet.png')
    manifest_path=DIR/'manifest.json'
    manifest=load_json(manifest_path)
    manifest.update(generator='pipeline/scripts/tab02_semantic_review.py',generated_at=audit['generated_at'],views=views,
        semantic_segmentation=dict(path=relative(audit_path),sha256=sha256(audit_path),classified_component_count=3),
        candidate_representations_sha256=sha256(candidate_path),review_status='visual_review_pending',
        approved_for_drawing_ifc=False,derived_ifc_write_allowed=False,formal_ifc_write='not performed',
        semantic_contact_sheet=relative(DIR/'semantic-contact-sheet.png'))
    manifest['review_contact_sheet']=relative(DIR/'review-contact-sheet.png')
    manifest['review_contact_sheet_sha256']=sha256(DIR/'review-contact-sheet.png')
    manifest['project_context']['manifest_sha256']=sha256(DIR/'project-context-manifest.json')
    write_json(manifest_path,manifest)
    approval_path=ROOT/'pipeline/decisions/tab02-drawing-approval.json'
    approval=load_json(approval_path)
    approval.update(candidate_manifest_sha256=sha256(manifest_path),status='pending',approved_views=[],
        derived_ifc_write_allowed=False,formal_authoritative_ifc_write_allowed=False,
        note='Semantic candidate uses 3 actual Body components, with per-view occlusion and real contact boundaries. Pending individual human approval; no IFC write.')
    write_json(approval_path,approval)
    import tab02_review
    tab02_review.write_index(manifest)
    index=DIR/'index.html'
    extra='<h2>Semantic single-product review</h2><p>基于原始高模几何生成的简化图纸表达。三组件：3 mm 圆角桌面、Ø50 立柱、Ø300×40 石材底座。Plan 遮挡隐藏底座和立柱。</p><a href="tab02-semantic-segmentation.json">Semantic audit</a><a href="semantic-contact-sheet.png">Latest semantic contact sheet</a>'
    extra+=''.join(f'<a href="{v}-body-comparison.svg">{v.title()} grey Body comparison</a>' for v in AXES)
    index.write_text(index.read_text().replace('<h2>Three-view review</h2>',extra+'<h2>Three-view review</h2>'),encoding='utf8')
    assert sha256(FORMAL)==BASELINE
    print({v:len(es) for v,es in semantics.items()})


if __name__=='__main__':
    main()
