#!/usr/bin/env python3
"""WD02 semantic projections of actual Body components, review artifacts only."""
import json
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone

import ifcopenshell
import ifcopenshell.geom
from shapely.geometry import Polygon, GeometryCollection
from shapely.ops import unary_union

import wd03_semantic_review as graphics
from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from int1_highpoly_type_review import display_edge_sample, projected_raw_edges
from render_sis04_project_context import write_uncached_png_preview

PRODUCT = ROOT / 'output/review/highpoly-types/wd02'
IFC = PRODUCT / 'Poliform-Senzafine-WD02-bonsai-isolated.ifc'
GID = '3yuXF4PtnBHgIXDlmXFJ$7'
SOURCE = 'geometry_derived_simplified_proxy'
LABEL = '基于原始高模几何生成的简化图纸表达'
AXES = {'plan': (0, 1), 'front': (0, 2), 'side': (1, 2)}
# IDs refer to immutable, hash-checked isolated Body representation items.
SEMANTICS = {
91:'bottom_plinth',123:'bottom_deck',195:'top_deck',267:'upper_shelf',275:'back_panel',
291:'left_side_panel',307:'right_side_panel',581:'clothes_rail',693:'left_rail_bracket',805:'right_rail_bracket',
855:'drawer_unit_top',901:'middle_drawer_front_frame',909:'middle_drawer_front_infill',1083:'middle_drawer_pull',1163:'middle_drawer_bottom',
1209:'lower_drawer_front_frame',1217:'lower_drawer_front_infill',1391:'lower_drawer_pull',1471:'lower_drawer_bottom',1479:'drawer_top_inset',
1559:'drawer_unit_right_side',1639:'drawer_unit_left_side',1719:'drawer_unit_back',1799:'drawer_right_front_edge',1879:'drawer_left_front_edge',
1959:'upper_drawer_bottom',1979:'upper_drawer_front',2006:'door_bottom_frame',2031:'door_top_frame',2039:'door_glass',
2058:'bottom_hinge',2077:'top_hinge',2099:'door_right_frame',2121:'door_left_frame',2129:'door_handle',
}


def load_components():
    assert sha256(IFC) == '87d1ca3990fc4dd77d54ab709518445fd87bfbc5837a4b5b63762517baaef1d4'
    model = ifcopenshell.open(IFC)
    shape = ifcopenshell.geom.create_shape(ifcopenshell.geom.settings(), model.by_guid(GID))
    mesh = shape.geometry
    vertices = [tuple(round(value * 1000, 6) for value in mesh.verts[i:i+3]) for i in range(0,len(mesh.verts),3)]
    faces = [tuple(mesh.faces[i:i+3]) for i in range(0,len(mesh.faces),3)]
    grouped = defaultdict(list)
    for face,item in zip(faces,mesh.item_ids): grouped[int(item)].append(face)
    assert set(grouped) == set(SEMANTICS), 'Body component identity changed'
    records=[]
    for item, triangles in grouped.items():
        ids=set(i for face in triangles for i in face)
        lo=[min(vertices[i][a] for i in ids) for a in range(3)]
        hi=[max(vertices[i][a] for i in ids) for a in range(3)]
        records.append(dict(item_id=item,semantic=SEMANTICS[item],minimum_mm=lo,maximum_mm=hi,
                            size_mm=[b-a for a,b in zip(lo,hi)],triangle_count=len(triangles),triangles=triangles))
    return vertices,faces,records


def polygons(geom):
    if geom.is_empty:return []
    if geom.geom_type == 'Polygon':return [geom]
    return [p for g in geom.geoms for p in polygons(g)]


def semantic_projection(view,vertices,records):
    axes=AXES[view]
    depth_axis={'plan':2,'front':1,'side':0}[view]
    nearest=max if view=='plan' else min
    # Actual projected face coverage; a part's bounding box is never line geometry.
    coverage={}
    for r in records:
        triangles=[Polygon([(vertices[i][axes[0]],vertices[i][axes[1]]) for i in face]) for face in r['triangles']]
        coverage[r['item_id']]=unary_union([p for p in triangles if p.area > 1e-7]).buffer(0)
    order=sorted(records,key=lambda r:nearest(r['maximum_mm'][depth_axis],r['minimum_mm'][depth_axis]),reverse=view=='plan')
    opaque=GeometryCollection()
    paths=[];visible_records=[]
    for record in order:
        item=record['item_id'];semantic=record['semantic'];full=coverage[item]
        visible=full.difference(opaque)
        # Closed glass is transparent only when looking through the door face.
        transparent=view=='front' and semantic=='door_glass'
        for idx,poly in enumerate(polygons(visible)):
            if poly.area < 1.0:continue
            simple=poly.simplify(0.3,preserve_topology=True)
            for ring_idx,ring in enumerate([simple.exterior,*simple.interiors]):
                points=[[round(x,6),round(y,6)] for x,y in ring.coords]
                paths.append(dict(id=f'{view}.{semantic}.{idx}.{ring_idx}',path=points,
                                  kind='visible_component_boundary',item_id=item,
                                  sides=[semantic,'adjacent_component_or_exterior'],
                                  transparent_drawing_convention=transparent))
        visible_records.append(dict(item_id=item,semantic=semantic,projected_area_mm2=full.area,
                                    visible_area_mm2=visible.area,occluded_area_mm2=full.area-visible.area,
                                    transparent=transparent))
        if not transparent:opaque=unary_union([opaque,full])
    return paths,visible_records


def main():
    assert sha256(graphics.FORMAL_IFC)==graphics.FORMAL_SHA256
    vertices,faces,records=load_components()
    graphics.GLOBAL_ID=GID
    graphics.PRODUCT_DIR=PRODUCT
    all_paths={};visibility={};views=[];previews=[];diagnostics=[]
    for view,axes in AXES.items():
        semantics,visibility[view]=semantic_projection(view,vertices,records)
        all_paths[view]=semantics
        raw=display_edge_sample(projected_raw_edges(vertices,faces,axes),maximum=975)
        svg=PRODUCT/f'{view}.svg'
        content=graphics.render_review_svg(view,raw,semantics).replace('WD03','WD02').replace('semantically validated drawing boundaries','visible component boundaries (glass transparent in Front)')
        # Many components are audited in JSON; keep the small drawing legend legible.
        import re
        content=re.sub(r'<text x="1080"[^>]*>•.*?</text>','',content)
        content=content.replace('Every internal line separates components or depth.','35 actual Body components; face-union boundaries.')
        svg.write_text(content,encoding='utf8')
        png=PRODUCT/f'{view}-preview.png';write_uncached_png_preview(svg,png);previews.append(png)
        diag=PRODUCT/f'wd02-semantic-diagnostic-{view}.svg'
        black=graphics.render_review_svg(view,[],semantics).replace('WD03','WD02').replace('Grey = actual IFC Body mesh · Black = semantically validated drawing boundaries · no official CAD / no blue line','Black-only semantic single-product drawing · glass transparent in Front').replace('stroke-width="3"','stroke-width="1.5"')
        black=re.sub(r'<text x="1080"[^>]*>•.*?</text>','',black)
        diag.write_text(black,encoding='utf8')
        diag_png=diag.with_suffix('.png');write_uncached_png_preview(diag,diag_png);diagnostics.append(diag_png)
        views.append(dict(view=view,svg=relative(svg),svg_sha256=sha256(svg),preview=relative(png),preview_sha256=sha256(png),projection_axes=list(axes),
                          raw_edge_count=len(projected_raw_edges(vertices,faces,axes)),displayed_raw_edge_count=len(raw),silhouette_path_count=len(semantics),semantic_path_count=len(semantics),
                          drawing_line_source_kind=SOURCE,drawing_line_source_label_zh=LABEL,official_cad_path_count=0,blue_line_present=False))
    audit_path=PRODUCT/'wd02-semantic-segmentation.json'
    audit=dict(schema_version=1,generator='pipeline/scripts/wd02_semantic_review.py',representative_global_id=GID,source_kind=SOURCE,source_label_zh=LABEL,
               isolated_ifc=relative(IFC),isolated_ifc_sha256=sha256(IFC),mesh_vertex_count=len(vertices),mesh_face_count=len(faces),representation_item_count=len(records),
               components=[{k:v for k,v in r.items() if k!='triangles'} for r in records],
               projection_method='union of actual projected triangles per semantic component, front-to-back component occlusion, 0.3 mm topology-preserving simplification',
               coordinate_basis='IFC product local XYZ in mm; no inferred translation or world-axis swapping',
               view_directions={'plan':'from local +Z looking -Z','front':'from local -Y looking +Y through glass door','side':'from local -X looking +X'},
               visibility=visibility,views={v:{'paths':p,'path_count':len(p)} for v,p in all_paths.items()},
               limitations=['Opaque-component depth ordering uses the nearest component depth; this is a component-level orthographic simplification, not a per-pixel raytracer.',
                            'Glass is intentionally transparent in Front so shelf, hanging rail and drawers are readable; the saved shaded render uses an opaque-looking glass material.',
                            'Shared semantic boundaries may occur in adjacent closed component loops; they represent one physical interface.'],
               simplification_tolerance_mm=0.3,official_cad_used=False,pass_checks=True)
    write_json(audit_path,audit)
    candidate_path=PRODUCT/'candidate-representations.json';candidate=load_json(candidate_path)
    candidate.update(source_kind=SOURCE,source_label_zh=LABEL,semantic_segmentation=relative(audit_path),semantic_segmentation_sha256=sha256(audit_path),review_status='visual_review_pending',formal_ifc_write_allowed=False)
    candidate['views']={v:dict(projection_axes=list(AXES[v]),proxy_paths_mm=[p['path'] for p in paths],semantic_paths=[{k:x for k,x in p.items() if k!='path'} for p in paths],source_kind=SOURCE,official_cad_paths_mm=[]) for v,paths in all_paths.items()}
    write_json(candidate_path,candidate)
    for name,title,items in [('review-contact-sheet','WD02 semantic single-product review',previews),('wd02-semantic-diagnostic-contact-sheet','WD02 black-only semantic drawings',diagnostics)]:
        target=PRODUCT/f'{name}.svg';graphics.render_contact_svg(target,title,items);write_uncached_png_preview(target,target.with_suffix('.png'))
    subprocess.run([sys.executable,str(ROOT/'pipeline/scripts/render_wd02_project_context.py')],cwd=ROOT,check=True)
    profile_path=PRODUCT/'profile.json';profile=load_json(profile_path)
    profile['profiles']['wd02'].update(semantic_audit=relative(audit_path),semantic_sections={v:[p['id'] for p in paths] for v,paths in all_paths.items()},
                                       profile_source='actual_highpoly_component_semantics; manufacturer_family_identity_only')
    write_json(profile_path,profile)
    manifest_path=PRODUCT/'manifest.json';manifest=load_json(manifest_path)
    manifest.update(generated_at=datetime.now(timezone.utc).isoformat(),generator='pipeline/scripts/wd02_semantic_review.py',views=views,
                    semantic_segmentation={'path':relative(audit_path),'sha256':sha256(audit_path),'representation_item_count':35,'classified_component_count':35,'actual_projected_triangles_used':True,'bbox_geometry_used':False},
                    profile_register_sha256=sha256(profile_path),candidate_representations_sha256=sha256(candidate_path),review_status='visual_review_pending',approved_for_drawing_ifc=False,derived_ifc_write_allowed=False,
                    review_contact_sheet_sha256=sha256(PRODUCT/'review-contact-sheet.png'),formal_ifc_bytes_unchanged=sha256(graphics.FORMAL_IFC)==graphics.FORMAL_SHA256)
    manifest['project_context']['manifest_sha256']=sha256(PRODUCT/'project-context-manifest.json')
    write_json(manifest_path,manifest)
    approval_path=ROOT/'pipeline/decisions/wd02-drawing-approval.json';approval=load_json(approval_path)
    approval.update(status='pending',candidate_manifest_sha256=sha256(manifest_path),derived_ifc_write_allowed=False,formal_authoritative_ifc_write_allowed=False,
                    note='User requested semantic geometry-derived WD02 review. New candidate awaits explicit human approval; no IFC written.')
    write_json(approval_path,approval)
    graphics.write_index(manifest)
    index=PRODUCT/'index.html';index.write_text(index.read_text().replace('WD03','WD02').replace('wd03','wd02'),encoding='utf8')
    print(json.dumps({'paths':{v:len(p) for v,p in all_paths.items()},'components':len(records),'formal_unchanged':manifest['formal_ifc_bytes_unchanged']}))


if __name__=='__main__':main()
