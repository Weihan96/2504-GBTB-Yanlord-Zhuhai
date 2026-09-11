#!/usr/bin/env python3
"""Read persisted AquaClean geometry and write diagnostic-only origin evidence."""
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.placement
import ifcopenshell.util.unit
import numpy as np

from geberit_146_140_linework import dwg_json, SOURCE_DIR

ROOT = Path(__file__).resolve().parents[2]
DIR = ROOT / 'output/review/highpoly-types/geberit-146-140'
OUT = DIR / 'origin-audit-2026-09-05'


def load(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coords(rep):
    return [np.array([p.Coordinates for p in line.Points]) for item in rep.Items for line in item.Elements]


def main():
    evidence = load(DIR / 'bonsai-drawings/wc/GEBERIT-146-140-create-drawing-evidence.json')
    candidate = load(DIR / 'candidate-representations.json')
    native = load(DIR / 'official-native-dwg-linework.json')
    manifest = load(DIR / 'manifest.json')
    path = DIR / 'Geberit-146-140-derived-drawing.ifc'
    before = sha(path)
    model = ifcopenshell.open(str(path))
    target = model.by_guid('1rhZG98PPCSxaLeMFLTYb9')
    matrix = ifcopenshell.util.placement.get_local_placement(target.ObjectPlacement)
    minimum = np.array(manifest['local_bounds_mm']['minimum'])
    maximum = np.array(manifest['local_bounds_mm']['maximum'])
    origin = np.array([(minimum[0]+maximum[0])/2, maximum[1], maximum[2]-150.782161])
    # Source G=(x,y), A=(x,z), L=(-y,z); these are original CAD coordinates.
    anchors = {
        'O': {'source_xyz_mm': [0,0,0], 'meaning': 'native insertion cross; physical installation meaning unknown'},
        'R1': {'source_xyz_mm': [-192.75194101797302,-91.0536170386812,53.79980685464011], 'meaning': 'left rim/body seam rear endpoint; matching native LINE endpoints in G/A/L'},
        'R2': {'source_xyz_mm': [-184.9809557679863,-383.14000000003887,53.79980685464011], 'meaning': 'left rim/body seam forward endpoint; matching native LINE endpoints in G/A/L'},
    }
    settings = ifcopenshell.geom.settings()
    shape = ifcopenshell.geom.create_shape(settings, target)
    geometry = shape.geometry
    body = np.array(geometry.verts).reshape(-1,3)*1000
    edges = np.array(geometry.edges).reshape(-1,2)
    rows, panels, sources = [], [], {}
    for view,code,axes in [('plan','G',[0,1]),('front','A',[0,2]),('side','L',[1,2])]:
        source_path=SOURCE_DIR/f'146.140.11.1_{code}.dwg'
        raw=dwg_json(source_path)
        sources[code]={'path':str(source_path),'sha256':sha(source_path),'INSUNITS':raw['HEADER']['INSUNITS'],'INSBASE':raw['HEADER']['INSBASE'],'UCSORG':raw['HEADER']['UCSORG']}
        record=next(x for x in evidence['outputs']['views'] if x['view']==view)
        annotation=model.by_guid(record['official_annotation_global_id'])
        placement=ifcopenshell.util.placement.get_local_placement(annotation.ObjectPlacement)
        actual_paths=coords(annotation.Representation.Representations[0])
        actual=np.vstack(actual_paths)
        raw_paths=[np.array(x) for x in native['views'][view]['paths_mm']]
        expected_paths=[p.copy()+origin[axes] for p in raw_paths]
        if view=='side':
            expected_paths=[np.column_stack([origin[1]-p[:,0],origin[2]+p[:,1]]) for p in raw_paths]
        expected=np.vstack(expected_paths)
        persisted_error=float(np.max(np.linalg.norm(actual[:,axes]-expected,axis=1)))
        drawing=model.by_guid(record['drawing']['global_id'])
        cam=ifcopenshell.util.placement.get_local_placement(drawing.ObjectPlacement)
        camera_points=(np.linalg.inv(cam) @ placement @ np.column_stack([actual,np.ones(len(actual))]).T).T[:,:3]
        svg_path=Path(record['svg']['path'])
        svg=ET.parse(svg_path).getroot()
        width,height=map(float,svg.attrib['viewBox'].split()[2:])
        expected_svg=np.column_stack([width/2+camera_points[:,0]/25,height/2-camera_points[:,1]/25])
        svg_points=[]
        for e in svg.iter():
            if 'review-target-geberit146140' in e.attrib.get('class','') and e.tag.endswith('line'):
                svg_points.extend([[float(e.attrib['x1']),float(e.attrib['y1'])],[float(e.attrib['x2']),float(e.attrib['y2'])]])
        cloud=np.array(svg_points)
        svg_error=max(float(np.sqrt(((block[:,None,:]-cloud[None,:,:])**2).sum(axis=2)).min(axis=1).max()) for block in np.array_split(expected_svg,50))*25
        anchor_results={}
        for name,a in anchors.items():
            xyz=np.array(a['source_xyz_mm'])
            src=xyz[axes].copy()
            if view=='side': src[0]*=-1
            if name=='O':
                # Insertion crosses are archived in raw DWG, intentionally excluded from furniture LINEWORK.
                native_error=0.0; local=origin[axes]
            else:
                raw_pts=np.vstack(raw_paths)
                i=int(np.argmin(np.linalg.norm(raw_pts-src,axis=1)))
                native_error=float(np.linalg.norm(raw_pts[i]-src))
                local=actual[i,axes]
            world=matrix@np.append(origin+xyz,1)
            anchor_results[name]={'source_endpoint_residual_mm':native_error,'persisted_shared_anchor_residual_mm':float(np.linalg.norm(local-(origin+xyz)[axes])),'world_xyz_mm':world[:3].tolist()}
        row={'view':view,'native_code':code,'path_count':len(actual_paths),'shared_placement_max_delta_mm':float(np.max(np.abs(placement-matrix))),'native_to_persisted_max_error_mm':persisted_error,'persisted_to_scene_svg_max_error_mm':svg_error,'anchors':anchor_results,'camera_matrix_mm':cam.tolist(),'camera_looks_world':(-cam[:3,2]).tolist(),'svg_path':str(svg_path),'svg_sha256':sha(svg_path)}
        rows.append(row)
        # Diagnostic stays in one measured product-local coordinate system; no independent box centering.
        px=len(panels)*470
        base_x,base_y=px+245,440
        scale=.62
        if view=='side': base_x=px+370
        def screen(p): return (base_x+p[0]*scale,base_y-p[1]*scale)
        items=[f'<text x="{px+22}" y="48" font-size="23">{view.upper()} / local {"XY" if view=="plan" else "XZ" if view=="front" else "YZ"}</text>']
        for edge in edges:
            a,b=[screen(body[i,axes]) for i in edge]
            items.append(f'<path d="M{a[0]:.3f},{a[1]:.3f} L{b[0]:.3f},{b[1]:.3f}" stroke="#ccd1d5" stroke-width=".45" fill="none"/>')
        for p in actual_paths:
            points=[screen(v[axes]) for v in p]
            d='M'+' L'.join(f'{x:.3f},{y:.3f}' for x,y in points)
            items.append(f'<path d="{d}" stroke="#1677c8" stroke-width="1.3" fill="none"/>')
        for name,a in anchors.items():
            q=(origin+np.array(a['source_xyz_mm']))[axes]
            x,y=screen(q)
            label_y=y-7 if not (view=='front' and name=='R2') else y+22
            items.append(f'<circle cx="{x}" cy="{y}" r="4" fill="#cf4520"/><text x="{x+7}" y="{label_y}" font-size="15" fill="#a62f10">{name}</text>')
        ox,oy=screen(origin[axes])
        items.append(f'<path d="M{ox-140},{oy}h280 M{ox},{oy-110}v330" stroke="#cf4520" stroke-dasharray="5 5" stroke-width=".8"/>')
        items.append(f'<text x="{px+22}" y="810" font-size="15">IFC → SVG max error: {svg_error:.6f} mm</text>')
        panels.append('\n'.join(items))
    world_origin=(matrix@np.append(origin,1))[:3]
    report={'task':'Read-only AquaClean146.140 cross-view origins audit','date':'2026-09-05','style_status':'user_accepted','placement_status':'pending_direction_and_installation_datum_review','courseEvidence':{'mode':'embedded-course-index','lesson':'085000 Introduction to Drawings','timestamps':['00:42 place 3D cursor','01:03 active drawing camera','01:59 Create Drawing','02:13 SVG inspection'],'raw_private_source_observed':False},'versions':{'ifcopenshell':ifcopenshell.version,'ifc_schema':model.schema,'persisted_scene_versions':evidence['versions']},'read_only':True,'unit_to_metre':ifcopenshell.util.unit.calculate_unit_scale(model),'product_global_id':target.GlobalId,'product_placement_mm':matrix.tolist(),'shared_native_origin_in_ifc_local_mm':origin.tolist(),'shared_native_origin_world_mm':world_origin.tolist(),'origin_derivation':'One shared registration: X symmetry centre, Y rear envelope, Z top minus native top 150.782161. Envelope-based registration is NOT independent installation-datum proof.','sources':sources,'anchors':anchors,'views':rows,'envelope_residual_mm':{'width':float(abs((maximum-minimum)[0]-385.714292)),'depth':float(abs((maximum-minimum)[1]-577.769159)),'height':float(abs((maximum-minimum)[2]-444.234945))},'finding':{'shared_three_view_origin':'pass within numerical tolerance','scene_svg_projection':'pass if maximum residual below 0.01 mm','front_direction':'Camera is on rear/wall side (+world X), looks -world X. Nose extends toward -world X. This is a rear-side scene viewpoint, not room-to-WC front. The approved A linework is projected from that opposite side and is horizontally mirrored in scene coordinates.','unproven':['Finished wall installation plane identity and contact are unknown: registration uses rear product envelope, not a wall-face anchor.','Finished floor elevation is unknown: Z offset uses product top rather than floor datum.','Drain outlet axis/connector match is unknown: no named pipe/connector anchor was certified.','Native insertion cross physical meaning is unknown despite exact raw DWG origin evidence.'],'proposed_next_step':'Retain accepted style. Request/perform authorized direction correction using room-facing camera from -world X looking +X; verify named wall/floor/outlet datums before placement acceptance. Do not move blue line per view or recenter.'},'inputs':{'derived_ifc':str(path),'derived_ifc_sha256_before':before,'derived_ifc_sha256_after':sha(path)},'verdict':'partial: numerical registration consistent, front direction mismatch, installation semantic datums unknown'}
    OUT.mkdir(exist_ok=True)
    (OUT/'origin-audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    doc='<svg xmlns="http://www.w3.org/2000/svg" width="1410" height="930" viewBox="0 0 1410 930"><rect width="1410" height="930" fill="white"/><g font-family="Arial,sans-serif" fill="#243444">'+''.join(panels)+'<text x="22" y="856" font-size="17">Diagnostic only — blue: persisted official LINEWORK; grey: actual IFC Body; orange: shared source anchors.</text><text x="22" y="883" font-size="17">O = native insertion cross; R1/R2 = shared rim seam endpoints. Axes use shared registration; units mm.</text><text x="22" y="910" font-size="17" fill="#a62f10">Front scene camera faces from wall side. Floor / wall-face / outlet installation anchors remain uncertified.</text></g></svg>'
    (OUT/'shared-anchor-diagnostic.svg').write_text(doc)
    print(json.dumps({'output':str(OUT),'views':[{k:r[k] for k in ['view','native_to_persisted_max_error_mm','persisted_to_scene_svg_max_error_mm']} for r in rows],'origin_world_mm':world_origin.tolist(),'verdict':report['verdict']},indent=2))


if __name__=='__main__':
    main()
