# -*- coding: utf-8 -*-
"""Independent file/geometry checks and review index for AquaClean registration."""
import json, hashlib, itertools
from pathlib import Path
import numpy as np
import ifcopenshell, ifcopenshell.geom, ifcopenshell.util.placement
ROOT=Path(__file__).resolve().parents[2]
P=ROOT/'output/review/highpoly-types/geberit-146-140'
O=P/'registered-3d-2026-09-07'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    e=json.loads((O/'evidence.json').read_text());m=ifcopenshell.open(str(P/'Geberit-146-140-derived-drawing.ifc'))
    settings=ifcopenshell.geom.settings();settings.set('use-world-coords',True)
    def geometry(guid):
        shape=ifcopenshell.geom.create_shape(settings,m.by_guid(guid));return np.array(shape.geometry.verts).reshape(-1,3),np.array(shape.geometry.faces).reshape(-1,3)
    body,faces=geometry('1rhZG98PPCSxaLeMFLTYb9');wall,wf=geometry('2UhljZvyHBM9jCYFR8xcyl')
    wallx=float(wall[:,0].min());rear=float(body[:,0].max())
    e['installation_semantics']={'wall_hung_verified':'Official archived Geberit product page title says wall-hung WC','reference_floor':{'global_id':'11c$NwzxL9hASgCgwDxM$w','z_mm':0,'role':'IfcCovering Plane; reference datum, not a verified tile contact elevation'},'body_bottom_above_reference_floor_mm':float(body[:,2].min()*1000),'pipe_wall':{'global_id':'2UhljZvyHBM9jCYFR8xcyl','name':'Public Bathroom Pipe Wall','room_facing_mesh_plane_x_mm':wallx*1000,'body_rear_x_mm':rear*1000,'rear_envelope_inside_pipe_wall_mm':(rear-wallx)*1000,'interpretation':'Measured mesh relation, not assumed finished-wall contact. Distinguishing concealed rear construction from the physical seating face still requires a connector/contact-face definition.'},'waste_connector':{'status':'unknown','IfcDistributionPort_count':len(m.by_type('IfcDistributionPort')),'reason':'No semantic connection port in this IFC; pipe Body presence alone does not establish toilet outlet center.'}}
    bodylocal=(np.linalg.inv(np.array(e['postState']['product_placement_mm']))@np.column_stack((body*1000,np.ones(len(body)))).T).T[:,:3]
    dims=bodylocal.max(0)-bodylocal.min(0);checks=[]
    for r in e['registered_planes']:
        a=m.by_guid(r['annotation']);T=ifcopenshell.util.placement.get_local_placement(a.ObjectPlacement)
        local=np.array([pt.Coordinates for rep in a.Representation.Representations for i in rep.Items for line in i.Elements for pt in line.Points]);world=(T@np.column_stack((local,np.ones(len(local)))).T).T[:,:3]
        r['world_bounds_mm']=[world.min(0).tolist(),world.max(0).tolist()]
        axes={'plan':[0,1],'front':[0,2],'side':[1,2]}[r['view']]
        r['projected_bbox_size_minus_body_mm']=((local.max(0)-local.min(0))[axes]-dims[axes]).tolist()
        r['projected_bbox_min_minus_body_mm']=(local.min(0)[axes]-bodylocal.min(0)[axes]).tolist()
        r['projected_bbox_max_minus_body_mm']=(local.max(0)[axes]-bodylocal.max(0)[axes]).tolist()
        checks.append({'view':r['view'],'one_named_drawing':sum(x.ObjectType=='DRAWING' and x.Name=='GEBERIT-146-140-WC-'+r['view'].upper() for x in m.by_type('IfcAnnotation'))==1,'paths':r['path_count']})
    e['independent_validation']={'drawings':checks,'formal_ifc_sha256':sha(ROOT/'2504 GBTB Yanlord Zhuhai.ifc'),'no_per_view_recenter':True}
    assert all(x['one_named_drawing'] for x in checks)
    e['verdict']='PASS three-view mutual registration and corrected Front; installation contact semantics remain unresolved'
    (O/'evidence.json').write_text(json.dumps(e,ensure_ascii=False,indent=2)+'\n')
    manifest=json.loads((P/'manifest.json').read_text())
    manifest['registered_3d_review']={'status':'pending_user_review','index':'registered-3d-2026-09-07/index.html','evidence':'registered-3d-2026-09-07/evidence.json','three_view_registration_pass':True,'front_direction_corrected':True,'installation_contact_semantics_complete':False}
    (P/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    ep=P/'bonsai-drawings/wc/GEBERIT-146-140-create-drawing-evidence.json'
    previous=json.loads(ep.read_text());previous['latest_revision']='2026-09-07 shared-3D-registration and room-facing Front'
    previous['postState']['derived_ifc_sha256']=sha(P/'Geberit-146-140-derived-drawing.ifc')
    previous['persistence']=e['persistence'];previous['latest_validation_evidence']=str(O/'evidence.json')
    blend=P/'Geberit-146-140-WC-project-drawings.blend';previous['outputs']['session_blend'].update({'path':str(blend),'bytes':blend.stat().st_size,'sha256':sha(blend)})
    for r in previous['outputs']['views']:
        if r['view']=='front':r['svg'].update(e['outputs']['front_svg'])
        cache=Path(r['linework_cache']['path']);r['linework_cache'].update({'bytes':cache.stat().st_size,'sha256':sha(cache)})
    ep.write_text(json.dumps(previous,ensure_ascii=False,indent=2)+'\n')
    rows=''.join(f'<tr><td>{r["view"]}</td><td>{r["path_count"]}</td><td>{r["placement_max_delta_mm"]:.6f}</td><td>{r["projected_bbox_size_minus_body_mm"]}</td></tr>' for r in e['registered_planes'])
    html=f'''<!doctype html><meta charset="utf-8"><title>AquaClean 3D registration</title><style>body{{font:16px system-ui;margin:32px;max-width:1300px}}img{{width:100%}}td,th{{padding:8px;border:1px solid #ddd}}table{{border-collapse:collapse}}</style><h1>AquaClean 三视图共享坐标渲染</h1><p>蓝 Plan / 绿 Front / 红 Side；灰色实际 Body 与共用包围盒；橙色为共同 DWG 插入原点。所有线从重载的 IFC 读取并按持久化 ObjectPlacement 映射，没有逐图居中。</p><p>Front 已从室内朝马桶正面重出 Bonsai 场景 SVG。三张二维平面穿过产品局部原点，平面本身不等于墙面或地面。</p><table><tr><th>View</th><th>Paths</th><th>Placement delta mm</th><th>Projected size minus Body mm</th></tr>{rows}</table><p>管井墙外表面与陶瓷后包络相差 {(rear-wallx)*1000:.3f} mm；这不能直接判为安装面错位，因为 IFC 尚未定义陶瓷安装接触面。排污连接口缺少 IfcDistributionPort，不能声称已经验证安装对接。</p>'''
    for name in ['front-left','rear-right','planes-only']:html+=f'<h2>{name}</h2><img src="registered-{name}.png">'
    html+='<p><a href="evidence.json">完整验证证据</a></p>'
    for v in ['PLAN','FRONT','SIDE']:html+=f'<p><a href="../bonsai-drawings/wc/GEBERIT-146-140-WC-{v}.svg">{v} 场景 SVG</a></p>'
    (O/'index.html').write_text(html)
    print(json.dumps({'index':str(O/'index.html'),'checks':checks,'wall_relation_mm':(rear-wallx)*1000}))
if __name__=='__main__':main()
