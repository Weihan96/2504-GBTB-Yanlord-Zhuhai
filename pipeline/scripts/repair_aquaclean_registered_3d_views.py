"""Task-owned short-lived Bonsai worker: camera repair and true-world line rendering."""
import sys, json, hashlib, traceback, contextlib
from pathlib import Path
import bpy, addon_utils, numpy as np
import ifcopenshell, ifcopenshell.util.placement, ifcopenshell.util.element
import ifcopenshell.api.geometry
import ifcopenshell.api.root
from mathutils import Matrix, Vector

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'pipeline/scripts'))
import types
# Legacy generator invokes main at module import; load definitions only.
wc=types.ModuleType('aquaclean_drawing_helpers')
wc.__file__=str(ROOT/'pipeline/scripts/create_geberit_146_140_wc_drawings.py')
exec(compile(Path(wc.__file__).read_text().split('\ntry:\n    main()')[0],wc.__file__,'exec'),wc.__dict__)
from bonsai import tool
OUT=wc.PRODUCT_DIR/'registered-3d-2026-09-07'
OUT.mkdir(exist_ok=True)
def record(p):
    return {'path':str(p),'sha256':wc.sha256(p),'bytes':p.stat().st_size}
def lines(annotation):
    return [np.array([p.Coordinates for p in l.Points]) for r in annotation.Representation.Representations for i in r.Items if i.is_a('IfcGeometricCurveSet') for l in i.Elements]
def material(name,color,alpha=1):
    m=bpy.data.materials.new(name); m.diffuse_color=(*color,alpha);m.use_nodes=True
    p=m.node_tree.nodes.get('Principled BSDF');p.inputs['Base Color'].default_value=(*color,1);p.inputs['Roughness'].default_value=.65;p.inputs['Alpha'].default_value=alpha
    m.surface_render_method='DITHERED'
    return m
def curve(scene,name,paths,mat,radius=.001):
    c=bpy.data.curves.new(name,'CURVE');c.dimensions='3D';c.bevel_depth=radius;c.bevel_resolution=1
    for path in paths:
        s=c.splines.new('POLY');s.points.add(len(path)-1)
        for p,co in zip(s.points,path):p.co=(*co,1)
    o=bpy.data.objects.new(name,c);scene.collection.objects.link(o);c.materials.append(mat);return o
def main():
    assert wc.sha256(wc.FORMAL_IFC)==wc.FORMAL_SHA256
    before=record(wc.DERIVED_IFC)
    frozen={str(wc.PRODUCT_DIR/f'{v}.svg'):wc.sha256(wc.PRODUCT_DIR/f'{v}.svg') for v in ['plan','front','side']}
    addon_utils.enable('bonsai_bridge',default_set=False,persistent=False)
    import bonsai_bridge
    with contextlib.suppress(Exception):addon_utils.disable('bl_ext.user_default.project_control',default_set=False,handle_error=None)
    assert bpy.ops.bim.load_project(filepath=str(wc.DERIVED_IFC),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
    model=tool.Ifc.get();target=model.by_guid(wc.TARGET_GLOBAL_ID);body=tool.Ifc.get_object(target)
    assert body
    T=ifcopenshell.util.placement.get_local_placement(target.ObjectPlacement)
    evidence=json.loads(wc.EVIDENCE.read_text());old=json.loads((wc.PRODUCT_DIR/'origin-audit-2026-09-05/origin-audit.json').read_text())
    views=evidence['outputs']['views'];removed=[]
    for rec in views:
        v=rec['view']
        for role,name in [('drawing','GEBERIT-146-140-WC-'+v.upper()),('annotation','Geberit 146.140.11.1 approved official linework / '+v)]:
            matches=sorted([e for e in model.by_type('IfcAnnotation') if e.Name==name],key=lambda e:e.id())
            keep=matches[0]
            if role=='drawing':rec['drawing']['global_id']=keep.GlobalId
            else:rec['official_annotation_global_id']=keep.GlobalId
            for e in matches[1:]:
                removed.append({'global_id':e.GlobalId,'name':e.Name})
                obj=tool.Ifc.get_object(e)
                if obj:bpy.data.objects.remove(obj,do_unlink=True)
                ifcopenshell.api.root.remove_product(model,product=e)
    r=next(v for v in views if v['view']=='front')
    drawing=model.by_guid(r['drawing']['global_id']);cam=tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
    before_cam=[list(row) for row in cam.matrix_world]
    bbox=wc.shared.world_bbox(body);center=Vector([(a+b)/2 for a,b in zip(*bbox)])
    # Product protrudes toward world -X; camera is on that room side, looking +X.
    M=Matrix(((0,0,-1,center.x-1.8),(-1,0,0,center.y),(0,1,0,cam.location.z),(0,0,0,1)))
    cam.matrix_world=M;cam.data.clip_end=4.5;bpy.context.view_layer.update()
    ifcopenshell.api.geometry.edit_object_placement(model,product=drawing,matrix=np.array(M),is_si=True,should_transform_children=False)
    override=wc.view3d_override()
    with bpy.context.temp_override(**override):assert bpy.ops.bim.activate_drawing(drawing=drawing.id(),should_view_from_camera=False)=={'FINISHED'}
    dp=tool.Drawing.get_document_props();dp.should_use_underlay_cache=False;dp.should_use_linework_cache=False;dp.should_use_annotation_cache=False
    with bpy.context.temp_override(**override):result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
    assert result=={'FINISHED'}
    svg=Path(r['svg']['path']);style=wc.style_svg(svg,r['official_annotation_global_id']);svg_check=wc.inspect_svg(svg,r['official_annotation_global_id']);assert svg_check['no_body_annotation_or_instance_duplicate']
    for rec in views:
        if rec['view']=='front':continue
        d=model.by_guid(rec['drawing']['global_id']);tool.Ifc.get_object(d) or tool.Drawing.import_drawing(d)
        with bpy.context.temp_override(**override):
            bpy.ops.bim.activate_drawing(drawing=d.id(),should_view_from_camera=False)
            assert bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)=={'FINISHED'}
        wc.style_svg(Path(rec['svg']['path']),rec['official_annotation_global_id'])
    persistence=bonsai_bridge._h_save_ifc_file({'output_path':str(wc.DERIVED_IFC),'overwrite':True,'reload':True})
    model=tool.Ifc.get();body=tool.Ifc.get_object(model.by_guid(wc.TARGET_GLOBAL_ID))
    actual=ifcopenshell.util.placement.get_local_placement(model.by_guid(drawing.GlobalId).ObjectPlacement)
    assert np.allclose(actual[:3,:3],np.array(M)[:3,:3],atol=1e-7)
    assert np.allclose(actual[:3,3],np.array(M)[:3,3]*1000,atol=.001)
    bpy.ops.wm.save_as_mainfile(filepath=str(wc.SESSION_BLEND),check_existing=False)
    refs=[]
    for elem,obj,bbmin,bbmax in wc.shared.room_elements(((center.x-1,center.y-1,-.3),(center.x+1,center.y+1,1.1))):
        if elem.is_a('IfcWall') or elem.is_a('IfcCovering') or elem.is_a('IfcSlab') or elem.is_a('IfcPipeSegment') or elem.is_a('IfcPipeFitting'):
            refs.append({'global_id':elem.GlobalId,'class':elem.is_a(),'name':elem.Name,'bbox_m':[bbmin,bbmax]})
    # Separate render Scene: exact persisted line coordinates transformed by placement, with no per-view recenter.
    scene=bpy.data.scenes.new('AquaClean shared WORLD coordinates');bpy.context.window.scene=scene
    colors={'plan':(.04,.35,.95),'front':(.02,.58,.23),'side':(.86,.14,.04)}
    mats={v:material(v,c) for v,c in colors.items()};gray=material('Shared Body bounding box',(.24,.28,.33));bodmat=material('Actual Body translucent',(.56,.61,.68),.13)
    bodycopy=bpy.data.objects.new('Actual persisted IFC Body',body.data.copy());scene.collection.objects.link(bodycopy);bodycopy.matrix_world=body.matrix_world.copy();bodycopy.data.materials.clear();bodycopy.data.materials.append(bodmat)
    registered=[];allpts=[]
    for rec in views:
        v=rec['view'];ann=model.by_guid(rec['official_annotation_global_id']);A=ifcopenshell.util.placement.get_local_placement(ann.ObjectPlacement);paths=lines(ann)
        assert np.max(np.abs(A-T))<1e-7
        world=[((A@np.column_stack((p,np.ones(len(p)))).T).T[:,:3]/1000) for p in paths]
        curve(scene,v.upper()+' persisted IFC LINEWORK',world,mats[v]);allpts.extend(world)
        flat=np.vstack(paths);axes={'plan':[0,1],'front':[0,2],'side':[1,2]}[v];fixed=({0,1,2}-set(axes)).pop()
        lo=flat.min(axis=0);hi=flat.max(axis=0);p=[]
        for a,b in [(0,0),(1,0),(1,1),(0,1),(0,0)]:
            q=lo.copy();q[axes[0]]=[lo[axes[0]],hi[axes[0]]][a];q[axes[1]]=[lo[axes[1]],hi[axes[1]]][b];p.append((A@np.append(q,1))[:3]/1000)
        curve(scene,v.upper()+' true plane bounds',[p],mats[v],.00055)
        registered.append({'view':v,'annotation':ann.GlobalId,'path_count':len(paths),'local_plane_axis':fixed,'local_plane_coordinate_mm':float(lo[fixed]),'local_plane_spread_mm':float(hi[fixed]-lo[fixed]),'placement_max_delta_mm':float(np.max(np.abs(A-T))),'local_min_mm':lo.tolist(),'local_max_mm':hi.tolist(),'world_matrix_mm':A.tolist(),'world_bounds_mm':[np.vstack(world).min(axis=0).tolist(),np.vstack(world).max(axis=0).tolist()]})
    corners=[Vector((x,y,z)) for x in [bbox[0][0],bbox[1][0]] for y in [bbox[0][1],bbox[1][1]] for z in [bbox[0][2],bbox[1][2]]]
    curve(scene,'Shared actual Body bbox',[[a,b] for i,a in enumerate(corners) for j,b in enumerate(corners) if j>i and sum(abs(a[k]-b[k])>1e-6 for k in range(3))==1],gray,.0007)
    origin=T[:3,3]/1000
    for k,v in enumerate(['plan','front','side']):
        end=origin.copy();end[k]+=.16;curve(scene,'WORLD '+ 'XYZ'[k],[[origin,end]],mats[v],.0014)
    # Native common insertion point is mapped once, not used to re-center any drawing.
    native=np.array(old['views'][0]['anchors']['O']['world_xyz_mm'])/1000
    orange=material('Native shared insertion',(.95,.55,.02))
    curve(scene,'Common insertion marker',[[native-np.eye(3)[i]*.012,native+np.eye(3)[i]*.012] for i in range(3)],orange,.0014)
    scene.render.engine='CYCLES';scene.cycles.samples=24;scene.cycles.use_denoising=True
    scene.render.resolution_x=1800;scene.render.resolution_y=1500;scene.render.resolution_percentage=100
    scene.world=bpy.data.worlds.new('White World');scene.world.use_nodes=True;scene.world.node_tree.nodes['Background'].inputs[0].default_value=(.8,.8,.8,1)
    scene.view_settings.view_transform='Standard'
    camera=bpy.data.objects.new('Registration evidence camera',bpy.data.cameras.new('Registration evidence camera'));scene.collection.objects.link(camera);scene.camera=camera;camera.data.type='ORTHO';camera.data.ortho_scale=1.03
    for i,(label,mat) in enumerate([('AQUACLEAN 146.140 | persisted IFC planes',gray),('BLUE  Plan XY   |   GREEN  Front XZ   |   RED  Side YZ',gray),('Grey: actual Body + shared bounding box',gray),('Orange cross: common native DWG insertion',gray)]):
        font=bpy.data.curves.new('legend','FONT');font.body=label;font.size=.018 if i==0 else .014
        txt=bpy.data.objects.new('legend',font);scene.collection.objects.link(txt);txt.parent=camera;txt.location=(-.47,.40-i*.026,-2);font.materials.append(mat)
    renders=[]
    for name,offset,show_body in [('front-left',(-1.3,-1.0,.9),True),('rear-right',(1.25,1.1,.8),True),('planes-only',(-1.3,-1.0,.9),False)]:
        bodycopy.hide_render=not show_body;camera.location=center+Vector(offset);camera.rotation_euler=(center-camera.location).to_track_quat('-Z','Y').to_euler()
        scene.render.filepath=str(OUT/f'registered-{name}.png');bpy.ops.render.render(write_still=True);renders.append(record(Path(scene.render.filepath)))
    bpy.ops.wm.save_as_mainfile(filepath=str(OUT/'AquaClean-three-registered-planes.blend'),check_existing=False)
    assert all(wc.sha256(Path(p))==h for p,h in frozen.items())
    assert wc.sha256(wc.FORMAL_IFC)==wc.FORMAL_SHA256
    report={'task':'AquaClean true-world Plan Front Side registration','versions':{'blender':bpy.app.version_string,'ifcopenshell':ifcopenshell.version,'bridge':bonsai_bridge.bl_info.get('version')},'courseEvidence':{'mode':'embedded-course-index','lesson':'085000','timestamps':['01:03','01:59','02:13']},'preState':{'ifc':before,'front_camera_m':before_cam},'execution':{'removed_duplicate_drawings_from_legacy_import':removed,'create_drawing':sorted(result),'adapter':'isolated short-lived Bonsai worker; installed public bridge save handler','render':'Cycles actual 3D curves read from reloaded IFC; no per-view centering'},'persistence':persistence,'postState':{'ifc':record(wc.DERIVED_IFC),'front_camera_m':[list(x) for x in M],'front_looks_world':[1,0,0],'product_placement_mm':T.tolist()},'registered_planes':registered,'body_bbox_m':bbox,'native_shared_insertion_world_mm':(native*1000).tolist(),'nearby_installation_candidates':refs,'installation_semantics':{'finished_floor':'pending geometric/semantic identification','finished_wall':'pending geometric/semantic identification','waste_connector':'unknown; not asserted from drawing origin'},'outputs':{'renders':renders,'front_svg':{**record(svg),**svg_check,**style},'blend':record(OUT/'AquaClean-three-registered-planes.blend')},'formal_ifc_sha256':wc.FORMAL_SHA256,'verdict':'registered geometry and Front repair verified; installation semantics pending'}
    for rec in views:
        rec['svg'].update(record(Path(rec['svg']['path'])))
        if rec['view']=='front':rec['camera']['matrix_world']=[list(x) for x in M]
    evidence['outputs']['views']=views
    wc.EVIDENCE.write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n')
    (OUT/'evidence.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,default=str)+'\n')
    print('AQUACLEAN_REGISTERED_RENDER_READY',str(OUT),flush=True)
if __name__=='__main__':
    try:main()
    except Exception:
        (OUT/'error.log').write_text(traceback.format_exc());traceback.print_exc()
    finally:bpy.app.timers.register(lambda:bpy.ops.wm.quit_blender() and None,first_interval=1)
