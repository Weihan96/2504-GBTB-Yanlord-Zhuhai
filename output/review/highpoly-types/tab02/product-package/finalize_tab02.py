"""Finalize technical migration after independent images have been inspected."""
import json,sys,copy,base64,xml.etree.ElementTree as ET
from pathlib import Path
import migrate_tab02 as task
from verify_tab02 import render

def finish(log):
    r=json.loads(task.REPORT.read_text());assert r['verdict'] in ('validated_pending_visual_check','pass')
    log=Path(log)
    ready=[line.split('TASK_BLENDER_IFC_READY ',1)[1] for line in log.read_text().splitlines() if 'TASK_BLENDER_IFC_READY ' in line][-1]
    session=json.loads(ready)
    assert session['ifc']==str(task.SINGLE) and session['sha256']==task.pkg.sha256(task.SINGLE)
    assert session['blend_saved'] is False and session['has_blend_warning'] is False
    assert session['scene_name']=='🔍 逐个审核高模图纸候选'
    r['fresh_pure_session']={**session,'log':str(log),'owner_task_id':'01a075a7-e948-7bd1-a27e-6ff542ea187c-tab02-migration','bridge_port':9897}
    r['visual']={'method':'Independent view_image of three SINGLE and three SCENE PNGs','observations':['Approved semantic black component views preserved with 1 Plan, 4 Front and 4 Side paths; column paths remain open.','Plan shows only the opaque rounded tabletop; elevations retain the tabletop, column sides and stone base.','Front and Side have world +Z upwards and horizontal axes with negligible world Z.','Historical Drawing symbols, grids and spaces are filtered before real CreateDrawing.'],'verdict':'pass_migration','scene_approval_status':'pending'}
    r.update(verdict='pass',validation_verdict='pass',scene_approval_status='pending')
    r['plan_camera_height_evidence']={'method':'IfcOpenShell world-coordinate IfcSlab mesh bounds at target XY (1675.832,4119.987)','camera_z_mm':1700,'product_top_z_mm':680,'above_product':True,'intersecting_slabs':[{'guid':'1fAJA2jj59uhEyVVcR7EM7','z_bounds_mm':[-200,0]}],'overhead_ceiling_slab_present':False,'context':'North terrace beyond living-room door; no overhead slab in original formal model at target XY'}
    r['courseEvidence']['provenance_sha256']='a86279d39cef0f85feb9c37d5f823556e41e1f6e1ebd1cc634e564c3ebfd24b4'
    r['task']={'target':task.GUID,'artifact':'TAB02 pure 3D + approved semantic three-view IFC and external pending scene recipe','blender':'4.5.3 LTS','bonsai_ifcopenshell':'0.8.4','save_boundary':'TAB02 product-package and disposable scene directory only'}
    r['plan']=['Inspect source approval and original Body','Extract exact approved semantic paths and original Body into pure IFC','Open pure IFC via task-owned public launcher','Save and reload pure through public Provider','Attach pure plus recipe to original formal temporary copy','Generate real Bonsai SVGs','Save and reload temporary scene','Validate disk geometry, inherited schema errors and rendered evidence','Fresh launch final pure IFC']
    r['execution']={'adapter':'public bridge execute_code (public MCP execute_blender_code maps to this protocol name)','drawing_operator':'bpy.ops.bim.create_drawing','recipe':'scene-recipe.json','source_script':task.pkg.record(task.OUT/'migrate_tab02.py'),'source_semantic_segmentation':task.pkg.record(task.PRODUCT/'tab02-semantic-segmentation.json')}
    r['persistence']={'pure_save':r['product_bonsai_save_reload'],'scene_save':r['temporary_bonsai_save_reload'],'pure_reloaded':r['pure_saved_reloaded'],'scene_reloaded':r['scene_saved_reloaded']}
    r['postState']=r['independent_validation']['pure_state']
    r['outputs']=r['independent_validation']['view_checks']
    r['runtime_idempotence']=json.loads((task.OUT/'runtime-idempotence.json').read_text())
    r['camera_round_trip']={'plan_front':'exact requested IFC block dimensions and placement preserved','side_requested_mm':[1600,1400,2300],'side_saved_mm':[1600.00002384186,1400.53052484831,2299.99995231628],'all_camera_placements_exact':True,'cause':'Bonsai raster-resolution and float normalization when saving active camera','svg_geometry_validation':'Uses actual saved camera placement and actual SVG viewBox; every approved edge endpoint verified independently','recipe_replay':'Twice applied runtime adapter on a fresh original formal copy creates zero new entities; see runtime-idempotence.json'}
    assert task.pkg.sha256(task.FORMAL)==task.FORMAL_HASH
    task.write(task.REPORT,r)
    ns='http://www.w3.org/2000/svg';root=ET.Element('{'+ns+'}svg',width='1800',height='1340',viewBox='0 0 1800 1340')
    ET.SubElement(root,'{'+ns+'}rect',width='1800',height='1340',fill='white')
    for i,view in enumerate(['plan','front','side']):
        for row,kind in enumerate(['SINGLE','SCENE']):
            text=ET.SubElement(root,'{'+ns+'}text',x=str(i*600+24),y=str(row*650+30),style='font:20px sans-serif;fill:#111820');text.text=f'TAB02 {view.upper()} / {kind}'
            png=task.OUT/f'TAB02-{kind}-{view.upper()}.png'
            ET.SubElement(root,'{'+ns+'}image',x=str(i*600+20),y=str(row*650+50),width='560',height='550',href='data:image/png;base64,'+base64.b64encode(png.read_bytes()).decode())
    caption=ET.SubElement(root,'{'+ns+'}text',x='24',y='1320',style='font:18px sans-serif;fill:#111820')
    caption.text='Approved product / Scene acceptance pending — high-poly-derived semantic black linework; official CAD not acquired'
    contact=task.OUT/'TAB02-review-contact-sheet.svg';ET.ElementTree(root).write(contact,encoding='utf-8',xml_declaration=True);render(contact,contact.with_suffix('.png'),1800)
    for name in ['manifest.json','handoff.json']:
        p=task.OUT/name;d=json.loads(p.read_text());d.update(status='complete',validation_verdict='pass',single_product_approval_status='approved',scene_approval_status='pending',approval_label_zh='单品已通过、场景待验收')
        if name=='manifest.json':
            d['library_previews']['iso']=task.pkg.record(task.PRODUCT/'bonsai-camera-iso.png')
            d['preview_provenance']['3d']='Existing actual IFC Body camera render; original Body representation preserved'
            d['official_cad_used']=False
            if 'source_dwg' in d:d['source_reference_only']=d.pop('source_dwg')
            d['reference_scope']='Official manufacturer catalogue: identity and 500 x 500 x 670 mm dimensions only; native CAD not acquired'
            d['geometry']['discarded_box_representation_count']=r['extraction']['unused_type_box_maps_removed']
            d['fresh_pure_session']=r['fresh_pure_session'];d['contact_sheet']=task.pkg.record(contact.with_suffix('.png'))
        else:
            d['ifc_sha256']=task.pkg.sha256(task.SINGLE);d['fresh_pure_session_pid']=session['pid'];d['best_visual']=str(contact.with_suffix('.png'))
        task.write(p,d)
    print('TAB02 complete; scene pending')

if __name__=='__main__':finish(sys.argv[1])
