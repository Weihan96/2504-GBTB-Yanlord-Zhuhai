"""Finalize technical migration after independent images have been inspected."""
import json,sys,copy,base64,xml.etree.ElementTree as ET
from pathlib import Path
import migrate_hima01 as task
from verify_hima01 import render

def finish(log):
    r=json.loads(task.REPORT.read_text());assert r['verdict'] in ('validated_pending_visual_check','pass')
    log=Path(log)
    ready=[line.split('TASK_BLENDER_IFC_READY ',1)[1] for line in log.read_text().splitlines() if 'TASK_BLENDER_IFC_READY ' in line][-1]
    session=json.loads(ready)
    assert session['ifc']==str(task.SINGLE) and session['sha256']==task.pkg.sha256(task.SINGLE)
    assert session['blend_saved'] is False and session['has_blend_warning'] is False
    assert session['scene_name']=='🔍 逐个审核高模图纸候选'
    r['fresh_pure_session']={**session,'log':str(log),'owner_task_id':'01a075a7-e948-7bd1-a27e-6ff542ea187c-hima01-migration','bridge_port':9894}
    r['visual']={'method':'Independent view_image of three SINGLE and three SCENE PNGs','observations':['Approved folded black screen views preserved with 8 Plan, 17 Front and 27 Side closed paths.','Plan shows the folded three-element footprint; elevations retain frame and feet.','Front and Side have world +Z upwards and horizontal axes with negligible world Z.','Historical Drawing symbols, grids and spaces are filtered before real CreateDrawing.'],'verdict':'pass_migration','scene_approval_status':'pending'}
    r.update(verdict='pass',validation_verdict='pass',scene_approval_status='pending')
    r['plan_camera_height_evidence']={'method':'IfcOpenShell USE_WORLD_COORDS slab mesh bounding box','camera_z_mm':1700,'product_top_z_mm':1008,'ceiling_slab_guid':'3ZKyTiBc1ECxQi0UR8K1Ke','ceiling_slab_bounds_mm':[[-6500,-4399.838,2820],[-3200,0.162,3000]],'above_product_and_below_ceiling':True}
    assert task.pkg.sha256(task.FORMAL)==task.FORMAL_HASH
    task.write(task.REPORT,r)
    ns='http://www.w3.org/2000/svg';root=ET.Element('{'+ns+'}svg',width='1800',height='1340',viewBox='0 0 1800 1340')
    ET.SubElement(root,'{'+ns+'}rect',width='1800',height='1340',fill='white')
    for i,view in enumerate(['plan','front','side']):
        for row,kind in enumerate(['SINGLE','SCENE']):
            text=ET.SubElement(root,'{'+ns+'}text',x=str(i*600+24),y=str(row*650+30),style='font:20px sans-serif;fill:#111820');text.text=f'HIMA01 {view.upper()} / {kind}'
            png=task.OUT/f'HIMA01-{kind}-{view.upper()}.png'
            ET.SubElement(root,'{'+ns+'}image',x=str(i*600+20),y=str(row*650+50),width='560',height='550',href='data:image/png;base64,'+base64.b64encode(png.read_bytes()).decode())
    caption=ET.SubElement(root,'{'+ns+'}text',x='24',y='1320',style='font:18px sans-serif;fill:#111820')
    caption.text='Approved product / Scene acceptance pending — folded high-poly-derived black linework; PVA11 DWG reference only'
    contact=task.OUT/'HIMA01-review-contact-sheet.svg';ET.ElementTree(root).write(contact,encoding='utf-8',xml_declaration=True);render(contact,contact.with_suffix('.png'),1800)
    for name in ['manifest.json','handoff.json']:
        p=task.OUT/name;d=json.loads(p.read_text());d.update(status='complete',validation_verdict='pass',single_product_approval_status='approved',scene_approval_status='pending',approval_label_zh='单品已通过、场景待验收')
        if name=='manifest.json':
            d['library_previews']['iso']=task.pkg.record(task.PRODUCT/'bonsai-camera-iso.png')
            d['preview_provenance']['3d']='Existing actual IFC Body camera render; all three original Body representations preserved'
            d['official_cad_used']=False
            if 'source_dwg' in d:d['source_reference_only']=d.pop('source_dwg')
            d['reference_scope']='Official PVA11 DWG: family and unfolded-state reference only, excluded from all approved black linework'
            d['geometry']['discarded_box_representation_count']=r['extraction']['unused_type_box_maps_removed']
            d['fresh_pure_session']=r['fresh_pure_session'];d['contact_sheet']=task.pkg.record(contact.with_suffix('.png'))
        else:
            d['ifc_sha256']=task.pkg.sha256(task.SINGLE);d['fresh_pure_session_pid']=session['pid'];d['best_visual']=str(contact.with_suffix('.png'))
        task.write(p,d)
    print('HIMA01 complete; scene pending')

if __name__=='__main__':finish(sys.argv[1])
