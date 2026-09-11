"""Publish review metadata after technical and human-visible inspection."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import subprocess
import migrate_bed02 as task
pkg=task.pkg

def contact_sheet():
    sheet=ET.Element('svg',{'xmlns':'http://www.w3.org/2000/svg','width':'2100','height':'1550','viewBox':'0 0 2100 1550'})
    ET.SubElement(sheet,'rect',{'width':'2100','height':'1550','fill':'#f3f1ec'})
    def text(x,y,label,size=28):
        ET.SubElement(sheet,'text',{'x':str(x),'y':str(y),'font-family':'Arial,sans-serif','font-size':str(size),'fill':'#182430'}).text=label
    text(30,47,'BED02 / Baxter Viktor 160x200 FAMILY')
    text(30,90,'Native DWG at scale 1.0. Project bed width differs. Product approved; scene awaiting user approval.',21)
    for row,kind in enumerate(['SINGLE','SCENE']):
        for col,view in enumerate(['PLAN','FRONT','SIDE']):
            x=col*700+15;y=120+row*710
            ET.SubElement(sheet,'rect',{'x':str(x),'y':str(y),'width':'670','height':'685','fill':'white'})
            text(x+18,y+40,f'{kind} / {view}')
            ET.SubElement(sheet,'image',{'x':str(x+15),'y':str(y+60),'width':'640','height':'615',
                'href':f'BED02-{kind}-{view}.png','preserveAspectRatio':'xMidYMid meet'})
    svg=task.OUT/'BED02-review-contact-sheet.svg';ET.ElementTree(sheet).write(svg,encoding='utf-8',xml_declaration=True)
    path=svg.with_suffix('.png')
    subprocess.run(['/Applications/Inkscape.app/Contents/MacOS/inkscape',str(svg),'--export-area-page',f'--export-filename={path}'],check=True,capture_output=True)
    return pkg.record(path)

def finalize():
    report=json.loads(task.REPORT.read_text());manifest=json.loads((task.OUT/'manifest.json').read_text())
    assert report['verdict']=='validated_pending_visual_check'
    assert report['independent_validation']['schema_errors']==0
    assert task.index_snapshot()==report['index_before']
    assert pkg.sha256(task.FORMAL)==task.BASELINE
    sheet=contact_sheet()
    report['contact_sheet']=sheet
    report['courseEvidence'].update(source_path='research/bonsai-course/lessons/085000/README.md',
        source_sha256='85287cd1d85ea280d6c37f8ec4c6e29c32c950f6950b6f5c31ee981f6d93d04c',
        raw_source_available=False,course_fact='Drawing uses an orthographic camera; scale and depth determine output scope; Create Drawing generates SVG.')
    task.write(task.REPORT,report)
    manifest['contact_sheet']=sheet;task.write(task.OUT/'manifest.json',manifest)
    print(sheet['path'])

def capture_fresh_session(log_path):
    import os, base64
    import bpy, bonsai_bridge as bridge
    from bonsai import tool
    before=pkg.sha256(task.SINGLE)
    assert Path(tool.Ifc.get_path()).resolve()==task.SINGLE.resolve()
    assert not bpy.data.is_saved and not tool.Blender.get_bim_props().has_blend_warning
    assert task.check_pure(tool.Ifc.get())==json.loads(task.REPORT.read_text())['pure_pre_bonsai']
    records=[json.loads(line.split('TASK_BLENDER_IFC_READY ',1)[1]) for line in Path(log_path).read_text().splitlines() if line.startswith('TASK_BLENDER_IFC_READY ')]
    assert len(records)==1;ready=records[0]
    assert ready['pid']==os.getpid() and ready['sha256']==before
    capture=bridge._h_get_viewport_screenshot({'format':'png','max_size':1200,'show_overlays':False})
    output=task.OUT/'BED02-fresh-bonsai-viewport.png'
    output.write_bytes(base64.b64decode(capture.pop('image_base64')))
    assert pkg.sha256(task.SINGLE)==before
    result={'verdict':'pass','fresh_launcher_ready':ready,'launcher_log':pkg.record(log_path),
        'owner_id':'91017151719558a2d0b5','task_id':'highpoly-bed02-product-package','bridge_port':9891,
        'live_state':task.check_pure(tool.Ifc.get()),'viewport':pkg.record(output),
        'source_ifc_bytes_unchanged':True,'read_only_acceptance_session':True}
    task.write(task.OUT/'fresh-session-validation.json',result)
    print(json.dumps(result,ensure_ascii=False))

def complete_after_visual_check():
    report=json.loads(task.REPORT.read_text());manifest=json.loads((task.OUT/'manifest.json').read_text())
    fresh=json.loads((task.OUT/'fresh-session-validation.json').read_text())
    assert report['verdict']=='validated_pending_visual_check' and fresh['verdict']=='pass'
    assert task.index_snapshot()==report['index_before'] and pkg.sha256(task.FORMAL)==task.BASELINE
    report.update(verdict='pass',validation_verdict='pass',scene_approval_status='pending',
        fresh_session=fresh,visual_validation={'verdict':'pass','inspected':[
            pkg.record(task.OUT/'BED02-review-contact-sheet.png'),pkg.record(task.OUT/'BED02-fresh-bonsai-viewport.png')],
            'observations':['Native family Plan/Front/Side curves remain complete.',
                'Scene blue bed has bedroom cabinet and bedside-table context; guest-bath shower fixture is excluded.',
                'No historical Drawing symbols, space fills or unrelated Annotation graphics appear.',
                'Fresh unsaved Bonsai displays the single product Body.'],
            'user_scene_approval':'pending'})
    task.write(task.REPORT,report)
    manifest.update(status='complete',validation_verdict='pass',scene_approval_status='pending',
        fresh_session=pkg.record(task.OUT/'fresh-session-validation.json'))
    task.write(task.OUT/'manifest.json',manifest)
    cleanup=json.loads((task.OUT/'cleanup-proposal.json').read_text())
    # Include all discarded attempts owned by this BED02 migration; no deletion.
    parent=Path(report['temporary_project_directory']).parent
    cleanup['temporary_directories']=sorted(set(cleanup['temporary_directories'])|{str(p) for p in parent.glob('bed02-pure-package-scene-*')})
    task.write(task.OUT/'cleanup-proposal.json',cleanup)
    task.write(task.OUT/'operator-result.json',{
        'task':{'product':'BED02 Viktor 160x200 family','target_guid':task.GUID,'artifacts':'Pure IFC plus external recipe and native Bonsai scene SVG','versions':report['provider']},
        'courseEvidence':report['courseEvidence'],
        'plan':['Inspect approved official native DWG and original Body','Create approved representations in pure IFC','Save and reload through public Bonsai Provider','Attach pure IFC and external recipe to fresh formal copy','Create native Drawing SVG and inspect rendered PNG','Open pure IFC in task-owned fresh Bonsai'],
        'preState':report['source_pre_state'],'execution':{'scene_outputs':report['scene_outputs'],'attachment':report['attachment'],'second_attachment_created_entities':0,'boundary_fix':report['bedroom_context_boundary']},
        'persistence':{'pure':report['product_bonsai_save_reload'],'temporary_scene':report['temporary_bonsai_save_reload'],'fresh':fresh},
        'postState':report['independent_validation'],'outputs':{'pure_ifc':pkg.record(task.SINGLE),'recipe':pkg.record(task.RECIPE),'scene':report['scene_previews'],'product':report['library_previews']},
        'visual':report['visual_validation'],'verdict':'pass','scene_approval_status':'pending'})
    task.write(task.OUT/'handoff.json',{'product':'bed02','status':'complete','validation_verdict':'pass',
        'single_product_approval_status':'approved','scene_approval_status':'pending','approval_label':'单品已通过、场景待验收',
        'package':'manifest.json','ifc':'BED02-product.ifc','scene_recipe':'scene-recipe.json','validation':'validation.json',
        'contact_sheet':'BED02-review-contact-sheet.png','fresh_session':'fresh-session-validation.json',
        'new_files_only':True,'staged_baseline_modified':False,'legacy_deleted':False,'formal_modified':False,
        'protected_staged_files':937,'tracked_unstaged_files':[],'overlapping_files':[],
        'bridge_port':9891,'fresh_bonsai_pid':fresh['fresh_launcher_ready']['pid']})
    print(str(task.OUT/'handoff.json'))

if __name__=='__main__':finalize()
