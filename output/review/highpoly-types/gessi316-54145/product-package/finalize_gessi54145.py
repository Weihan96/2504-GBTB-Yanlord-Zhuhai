"""Finalize technical migration after explicit independent image inspection."""
import json
from pathlib import Path
import migrate_gessi54145 as task

def finish(log):
    r=json.loads(task.REPORT.read_text());assert r['verdict']=='validated_pending_visual_check'
    log=Path(log);ready=[line.split('TASK_BLENDER_IFC_READY ',1)[1] for line in log.read_text().splitlines() if 'TASK_BLENDER_IFC_READY ' in line][-1]
    session=json.loads(ready)
    assert session['ifc']==str(task.SINGLE) and session['sha256']==task.pkg.sha256(task.SINGLE)
    assert session['blend_saved'] is False and session['has_blend_warning'] is False
    assert session['scene_name']=='🔍 逐个审核高模图纸候选'
    r['fresh_pure_session']={**session,'log':str(log),'owner_task_id':'01a075a7-e948-7bd1-a27e-6ff542ea187c-gessi54145-migration','bridge_port':9892}
    r['visual']={'method':'Independent view_image of 3 SINGLE and 3 SCENE PNGs','observations':['Plan preserves horizontal wall arm and circular spray head.','Front preserves wall arm above joint and spray plate; original local approval orientation retained in single preview.','Side places circular wall anchor and joint above spray plate; world +Z projects upwards.','Scene references use only current product annotation; historical Drawing/grid symbols and spatial filling excluded in external recipe.'],'verdict':'pass_migration','scene_approval_status':'pending'}
    r['visual']['parent_agent_cross_check']='Parent independently viewed all three regenerated scene PNGs: historical markers and Plan black fill absent, wall context readable, blue lines complete, Side upright. Technical visual review only; not user scene approval.'
    r['verdict']='pass';r['validation_verdict']='pass';r['scene_approval_status']='pending'
    assert task.pkg.sha256(task.FORMAL)==task.FORMAL_HASH
    task.write(task.REPORT,r)
    for name in ['manifest.json','handoff.json']:
        path=task.OUT/name;data=json.loads(path.read_text());data.update(status='complete',validation_verdict='pass',single_product_approval_status='approved',scene_approval_status='pending',approval_label_zh='单品已通过、场景待验收')
        if name=='manifest.json':
            data['library_previews']['iso']=task.pkg.record(task.PRODUCT/'bonsai-camera-iso.png')
            data['preview_provenance']['3d']='Existing actual IFC Body camera render; both Body representations preserved'
            data['fresh_pure_session']=r['fresh_pure_session']
        else:
            data['ifc_sha256']=task.pkg.sha256(task.SINGLE)
            data['fresh_pure_session_pid']=session['pid']
        task.write(path,data)
    print('GESSI54145 migration complete; scene pending')

if __name__=='__main__':
    import sys
    finish(sys.argv[1])
