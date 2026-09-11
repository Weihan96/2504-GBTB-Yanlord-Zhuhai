"""Record completed technical checks without granting scene approval."""
import json
from pathlib import Path
import migrate_falper as task

report=json.loads(task.REPORT.read_text())
fresh=json.loads((task.OUT/'fresh-session.json').read_text())
assert report['verdict']=='validated_pending_visual_check'
assert fresh['pure_state']==report['pure_pre_bonsai']
assert fresh['scene_name']=='🔍 逐个审核高模图纸候选'
assert not fresh['blend_saved'] and not fresh['has_blend_warning']
assert report['all_formal_bodies_unchanged'] and report['all_formal_placements_unchanged']
assert report['formal_physical_elements']==714 and report['second_attachment_created_entities']==0
assert task.pkg.sha256(task.FORMAL)==task.FORMAL_HASH
assert all(task.pkg.sha256(r['path'])==r['sha256'] for r in report['protected_files'])
report.update({'status':'complete','verdict':'pass','validation_verdict':'pass','scene_approval_status':'pending',
    'visual':{'status':'pass','checked_by':['migration_agent','main_agent'],
              'checks':['Three PNGs show complete official blue product lines and readable wall/floor/adjacent component context.',
                        'No historical Drawing symbols, grid references, IfcSpace shapes or embedded raster images remain.',
                        'Technical review does not constitute user scene acceptance.'],
              'images':[task.pkg.record(task.OUT/f'FALPER-WFB-SCENE-{v.upper()}.png') for v in ('plan','front','side')]},
    'fresh_session':fresh,'formal_sha256_after':task.pkg.sha256(task.FORMAL),
    'worktree_split':{'protected_staged_baseline_modified':False,'new_unstaged_scope':str(task.OUT),'overlapping_tracked_files':[],
                      'stage_commit_stash_performed':False},
    'reference_filter_policy':json.loads(task.RECIPE.read_text())['reference_filter_policy']})
task.write(task.REPORT,report)
for filename in ('manifest.json','handoff.json'):
    path=task.OUT/filename
    data=json.loads(path.read_text())
    data.update({'status':'complete','validation_verdict':'pass','scene_approval_status':'pending','label_zh':'单品已通过、场景待验收'})
    if filename=='handoff.json':
        data.update({'fresh_session':'fresh-session.json','pid':fresh['pid'],'schema_errors':0,'formal_physical_elements':714,
                     'all_formal_bodies_and_placements_unchanged':True,'second_attachment_created_entities':0})
    else:
        data['packaged_source_evidence']=report['packaged_source_evidence']
        data['scene_recipe']=task.pkg.record(task.RECIPE)
    task.write(path,data)
cleanup_path=task.OUT/'cleanup-proposal.json'
cleanup=json.loads(cleanup_path.read_text())
temporary_parent=Path(report['temporary_project_directory']).parent
cleanup['temporary_directories']=sorted(str(p) for p in temporary_parent.glob('falper-pure-package-scene-*') if p.is_dir())
task.write(cleanup_path,cleanup)
print(json.dumps({'status':'complete','validation_verdict':'pass','scene_approval_status':'pending','pid':fresh['pid'],
                  'pure_ifc':task.pkg.record(task.SINGLE),'handoff':str(task.OUT/'handoff.json')},ensure_ascii=False))
