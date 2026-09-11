"""Record observed PNG review and the final fresh, task-owned IFC startup."""
import json,sys,subprocess
from pathlib import Path
import migrate_duofix as task
pkg=task.pkg

def main():
 log=Path(sys.argv[1]);ready=next(json.loads(line.split('TASK_BLENDER_IFC_READY ',1)[1]) for line in log.read_text().splitlines() if line.startswith('TASK_BLENDER_IFC_READY '))
 assert ready['ifc']==str(task.SINGLE) and ready['sha256']==pkg.sha256(task.SINGLE)
 assert not ready['blend_saved'] and not ready['has_blend_warning']
 report=json.loads(task.REPORT.read_text());assert report['verdict']=='validated_pending_visual_check'
 assert report['independent_validation']['schema_errors']==0
 report['visual']={'verdict':'pass','evidence_kind':'runtime-output',
  'observations':['Plan: one blue cistern/frame behind the WC; no duplicate black target Body projection.',
   'Front: full frame and support feet retained at the approved project height.',
   'Side: approved left-side direction, tank/frame/feet retained; original scene context unchanged.',
   'Single-product previews use only the approved blue geometry emitted from the pure IFC.'],
  'images':[pkg.record(task.OUT/f'DUOFIX-SCENE-{v}.png') for v in ['PLAN','FRONT','SIDE']]}
 report['final_ifc_acceptance_startup']={'log':pkg.record(log),'ready':ready,'prior_temporary_scene_discarded':True,'fresh_process':True}
 report['formal_sha256_after']=pkg.sha256(task.FORMAL);assert report['formal_sha256_after']==task.FORMAL_HASH
 staged=subprocess.check_output(['git','diff','--cached','--name-only','-z'],cwd=task.ROOT).split(b'\0')
 unstaged=subprocess.check_output(['git','diff','--name-only','-z'],cwd=task.ROOT).split(b'\0')
 report['git_split']={'protected_staged_files':sum(bool(s) for s in staged),
  'tracked_unstaged_files':[s.decode() for s in unstaged if s],
  'overlap':[s.decode() for s in unstaged if s and s in staged],
  'agent_new_files_directory':str(task.OUT),'agent_staged_anything':False}
 assert report['git_split']['protected_staged_files']==937
 report['verdict']='pass';task.write(task.REPORT,report)
 for name in ['manifest.json','handoff.json']:
  p=task.OUT/name;d=json.loads(p.read_text());d['status']='pass';d['validation_verdict']='pass'
  d['final_ifc_acceptance_startup']=report['final_ifc_acceptance_startup']
  if name=='manifest.json':d['persistent_ifc']=pkg.record(task.SINGLE)
  task.write(p,d)
 cleanup=task.OUT/'cleanup-proposal.json';d=json.loads(cleanup.read_text())
 previous=Path('/var/folders/rz/d8p6s4y50ws53rd150nm2s3r0000gn/T/duofix-pure-package-scene-gv45ethy')
 if previous.is_dir() and str(previous) not in d['temporary_directories']:d['temporary_directories'].append(str(previous))
 d['legacy_files'].append(pkg.record(task.PRODUCT/'Geberit-Duofix-Sigma-224-212-derived-drawing.ifc'))
 task.write(cleanup,d)
 print(json.dumps({'verdict':'pass','pure_ifc':pkg.record(task.SINGLE),'final_pid':ready['pid'],'git':report['git_split']}))

if __name__=='__main__':main()
