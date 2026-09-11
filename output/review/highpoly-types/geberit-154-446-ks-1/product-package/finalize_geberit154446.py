"""Record independent image review and verified fresh pure-IFC session."""
import json,sys
from pathlib import Path
import migrate_geberit154446 as task

def finish(log):
    r=json.loads(task.REPORT.read_text())
    assert r['verdict']=='validated_pending_visual_check'
    log=Path(log)
    ready=[x.split('TASK_BLENDER_IFC_READY ',1)[1] for x in log.read_text().splitlines() if 'TASK_BLENDER_IFC_READY ' in x][-1]
    session=json.loads(ready)
    assert session['ifc']==str(task.SINGLE) and session['sha256']==task.pkg.sha256(task.SINGLE)
    assert not session['blend_saved'] and not session['has_blend_warning']
    assert session['scene_name']=='🔍 逐个审核高模图纸候选'
    session.update(log=str(log),owner_task_id='01a075a7-e948-7bd1-a27e-6ff542ea187c-geberit154446-migration',bridge_port=9895)
    r['fresh_pure_session']=session
    r['visual']={'method':'Independent image inspection of pure-view and actual Bonsai scene PNGs',
      'observations':['Plan and Front preserve full approved 900 mm channel length.',
       'Side retains the approved 12.241154 mm datum translation and 53.4 mm visible width.',
       'Scene Plan is below the ceiling; grey room geometry is readable and historical annotation symbols are absent.',
       'Front and Side retain world Z upright; drawing framing remains pending user scene acceptance.'],
      'verdict':'pass_migration','scene_approval_status':'pending'}
    r.update(verdict='pass',validation_verdict='pass',scene_approval_status='pending')
    for item in r['official_sources'].values():
        assert task.pkg.sha256(task.ROOT/item['path'])==item['sha256']
    r['all_archived_dwg_hashes_match']=True
    r['plan_camera_z_evidence']={'camera_z_mm':2350,'target_top_z_mm':20,'slab_underside_z_mm':2450,'slab_guid':'1Yf7cK9MTAoBRvfK3BluSI','slab_world_bounds_m':[[-6.50076208114624,.1501579284668002,2.45],[-4.90076208114624,2.6001579284668006,3.0]],'source':'IfcOpenShell world-coordinate physical geometry bounding box','product_translation_mm':[0,0,0]}
    r['native_generator_degenerate_path_exception']={'view':'side','path_index_zero_based':107,'length_mm':0.00004365776002709553,'outcome':'Native Bonsai omits the near-zero-length two-point path; its full exact coordinates remain in ApprovedSide. The 0.05 mm projection threshold is retained for every non-degenerate contour.','pure_ifc_all_176_paths_retained':True,'official_dwg_retained':True}
    r['scene_installation_observation']={'status':'requires_user_scene_judgment','observation':'In Plan the blue drain lies within the grey bathtub projection. This is retained from the formal model, not a migration translation.','source_body_and_placement_unchanged':True,'installation_position_approved':False,'parent_agent_visual_cross_check':'Three SCENE PNGs independently viewed: blue lines complete; no historical markers; bathtub/drain overlap retained for user judgment.'}
    assert task.pkg.sha256(task.FORMAL)==task.FORMAL_HASH
    task.write(task.REPORT,r)
    for name in ['manifest.json','handoff.json']:
        path=task.OUT/name;data=json.loads(path.read_text())
        data.update(status='complete',validation_verdict='pass',single_product_approval_status='approved',scene_approval_status='pending',approval_label_zh='单品已通过、场景待验收')
        data['fresh_pure_session']=session
        data['scene_installation_observation']=r['scene_installation_observation']
        data['native_generator_degenerate_path_exception']=r['native_generator_degenerate_path_exception']
        if name=='manifest.json':
            data['library_previews']['iso']=task.pkg.record(task.PRODUCT/'bonsai-camera-iso.png')
            data['preview_provenance']['3d']='Existing actual IFC Body camera render; all three original Body representations preserved.'
        else:data['ifc_sha256']=task.pkg.sha256(task.SINGLE)
        task.write(path,data)
    print('GEBERIT154446 migration complete; scene pending')

if __name__=='__main__':finish(sys.argv[1])
