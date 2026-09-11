"""Remove only newly introduced duplicate metadata; preserve formal baseline defects."""
import json,sys
from collections import Counter
import ifcopenshell,ifcopenshell.util.element as eu
import migrate_tab02 as task

def repair(model,formal):
    changes=[]
    for cls,owner,related in [('IfcRelAggregates','RelatingObject','RelatedObjects'),('IfcRelDefinesByType','RelatingType','RelatedObjects'),('IfcRelContainedInSpatialStructure','RelatingStructure','RelatedElements')]:
        known={r.GlobalId for r in formal.by_type(cls)}
        for rel in list(model.by_type(cls)):
            if rel.GlobalId in known:continue
            members={x.GlobalId for x in getattr(rel,related)}
            matches=[r for r in formal.by_type(cls) if getattr(r,owner).GlobalId==getattr(rel,owner).GlobalId and members.issubset({x.GlobalId for x in getattr(r,related)})]
            if matches:
                assert len(matches)==1
                original=model.by_guid(matches[0].GlobalId)
                assert members.issubset({x.GlobalId for x in getattr(original,related)})
                changes.append({'class':cls,'discarded_new_duplicate':rel.GlobalId,'preserved_formal':original.GlobalId})
                model.remove(rel)
    counts=Counter(task.pkg.fingerprint(a) for a in formal.by_type('IfcApplication'))
    remaining=counts.copy();retained=[];extra=[]
    for app in model.by_type('IfcApplication'):
        h=task.pkg.fingerprint(app)
        if remaining[h]>0:remaining[h]-=1;retained.append(app)
        else:extra.append(app)
    assert not +remaining
    for app in extra:
        matches=[a for a in retained if (a.ApplicationIdentifier,a.Version,a.ApplicationFullName)==(app.ApplicationIdentifier,app.Version,app.ApplicationFullName)]
        assert matches,app
        keep=matches[0]
        for inverse in model.get_inverse(app):eu.replace_attribute(inverse,app,keep)
        changes.append({'class':'IfcApplication','removed_new_step':app.id(),'reused_baseline_step':keep.id(),'identity':[keep.ApplicationIdentifier,keep.Version,keep.ApplicationFullName]})
        model.remove(app)
    assert counts==Counter(task.pkg.fingerprint(a) for a in model.by_type('IfcApplication'))
    return changes

def finish():
    import bpy,bonsai_bridge as bridge
    from bonsai import tool
    import create_wd03_wardrobe_scene_drawings as context
    r=json.loads(task.REPORT.read_text());path=r['temporary_project']['path']
    assert tool.Ifc.get_path()==str(task.SINGLE)
    with bpy.context.temp_override(**context.view3d_override()):
        r['product_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':str(task.SINGLE),'overwrite':True,'reload':True})
    assert task.check_pure(ifcopenshell.open(str(task.SINGLE)))==r['pure_pre_bonsai']
    r['single_product']=task.pkg.record(task.SINGLE)
    assert bpy.ops.bim.load_project(filepath=path,should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}
    f=tool.Ifc.get();formal=ifcopenshell.open(str(task.FORMAL))
    before={e.GlobalId:(task.pkg.body_fingerprint(e) if e.Representation else None,task.pkg.placement(e).tolist()) for e in f.by_type('IfcElement')}
    r['runtime_metadata_reuse']=repair(f,formal)
    with bpy.context.temp_override(**context.view3d_override()):
        r['temporary_bonsai_save_reload']=bridge._h_save_ifc_file({'output_path':path,'overwrite':True,'reload':True})
    reload=ifcopenshell.open(path)
    assert before=={e.GlobalId:(task.pkg.body_fingerprint(e) if e.Representation else None,task.pkg.placement(e).tolist()) for e in reload.by_type('IfcElement')}
    assert Counter(task.pkg.fingerprint(a) for a in formal.by_type('IfcApplication'))==Counter(task.pkg.fingerprint(a) for a in reload.by_type('IfcApplication'))
    r['temporary_project']=task.pkg.record(path);r['metadata_repair_saved_reloaded']=True
    task.write(task.REPORT,r)
    assert bpy.ops.bim.load_project(filepath=str(task.SINGLE),should_start_fresh_session=False,use_relative_path=False)=={'FINISHED'}

if __name__=='__main__':finish()
