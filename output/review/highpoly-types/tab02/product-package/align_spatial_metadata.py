"""Reuse the formal ancestor relationship identities in pure and recipe metadata."""
import json,ifcopenshell
import migrate_tab02 as task
from pure_product_package import graph_from_json,graph_to_json

def align(model,formal):
    changes=[]
    for rel in model.by_type('IfcRelAggregates'):
        related={x.GlobalId for x in rel.RelatedObjects}
        matches=[r for r in formal.by_type('IfcRelAggregates') if r.RelatingObject.GlobalId==rel.RelatingObject.GlobalId and related.issubset({x.GlobalId for x in r.RelatedObjects})]
        assert len(matches)==1
        if rel.GlobalId!=matches[0].GlobalId:
            changes.append({'old':rel.GlobalId,'formal':matches[0].GlobalId,'owner':rel.RelatingObject.GlobalId,'related':sorted(related)})
            rel.GlobalId=matches[0].GlobalId
    for cls,owner,related in [('IfcRelDefinesByType','RelatingType','RelatedObjects'),('IfcRelContainedInSpatialStructure','RelatingStructure','RelatedElements')]:
        for rel in model.by_type(cls):
            members={x.GlobalId for x in getattr(rel,related)}
            matches=[r for r in formal.by_type(cls) if getattr(r,owner).GlobalId==getattr(rel,owner).GlobalId and members.issubset({x.GlobalId for x in getattr(r,related)})]
            assert len(matches)==1
            if rel.GlobalId!=matches[0].GlobalId:
                changes.append({'old':rel.GlobalId,'formal':matches[0].GlobalId,'owner':getattr(rel,owner).GlobalId,'related':sorted(members)})
                rel.GlobalId=matches[0].GlobalId
    return changes

if __name__=='__main__':
    formal=ifcopenshell.open(str(task.FORMAL));pure=ifcopenshell.open(str(task.SINGLE))
    before=task.check_pure(pure);changes=align(pure,formal);assert task.check_pure(pure)==before
    pure.write(str(task.SINGLE))
    recipe=json.loads(task.RECIPE.read_text());graph=graph_from_json(recipe['graph']);recipe_changes=align(graph,formal)
    recipe['graph']=graph_to_json(graph);recipe['spatial_relationship_policy']='Reuse original formal ancestor relationship GUIDs; prevents attaching duplicate IfcRelAggregates for existing site/building.'
    task.write(task.RECIPE,recipe)
    r=json.loads(task.REPORT.read_text());r['spatial_relationship_alignment']={'pure':changes,'recipe':recipe_changes,'body_placement_and_representations_unchanged':True};task.write(task.REPORT,r)
