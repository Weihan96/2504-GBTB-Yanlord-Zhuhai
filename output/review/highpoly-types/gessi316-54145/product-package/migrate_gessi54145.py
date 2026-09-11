"""Preserve approved 54145 geometry; reconstruct pending R17 scene externally."""
from pathlib import Path
import sys, json, shutil, tempfile, traceback
import numpy as np
import ifcopenshell
import ifcopenshell.api.pset
import ifcopenshell.util.element as eu
import ifcopenshell.util.placement as pu

OUT = Path(__file__).resolve().parent
PRODUCT = OUT.parent
ROOT = OUT.parents[4]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
import review_product_package as pkg
from pure_product_package import build_pure_package, attach_recipe
FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_HASH = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OLD = PRODUCT / "Gessi316-54145-bonsai-isolated.ifc"
APPROVAL = ROOT / "pipeline/decisions/gessi316-54145-drawing-approval.json"
SINGLE = OUT / "GESSI54145-product.ifc"
RECIPE = OUT / "scene-recipe.json"
REPORT = OUT / "validation.json"
GUID = "3jT4sCgpHC98VSIUdGUNYH"
SOURCE_LABEL = "基于官方54145 G000原生DWG轮廓的简化蓝线审核表达"
DWG_HASH = "9978b68468a61875acb0736aab45a62a98efadfd6be61d347e7ecaa94e08fd09"

def write(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")

def snapshot(model):
    t = model.by_guid(GUID)
    return {"body":pkg.body_fingerprint(t), "representations":pkg.fingerprint(t.Representation),
            "placement":pkg.placement(t).tolist(), "physical_elements":sorted(e.GlobalId for e in model.by_type("IfcElement")),
            "annotation_count":len(model.by_type("IfcAnnotation")), "drawing_pset_count":len([p for p in model.by_type("IfcPropertySet") if p.Name=="EPset_Drawing"]),
            "unit_scale_to_m":pkg.unit_util.calculate_unit_scale(model)}

def check_pure(model):
    s = snapshot(model)
    assert s["physical_elements"] == [GUID]
    assert s["annotation_count"] == s["drawing_pset_count"] == 0
    assert not model.by_type("IfcGroup")
    assert not [x for x in model.by_type("IfcSpatialElement") if x.Representation]
    assert not [x for x in model.by_type("IfcPropertySingleValue") if x.Name in ("Include", "Exclude")]
    assert not [x for x in model.by_type("IfcShapeRepresentation") if x.RepresentationIdentifier == "Box"]
    assert len([r for r in model.by_guid(GUID).Representation.Representations if r.RepresentationIdentifier=="Body"]) == 2
    assert len(model.by_type("IfcRoot")) == len({x.GlobalId for x in model.by_type("IfcRoot")})
    return s

def discard_type_box(model):
    target=model.by_guid(GUID);typ=eu.get_type(target)
    boxes=[m for m in typ.RepresentationMaps if m.MappedRepresentation.RepresentationIdentifier=="Box"]
    typ.RepresentationMaps=tuple(m for m in typ.RepresentationMaps if m not in boxes)
    for box in boxes:eu.remove_deep2(model,box)
    return len(boxes)

def pset(model, target, name, values):
    p = ifcopenshell.api.pset.add_pset(model, product=target, name=name)
    ifcopenshell.api.pset.edit_pset(model, pset=p, properties=values)

def prepare():
    import gessi316_54145_drawing_ifc as approved
    assert pkg.sha256(FORMAL)==FORMAL_HASH
    assert json.loads(APPROVAL.read_text())["status"] == "approved"
    candidate_path=PRODUCT/"candidate-representations.json"
    candidate=json.loads(candidate_path.read_text())
    paths=approved.candidate_paths(candidate, json.loads((PRODUCT/"official-source/source-access-record.json").read_text()))
    assert pkg.sha256(PRODUCT/"official-source/GPF5414500000G000_3.dwg")==DWG_HASH
    source=ifcopenshell.open(str(OLD)); formal=ifcopenshell.open(str(FORMAL)); target=source.by_guid(GUID)
    assert pkg.body_fingerprint(target)==pkg.body_fingerprint(formal.by_guid(GUID))
    assert np.array_equal(pkg.placement(target),pkg.placement(formal.by_guid(GUID)))
    pre=snapshot(source)
    protected=[pkg.record(p) for p in [FORMAL, OLD, APPROVAL, PRODUCT/"manifest.json", candidate_path,
               *[p for p in (PRODUCT/"official-source").rglob("*") if p.is_file()],
               *[PRODUCT/f"{v}.svg" for v in paths]]]
    definitions={"plan":"FFL PLAN", "front":"EL-08-31-R17-NY", "side":"EL-08-32-R17-NX"}
    includes=eu.get_pset(next(d for d in formal.by_type("IfcAnnotation") if d.Name==definitions["front"]),"EPset_Drawing")["Include"].split(",")
    specs={}; evidence=[]
    for view,name in definitions.items():
        original=next(d for d in formal.by_type("IfcAnnotation") if d.Name==name)
        # Copy only the camera's forward geometry/placement. New root identities
        # keep formal Drawings and their groups outside this pending recipe.
        copier=pkg.ScopedCopy(formal,source,{original},skip_inverse_ids=[original.id()])
        drawing=copier.copy(original)
        drawing.GlobalId=ifcopenshell.guid.new();drawing.Name=f"GESSI54145-SCENE-{view.upper()}"
        drawing.Description="Pending scene validation; approved single-product outline only"
        values={k:v for k,v in eu.get_pset(original,"EPset_Drawing").items() if k not in ("id","Exclude","Include")}
        values.update(HasAnnotation=True, HasUnderlay=False, GlobalReferencing=False,
                      Include=",".join(g for g in includes if g!=GUID),
                      Exclude=",".join(['IfcSpace','IfcGrid','IfcBuildingStorey',*[a.GlobalId for a in formal.by_type('IfcAnnotation')]]))
        for key in ("Stylesheet","Markers","Symbols","Patterns","ShadingStyles"):
            values[key]=str(ROOT/values[key])
        if view=="plan":
            # Same R17 context scope, viewed from above the 2210 mm wall head.
            loc=drawing.ObjectPlacement.RelativePlacement.Location
            loc.Coordinates=(-1329.32305335999, 2689.7349357605, 2350.0)
            block=next(x for x in source.traverse(drawing.Representation) if x.is_a("IfcBlock"))
            block.XLength=2000.62797881745;block.YLength=2280.37680636395;block.ZLength=3200.0
            block.Position.Location.Coordinates=(-block.XLength/2,-block.YLength/2,-block.ZLength)
        pset(source,drawing,"EPset_Drawing",values)
        document=source.create_entity("IfcDocumentReference",Location=str(OUT/f"{drawing.Name}.svg"),Identification=drawing.Name,Name=drawing.Name)
        source.create_entity("IfcRelAssociatesDocument",GlobalId=ifcopenshell.guid.new(),RelatedObjects=[drawing],RelatingDocument=document)
        context=approved.shared.representation_context(source,"Annotation",values["TargetView"])
        rep=approved.shared.curve_representation(source,context,"Annotation",view,paths[view])
        color=source.create_entity("IfcColourRgb",Red=22/255,Green=119/255,Blue=200/255)
        style=source.create_entity("IfcCurveStyle",Name="Gessi 54145 approved simplified blue",CurveWidth=source.create_entity("IfcPositiveLengthMeasure",0.35),CurveColour=color,ModelOrDraughting=True)
        source.create_entity("IfcStyledItem",Item=rep.Items[0],Styles=[style])
        annotation=source.create_entity("IfcAnnotation",GlobalId=ifcopenshell.guid.new(),Name=f"Gessi 54145 approved {view}",ObjectType="LINEWORK",ObjectPlacement=target.ObjectPlacement,Representation=source.create_entity("IfcProductDefinitionShape",Representations=[rep]))
        pset(source,annotation,"EPset_Annotation",{"Classes":"review-target-gessi54145 native-dwg-review-simplification","TargetGlobalId":GUID,"SourceKind":"native_dwg_review_simplification","SourceDwgSha256":DWG_HASH})
        group=source.create_entity("IfcGroup",GlobalId=ifcopenshell.guid.new(),Name=drawing.Name,ObjectType="DRAWING")
        source.create_entity("IfcRelAssignsToGroup",GlobalId=ifcopenshell.guid.new(),RelatedObjects=[drawing,annotation],RelatingGroup=group)
        specs[view]={"annotation_guid":annotation.GlobalId,"drawing_guid":drawing.GlobalId}
        evidence.append({"view":view,**specs[view],"source_camera_guid":original.GlobalId,"source_camera_name":name,"camera_matrix":pkg.placement(drawing).tolist(),"approved_path_count":len(paths[view])})
    pset(source,target,"Pset_Gessi54145ApprovedSource",{"SourceKind":"native_dwg_review_simplification","SourceLabelZh":SOURCE_LABEL,"SourceDwgSha256":DWG_HASH,"OfficialDownloadUrl":approved.DWG_ZIP,"ArticleNumber":"54145","Configuration":"G000","UnalteredOfficialDwg":False,"SingleProductApprovalStatus":"approved","SceneApprovalStatus":"pending","FormalIfcWriteAllowed":False})
    pure,recipe,audit=build_pure_package(source,GUID,specs)
    audit['unused_type_box_maps_removed']=discard_type_box(pure)
    pure.write(str(SINGLE));write(RECIPE,recipe)
    state=check_pure(ifcopenshell.open(str(SINGLE)))
    deps=[]
    for location in sorted({x.Location for x in pure.by_type("IfcExternallyDefinedSurfaceStyle") if x.Location}):
        rel=Path(location);assert not rel.is_absolute() and ".." not in rel.parts
        dest=OUT/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,dest)
        deps.append({"source":pkg.record(ROOT/rel),"packaged":pkg.record(dest),"location":location})
    write(REPORT,{"product":"GESSI54145","target_guid":GUID,"verdict":"pending_bonsai_scene","formal_sha256_before":FORMAL_HASH,
        "protected_files":protected,"preState":pre,"pure_pre_bonsai":state,"extraction":audit,"views":evidence,
        "single_product":pkg.record(SINGLE),"scene_recipe":pkg.record(RECIPE),"necessary_external_style_dependencies":deps,
        "single_product_approval_status":"approved","scene_approval_status":"pending","formal_write_allowed":False,"cleanup_performed":False,
        "source_kind":"native_dwg_review_simplification","source_label_zh":SOURCE_LABEL,
        "scene_camera_policy":"R17 original NY/NX camera metadata; plan at 2350 mm, above shower top 2241 and below slab underside 2450, with R17 crop. Candidate scene only.",
        "courseEvidence":{"mode":"embedded-course-index","lesson":"085000","timestamp":"01:59 Create Drawing","screenshot":"research/bonsai-course/lessons/085000/screenshots/085000-01m59s-create-drawing-button.png","provenance_sha256":"6997420afc8e3d91d0016b8c68c7db904cc782781317a0de7ff50e9d10f945dc","label":"course_fact"}})
    print(json.dumps({"pure":pkg.record(SINGLE),"state":state}))

def tighten_existing_pure():
    model=ifcopenshell.open(str(SINGLE));before=snapshot(model)
    removed=discard_type_box(model)
    assert snapshot(model)==before
    model.write(str(SINGLE));check_pure(ifcopenshell.open(str(SINGLE)))
    report=json.loads(REPORT.read_text());report['extraction']['unused_type_box_maps_removed']=removed
    report['single_product']=pkg.record(SINGLE);report['pure_saved_reloaded']=False
    write(REPORT,report)

def filter_recipe_references():
    from pure_product_package import graph_from_json,graph_to_json
    recipe=json.loads(RECIPE.read_text());graph=graph_from_json(recipe['graph'])
    formal=ifcopenshell.open(str(FORMAL))
    excluded=['IfcSpace','IfcGrid','IfcBuildingStorey',*[a.GlobalId for a in formal.by_type('IfcAnnotation')]]
    for spec in recipe['views'].values():
        drawing=graph.by_guid(spec['drawing_guid']);data=eu.get_pset(drawing,'EPset_Drawing')
        ifcopenshell.api.pset.edit_pset(graph,pset=graph.by_id(data['id']),properties={'Exclude':','.join(excluded)})
    recipe['graph']=graph_to_json(graph);write(RECIPE,recipe)
    report=json.loads(REPORT.read_text());report['scene_recipe']=pkg.record(RECIPE)
    report['scene_reference_filter']={'excluded_historical_annotation_count':len(excluded)-3,'excluded_classes':excluded[:3],'source':'Installed Bonsai get_potential_reference_elements and get_drawing_spaces each consult Exclude independently of Include','product_geometry_unchanged':True}
    write(REPORT,report)

def lower_plan_camera_below_ceiling():
    from pure_product_package import graph_from_json,graph_to_json
    recipe=json.loads(RECIPE.read_text());graph=graph_from_json(recipe['graph'])
    camera=graph.by_guid(recipe['views']['plan']['drawing_guid'])
    location=camera.ObjectPlacement.RelativePlacement.Location
    location.Coordinates=(*location.Coordinates[:2],2350.0)
    recipe['graph']=graph_to_json(graph);write(RECIPE,recipe)
    report=json.loads(REPORT.read_text());report['scene_recipe']=pkg.record(RECIPE)
    report['scene_camera_policy']='Plan camera 2350 mm: above original shower top 2241 mm; below original slab underside 2450 mm and beam underside 2600 mm. Physical geometry unchanged.'
    report['plan_camera_z_evidence']={'camera_z_mm':2350,'target_top_z_mm':2241,'slab_underside_z_mm':2450,'slab_guid':'3qo2s5gKPA$8hNQznuTTCf','beam_underside_z_mm':2600,'beam_guid':'3gx_uavED0yvBMYfw86Rw4','measurement':'IfcOpenShell world-coordinate geometry bounding boxes','blue_line_translation_mm':[0,0,0]}
    for v in report['views']:
        if v['view']=='plan':v['camera_matrix']=pkg.placement(camera).tolist()
    write(REPORT,report)

def style_svg(svg, annotation_guid, view):
    import xml.etree.ElementTree as ET
    raw=pkg.record(svg);tree=ET.parse(svg);root=tree.getroot();count=0
    for e in root.iter():
        if annotation_guid in e.get("class","") and e.tag.rsplit("}",1)[-1] in ("line","polyline","path"):
            e.set("style","stroke:#1677c8;stroke-width:0.18;fill:none");e.set("data-target-global-id",GUID);count+=1
    assert count>0,"No 54145 generated annotation geometry"
    root.set("data-create-drawing-result","FINISHED");root.set("data-scene-approval-status","pending")
    tree.write(svg,encoding="utf-8",xml_declaration=True)
    return {"raw":raw,"blue_geometry_count":count,"post_style_only":True,'geometry_moved_removed_or_redrawn':False}

def scene():
    import bpy,bonsai_bridge as bridge
    from bonsai import tool
    import create_wd03_wardrobe_scene_drawings as context
    report=json.loads(REPORT.read_text())
    if report.get('temporary_project_directory'):
        report.setdefault('previous_temporary_attempts',[]).append(report['temporary_project_directory'])
    assert Path(tool.Ifc.get_path()).resolve()==SINGLE.resolve() and not bpy.data.is_saved
    temporary=Path(tempfile.mkdtemp(prefix="gessi54145-pure-package-scene-"));temp_ifc=temporary/"scene.ifc"
    report.update(temporary_project_directory=str(temporary),runtime_inputs=[str(SINGLE),str(RECIPE),str(FORMAL)],legacy_ifc_used_for_runtime=False)
    write(REPORT,report)
    try:
        report["provider"]={"status":"supported","version":list(bridge.bl_info["version"]),"port":9892,"pid":__import__('os').getpid(),"blender":bpy.app.version_string,"ifcopenshell":ifcopenshell.version,"bridge_source":pkg.record(bridge.__file__)}
        with bpy.context.temp_override(**context.view3d_override()):
            report["product_bonsai_save_reload"]=bridge._h_save_ifc_file({"output_path":str(SINGLE),"overwrite":True,"reload":True})
        pure=ifcopenshell.open(str(SINGLE));assert check_pure(pure)==report["pure_pre_bonsai"]
        report["pure_saved_reloaded"]=True
        assert pkg.sha256(FORMAL)==FORMAL_HASH
        shutil.copy2(FORMAL,temp_ifc);model=ifcopenshell.open(str(temp_ifc))
        baseline={e.GlobalId:(pkg.body_fingerprint(e) if e.Representation else None,pkg.placement(e).tolist()) for e in model.by_type("IfcElement")}
        report["formal_physical_elements"]=len(baseline)
        for loc in sorted({s.Location for s in model.by_type("IfcExternallyDefinedSurfaceStyle") if s.Location}):
            rel=Path(loc);assert not rel.is_absolute() and ".." not in rel.parts
            if (ROOT/rel).is_file():
                dest=temporary/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,dest)
        recipe=json.loads(RECIPE.read_text());report["attachment"]=attach_recipe(model,pure,recipe)
        count=len(list(model));attach_recipe(model,pure,recipe);assert len(list(model))==count
        report["second_attachment_created_entities"]=0
        model.write(str(temp_ifc))
        assert bpy.ops.bim.load_project(filepath=str(temp_ifc),should_start_fresh_session=False,use_relative_path=False)=={"FINISHED"}
        outputs=[]
        for v in report["views"]:
            entity=tool.Ifc.get().by_guid(v["drawing_guid"])
            tool.Ifc.get_object(entity) or tool.Drawing.import_drawing(entity)
            with bpy.context.temp_override(**context.view3d_override()):
                assert bpy.ops.bim.activate_drawing(drawing=entity.id(),should_view_from_camera=False)=={"FINISHED"}
                props=tool.Drawing.get_document_props();props.should_use_underlay_cache=False;props.should_use_linework_cache=False;props.should_use_annotation_cache=False
                result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
            assert result=={"FINISHED"}
            svg=OUT/f"GESSI54145-SCENE-{v['view'].upper()}.svg"
            styles=style_svg(svg,v["annotation_guid"],v["view"])
            outputs.append({"view":v["view"],"svg":pkg.record(svg),"styles":styles,"operator":"bpy.ops.bim.create_drawing","result":sorted(result)})
            report["scene_outputs"]=outputs;write(REPORT,report)
        with bpy.context.temp_override(**context.view3d_override()):
            report["temporary_bonsai_save_reload"]=bridge._h_save_ifc_file({"output_path":str(temp_ifc),"overwrite":True,"reload":True})
        reloaded=ifcopenshell.open(str(temp_ifc))
        assert baseline=={e.GlobalId:(pkg.body_fingerprint(e) if e.Representation else None,pkg.placement(e).tolist()) for e in reloaded.by_type("IfcElement")}
        report.update(scene_saved_reloaded=True,all_formal_bodies_unchanged=True,all_formal_placements_unchanged=True,target_instances=1,temporary_project=pkg.record(temp_ifc),verdict="pending_independent_visual_validation")
        assert bpy.ops.bim.load_project(filepath=str(SINGLE),should_start_fresh_session=False,use_relative_path=False)=={"FINISHED"}
    except Exception:
        report["error"]=traceback.format_exc();report["verdict"]="fail";raise
    finally:
        report["formal_sha256_after"]=pkg.sha256(FORMAL);report["single_product"]=pkg.record(SINGLE);write(REPORT,report)
        assert report["formal_sha256_after"]==FORMAL_HASH

if __name__=="__main__":
    prepare()
