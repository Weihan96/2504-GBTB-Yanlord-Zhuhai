"""Preserve approved CleanLine50 geometry; reconstruct pending R12 scene externally."""
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
OLD = PRODUCT / "Geberit-154-446-KS-1-bonsai-isolated.ifc"
APPROVAL = ROOT / "pipeline/decisions/geberit-154-446-ks-1-drawing-approval.json"
SINGLE = OUT / "GEBERIT154446-product.ifc"
RECIPE = OUT / "scene-recipe.json"
REPORT = OUT / "validation.json"
GUID = "2S2c498tb7$gzdukjhCGVQ"
SOURCE_LABEL = "Geberit 官方精确归档154.446.KS.1原生DWG蓝线（批准Side基准修正版）"
DWG_HASH = "dedb47964cc69310bd9af3982379c428c8b20494b7588296aedc5d07380e3421"

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
    assert len([r for r in model.by_guid(GUID).Representation.Representations if r.RepresentationIdentifier=="Body"]) == 3
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
    import geberit_154_446_ks_1_drawing_ifc as approved
    assert pkg.sha256(FORMAL)==FORMAL_HASH
    assert json.loads(APPROVAL.read_text())["status"] == "approved"
    candidate_path=PRODUCT/"candidate-representations.json"
    candidate=json.loads(candidate_path.read_text())
    paths=approved.official_paths(candidate, json.loads((PRODUCT/"official-native-dwg-linework.json").read_text()))
    assert pkg.sha256(PRODUCT/"manifest.json")==json.loads(APPROVAL.read_text())["candidate_manifest_sha256"]
    official_manifest=json.loads((PRODUCT/"manifest.json").read_text())
    alignment=official_manifest["official_reference"]["mechanical_cross_checks"]["side"]["horizontal_alignment"]
    assert alignment["pass"] and alignment["applied_translation_mm"]==12.241154
    assert pkg.sha256(PRODUCT/"official-source/154.446.KS.1_L.dwg")==DWG_HASH
    source=ifcopenshell.open(str(OLD)); formal=ifcopenshell.open(str(FORMAL)); target=source.by_guid(GUID)
    assert pkg.body_fingerprint(target)==pkg.body_fingerprint(formal.by_guid(GUID))
    assert np.array_equal(pkg.placement(target),pkg.placement(formal.by_guid(GUID)))
    pre=snapshot(source)
    protected=[pkg.record(p) for p in [FORMAL, OLD, APPROVAL, PRODUCT/"manifest.json", candidate_path,
               *[p for p in (PRODUCT/"official-source").rglob("*") if p.is_file()],
               *[PRODUCT/f"{v}.svg" for v in paths]]]
    definitions={"plan":"FFL PLAN", "front":"EL-06-21-R12-PX", "side":"EL-06-20-R12-NY"}
    includes=eu.get_pset(next(d for d in formal.by_type("IfcAnnotation") if d.Name==definitions["front"]),"EPset_Drawing")["Include"].split(",")
    # Restrict movable fixtures to the actual wet-room interior; no bbox tolerance.
    from ifcopenshell import geom
    import ifcopenshell.util.shape as su
    settings=geom.settings();settings.set("use-world-coords",True)
    filtered=[];context_audit=[]
    for guid in includes:
        e=formal.by_guid(guid)
        keep=e.is_a("IfcWall") or e.is_a("IfcSlab") or e.is_a("IfcBeam") or e.is_a("IfcCovering")
        if not keep and e.Representation:
            try:
                shape=geom.create_shape(settings,e);vertices=su.get_vertices(shape.geometry)
                center=(vertices.min(0)+vertices.max(0))/2
                keep=-6.6<center[0]<-4.8 and .1<center[1]<2.7
                context_audit.append({"guid":guid,"center_m":center.tolist(),"included":bool(keep)})
            except RuntimeError:keep=False
        if keep:filtered.append(guid)
    includes=filtered
    specs={}; evidence=[]
    for view,name in definitions.items():
        original=next(d for d in formal.by_type("IfcAnnotation") if d.Name==name)
        # Copy only the camera's forward geometry/placement. New root identities
        # keep formal Drawings and their groups outside this pending recipe.
        copier=pkg.ScopedCopy(formal,source,{original},skip_inverse_ids=[original.id()])
        drawing=copier.copy(original)
        drawing.GlobalId=ifcopenshell.guid.new();drawing.Name=f"GEBERIT154446-SCENE-{view.upper()}"
        drawing.Description="Pending scene validation; approved single-product outline only"
        values={k:v for k,v in eu.get_pset(original,"EPset_Drawing").items() if k not in ("id","Exclude","Include")}
        values.update(HasAnnotation=True, HasUnderlay=False, GlobalReferencing=False,
                      Include=",".join(g for g in includes if g!=GUID),
                      Exclude=",".join(['IfcSpace','IfcGrid','IfcBuildingStorey',*[a.GlobalId for a in formal.by_type('IfcAnnotation')]]))
        for key in ("Stylesheet","Markers","Symbols","Patterns","ShadingStyles"):
            values[key]=str(ROOT/values[key])
        if view=="plan":
            # Same R12 context scope, viewed from above the 2210 mm wall head.
            loc=drawing.ObjectPlacement.RelativePlacement.Location
            loc.Coordinates=(-5700.0,1400.0,2350.0)
            block=next(x for x in source.traverse(drawing.Representation) if x.is_a("IfcBlock"))
            block.XLength=1800.0;block.YLength=2600.0;block.ZLength=3000.0
            block.Position.Location.Coordinates=(-block.XLength/2,-block.YLength/2,-block.ZLength)
        pset(source,drawing,"EPset_Drawing",values)
        document=source.create_entity("IfcDocumentReference",Location=str(OUT/f"{drawing.Name}.svg"),Identification=drawing.Name,Name=drawing.Name)
        source.create_entity("IfcRelAssociatesDocument",GlobalId=ifcopenshell.guid.new(),RelatedObjects=[drawing],RelatingDocument=document)
        context=approved.representation_context(source,"Annotation",values["TargetView"])
        rep=approved.curve_representation(source,context,"Annotation",view,paths[view])
        color=source.create_entity("IfcColourRgb",Red=22/255,Green=119/255,Blue=200/255)
        style=source.create_entity("IfcCurveStyle",Name="Geberit 154.446.KS.1 approved official DWG blue",CurveWidth=source.create_entity("IfcPositiveLengthMeasure",0.35),CurveColour=color,ModelOrDraughting=True)
        source.create_entity("IfcStyledItem",Item=rep.Items[0],Styles=[style])
        annotation=source.create_entity("IfcAnnotation",GlobalId=ifcopenshell.guid.new(),Name=f"Geberit 154.446.KS.1 approved {view}",ObjectType="LINEWORK",ObjectPlacement=target.ObjectPlacement,Representation=source.create_entity("IfcProductDefinitionShape",Representations=[rep]))
        view_source=official_manifest["official_reference"]["official_sources"][{"plan":"G","front":"A","side":"L"}[view]]
        pset(source,annotation,"EPset_Annotation",{"Classes":"review-target-geberit154446 official-native-dwg","TargetGlobalId":GUID,"SourceKind":"official_native_dwg_paths_mm","SourceDwgSha256":view_source["sha256"],"OfficialDownloadUrl":view_source["url"]})
        group=source.create_entity("IfcGroup",GlobalId=ifcopenshell.guid.new(),Name=drawing.Name,ObjectType="DRAWING")
        source.create_entity("IfcRelAssignsToGroup",GlobalId=ifcopenshell.guid.new(),RelatedObjects=[drawing,annotation],RelatingGroup=group)
        specs[view]={"annotation_guid":annotation.GlobalId,"drawing_guid":drawing.GlobalId}
        evidence.append({"view":view,**specs[view],"source_camera_guid":original.GlobalId,"source_camera_name":name,"camera_matrix":pkg.placement(drawing).tolist(),"approved_path_count":len(paths[view])})
    pset(source,target,"Pset_Geberit154446ApprovedSource",{"SourceKind":"official_native_dwg_paths_mm","SourceLabelZh":SOURCE_LABEL,"SourceDwgSha256":DWG_HASH,"OfficialDownloadUrl":"https://cdn.data.geberit.com/cad/154.446.KS.1_L.dwg","ArticleNumber":"154.446.KS.1","UnalteredOfficialDwg":True,"SingleProductApprovalStatus":"approved","SceneApprovalStatus":"pending","FormalIfcWriteAllowed":False})
    pure,recipe,audit=build_pure_package(source,GUID,specs)
    audit['unused_type_box_maps_removed']=discard_type_box(pure)
    pure.write(str(SINGLE));write(RECIPE,recipe)
    state=check_pure(ifcopenshell.open(str(SINGLE)))
    deps=[]
    for location in sorted({x.Location for x in pure.by_type("IfcExternallyDefinedSurfaceStyle") if x.Location}):
        rel=Path(location);assert not rel.is_absolute() and ".." not in rel.parts
        dest=OUT/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,dest)
        deps.append({"source":pkg.record(ROOT/rel),"packaged":pkg.record(dest),"location":location})
    write(REPORT,{"product":"GEBERIT154446","target_guid":GUID,"verdict":"pending_bonsai_scene","formal_sha256_before":FORMAL_HASH,
        "protected_files":protected,"preState":pre,"pure_pre_bonsai":state,"extraction":audit,"views":evidence,"side_alignment_evidence":alignment,"official_sources":official_manifest["official_reference"]["official_sources"],"scene_context_selection":context_audit,
        "single_product":pkg.record(SINGLE),"scene_recipe":pkg.record(RECIPE),"necessary_external_style_dependencies":deps,
        "single_product_approval_status":"approved","scene_approval_status":"pending","formal_write_allowed":False,"cleanup_performed":False,
        "source_kind":"official_native_dwg_paths_mm","source_label_zh":SOURCE_LABEL,
        "scene_camera_policy":"R12 original PX/NY cameras inside wet room; Plan at2350mm above target20mm and below slab2450mm. Candidate scene only.",
        "courseEvidence":{"mode":"embedded-course-index","lesson":"085000","timestamp":"01:59 Create Drawing","screenshot":"research/bonsai-course/lessons/085000/screenshots/085000-01m59s-create-drawing-button.png","provenance_sha256":"6997420afc8e3d91d0016b8c68c7db904cc782781317a0de7ff50e9d10f945dc","label":"course_fact"}})
    print(json.dumps({"pure":pkg.record(SINGLE),"state":state}))

def style_svg(svg, annotation_guid, view):
    import xml.etree.ElementTree as ET
    raw=pkg.record(svg);tree=ET.parse(svg);root=tree.getroot();count=0
    for e in root.iter():
        if e.tag.rsplit("}",1)[-1] in ("line","polyline","path","polygon"):
            e.set("style","stroke:#8b949e;stroke-width:0.18;fill:none")
        if annotation_guid in e.get("class","") and e.tag.rsplit("}",1)[-1] in ("line","polyline","path"):
            e.set("style","stroke:#1677c8;stroke-width:0.18;fill:none");e.set("data-target-global-id",GUID);e.set("class",e.get("class","")+" targetGlobalId-"+GUID);count+=1
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
    temporary=Path(tempfile.mkdtemp(prefix="geberit154446-pure-package-scene-"));temp_ifc=temporary/"scene.ifc"
    report.update(temporary_project_directory=str(temporary),runtime_inputs=[str(SINGLE),str(RECIPE),str(FORMAL)],legacy_ifc_used_for_runtime=False)
    write(REPORT,report)
    try:
        report["provider"]={"status":"supported","version":list(bridge.bl_info["version"]),"port":9895,"pid":__import__('os').getpid(),"blender":bpy.app.version_string,"ifcopenshell":ifcopenshell.version,"bridge_source":pkg.record(bridge.__file__)}
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
            svg=OUT/f"GEBERIT154446-SCENE-{v['view'].upper()}.svg"
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
