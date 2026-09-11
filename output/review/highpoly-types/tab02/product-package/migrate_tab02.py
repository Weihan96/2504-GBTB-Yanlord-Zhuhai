"""Migrate approved TAB02 component boundaries into pure IFC and pending living-room views."""
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
OLD = PRODUCT / "RODA-Bernardo-367-TAB02-bonsai-isolated.ifc"
APPROVAL = ROOT / "pipeline/decisions/tab02-drawing-approval.json"
SINGLE = OUT / "TAB02-product.ifc"
RECIPE = OUT / "scene-recipe.json"
REPORT = OUT / "validation.json"
GUID = "2eX84IyLr8_e34nuOoRPLQ"
SOURCE_LABEL = "基于原始高模几何生成的简化图纸表达"

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
    assert len([r for r in model.by_guid(GUID).Representation.Representations if r.RepresentationIdentifier=="Body"]) == 1
    assert len(model.by_type("IfcRoot")) == len({x.GlobalId for x in model.by_type("IfcRoot")})
    return s

def discard_type_box(model):
    target=model.by_guid(GUID);typ=eu.get_type(target)
    boxes=[m for m in (typ.RepresentationMaps or ()) if m.MappedRepresentation.RepresentationIdentifier=="Box"]
    typ.RepresentationMaps=tuple(m for m in (typ.RepresentationMaps or ()) if m not in boxes)
    for box in boxes:eu.remove_deep2(model,box)
    return len(boxes)

def pset(model, target, name, values):
    p = ifcopenshell.api.pset.add_pset(model, product=target, name=name)
    ifcopenshell.api.pset.edit_pset(model, pset=p, properties=values)

def curve_representation(model,context,view,paths):
    axes={"plan":(0,1),"front":(0,2),"side":(1,2)}[view]
    lines=[]
    for path in paths:
        points=[]
        for p in path:
            xyz=[0.0,0.0,0.0];xyz[axes[0]]=float(p[0]);xyz[axes[1]]=float(p[1])
            points.append(model.create_entity("IfcCartesianPoint",Coordinates=xyz))
        lines.append(model.create_entity("IfcPolyline",Points=points))
    return model.create_entity("IfcShapeRepresentation",ContextOfItems=context,RepresentationIdentifier="Annotation",RepresentationType="GeometricCurveSet",Items=[model.create_entity("IfcGeometricCurveSet",Elements=lines)])

def prepare():
    import hima01_drawing_ifc as approved
    assert pkg.sha256(FORMAL)==FORMAL_HASH
    assert json.loads(APPROVAL.read_text())["status"] == "approved"
    assert json.loads(APPROVAL.read_text())["candidate_manifest_sha256"]==pkg.sha256(PRODUCT/"manifest.json")
    authorization=json.loads((ROOT/"output/review/approved-product-library/migration-authorization.json").read_text())
    assert authorization["pure_product_ifc_write_authorized"] and "tab02" in authorization["products"]
    candidate_path=PRODUCT/"candidate-representations.json"
    candidate=json.loads(candidate_path.read_text())
    semantics=json.loads((PRODUCT/"tab02-semantic-segmentation.json").read_text())
    paths={v:[e["path_mm"] for e in entries] for v,entries in semantics["views"].items()}
    assert {v:len(p) for v,p in paths.items()}=={"plan":1,"front":4,"side":4}
    assert not candidate["official_cad_used"]
    assert all(paths[v]==candidate["views"][v]["proxy_paths_mm"] for v in paths)
    source=ifcopenshell.open(str(OLD)); formal=ifcopenshell.open(str(FORMAL)); target=source.by_guid(GUID)
    assert pkg.body_fingerprint(target)==pkg.body_fingerprint(formal.by_guid(GUID))
    assert np.array_equal(pkg.placement(target),pkg.placement(formal.by_guid(GUID)))
    from align_spatial_metadata import align
    align(source,formal)
    pre=snapshot(source)
    protected=[pkg.record(p) for p in [FORMAL, OLD, APPROVAL, PRODUCT/"manifest.json", candidate_path,
               *[p for p in (PRODUCT/"official-source").rglob("*") if p.is_file()],
               *[PRODUCT/f"{v}.svg" for v in paths]]]
    definitions={"plan":"FFL PLAN", "front":"EL-01-02-R20-PY", "side":"EL-01-03-R20-PX"}
    includes=eu.get_pset(next(d for d in formal.by_type("IfcAnnotation") if d.Representation and d.Name==definitions["front"]),"EPset_Drawing")["Include"].split(",")
    specs={}; evidence=[]
    for view,name in definitions.items():
        original=next(d for d in formal.by_type("IfcAnnotation") if d.Representation and d.Name==name)
        # Copy only the camera's forward geometry/placement. New root identities
        # keep formal Drawings and their groups outside this pending recipe.
        copier=pkg.ScopedCopy(formal,source,{original},skip_inverse_ids=[original.id()])
        drawing=copier.copy(original)
        drawing.GlobalId=ifcopenshell.guid.new();drawing.Name=f"TAB02-SCENE-{view.upper()}"
        drawing.Description="Pending scene validation; approved single-product outline only"
        values={k:v for k,v in eu.get_pset(original,"EPset_Drawing").items() if k not in ("id","Exclude","Include")}
        values.update(HasAnnotation=True, HasUnderlay=False, GlobalReferencing=False,
                      Include=",".join(g for g in includes if g!=GUID),
                      Exclude=",".join(['IfcSpace','IfcGrid','IfcBuildingStorey',*[a.GlobalId for a in formal.by_type('IfcAnnotation')]]))
        for key in ("Stylesheet","Markers","Symbols","Patterns","ShadingStyles"):
            values[key]=str(ROOT/values[key])
        block=next(x for x in source.traverse(drawing.Representation) if x.is_a("IfcBlock"))
        if view=="plan":
            drawing.ObjectPlacement.RelativePlacement.Location.Coordinates=(1750.0,3900.0,1700.0)
            block.XLength=3200.0;block.YLength=2400.0;block.ZLength=1900.0
        elif view=="front":
            # Look along +Y from ahead of the table, with world Z up.
            drawing.ObjectPlacement.RelativePlacement.Location.Coordinates=(1750.0,3600.0,650.0)
            block.XLength=2200.0;block.YLength=1400.0;block.ZLength=1600.0
        else:
            drawing.ObjectPlacement.RelativePlacement.Location.Coordinates=(900.0,3950.0,650.0)
            block.XLength=1600.0;block.YLength=1400.0;block.ZLength=2300.0
        block.Position.Location.Coordinates=(-block.XLength/2,-block.YLength/2,-block.ZLength)
        pset(source,drawing,"EPset_Drawing",values)
        document=source.create_entity("IfcDocumentReference",Location=str(OUT/f"{drawing.Name}.svg"),Identification=drawing.Name,Name=drawing.Name)
        source.create_entity("IfcRelAssociatesDocument",GlobalId=ifcopenshell.guid.new(),RelatedObjects=[drawing],RelatingDocument=document)
        context=approved.representation_context(source,"Annotation",values["TargetView"])
        rep=curve_representation(source,context,view,paths[view])
        color=source.create_entity("IfcColourRgb",Red=17/255,Green=24/255,Blue=32/255)
        style=source.create_entity("IfcCurveStyle",Name="TAB02 approved semantic black",CurveWidth=source.create_entity("IfcPositiveLengthMeasure",0.35),CurveColour=color,ModelOrDraughting=True)
        source.create_entity("IfcStyledItem",Item=rep.Items[0],Styles=[style])
        annotation=source.create_entity("IfcAnnotation",GlobalId=ifcopenshell.guid.new(),Name=f"TAB02 approved {view}",ObjectType="LINEWORK",ObjectPlacement=target.ObjectPlacement,Representation=source.create_entity("IfcProductDefinitionShape",Representations=[rep]))
        pset(source,annotation,"EPset_Annotation",{"Classes":"review-target-tab02 geometry-derived-simplified-proxy","TargetGlobalId":GUID,"SourceKind":"geometry_derived_simplified_proxy"})
        group=source.create_entity("IfcGroup",GlobalId=ifcopenshell.guid.new(),Name=drawing.Name,ObjectType="DRAWING")
        source.create_entity("IfcRelAssignsToGroup",GlobalId=ifcopenshell.guid.new(),RelatedObjects=[drawing,annotation],RelatingGroup=group)
        specs[view]={"annotation_guid":annotation.GlobalId,"drawing_guid":drawing.GlobalId}
        evidence.append({"view":view,**specs[view],"source_camera_guid":original.GlobalId,"source_camera_name":name,"camera_matrix":pkg.placement(drawing).tolist(),"approved_path_count":len(paths[view])})
    pset(source,target,"Pset_Tab02ApprovedSource",{"SourceKind":"geometry_derived_simplified_proxy","SourceLabelZh":SOURCE_LABEL,"OfficialCadStatus":"reserved_area_authentication_required_not_acquired","OfficialReferenceUse":"manufacturer identity and catalogue dimensions only","OfficialCadUsed":False,"SingleProductApprovalStatus":"approved","SceneApprovalStatus":"pending","FormalIfcWriteAllowed":False})
    pure,recipe,audit=build_pure_package(source,GUID,specs)
    recipe['annotation_import_adapter']={'script':'migrate_tab02.py','function':'exact_annotation_mesh','run_before':'each bpy.ops.bim.create_drawing after activate_drawing','geometry_source':'persisted pure IFC ApprovedPlan/ApprovedFront/ApprovedSide polylines only','strategy':'create distinct Blender mesh vertices for every stored polyline point and one edge for every successive pair; no merging or simplification','scope':'temporary Blender annotation meshes only; IFC representations and placement unchanged','installed_api':'Blender 4.5.3 LTS; Bonsai and IfcOpenShell 0.8.4; SvgWriter.draw_misc_annotation and draw_edge_annotation emit LINEWORK mesh edges','historical_cache_used':False}
    audit['unused_type_box_maps_removed']=discard_type_box(pure)
    pure.write(str(SINGLE));write(RECIPE,recipe)
    state=check_pure(ifcopenshell.open(str(SINGLE)))
    deps=[]
    for location in sorted({x.Location for x in pure.by_type("IfcExternallyDefinedSurfaceStyle") if x.Location}):
        rel=Path(location);assert not rel.is_absolute() and ".." not in rel.parts
        dest=OUT/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,dest)
        deps.append({"source":pkg.record(ROOT/rel),"packaged":pkg.record(dest),"location":location})
    write(REPORT,{"product":"TAB02","target_guid":GUID,"verdict":"pending_bonsai_scene","formal_sha256_before":FORMAL_HASH,
        "protected_files":protected,"preState":pre,"pure_pre_bonsai":state,"extraction":audit,"views":evidence,
        "single_product":pkg.record(SINGLE),"scene_recipe":pkg.record(RECIPE),"necessary_external_style_dependencies":deps,
        "single_product_approval_status":"approved","scene_approval_status":"pending","formal_write_allowed":False,"cleanup_performed":False,
        "source_kind":"geometry_derived_simplified_proxy","source_label_zh":SOURCE_LABEL,
        "scene_camera_policy":"Pending living-room context cameras; Plan 1700 mm above 680 mm product and below ceiling; Front from -Y, Side from -X. World Z remains up.",
        "courseEvidence":{"mode":"embedded-course-index","lesson":"085000","timestamp":"01:59 Create Drawing","screenshot":"research/bonsai-course/lessons/085000/screenshots/085000-01m59s-create-drawing-button.png","provenance_sha256":"6997420afc8e3d91d0016b8c68c7db904cc782781317a0de7ff50e9d10f945dc","label":"course_fact"}})
    print(json.dumps({"pure":pkg.record(SINGLE),"state":state}))

def style_svg(svg, annotation_guid, view):
    import xml.etree.ElementTree as ET
    raw=pkg.record(svg);tree=ET.parse(svg);root=tree.getroot();count=0
    css=ET.SubElement(root,'{http://www.w3.org/2000/svg}style',type='text/css')
    css.text='.projection{fill:none;stroke:#aaa;stroke-width:0.06}.cut{fill:#eee;stroke:#999;stroke-width:0.10}'
    for e in root.iter():
        if annotation_guid in e.get("class","") and e.tag.rsplit("}",1)[-1] in ("line","polyline","path"):
            e.set("style","stroke:#111820;stroke-width:0.18;fill:none");e.set("data-target-global-id",GUID);count+=1
    assert count>0,"No TAB02 generated annotation geometry"
    root.set("data-create-drawing-result","FINISHED");root.set("data-scene-approval-status","pending")
    tree.write(svg,encoding="utf-8",xml_declaration=True)
    return {"raw":raw,"approved_black_geometry_count":count,"post_style_only":True,'geometry_moved_removed_or_redrawn':False}

def exact_annotation_mesh(model, pure, spec, view):
    """Avoid tessellator vertex merging: the saved pure IFC owns every point.

    Bonsai SVG writer draws LINEWORK from its Blender mesh. Rebuild only that
    transient mesh from exact persisted polylines; no IFC mutation or redraw.
    """
    import bpy
    from bonsai import tool
    annotation=model.by_guid(spec['annotation_guid'])
    obj=tool.Ifc.get_object(annotation)
    assert obj is not None
    rep=next(x for x in pure.by_guid(GUID).Representation.Representations if x.RepresentationIdentifier=='Approved'+view.title())
    vertices=[];edges=[];scale=pkg.unit_util.calculate_unit_scale(pure)
    for line in rep.Items[0].Elements:
        offset=len(vertices)
        vertices.extend([tuple(float(c)*scale for c in point.Coordinates) for point in line.Points])
        edges.extend((offset+i,offset+i+1) for i in range(len(line.Points)-1))
    mesh=bpy.data.meshes.new('TAB02 exact persisted '+view)
    mesh.from_pydata(vertices,edges,[]);mesh.update();obj.data=mesh
    return {'view':view,'vertices':len(vertices),'edges':len(edges),'source':'pure IFC Approved polylines','ifc_modified':False,'reason':'Avoid default IfcOpenShell tessellation merging near-coincident annotation vertices'}

def scene():
    import bpy,bonsai_bridge as bridge
    from bonsai import tool
    import create_wd03_wardrobe_scene_drawings as context
    report=json.loads(REPORT.read_text())
    if report.get('temporary_project_directory'):
        report.setdefault('previous_temporary_attempts',[]).append(report['temporary_project_directory'])
    assert Path(tool.Ifc.get_path()).resolve()==SINGLE.resolve() and not bpy.data.is_saved
    temporary=Path(tempfile.mkdtemp(prefix="tab02-pure-package-scene-"));temp_ifc=temporary/"scene.ifc"
    report.update(temporary_project_directory=str(temporary),runtime_inputs=[str(SINGLE),str(RECIPE),str(FORMAL)],legacy_ifc_used_for_runtime=False)
    write(REPORT,report)
    try:
        report["provider"]={"status":"supported","version":list(bridge.bl_info["version"]),"port":9897,"pid":__import__('os').getpid(),"blender":bpy.app.version_string,"ifcopenshell":ifcopenshell.version,"bridge_source":pkg.record(bridge.__file__)}
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
        from repair_runtime_metadata import repair
        report['runtime_metadata_reuse']=repair(model,ifcopenshell.open(str(FORMAL)))
        model.write(str(temp_ifc))
        assert bpy.ops.bim.load_project(filepath=str(temp_ifc),should_start_fresh_session=False,use_relative_path=False)=={"FINISHED"}
        outputs=[]
        for v in report["views"]:
            entity=tool.Ifc.get().by_guid(v["drawing_guid"])
            tool.Ifc.get_object(entity) or tool.Drawing.import_drawing(entity)
            with bpy.context.temp_override(**context.view3d_override()):
                assert bpy.ops.bim.activate_drawing(drawing=entity.id(),should_view_from_camera=False)=={"FINISHED"}
                report.setdefault('exact_annotation_mesh',[]).append(exact_annotation_mesh(tool.Ifc.get(),pure,v,v['view']))
                props=tool.Drawing.get_document_props();props.should_use_underlay_cache=False;props.should_use_linework_cache=False;props.should_use_annotation_cache=False
                result=bpy.ops.bim.create_drawing(print_all=False,open_viewer=False,sync=False)
            assert result=={"FINISHED"}
            svg=OUT/f"TAB02-SCENE-{v['view'].upper()}.svg"
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
