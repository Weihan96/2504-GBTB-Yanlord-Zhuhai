"""TRAP01 approved v4 -> pure product package; no formal/legacy writes.

Preparation may read the approved legacy source. Scene generation must only
read the new product, external recipe, and a temporary copy of the formal IFC.
"""
from pathlib import Path
import json
import sys
import shutil
import tempfile
import traceback
import xml.etree.ElementTree as ET
import numpy as np
import ifcopenshell
import ifcopenshell.util.element as eu

OUT = Path(__file__).resolve().parent
PRODUCT = OUT.parent
ROOT = OUT.parents[4]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
import review_product_package as pkg
import pure_product_package as pure_pkg

GUID = "2Ak2ma0lvBEA49UpplzUqi"
SINGLE = OUT / "TRAP01-product.ifc"
RECIPE = OUT / "scene-recipe.json"
REPORT = OUT / "validation.json"
FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_SHA = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OLD = PRODUCT / "Geberit-151.116.11.1-TRAP01-official-detail-v4.ifc"
OLD_SHA = "e0f5955d96f7d4bf73dd9e34f9da05b634262def0f32ca5c76b3976f3eb99dbb"
APPROVAL = ROOT / "pipeline/decisions/trap01-drawing-approval.json"
OLD_MANIFEST = PRODUCT / "TRAP01-official-detail-v4-manifest.json"
OLD_EVIDENCE = PRODUCT / "bonsai-drawings/cabinet-official-detail-v4/TRAP01-official-detail-v4-evidence.json"
VIEWS = {
    "plan": {"annotation_guid": "2wMANTAEP4WRRf8Kz4nQhG", "drawing_guid": "3EZt1DKnH8sQoOeH8BbH5y"},
    "front": {"annotation_guid": "1eBrRf2_5BN9iGp9_AkRa2", "drawing_guid": "0oZNWYRa1BxwCLs0yWOQAj"},
    "side": {"annotation_guid": "3KI_E9ncXERQwoWAJcOxeI", "drawing_guid": "1YeH3d8DzB4hLiD2Rj4Zua"},
}


def write(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def source_audit():
    assert pkg.sha256(FORMAL) == FORMAL_SHA
    assert pkg.sha256(OLD) == OLD_SHA
    approval = json.loads(APPROVAL.read_text())
    assert approval["latest_scene_review"]["approved_version"] == "official-detail-v4"
    assert approval["scene_svg_approved"] and approval["derived_ifc_write_allowed"]
    assert not approval["formal_authoritative_ifc_write_allowed"]
    evidence = json.loads(OLD_EVIDENCE.read_text())
    source = ifcopenshell.open(str(OLD))
    target = source.by_guid(GUID)
    scope = [FORMAL, OLD, APPROVAL, OLD_MANIFEST, OLD_EVIDENCE]
    scope += [p for p in (PRODUCT / "official-source").rglob("*") if p.is_file()]
    audit = {
        "schema_version": 1,
        "source": pkg.record(OLD),
        "approval": pkg.record(APPROVAL),
        "formal_before": pkg.record(FORMAL),
        "protected_files": [pkg.record(p) for p in scope],
        "target_guid": GUID,
        "target_body_fingerprint": pkg.body_fingerprint(target),
        "target_placement": pkg.placement(target).tolist(),
        "body_representation_count": sum(r.RepresentationIdentifier == "Body" for r in target.Representation.Representations),
        "superseded_target_proxy_representation_identifiers": [r.RepresentationIdentifier for r in target.Representation.Representations if r.RepresentationIdentifier != "Body"],
        "legacy_target_psets": eu.get_psets(target),
        "source_view_mapping": evidence["view_mapping"],
        "approved_scene_views": json.loads(OLD_MANIFEST.read_text())["views"],
        "adjustable_components": evidence["adjustable_components"],
        "adjustable_official_source_audit": evidence["official_linework_audit"],
        "line_semantics": evidence["line_semantics"],
        "blue_dashed_reference_policy": "Source v4 IFC contains no dashed reference annotation. Reference excess paths and adjustment parameters stay in this external audit; approved installed scenes omit them.",
        "legacy_property_discrepancy": "Target Pset_Trap01DrawingSource describes superseded proxy reps and OfficialCadUsed=false. Approved v4 LINEWORK psets describe actual configured official DWG. Only the new pure package metadata is corrected; old source and approval stay unchanged.",
        "courseEvidence": {"mode": "embedded-course-index", "lesson": "085000", "timestamps": ["01:59 Create Drawing", "02:13 SVG"], "private_screenshot_observed": False},
        "cleanup_performed": False,
    }
    write(OUT / "source-audit.json", audit)
    print(json.dumps({"source_verified": True, "body_representations": audit["body_representation_count"], "official_path_counts": {v["view"]: v["persisted_path_count_expected"] for v in evidence["outputs"]["views"]}}))
    return source, audit


def prepare():
    source, audit = source_audit()
    if SINGLE.exists():
        assert "--rebuild-owned" in sys.argv, "Do not overwrite a previously generated pure package"
        prior = json.loads(REPORT.read_text())
        assert pkg.sha256(SINGLE) == prior["pure_product"]["sha256"], "Existing product changed outside this adapter"
    pure, recipe, extraction = pure_pkg.build_pure_package(source, GUID, VIEWS)
    import ifcopenshell.api.pset
    target = pure.by_guid(GUID)
    old_psets = eu.get_psets(target)
    if "Pset_Trap01DrawingSource" in old_psets:
        pset = pure.by_id(old_psets["Pset_Trap01DrawingSource"]["id"])
        ifcopenshell.api.pset.edit_pset(pure, pset=pset, properties={
            "SourceKind": "native_dwg_configured_installation",
            "SourceLabelZh": "基于官方原生DWG固定接头刚性移位及直管长度调整的二维安装细节表达",
            "OfficialCadUsed": "true", "OfficialCadGeometryIncluded": "true",
            "OfficialCadRole": "Approved v4 configured full 2D detail; fixed joints rigidly translated, straight adjustable pipes shortened; not a project shop drawing",
            "RepresentationGeometrySource": "approved official-detail-v4 LINEWORK transferred without geometric edits",
            "ApprovedRepresentationIdentifiers": "ApprovedPlan;ApprovedFront;ApprovedSide",
            "ProjectConfiguration": "shortened installed configuration; v4 complete native fixed-joint details",
            "ReviewDate": "2026-09-07", "ApprovalEvidence": "这两个也都验收通过。 official-detail-v4 scene approval",
            "PlanPathCount": "95", "FrontPathCount": "117", "SidePathCount": "77",
        })
    # Configuration references to legacy DRAWING_COMPONENT entities become
    # external dependencies, not dangling in-file identities.
    if "Pset_Trap01AdjustableConfiguration" in old_psets:
        pset = pure.by_id(old_psets["Pset_Trap01AdjustableConfiguration"]["id"])
        props = {k: None for k in ("ComponentGroupGlobalId", "FixedBodyGlobalId", "HorizontalAdjustableGlobalId", "VerticalAdjustableGlobalId")}
        props["AdjustmentRecipe"] = "source-audit.json#adjustable_components"
        ifcopenshell.api.pset.edit_pset(pure, pset=pset, properties=props)
    repairs = pkg.repair_missing_metadata(pure)
    pure.write(str(SINGLE))
    reloaded = ifcopenshell.open(str(SINGLE))
    assert pkg.body_fingerprint(reloaded.by_guid(GUID)) == audit["target_body_fingerprint"]
    assert np.array_equal(pkg.placement(reloaded.by_guid(GUID)), np.array(audit["target_placement"]))
    assert not reloaded.by_type("IfcAnnotation")
    recipe["trap01_line_semantics"] = audit["line_semantics"]
    recipe["trap01_blue_dashed_reference_policy"] = audit["blue_dashed_reference_policy"]
    write(RECIPE, recipe)
    write(REPORT, {"stage": "pure_package_prepared", "verdict": "pending_bonsai_roundtrip_and_scene", "source_audit": "source-audit.json", "extraction": extraction, "metadata_repairs": repairs, "pure_product": pkg.record(SINGLE), "recipe": pkg.record(RECIPE), "formal_write_allowed": False, "legacy_copy_retained": True})
    print(json.dumps({"pure_product": pkg.record(SINGLE), "recipe": pkg.record(RECIPE)}))


def package_scene_assets():
    """Migration-time copying only; runtime never reads legacy drawing assets."""
    destination = OUT / "assets"
    destination.mkdir(exist_ok=True)
    records = []
    for source in sorted((PRODUCT / "drawings/assets").iterdir()):
        if source.is_file():
            dest = destination / source.name
            shutil.copy2(source, dest)
            assert pkg.sha256(dest) == pkg.sha256(source)
            records.append({"source": pkg.record(source), "package": pkg.record(dest)})
    write(OUT / "assets-manifest.json", records)


def view3d_override():
    import bpy
    for window in bpy.context.window_manager.windows:
        for area in window.screen.areas:
            if area.type == "VIEW_3D":
                region = next((r for r in area.regions if r.type == "WINDOW"), None)
                if region:
                    return {"window": window, "screen": window.screen, "area": area, "region": region, "scene": bpy.context.scene}
    raise RuntimeError("Bonsai Drawing requires a real VIEW_3D area")


def pure_state(model):
    target = model.by_guid(GUID)
    assert len(model.by_type("IfcElement")) == 1
    assert not model.by_type("IfcAnnotation") and not model.by_type("IfcGroup")
    assert not [p for p in model.by_type("IfcPropertySet") if p.Name == "EPset_Drawing"]
    assert not [p for p in model.by_type("IfcPropertySingleValue") if p.Name in ("Include", "Exclude")]
    assert all(not e.Representation for e in model.by_type("IfcSpatialElement"))
    reps = target.Representation.Representations
    assert [r.RepresentationIdentifier for r in reps] == ["Body", "Body", "ApprovedPlan", "ApprovedFront", "ApprovedSide"]
    return {"body_fingerprint": pkg.body_fingerprint(target), "placement": pkg.placement(target).tolist(),
            "view_content": {r.RepresentationIdentifier: pure_pkg.representation_content(r) for r in reps if r.RepresentationIdentifier != "Body"},
            "represented_elements": 1, "annotation_count": 0, "group_count": 0, "spatial_geometry_count": 0,
            "scene_camera_and_filter_count": 0}


def style_inspect_svg(path, annotation_guid):
    """Apply the approved v4 presentation, without moving any geometry."""
    namespace = "http://www.ifcopenshell.org/ns"
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("ifc", namespace)
    tree = ET.parse(path)
    root = tree.getroot()
    parent = {child: p for p in root.iter() for child in p}
    def identity(e, guid):
        return e.get(f"{{{namespace}}}guid") == guid or f"GlobalId-{guid}" in e.get("class", "").split()
    targets = [e for e in root.iter() if identity(e, annotation_guid)]
    assert targets
    raw = pkg.record(path)
    for e in targets:
        if e in parent:
            parent[e].remove(e)
            parent[e].append(e)
    counts = {"blue": 0, "grey": 0}
    geometry = {"path", "polyline", "polygon", "line", "circle", "ellipse", "rect"}
    for e in root.iter():
        if e.tag.split("}")[-1] in geometry:
            blue = identity(e, annotation_guid)
            style = "stroke:#1677c8;stroke-width:0.45;fill:none" if blue else "stroke:#a3abb3;stroke-width:0.20;fill:none;stroke-opacity:0.62"
            e.set("style", e.get("style", "").rstrip(";") + ";" + style)
            counts["blue" if blue else "grey"] += 1
    assert counts["blue"] > 0 and counts["grey"] > 0
    assert not any(identity(e, GUID) for e in root.iter()), "Duplicate target Body leaked through installed-detail filter"
    occluders = ["1FgLPMw$5B4wBH2ySMkXE1", "2YnGmmoYbF28Vy25$aNKT8", "3OVQygdDn17huGOgJJFTOY"]
    assert not any(identity(e, g) for e in root.iter() for g in occluders)
    assert not any(e.tag.endswith("image") for e in root.iter())
    root.set("data-blue-dashed-reference-displayed", "false")
    root.set("data-trap01-source", "pure-product-ifc-plus-external-scene-recipe")
    root.set("data-create-drawing-result", "FINISHED")
    tree.write(path, encoding="utf-8", xml_declaration=True)
    return {"raw_svg": raw, "presentation_only": True, "geometry_counts": counts,
            "target_body_absent": True, "all_three_cabinet_occluders_absent": True,
            "blue_dashed_reference_absent": True, "scale": root.get("data-scale"), "view_box": root.get("viewBox")}


def scene():
    import bpy
    import bonsai_bridge as bridge
    from bonsai import tool
    import ifcopenshell.api.pset
    report = json.loads(REPORT.read_text())
    assert Path(tool.Ifc.get_path()).resolve() == SINGLE.resolve()
    assert not bpy.data.is_saved
    assert pkg.sha256(FORMAL) == FORMAL_SHA
    assert bridge.bl_info["version"] == (1, 1, 0)
    recipe = json.loads(RECIPE.read_text())
    temporary = Path(tempfile.mkdtemp(prefix="trap01-pure-scene-"))
    temp_ifc = temporary / "scene.ifc"
    report.update({"stage": "scene_running", "temporary_directory": str(temporary), "temporary_ifc": str(temp_ifc),
                   "runtime_inputs": [str(SINGLE), str(RECIPE), str(FORMAL)],
                   "runtime_legacy_ifc_read": False, "formal_sha256_before_scene": pkg.sha256(FORMAL),
                   "provider": {"version": list(bridge.bl_info["version"]), "blender": bpy.app.version_string, "ifcopenshell": ifcopenshell.version, "bridge_file": bridge.__file__, "bridge_sha256": pkg.sha256(bridge.__file__), "port": 9884, "pid": __import__('os').getpid()}})
    write(REPORT, report)
    try:
        before = pure_state(ifcopenshell.open(str(SINGLE)))
        with bpy.context.temp_override(**view3d_override()):
            report["product_bonsai_save_reload"] = bridge._h_save_ifc_file({"output_path": str(SINGLE), "overwrite": True, "reload": True})
        pure = ifcopenshell.open(str(SINGLE))
        after = pure_state(pure)
        assert before == after, "Bonsai pure product roundtrip drift"
        assert tool.Ifc.get_object(tool.Ifc.get().by_guid(GUID)) is not None
        report["pure_roundtrip"] = after
        # Only a disposable project copy is saved or loaded for scene work.
        shutil.copy2(FORMAL, temp_ifc)
        project = ifcopenshell.open(str(temp_ifc))
        resources = []
        for location in sorted({s.Location for s in project.by_type("IfcExternallyDefinedSurfaceStyle") if s.Location}):
            path = Path(location)
            assert not path.is_absolute() and ".." not in path.parts
            src = ROOT / path
            if src.is_file():
                dst = temporary / path
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
                resources.append({"source": pkg.record(src), "temporary": str(dst)})
            else:
                resources.append({"path": location, "status": "already_missing_at_formal_source"})
        report["temporary_resources"] = resources
        body = {e.GlobalId: pkg.body_fingerprint(e) for e in project.by_type("IfcElement") if e.Representation}
        report["attachment"] = pure_pkg.attach_recipe(project, pure, recipe)
        count = len(list(project))
        pure_pkg.attach_recipe(project, pure, recipe)
        assert len(list(project)) == count, "Second attach duplicated entities"
        report["attachment"]["second_attachment_created_entities"] = 0
        expected_annotations = {}
        for view, spec in recipe["views"].items():
            drawing = project.by_guid(spec["drawing_guid"])
            ann = project.by_guid(spec["annotation_guid"])
            pure_rep = next(r for r in pure.by_guid(GUID).Representation.Representations if r.RepresentationIdentifier == f"Approved{view.title()}")
            assert pure_pkg.representation_content(ann.Representation.Representations[0]) == pure_pkg.representation_content(pure_rep)
            expected_annotations[view] = {"geometry": pure_pkg.representation_content(ann.Representation.Representations[0]), "placement": pkg.placement(ann).tolist()}
            refs = [r.RelatingDocument for r in drawing.HasAssociations if r.is_a("IfcRelAssociatesDocument")]
            assert len(refs) == 1
            refs[0].Location = str(OUT / f"TRAP01-SCENE-{view.upper()}.svg")
            pset = project.by_id(eu.get_psets(drawing)["EPset_Drawing"]["id"])
            properties = {k: str(OUT / "assets" / filename) for k, filename in {
                "Stylesheet": "default.css", "Markers": "markers.svg", "Symbols": "symbols.svg", "Patterns": "patterns.svg", "ShadingStyles": "shading_styles.json"}.items()}
            ifcopenshell.api.pset.edit_pset(project, pset=pset, properties=properties)
        assert body == {e.GlobalId: pkg.body_fingerprint(e) for e in project.by_type("IfcElement") if e.Representation}
        project.write(str(temp_ifc))
        bridge._reload_ifc_project(str(temp_ifc))
        assert Path(tool.Ifc.get_path()).resolve() == temp_ifc.resolve()
        outputs = []
        for view, spec in recipe["views"].items():
            drawing = tool.Ifc.get().by_guid(spec["drawing_guid"])
            tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
            with bpy.context.temp_override(**view3d_override()):
                assert bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False) == {"FINISHED"}
                props = tool.Drawing.get_document_props()
                props.should_use_underlay_cache = False
                props.should_use_linework_cache = False
                props.should_use_annotation_cache = False
                result = bpy.ops.bim.create_drawing(print_all=False, open_viewer=False, sync=False)
            assert result == {"FINISHED"}
            svg = OUT / f"TRAP01-SCENE-{view.upper()}.svg"
            inspection = style_inspect_svg(svg, spec["annotation_guid"])
            outputs.append({"view": view, **spec, "generator": "bpy.ops.bim.create_drawing", "create_result": sorted(result), "inspection": inspection, "svg": pkg.record(svg)})
        with bpy.context.temp_override(**view3d_override()):
            report["temporary_bonsai_save_reload"] = bridge._h_save_ifc_file({"output_path": str(temp_ifc), "overwrite": True, "reload": True})
        reloaded = tool.Ifc.get()
        assert body == {e.GlobalId: pkg.body_fingerprint(e) for e in reloaded.by_type("IfcElement") if e.Representation}
        assert len([e for e in reloaded.by_type("IfcElement") if e.GlobalId == GUID]) == 1
        assert np.array_equal(pkg.placement(reloaded.by_guid(GUID)), pkg.placement(pure.by_guid(GUID)))
        for view, spec in recipe["views"].items():
            ann = reloaded.by_guid(spec["annotation_guid"])
            assert expected_annotations[view] == {"geometry": pure_pkg.representation_content(ann.Representation.Representations[0]), "placement": pkg.placement(ann).tolist()}
        report.update({"scene_outputs": outputs, "temporary_project": pkg.record(temp_ifc), "scene_reload_verified": True, "all_formal_bodies_unchanged": True, "stage": "bonsai_roundtrip_and_scene_complete", "verdict": "pending_independent_svg_and_visual_verification"})
        bridge._reload_ifc_project(str(SINGLE))
        assert Path(tool.Ifc.get_path()).resolve() == SINGLE.resolve()
    except Exception:
        report["verdict"] = "fail"
        report["error"] = traceback.format_exc()
        raise
    finally:
        report["formal_sha256_after_scene"] = pkg.sha256(FORMAL)
        report["pure_product"] = pkg.record(SINGLE)
        write(REPORT, report)
        assert report["formal_sha256_after_scene"] == FORMAL_SHA


def capture_pure_viewport():
    import bpy
    import bonsai_bridge as bridge
    import base64
    from bonsai import tool
    assert Path(tool.Ifc.get_path()).resolve() == SINGLE.resolve()
    target = tool.Ifc.get_object(tool.Ifc.get().by_guid(GUID))
    with bpy.context.temp_override(**view3d_override()):
        bpy.ops.object.select_all(action="DESELECT")
        target.select_set(True)
        bpy.context.view_layer.objects.active = target
        # Use orthographic projection and a measured Body diagonal to avoid
        # bridge's perspective auto-fit clipping this small plumbing product.
        from mathutils import Vector
        bridge._h_get_viewport_screenshot({"view": "iso", "fit": "selected", "max_size": 160})
        region = bpy.context.area.spaces.active.region_3d
        region.view_perspective = "ORTHO"
        region.view_location = sum((target.matrix_world @ Vector(p) for p in target.bound_box), Vector()) / 8
        region.view_distance = max(target.dimensions) * 2.5
        result = bridge._h_get_viewport_screenshot({"max_size": 1400, "format": "png", "include_objects": True})
    image_data = result.pop("image_base64")
    suffix = ".png" if result.get("format") == "png" else ".jpg"
    path = OUT / ("TRAP01-pure-bonsai-3d" + suffix)
    path.write_bytes(base64.b64decode(image_data))
    result["image"] = pkg.record(path)
    result["active_ifc"] = str(SINGLE)
    result["ifc_sha256"] = pkg.sha256(SINGLE)
    result["blend_saved"] = bpy.data.is_saved
    write(OUT / "bonsai-product-viewport.json", result)
    print(json.dumps({"image": str(path), "ifc": str(SINGLE)}))


if __name__ == "__main__":
    if "--audit-only" in sys.argv:
        source_audit()
    else:
        prepare()
