"""GESSI54294 approved-source extraction and independent temporary-scene proof.

prepare reads the legacy source once. scene uses only the new IFC, its JSON
recipe and the formal baseline; it must never open the legacy IFC.
"""
from pathlib import Path
import json
import sys
import tempfile
import shutil
import traceback
import numpy as np
import ifcopenshell
import ifcopenshell.util.element as element_util

OUT = Path(__file__).resolve().parent
PRODUCT = OUT.parent
ROOT = OUT.parents[4]
sys.path.insert(0, str(ROOT / "pipeline/scripts"))
import review_product_package as pkg

FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_HASH = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
APPROVAL = ROOT / "pipeline/decisions/gessi316-54294-drawing-approval.json"
OLD_MANIFEST = PRODUCT / "main-bathroom-bonsai-drawing-manifest.json"
SINGLE = OUT / "GESSI54294-product.ifc"
RECIPE = OUT / "scene-recipe.json"
REPORT = OUT / "validation.json"
GUID = "2iKOL78$H0N9Yd9$ky3pW4"


def write(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def snapshot(model):
    target = model.by_guid(GUID)
    return {"body": pkg.body_fingerprint(target),
            "representations": pkg.fingerprint(target.Representation),
            "placement": pkg.placement(target).tolist(),
            "physical_elements": sorted(e.GlobalId for e in model.by_type("IfcElement")),
            "annotation_count": len(model.by_type("IfcAnnotation")),
            "drawing_pset_count": len([p for p in model.by_type("IfcPropertySet") if p.Name == "EPset_Drawing"]),
            "unit_scale_to_m": pkg.unit_util.calculate_unit_scale(model)}


def check_pure(model):
    result = snapshot(model)
    assert result["physical_elements"] == [GUID]
    assert result["annotation_count"] == result["drawing_pset_count"] == 0
    assert not [s for s in model.by_type("IfcSpatialElement") if s.Representation]
    roots = model.by_type("IfcRoot")
    assert len(roots) == len({e.GlobalId for e in roots})
    return result


def package_resources():
    """Keep referenced material library, not a saved complete-project session."""
    model = ifcopenshell.open(str(SINGLE))
    dependencies = []
    for location in sorted({s.Location for s in model.by_type("IfcExternallyDefinedSurfaceStyle") if s.Location}):
        relative = Path(location)
        assert not relative.is_absolute() and ".." not in relative.parts
        source = ROOT / relative
        if not source.is_file():
            raise ValueError(f"Pure product requires unavailable material dependency: {source}")
        dest = OUT / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            assert pkg.sha256(dest) == pkg.sha256(source)
        else:
            shutil.copy2(source, dest)
        dependencies.append({"location": location, "source": pkg.record(source), "packaged": pkg.record(dest),
                             "identifications": [s.Identification for s in model.by_type("IfcExternallyDefinedSurfaceStyle") if s.Location == location]})
    report = json.loads(REPORT.read_text())
    report["necessary_external_style_dependencies"] = dependencies
    write(REPORT, report)
    return dependencies


def prepare():
    from pure_product_package import build_pure_package, spatial_scope
    approval = json.loads(APPROVAL.read_text())
    assert approval["status"] == "approved"
    assert pkg.sha256(FORMAL) == FORMAL_HASH
    m = json.loads(OLD_MANIFEST.read_text())
    assert pkg.sha256(OLD_MANIFEST) == approval["approved_scene_svgs"]["manifest_sha256"]
    combined = None
    protected = [pkg.record(p) for p in [FORMAL, APPROVAL, OLD_MANIFEST,
        ROOT / "pipeline/decisions/highpoly-subtask-completion-2026-09-05.json"]]
    specs = {}
    for v in m["views"]:
        path = ROOT / v["drawing_session_ifc"]["path"]
        assert pkg.sha256(path) == v["drawing_session_ifc"]["sha256"]
        source = ifcopenshell.open(str(path))
        protected += [pkg.record(path), pkg.record(ROOT / v["svg"]["path"])]
        spec = {"annotation_guid": v["linework_annotation_global_id"], "drawing_guid": v["drawing_global_id"]}
        specs[v["view"]] = spec
        if combined is None:
            combined = source
        else:
            target = source.by_guid(GUID)
            assert pkg.body_fingerprint(target) == pkg.body_fingerprint(combined.by_guid(GUID))
            annotations = [source.by_guid(g) for g in spec.values()]
            scope = spatial_scope(source, target) | set(annotations)
            for e in annotations:
                scope.update(r.RelatingGroup for r in e.HasAssignments if r.is_a("IfcRelAssignsToGroup"))
            clone = pkg.ScopedCopy(source, combined, scope, reuse_roots=True, skip_inverse_ids=[target.id()])
            clone.memo[target.id()] = combined.by_guid(GUID)
            clone.memo[target.ObjectPlacement.id()] = combined.by_guid(GUID).ObjectPlacement
            for e in annotations:
                clone.copy(e)
            for e in scope:
                if e.is_a("IfcGroup"):
                    clone.copy(e)
            clone.inverses()
    for p in PRODUCT.rglob("*"):
        if p.is_file() and ("official" in str(p.relative_to(PRODUCT)) or p.suffix.lower() == ".dwg") and OUT not in p.parents:
            protected.append(pkg.record(p))
    # Independently saved sessions each carry an identical Bonsai application.
    # Reuse identical metadata in the in-memory merge, preserving IFC UR1/UR2.
    applications = {}
    merged_applications = 0
    for app in list(combined.by_type("IfcApplication")):
        key = pkg.fingerprint(app)
        if key in applications:
            for inverse in combined.get_inverse(app):
                element_util.replace_attribute(inverse, app, applications[key])
            combined.remove(app)
            merged_applications += 1
        else:
            applications[key] = app
    pure, recipe, audit = build_pure_package(combined, GUID, specs)
    audit["identical_application_records_deduplicated"] = merged_applications
    pure.write(str(SINGLE))
    write(RECIPE, recipe)
    state = check_pure(ifcopenshell.open(str(SINGLE)))
    assert state["body"] == pkg.body_fingerprint(combined.by_guid(GUID))
    write(REPORT, {"product":"GESSI54294","target_guid":GUID,"verdict":"pending_bonsai_scene",
        "formal_sha256_before":FORMAL_HASH,"protected_files":protected,
        "source_manifest":pkg.record(OLD_MANIFEST),"legacy_full_projects":[v["drawing_session_ifc"] for v in m["views"]],
        "views":[{"view":v["view"],"drawing_guid":v["drawing_global_id"],"annotation_guid":v["linework_annotation_global_id"],"old_svg":v["svg"]} for v in m["views"]],
        "extraction":audit,"pure_pre_bonsai":state,"single_product":pkg.record(SINGLE),"scene_recipe":pkg.record(RECIPE),
        "formal_write_allowed":False,"cleanup_performed":False,
        "approval_scope":"User authorized storage migration only; previously approved scene-only design preserved",
        "source_kind":m["source_kind"],"source_label_zh":m["source_label_zh"],
        "handle_texture_detail_path_count":0,"approved_side_wall_residual_mm":0.02276718446054815,
        "courseEvidence":{"mode":"embedded-course-index","lesson":"085000","timestamps":["01:59 Create Drawing","02:13 SVG"]}})
    package_resources()
    print(json.dumps({"product":pkg.record(SINGLE),"state":state}))


def scene():
    """No reads of OLD or OLD_MANIFEST: reconstruction is package-only."""
    import bpy
    import bonsai_bridge as bridge
    from bonsai import tool
    from pure_product_package import attach_recipe
    import create_gessi316_54294_main_bathroom_drawing as drawing
    import create_wd03_wardrobe_scene_drawings as context
    report = json.loads(REPORT.read_text())
    if report.get("error"):
        report.setdefault("resolved_attempt_errors", []).append(report.pop("error"))
    if report.get("temporary_project_directory"):
        report.setdefault("previous_temporary_attempts", []).append(report["temporary_project_directory"])
    assert Path(tool.Ifc.get_path()).resolve() == SINGLE.resolve()
    assert pkg.sha256(FORMAL) == FORMAL_HASH
    assert not bpy.data.is_saved
    temporary = Path(tempfile.mkdtemp(prefix="gessi54294-pure-package-scene-"))
    temp_ifc = temporary / "scene.ifc"
    report["temporary_project_directory"] = str(temporary)
    report["runtime_inputs"] = [str(SINGLE), str(RECIPE), str(FORMAL)]
    report["legacy_ifc_used_for_runtime"] = False
    write(REPORT, report)
    try:
        report["provider"] = {"version": list(bridge.bl_info["version"]), "port": 9889,
                              "blender": bpy.app.version_string, "ifcopenshell": ifcopenshell.version,
                              "bridge_source": pkg.record(bridge.__file__)}
        with bpy.context.temp_override(**context.view3d_override()):
            report["product_bonsai_save_reload"] = bridge._h_save_ifc_file(
                {"output_path": str(SINGLE), "overwrite": True, "reload": True})
        pure = ifcopenshell.open(str(SINGLE))
        assert snapshot(pure) == report["pure_pre_bonsai"], "Bonsai package roundtrip drift"
        check_pure(pure)
        report["pure_saved_reloaded"] = True
        shutil.copy2(FORMAL, temp_ifc)
        model = ifcopenshell.open(str(temp_ifc))
        report["formal_physical_elements"] = len(model.by_type("IfcElement"))
        resources = []
        for location in sorted({s.Location for s in model.by_type("IfcExternallyDefinedSurfaceStyle") if s.Location}):
            relative = Path(location)
            assert not relative.is_absolute() and ".." not in relative.parts
            source_resource = ROOT / relative
            if source_resource.is_file():
                dest = temporary / relative
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_resource, dest)
                resources.append(pkg.record(source_resource))
        report["external_resources_copied_to_temporary"] = resources
        baseline_bodies = {e.GlobalId: pkg.body_fingerprint(e) for e in model.by_type("IfcElement") if e.Representation}
        recipe = json.loads(RECIPE.read_text())
        report["attachment"] = attach_recipe(model, pure, recipe)
        first_count = len(list(model))
        attach_recipe(model, pure, recipe)
        assert len(list(model)) == first_count
        report["second_attachment_created_entities"] = 0
        for v in report["views"]:
            entity = model.by_guid(v["drawing_guid"])
            refs = [r.RelatingDocument for r in entity.HasAssociations if r.is_a("IfcRelAssociatesDocument")]
            assert len(refs) == 1
            refs[0].Location = str(OUT / f"GESSI54294-SCENE-{v['view'].upper()}.svg")
        model.write(str(temp_ifc))
        assert bpy.ops.bim.load_project(filepath=str(temp_ifc), should_start_fresh_session=False, use_relative_path=False) == {"FINISHED"}
        outputs = []
        for v in report["views"]:
            entity = tool.Ifc.get().by_guid(v["drawing_guid"])
            tool.Ifc.get_object(entity) or tool.Drawing.import_drawing(entity)
            with bpy.context.temp_override(**context.view3d_override()):
                assert bpy.ops.bim.activate_drawing(drawing=entity.id(), should_view_from_camera=False) == {"FINISHED"}
                props = tool.Drawing.get_document_props()
                props.should_use_underlay_cache = False
                props.should_use_linework_cache = False
                props.should_use_annotation_cache = False
                result = bpy.ops.bim.create_drawing(print_all=False, open_viewer=False, sync=False)
            assert result == {"FINISHED"}
            svg = OUT / f"GESSI54294-SCENE-{v['view'].upper()}.svg"
            styles = drawing.style_review_annotation(svg, v["annotation_guid"], v["view"])
            inspection = drawing.inspect_svg(svg, GUID, v["annotation_guid"])
            outputs.append({"view": v["view"], "svg": pkg.record(svg), "styles": styles,
                            "inspection": inspection, "operator": "bpy.ops.bim.create_drawing", "result": sorted(result)})
        with bpy.context.temp_override(**context.view3d_override()):
            report["temporary_bonsai_save_reload"] = bridge._h_save_ifc_file(
                {"output_path": str(temp_ifc), "overwrite": True, "reload": True})
        reloaded = tool.Ifc.get()
        assert len(reloaded.by_type("IfcElement")) == report["formal_physical_elements"]
        assert baseline_bodies == {e.GlobalId: pkg.body_fingerprint(e) for e in reloaded.by_type("IfcElement") if e.Representation}
        assert np.array_equal(pkg.placement(reloaded.by_guid(GUID)), pkg.placement(pure.by_guid(GUID)))
        assert len([e for e in reloaded.by_type("IfcElement") if e.GlobalId == GUID]) == 1
        report["scene_saved_reloaded"] = True
        report["all_formal_bodies_unchanged"] = True
        report["target_instances"] = 1
        report["scene_outputs"] = outputs
        report["temporary_project"] = pkg.record(temp_ifc)
        report["verdict"] = "pending_independent_visual_validation"
        assert bpy.ops.bim.load_project(filepath=str(SINGLE), should_start_fresh_session=False, use_relative_path=False) == {"FINISHED"}
    except Exception:
        report["error"] = traceback.format_exc()
        report["verdict"] = "fail"
        raise
    finally:
        report["formal_sha256_after"] = pkg.sha256(FORMAL)
        report["single_product"] = pkg.record(SINGLE)
        write(REPORT, report)
        assert report["formal_sha256_after"] == FORMAL_HASH



if __name__ == "__main__":
    prepare()
