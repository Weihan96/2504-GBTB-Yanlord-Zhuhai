"""WD02 pilot: one persistent product IFC, temporary whole-project validation.

prepare runs with IfcOpenShell Python; scene runs in the task-owned Bonsai UI.
No legacy artifact is overwritten and no cleanup is performed by this script.
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
import review_product_package as pkg

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = ROOT / "output/review/highpoly-types/wd02"
OUT = PRODUCT / "product-storage-pilot"
FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
FORMAL_HASH = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
OLD = PRODUCT / "Poliform-Senzafine-WD02-derived-drawing.ifc"
SINGLE = OUT / "WD02-product.ifc"
GUID = "3yuXF4PtnBHgIXDlmXFJ$7"
APPROVAL = ROOT / "pipeline/decisions/wd02-drawing-approval.json"
OLD_EVIDENCE = PRODUCT / "bonsai-drawings/wardrobe/WD02-WARDROBE-create-drawing-evidence.json"
REPORT = OUT / "validation.json"


def write(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def view_records():
    return json.loads(OLD_EVIDENCE.read_text())["outputs"]["views"]


def annotation_guids():
    return [v[k] for v in view_records() for k in ("drawing_global_id", "annotation_global_id")]


def protected_files():
    return [FORMAL, OLD, APPROVAL, OLD_EVIDENCE, PRODUCT / "candidate-representations.json",
            ROOT / "Materials.blend",
            *[PRODUCT / f"{v}.svg" for v in ("plan", "front", "side")],
            *[p for p in (PRODUCT / "official-source").rglob("*") if p.is_file()]]


def prepare():
    assert pkg.sha256(FORMAL) == FORMAL_HASH
    approval = json.loads(APPROVAL.read_text())
    assert approval["scene_svg_approved"] and approval["derived_ifc_write_allowed"]
    assert not approval["formal_authoritative_ifc_write_allowed"]
    for path, sha in approval["approved_artifact_sha256"].items():
        assert pkg.sha256(path) == sha, f"Approved input drift: {path}"
    assert not SINGLE.exists(), "Do not overwrite a previous pilot silently"
    OUT.mkdir(parents=True, exist_ok=True)
    protected = [pkg.record(p) for p in protected_files()]
    source = ifcopenshell.open(str(OLD))
    package = pkg.extract(source, GUID, annotation_guids())
    repairs = pkg.repair_missing_metadata(package)
    for v in view_records():
        for rel in package.by_guid(v["drawing_global_id"]).HasAssociations:
            if rel.is_a("IfcRelAssociatesDocument"):
                rel.RelatingDocument.Location = f"WD02-SCENE-{v['view'].upper()}.svg"
    package.write(str(SINGLE))
    reloaded = ifcopenshell.open(str(SINGLE))
    validation = pkg.validate(source, reloaded, GUID, annotation_guids())
    # Mapping recipes and paths may reference project elements absent in a
    # product-only file. They are resolved ONLY during temporary scene use.
    report = {"task": "WD02 persistent product IFC / temporary scene pilot", "stage": "extracted",
              "courseEvidence": {"mode": "embedded-course-index", "lesson": "085000", "timestamps": ["01:59 Create Drawing", "02:13 SVG"]},
              "preState": {"protected_files": protected, "formal_sha256": FORMAL_HASH},
              "extraction": validation, "single_product": pkg.record(SINGLE),
              "inherited_schema_repairs": repairs,
              "legacy_full_project": pkg.record(OLD), "formal_write_allowed": False,
              "batch_conversion_started": False, "cleanup_performed": False,
              "verdict": "pending_bonsai_roundtrip_and_scene"}
    write(REPORT, report)
    print(json.dumps({"single_product": report["single_product"], "validation": validation}, ensure_ascii=False))


def scene():
    import bpy
    import bonsai_bridge as bridge
    from bonsai import tool
    import ifcopenshell.api.document
    import create_wd03_wardrobe_scene_drawings as drawing_helpers
    assert bridge.bl_info["version"] == (1, 1, 0)
    report = json.loads(REPORT.read_text())
    if report.get("error"):
        report.setdefault("resolved_attempt_errors", []).append(report.pop("error"))
    report["verdict"] = "running"
    if report.get("temporary_project_ifc"):
        report.setdefault("previous_temporary_attempts", []).append({"path": report["temporary_project_ifc"], "retained_for_cleanup_confirmation": True})
    temporary = Path(tempfile.mkdtemp(prefix="wd02-review-scene-"))
    temp_ifc = temporary / "scene.ifc"
    report["temporary_project_directory"] = str(temporary)
    report["temporary_project_ifc"] = str(temp_ifc)
    write(REPORT, report)
    try:
        assert Path(tool.Ifc.get_path()).resolve() == SINGLE.resolve()
        assert pkg.sha256(FORMAL) == FORMAL_HASH
        assert not bpy.data.is_saved
        report["provider"] = {"adapter": "installed public bonsai-mcp handlers", "version": list(bridge.bl_info["version"]),
                              "blender": bpy.app.version_string, "ifcopenshell": ifcopenshell.version,
                              "bridge_file": bridge.__file__, "bridge_sha256": pkg.sha256(bridge.__file__)}
        # Public save handler proves Bonsai can persist and reload the package.
        with bpy.context.temp_override(**drawing_helpers.view3d_override()):
            report["product_bonsai_save_reload"] = bridge._h_save_ifc_file({"output_path": str(SINGLE), "overwrite": True, "reload": True})
        source = ifcopenshell.open(str(OLD))
        package = ifcopenshell.open(str(SINGLE))
        report["product_bonsai_roundtrip"] = pkg.validate(source, package, GUID, annotation_guids())
        assert tool.Ifc.get_object(tool.Ifc.get().by_guid(GUID)) is not None
        report["single_product"] = pkg.record(SINGLE)
        # Only the temporary path is ever loaded for project-scene work.
        shutil.copy2(FORMAL, temp_ifc)
        model = ifcopenshell.open(str(temp_ifc))
        resources = []
        for location in sorted({s.Location for s in model.by_type("IfcExternallyDefinedSurfaceStyle") if s.Location}):
            relative = Path(location)
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"Unreviewed external resource path: {location}")
            source_resource = ROOT / relative
            if source_resource.is_file():
                dest_resource = temporary / relative
                dest_resource.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source_resource, dest_resource)
                resources.append({"location": location, "source": pkg.record(source_resource), "temporary_copy": str(dest_resource)})
            else:
                resources.append({"location": location, "status": "already_missing_at_formal_source", "blocks_2d_linework": False})
        report["temporary_external_style_resources"] = resources
        formal_body = {e.GlobalId: pkg.body_fingerprint(e) for e in model.by_type("IfcElement") if e.Representation}
        report["attachment"] = pkg.attach(model, package, GUID, annotation_guids())
        first_count = len(list(model))
        pkg.attach(model, package, GUID, annotation_guids())
        assert len(list(model)) == first_count, "Repeated attachment must be idempotent"
        report["attachment"]["second_attachment_created_entities"] = 0
        assert formal_body == {e.GlobalId: pkg.body_fingerprint(e) for e in model.by_type("IfcElement") if e.Representation}
        report["attachment"]["all_formal_bodies_unchanged"] = True
        # Redirect all imported Drawing document locations away from legacy files.
        for v in view_records():
            drawing = model.by_guid(v["drawing_global_id"])
            refs = [r.RelatingDocument for r in drawing.HasAssociations if r.is_a("IfcRelAssociatesDocument")]
            assert len(refs) == 1
            refs[0].Location = str(OUT / f"WD02-SCENE-{v['view'].upper()}.svg")
        model.write(str(temp_ifc))
        assert bpy.ops.bim.load_project(filepath=str(temp_ifc), should_start_fresh_session=False, use_relative_path=False) == {"FINISHED"}
        tool.Ifc.get().by_guid(GUID)
        drawing_helpers.TARGET_GLOBAL_ID = GUID
        drawing_helpers.EXPECTED_PATH_COUNTS = {"plan": 15, "front": 41, "side": 7}
        outputs = []
        for v in view_records():
            drawing = tool.Ifc.get().by_guid(v["drawing_global_id"])
            tool.Ifc.get_object(drawing) or tool.Drawing.import_drawing(drawing)
            with bpy.context.temp_override(**drawing_helpers.view3d_override()):
                assert bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False) == {"FINISHED"}
                props = tool.Drawing.get_document_props()
                props.should_use_underlay_cache = False
                props.should_use_linework_cache = False
                props.should_use_annotation_cache = False
                result = bpy.ops.bim.create_drawing(print_all=False, open_viewer=False, sync=False)
            assert result == {"FINISHED"}
            svg = OUT / f"WD02-SCENE-{v['view'].upper()}.svg"
            inspection = drawing_helpers.style_and_inspect(svg, v["annotation_global_id"])
            outputs.append({"view": v["view"], "drawing_global_id": drawing.GlobalId,
                            "annotation_global_id": v["annotation_global_id"], "svg": pkg.record(svg),
                            "inspection": inspection, "generator": "bpy.ops.bim.create_drawing"})
        with bpy.context.temp_override(**drawing_helpers.view3d_override()):
            report["temporary_bonsai_save_reload"] = bridge._h_save_ifc_file({"output_path": str(temp_ifc), "overwrite": True, "reload": True})
        reloaded = tool.Ifc.get()
        assert np.array_equal(pkg.placement(reloaded.by_guid(GUID)), pkg.placement(package.by_guid(GUID)))
        assert formal_body == {e.GlobalId: pkg.body_fingerprint(e) for e in reloaded.by_type("IfcElement") if e.Representation}
        assert len([e for e in reloaded.by_type("IfcElement") if e.GlobalId == GUID]) == 1
        for v in view_records():
            ann = reloaded.by_guid(v["annotation_global_id"])
            assert pkg.fingerprint(ann.Representation) == pkg.fingerprint(package.by_guid(ann.GlobalId).Representation)
            assert np.array_equal(pkg.placement(ann), pkg.placement(package.by_guid(ann.GlobalId)))
        report["scene_outputs"] = outputs
        report["temporary_project"] = pkg.record(temp_ifc)
        report["scene_reload_verified"] = True
        report["verdict"] = "pending_visual_inspection"
        report["stage"] = "bonsai_roundtrip_and_scene_complete"
        # Leave the accepted product-only IFC open, never the formal project.
        assert bpy.ops.bim.load_project(filepath=str(SINGLE), should_start_fresh_session=False, use_relative_path=False) == {"FINISHED"}
    except Exception:
        report["verdict"] = "fail"
        report["error"] = traceback.format_exc()
        raise
    finally:
        report["formal_sha256_after"] = pkg.sha256(FORMAL)
        report["protected_files_unchanged"] = all(pkg.sha256(r["path"]) == r["sha256"] for r in report["preState"]["protected_files"])
        report["single_product"] = pkg.record(SINGLE)
        write(REPORT, report)
        assert report["formal_sha256_after"] == FORMAL_HASH
        assert report["protected_files_unchanged"]


if __name__ == "__main__":
    prepare()
