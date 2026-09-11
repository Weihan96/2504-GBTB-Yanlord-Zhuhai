"""Isolated Blender worker: native Create Drawing, provider save/reload.

Writes only a temporary IFC fixture and named review SVG/validation artifacts.
"""
from pathlib import Path
import importlib.util
import json
import sys
import tempfile
import shutil
import hashlib
import numpy as np
import bpy
import addon_utils
import ifcopenshell
import ifcopenshell.api.geometry
import ifcopenshell.api.pset
import ifcopenshell.util.element as eu
import ifcopenshell.util.unit as uu
import bonsai.tool as tool
from bonsai.core import drawing as core
from mathutils import Matrix, Vector
from lxml import etree

import ifcopenshell.api.root
import ifcopenshell.api.unit
import ifcopenshell.api.aggregate
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--runtime", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--slug", required=True)
parser.add_argument("--blocked-root", required=True)
args = parser.parse_args(sys.argv[sys.argv.index("--")+1:])
RUNTIME = Path(args.runtime).resolve()
OUT = Path(args.out).resolve()
OUT.mkdir(parents=True, exist_ok=True)
slug = args.slug
cut_mode = "BISECT"
digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
try:
    Path(args.blocked_root, "2504 GBTB Yanlord Zhuhai.ifc").read_bytes()
except PermissionError:
    pass
else:
    raise AssertionError("OS sandbox did not deny old repository access")
addon_utils.enable("bonsai_bridge", default_set=False, persistent=False)
import bonsai_bridge
assert bonsai_bridge.bl_info["version"] == (1, 1, 0)
# Native Provider is imported and version checked before IFC authoring.
spec = importlib.util.spec_from_file_location("highpoly_review_library", RUNTIME / "addon/highpoly_review_library/__init__.py")
addon = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = addon
spec.loader.exec_module(addon)
addon.register()
assert addon._catalog_path == RUNTIME / "catalog.json"
assert len(addon._catalog["products"]) == 39
from highpoly_review_library import body_views, drawing_integration, loading_feedback
from highpoly_review_library.native_assets import drop_matrix
entry = addon.entry(slug)
source = addon.resolve_path(entry["ifc_path"])
source_hash = digest(source)
assert source_hash == entry["ifc_sha256"]
fixture_dir = OUT
fixture = fixture_dir / "test.ifc"
f = ifcopenshell.file(schema="IFC4")
project = ifcopenshell.api.root.create_entity(f, ifc_class="IfcProject", name="Portable library isolation test")
ifcopenshell.api.unit.assign_unit(f)
site = ifcopenshell.api.root.create_entity(f, ifc_class="IfcSite", name="Test site")
storey = ifcopenshell.api.root.create_entity(f, ifc_class="IfcBuildingStorey", name="Test storey")
ifcopenshell.api.aggregate.assign_object(f, products=[site], relating_object=project)
ifcopenshell.api.aggregate.assign_object(f, products=[storey], relating_object=site)
f.write(str(fixture))
before_hash = digest(fixture)
assert bpy.ops.bim.load_project(filepath=str(fixture), should_start_fresh_session=False, use_relative_path=False) == {"FINISHED"}
tool.Spatial.set_default_container(tool.Ifc.get().by_guid(storey.GlobalId))
insertion_matrix = drop_matrix(entry, (2., 3., 0.))
class Observer:
    product_id = slug
    stages = []
    def advance(self, fraction, label):
        self.stages.append(fraction)
observer = Observer()
loading_feedback._active = observer
try:
    assert bpy.ops.review_library.insert_ifc(product_id=slug, project_matrix=[float(v) for row in insertion_matrix for v in row]) == {"FINISHED"}
finally:
    loading_feedback._active = None
assert observer.stages == [.05,.25,.5,.65,.85,.95]
assert digest(fixture) == before_hash
model = tool.Ifc.get()
product = next(e for e in model.by_type("IfcElement") if not e.is_a("IfcOpeningElement"))
guid = product.GlobalId
directional = body_views.validate(product)
obj = tool.Ifc.get_object(product)
assert obj
project_matrix = ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement)
points = [obj.matrix_world @ Vector(p) for p in obj.bound_box]
centre = sum(points, Vector()) / len(points)
diameter = max((p-centre).length for p in points) * 2 + .2
drawings = []
svg_dir = OUT / "svg"
svg_dir.mkdir(parents=True, exist_ok=True)
for role in body_views.ROLES:
    target_view = "PLAN_VIEW" if role == "Plan" else "ELEVATION_VIEW"
    before = {e.id() for e in model.by_type("IfcAnnotation")}
    core.add_drawing(tool.Ifc, tool.Collector, tool.Drawing, target_view=target_view, location_hint=0 if role == "Plan" else "SOUTH")
    drawing = next(e for e in model.by_type("IfcAnnotation") if e.id() not in before)
    core.update_drawing_name(tool.Ifc, tool.Drawing, drawing=drawing, name=f"BODY-{slug}-{role}")
    camera = tool.Ifc.get_object(drawing)
    normal = Vector(project_matrix[:3,:3] @ directional[role][1]).normalized()
    up = Vector((0,1,0)) if abs(normal.z) > .9 else Vector((0,0,1))
    x = up.cross(normal).normalized()
    y = normal.cross(x).normalized()
    rotation = Matrix((x, y, normal)).transposed()
    matrix = rotation.to_4x4()
    matrix.translation = centre + normal * diameter
    camera.matrix_world = matrix
    camera.data.type = "ORTHO"
    camera.data.clip_start = .002
    camera.data.clip_end = diameter * 3
    props = tool.Drawing.get_camera_props(camera)
    props.update_props = False
    props.camera_type = "ORTHO"
    props.target_view = target_view
    props.custom_scale_numerator = "1"
    props.custom_scale_denominator = "10"
    props.diagram_scale = "CUSTOM"
    props.has_underlay = False
    props.has_linework = True
    props.has_annotation = False
    props.linework_mode = "OPENCASCADE"
    props.fill_mode = "NONE"
    props.cut_mode = cut_mode
    projected = [rotation.inverted() @ (p-centre) for p in points]
    props.width = max(p.x for p in projected)-min(p.x for p in projected)+.2
    props.height = max(p.y for p in projected)-min(p.y for p in projected)+.2
    props.update_camera_resolution()
    props.update_props = True
    tool.Drawing.sync_object_representation(camera)
    ifcopenshell.api.geometry.edit_object_placement(model, product=drawing, matrix=np.asarray(matrix), is_si=True)
    pset = model.by_id(eu.get_pset(drawing, "EPset_Drawing")["id"])
    ifcopenshell.api.pset.edit_pset(model, pset=pset, properties={"TargetView":target_view,"Scale":"1/10","HumanScale":"1:10",
        "HasUnderlay":False,"HasLinework":True,"HasAnnotation":False,"Include":guid,"LineworkMode":"OPENCASCADE","FillMode":"NONE","CutMode":cut_mode})
    svg = svg_dir / (role.lower()+".svg")
    tool.Drawing.get_drawing_document(drawing).Location = str(svg)
    drawings.append({"role":role,"guid":drawing.GlobalId,"svg":str(svg)})

persistence = bonsai_bridge._h_save_ifc_file({"output_path":str(fixture),"overwrite":True,"reload":True})
model = tool.Ifc.get()
body_views.validate(model.by_guid(guid))
area = next(a for a in bpy.context.screen.areas if a.type == "VIEW_3D")
region = next(r for r in area.regions if r.type == "WINDOW")
results = []
for item in drawings:
    drawing = model.by_guid(item["guid"])
    with bpy.context.temp_override(area=area, region=region):
        assert bpy.ops.bim.activate_drawing(drawing=drawing.id(), should_view_from_camera=False) == {"FINISHED"}
        props = tool.Drawing.get_document_props()
        props.should_use_linework_cache = False
        props.should_use_underlay_cache = False
        props.should_use_annotation_cache = False
        drawing_integration.selection_log.clear()
        drawing_integration.bisect_log.clear()
        assert bpy.ops.bim.create_drawing(print_all=False, open_viewer=False, sync=False) == {"FINISHED"}
    svg = Path(item["svg"])
    root = etree.parse(str(svg)).getroot()
    groups = root.xpath('//*[@*[local-name()="guid"]=$guid]',guid=guid)
    paths = [p for g in groups for p in g.iter() if etree.QName(p).localname in ("path","polyline","line","polygon","circle","ellipse")]
    assert paths, (item,"empty target geometry")
    selected = [r for r in drawing_integration.selection_log if r["global_id"] == guid]
    assert len(selected) == 1 and selected[0]["role"] == item["role"], selected
    if cut_mode == "BISECT":
        assert any(r["excluded_mesh_objects"] for r in drawing_integration.bisect_log), drawing_integration.bisect_log
    results.append({**item,"sha256":digest(svg),"bytes":svg.stat().st_size,"target_geometry_nodes":len(paths),"selected":selected})
    print("BODY_DRAWING_PASS",slug,item["role"],len(paths),flush=True)
assert digest(source) == source_hash
product = tool.Ifc.get().by_guid(guid)
body_views.validate(product)
actual = ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement)
actual[:3,3] *= uu.calculate_unit_scale(tool.Ifc.get())
assert np.allclose(actual, insertion_matrix, atol=1e-6)
assert len([e for e in tool.Ifc.get().by_type("IfcElement") if not e.is_a("IfcOpeningElement")]) == 1
assert len([o for o in bpy.context.scene.objects if tool.Ifc.get_entity(o) == product]) == 1
assert eu.get_pset(product, "ReviewLibrarySource", "ProductId") == slug
assert guid != entry["global_id"]
# Texture files cannot silently fall back to a deleted external directory.
images = []
for material in bpy.data.materials:
    if not material.use_nodes:
        continue
    for node in material.node_tree.nodes:
        image = getattr(node, "image", None)
        if image and image.source in ("FILE","TILED","SEQUENCE","MOVIE"):
            path = Path(bpy.path.abspath(image.filepath)).resolve()
            assert image.packed_file or list(image.packed_files) or (path.is_relative_to(RUNTIME) and path.is_file()), (material.name,image.name,str(path))
            images.append({"name":image.name,"packed":bool(image.packed_file or list(image.packed_files))})
report = {"verdict":"structured_pass_pending_raster_check","product":slug,"source_sha256":source_hash,"drawings":results,
          "temporary_fixture":str(fixture),"provider_persistence":persistence,"generator":"Bonsai Create Drawing / OpenCASCADE + BISECT",
          "old_repository_read_denied_by_OS":True,"catalog_products_loaded":39,"native_insert":True,
          "disk_unchanged_until_save":True,"unique_product_and_object":True,"reloaded":True,
          "loading_stages":observer.stages,"texture_dependencies":images,"global_id_changed_on_insert":True,
          "placement_m":[list(r) for r in actual]}
(OUT / "result.json").write_text(json.dumps(report,indent=2,default=str)+"\n")
print("PORTABLE_PRODUCT_PASS",slug,flush=True)
