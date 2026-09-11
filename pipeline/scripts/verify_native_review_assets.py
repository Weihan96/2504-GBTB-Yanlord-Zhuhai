"""Blender integration worker. Saves only a disposable IFC fixture."""
from pathlib import Path
import importlib.util
import hashlib
import json
import sys
import tempfile
import numpy as np
import bpy
import ifcopenshell
import ifcopenshell.api.root
import ifcopenshell.api.unit
import ifcopenshell.api.aggregate
import ifcopenshell.util.placement
import ifcopenshell.util.unit
import ifcopenshell.util.element
import bonsai.tool as tool

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
formal = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
formal_before = digest(formal)
spec = importlib.util.spec_from_file_location("highpoly_review_library", ROOT / "pipeline/addons/highpoly_review_library/__init__.py")
addon = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = addon
spec.loader.exec_module(addon)
addon.register()
from highpoly_review_library.native_assets import drop_matrix

temp = Path(tempfile.mkdtemp(prefix="review-native-ifc-"))
fixture = temp / "asset-insertion.ifc"
f = ifcopenshell.file(schema="IFC4")
project = ifcopenshell.api.root.create_entity(f, ifc_class="IfcProject", name="Disposable asset integration test")
ifcopenshell.api.unit.assign_unit(f)
site = ifcopenshell.api.root.create_entity(f, ifc_class="IfcSite", name="Test site")
storey = ifcopenshell.api.root.create_entity(f, ifc_class="IfcBuildingStorey", name="Test storey")
ifcopenshell.api.aggregate.assign_object(f, products=[site], relating_object=project)
ifcopenshell.api.aggregate.assign_object(f, products=[storey], relating_object=site)
f.write(str(fixture))
fixture_before = digest(fixture)
bpy.ops.bim.load_project(filepath=str(fixture), should_start_fresh_session=False, use_relative_path=False)
tool.Spatial.set_default_container(tool.Ifc.get().by_guid(storey.GlobalId))
records = []
products = addon._catalog["products"]
args = sys.argv[sys.argv.index("--")+1:] if "--" in sys.argv else []
if "--only" in args:
    requested = args[args.index("--only")+1].split(",")
    products = [p for p in products if p["id"] in requested]
    assert len(products) == len(requested)
from highpoly_review_library import loading_feedback
class StageObserver:
    def __init__(self, product_id):
        self.product_id = product_id
        self.stages = []
    def advance(self, fraction, label):
        self.stages.append((fraction, label))
for index, product in enumerate(products + [products[0]]):
    matrix = drop_matrix(product, ((index % 6)*4., (index // 6)*4., 0.))
    prior = {p.GlobalId for p in tool.Ifc.get().by_type("IfcElement")}
    observer = StageObserver(product["id"])
    loading_feedback._active = observer
    try:
        result = bpy.ops.review_library.insert_ifc(product_id=product["id"],
            project_matrix=[v for row in matrix for v in row])
    finally:
        loading_feedback._active = None
    assert [v for v, _ in observer.stages] == [.05, .25, .5, .65, .85, .95]
    assert result == {"FINISHED"}
    current = tool.Ifc.get()
    created = [p for p in current.by_type("IfcElement") if p.GlobalId not in prior and not p.is_a("IfcOpeningElement")]
    assert len(created) == 1, (product["id"], len(created))
    element = created[0]
    if product.get("insertion_content") == "review_body_views":
        from highpoly_review_library.body_views import validate
        validate(element)
    assert element.GlobalId != product["global_id"]
    source = ifcopenshell.open(str(addon.resolve_path(product["ifc_path"])))
    original = source.by_guid(product["global_id"])
    assert len(element.HasOpenings) == len(original.HasOpenings)
    expected = [(r.RepresentationIdentifier, r.RepresentationType, len(r.Items)) for r in original.Representation.Representations]
    actual = [(r.RepresentationIdentifier, r.RepresentationType, len(r.Items)) for r in element.Representation.Representations]
    assert expected == actual, (product["id"], expected, actual)
    actual_matrix = ifcopenshell.util.placement.get_local_placement(element.ObjectPlacement)
    actual_matrix[:3, 3] *= ifcopenshell.util.unit.calculate_unit_scale(current)
    assert np.allclose(actual_matrix, matrix, atol=1e-6)
    assert tool.Ifc.get_entity(tool.Ifc.get_object(element)) == element
    linked_objects = [o for o in bpy.context.scene.objects if tool.Ifc.get_entity(o) == element]
    assert len(linked_objects) == 1, (product["id"], "duplicate Blender product objects", len(linked_objects))
    assert digest(fixture) == fixture_before, "Insertion saved IFC unexpectedly"
    assert digest(addon.resolve_path(product["ifc_path"])) == product["ifc_sha256"]
    print("SERIALIZE_CHECK", product["id"], flush=True)
    assert current.to_string().startswith("ISO-10303-21;")
    records.append({"product": product["id"], "instance_guid": element.GlobalId, "loading_stages": observer.stages,
        "representation_count": len(actual), "project_matrix_m": [list(r) for r in matrix],
        "source_unchanged": True, "unsaved_fixture_unchanged": True, "bonsai_link": True})
    print("INSERT_PASS", product["id"], flush=True)

# Undo/redo belongs to the interactive acceptance test, with event-loop
# boundaries between operations. Do not simulate a UI undo stack inside one
# long-running background Python operator containing the entire batch.

# The installed Bonsai keymap, not an add-on save handler, owns Ctrl+S.
save_keys = [k for km in bpy.context.window_manager.keyconfigs.addon.keymaps for k in km.keymap_items
    if k.idname == "bim.save_project" and k.type == "S" and k.ctrl and not k.properties.should_save_as]
assert save_keys
bpy.ops.bim.save_project(filepath=str(fixture))
assert digest(fixture) != fixture_before
reloaded = ifcopenshell.open(str(fixture))
assert len([p for p in reloaded.by_type("IfcElement") if not p.is_a("IfcOpeningElement")]) == len(products) + 1
for record in records:
    element = reloaded.by_guid(record["instance_guid"])
    if addon.entry(record["product"]).get("insertion_content") == "review_body_views":
        validate(element)
    assert len(element.Representation.Representations) == record["representation_count"]
    assert ifcopenshell.util.element.get_psets(element)["ReviewLibrarySource"]["ProductId"] == record["product"]
    assert ifcopenshell.util.element.get_container(element).GlobalId == storey.GlobalId
assert digest(formal) == formal_before
report = {"status": "pass", "fixture_path": str(fixture), "formal_sha256_unchanged": formal_before,
    "products": records, "native_ctrl_s_keymap": True, "native_bonsai_save_reload": True,
    "physical_drag_ui_test": "pending", "undo_redo_evidence": "native-interactive-validation.json"}
report_name = "native-loading-insertion-validation.json" if "--only" in args else "native-integration-validation.json"
(OUT / report_name).write_text(json.dumps(report, indent=2) + "\n")
print("NATIVE_IFC_INTEGRATION_PASS", flush=True)
