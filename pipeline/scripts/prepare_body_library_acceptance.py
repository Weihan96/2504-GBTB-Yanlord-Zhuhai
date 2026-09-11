"""Create a disposable 39-product acceptance IFC from the passing test fixture."""
from pathlib import Path
import json
import tempfile
import ifcopenshell
import ifcopenshell.api.root

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
report = json.loads((OUT / "native-integration-validation.json").read_text())
assert report["status"] == "pass" and len(report["products"]) == 40
model = ifcopenshell.open(report["fixture_path"])
# The fortieth item tests intentional repeated insertion, not another catalog
# product. Omit that test instance only from the separate acceptance fixture.
repeat = model.by_guid(report["products"][-1]["instance_guid"])
# This is a disposable multi-product display, not a published package. Detach
# the repeated instance's geometry before removal rather than recursively
# pruning shared display data across all 39 assets. Orphan data can remain in
# this temporary fixture; the validated single-product sources stay untouched.
repeat.Representation = None
ifcopenshell.api.root.remove_product(model, product=repeat)
assert len([e for e in model.by_type("IfcElement") if not e.is_a("IfcOpeningElement")]) == 39
target = Path(tempfile.mkdtemp(prefix="body-library-acceptance-")) / "library-review.ifc"
model.write(str(target))
(OUT / "body-library-acceptance-session.json").write_text(json.dumps({"ifc_path":str(target),"purpose":"39 product Body review; disposable, not formal IFC"},indent=2)+"\n")
print(target)
