"""Collect the explicit runtime closure; do not delete any historical input."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import ifcopenshell

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
DEST = OUT / "runtime"
FORMAL = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
EXPECTED = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert digest(FORMAL) == EXPECTED
index = subprocess.check_output(["git","ls-files","--stage","-z"])
assert hashlib.sha256(index).hexdigest() == "fa9b7329175207683069cb18cc2f5015e1e7043ac52e64f882e95720168c53ea"
catalog = json.loads((OUT / "all-review-catalog.json").read_text())
DEST.mkdir(exist_ok=True)
copies = {}
def collect(source, relative):
    source = Path(source).resolve()
    target = DEST / relative
    assert source.is_file() and target.resolve().is_relative_to(DEST.resolve())
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,target)
    assert digest(source) == digest(target)
    copies[str(relative)] = {"source":str(source.relative_to(ROOT)),"sha256":digest(target),"bytes":target.stat().st_size}
    return str(relative)

fields = ("id","name","global_id","source_label","display_front_normal_project","display_front_basis",
          "single_product_approval_status","scene_approval_status","approval_label","package_validation",
          "ifc_sha256","ifc_bytes","insertion_content","representation_contract","review_category",
          "category_label","formal_ifc_write_authorized","skipped")
products = []
for entry in catalog["products"]:
    product = {key:entry[key] for key in fields if key in entry}
    source = (OUT / entry["ifc_path"]).resolve()
    product["ifc_path"] = collect(source,Path("ifc")/source.name)
    model = ifcopenshell.open(str(source))
    # IFC-resident external assets must resolve inside this same runtime.
    for style in model.by_type("IfcExternallyDefinedSurfaceStyle"):
        if style.Location:
            assert "://" not in style.Location and not Path(style.Location).is_absolute(), style.Location
            collect(source.parent/style.Location,Path("ifc")/style.Location)
    for texture in model.by_type("IfcImageTexture"):
        assert not texture.URLReference, "Texture dependencies need explicit packaging before activation"
    product["previews"] = {view:collect((OUT/path).resolve(),Path("previews")/entry["id"]/(view+Path(path).suffix)) for view,path in entry["previews"].items()}
    record = entry["approval_record"]
    product["approval_record"] = {"path":collect((OUT/record["path"]).resolve(),Path("approvals")/(entry["id"]+".json")),"sha256":record["sha256"]}
    collect(OUT/"native-assets"/(entry["id"]+".blend"),Path("native-assets")/(entry["id"]+".blend"))
    products.append(product)
collect(OUT/"native-assets/placements.json",Path("native-assets/placements.json"))
for source in (ROOT/"pipeline/addons/highpoly_review_library").glob("*.py"):
    collect(source,Path("addon/highpoly_review_library")/source.name)
result = {"schema_version":4,"body_view_adapter_required":True,"portable":True,"products":products}
(DEST/"catalog.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
report = {"status":"collected_pending_isolation_test","products":len(products),"files":copies,
          "formal_sha256":EXPECTED,"protected_index_sha256":hashlib.sha256(index).hexdigest(),
          "save_boundary":"new self-contained runtime and disposable test IFCs only; historical inputs retained"}
(OUT/"portable-build.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
assert subprocess.check_output(["git","ls-files","--stage","-z"]) == index
assert digest(FORMAL) == EXPECTED
print(json.dumps({"runtime":str(DEST),"products":len(products),"copied_files":len(copies),"bytes":sum(v["bytes"] for v in copies.values())}))
