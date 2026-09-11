"""Cold-start OS-isolated Blender workers; no legacy data is readable."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/review/approved-product-library"
BLENDER = "/Applications/Blender.app/Contents/MacOS/Blender"
digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
parser = argparse.ArgumentParser()
parser.add_argument("--only")
args = parser.parse_args()
capsule = Path(tempfile.mkdtemp(prefix="portable-library-isolation-"))
shutil.copytree(OUT/"runtime",capsule/"runtime")
worker = capsule/"verify_product.py"
shutil.copy2(ROOT/"pipeline/scripts/verify_portable_review_product.py",worker)
catalog = json.loads((capsule/"runtime/catalog.json").read_text())
products = catalog["products"]
if args.only:
    products = [p for p in products if p["id"] in args.only.split(",")]
profile = '(version 1)(allow default)(deny file-read* file-write* (subpath '+json.dumps(str(ROOT))+'))(deny file-read* (subpath "/Users/jiaxinchen/Poliigon"))'
assert subprocess.run(["/usr/bin/sandbox-exec","-p",profile,"/bin/cat",str(ROOT/"2504 GBTB Yanlord Zhuhai.ifc")],capture_output=True).returncode != 0
session = {"capsule":str(capsule),"blocked_root":str(ROOT),"products":[p["id"] for p in products],
           "runtime_catalog_sha256":digest(capsule/"runtime/catalog.json"),"worker_sha256":digest(worker),"profile":profile}
(OUT/"portable-isolation-session.json").write_text(json.dumps(session,indent=2)+"\n")
def run(product):
    slug = product["id"]
    output = capsule/"results"/slug
    output.mkdir(parents=True)
    command = ["/usr/bin/sandbox-exec","-p",profile,BLENDER,"-b","--python",str(worker),"--",
               "--runtime",str(capsule/"runtime"),"--out",str(output),"--slug",slug,"--blocked-root",str(ROOT)]
    with (output/"worker.log").open("w") as log:
        process = subprocess.run(command,cwd=capsule,stdout=log,stderr=subprocess.STDOUT,timeout=300)
    result = {"id":slug,"exit_code":process.returncode,"log":str(output/"worker.log")}
    if process.returncode == 0 and (output/"result.json").is_file():
        evidence = json.loads((output/"result.json").read_text())
        assert evidence["source_sha256"] == product["ifc_sha256"]
        target = OUT/"portable-validation"/slug
        target.mkdir(parents=True,exist_ok=True)
        for drawing in evidence["drawings"]:
            source = Path(drawing["svg"])
            destination = target/source.name
            shutil.copy2(source,destination)
            assert digest(destination) == drawing["sha256"]
            drawing["svg"] = str(destination)
        shutil.copy2(output/"worker.log",target/"worker.log")
        evidence["worker_log"] = str(target/"worker.log")
        (target/"result.json").write_text(json.dumps(evidence,indent=2)+"\n")
        result.update(status="structured_pass",evidence=str(target/"result.json"))
    else:
        result["status"] = "fail"
    print(json.dumps(result),flush=True)
    return result
results = []
with ThreadPoolExecutor(max_workers=2) as pool:
    for future in as_completed([pool.submit(run,p) for p in products]):
        results.append(future.result())
summary = {**session,"status":"structured_pass" if all(r["status"]=="structured_pass" for r in results) else "fail",
           "results":results,"svg_count":3*sum(r["status"]=="structured_pass" for r in results),"visual_check":"pending"}
(OUT/"portable-isolation-validation.json").write_text(json.dumps(summary,indent=2)+"\n")
assert summary["status"] == "structured_pass"
