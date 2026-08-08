import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a105_floor_build_up_audit.py");

test("A105 floor build-up audit computes XY evidence overlap", () => {
  const script = `
import importlib.util, json, pathlib, sys
spec=importlib.util.spec_from_file_location("audit", ${JSON.stringify(modulePath)})
module=importlib.util.module_from_spec(spec); sys.modules[spec.name]=module; spec.loader.exec_module(module)
a={"min_mm":[0,0,0],"max_mm":[100,100,0]}
b={"min_mm":[50,-20,-50],"max_mm":[150,50,50]}
print(json.dumps(module.bbox_overlap_area_mm2(a,b)))
`;
  const result = Bun.spawnSync(["python3", "-c", script], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toBe(2500);
});

test("A105 floor build-up audit is read-only and tolerance-driven", async () => {
  const source = await Bun.file(modulePath).text();
  expect(source).toContain('parser.add_argument("--tolerance-mm"');
  expect(source).toContain('"automatic_ifc_write_allowed": False');
  expect(source).not.toContain("model.write(");
});
