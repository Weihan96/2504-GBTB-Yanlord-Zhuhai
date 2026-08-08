import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a105_semantics_candidate.py");

test("A105 semantic boundary contains three finishes and two depressed slabs", () => {
  const script = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec=importlib.util.spec_from_file_location("a105s", ${JSON.stringify(modulePath)})
module=importlib.util.module_from_spec(spec); sys.modules[spec.name]=module; spec.loader.exec_module(module)
print(json.dumps({"finishes":module.REFERENCE_MATERIALS,"slabs":module.DEPRESSED_SLABS,"guid":module.deterministic_guid("TEST")}))
`;
  const result = Bun.spawnSync(["python3", "-c", script], { cwd: root });
  expect(result.exitCode).toBe(0);
  const data = JSON.parse(result.stdout.toString());
  expect(Object.keys(data.finishes).length).toBe(3);
  expect(Object.keys(data.slabs).length).toBe(2);
  expect(data.guid.length).toBe(22);
});

test("A105 semantics candidate stays outside the formal IFC", async () => {
  const source = await Bun.file(modulePath).text();
  expect(source).toContain("candidate output must not overwrite the formal IFC");
  expect(source).toContain('"formal_ifc_write_allowed": False');
  expect(source).toContain("FINISH_REFERENCE_PLANE");
});
