import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a104_semantics_candidate.py");

test("A104 semantic groups have deterministic identities and exact members", () => {
  const script = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec=importlib.util.spec_from_file_location("a104s", ${JSON.stringify(modulePath)})
module=importlib.util.module_from_spec(spec); sys.modules[spec.name]=module; spec.loader.exec_module(module)
print(json.dumps({"ids":[module.deterministic_guid(k+":GROUP") for k in module.GROUPS],"members":[len(v["members"]) for v in module.GROUPS.values()]}))
`;
  const result = Bun.spawnSync(["python3", "-c", script], { cwd: root });
  expect(result.exitCode).toBe(0);
  const data = JSON.parse(result.stdout.toString());
  expect(new Set(data.ids).size).toBe(3);
  expect(data.members).toEqual([2, 2, 2]);
});

test("A104 semantics candidate never overwrites the formal IFC", async () => {
  const source = await Bun.file(modulePath).text();
  expect(source).toContain("candidate output must not overwrite the formal IFC");
  expect(source).toContain('"formal_ifc_write_allowed": False');
  expect(source).toContain('"IfcGroup"');
});
