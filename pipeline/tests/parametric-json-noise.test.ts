import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/parametric_json_noise_candidate.py");

test("parametric JSON cleanup snaps only near-integer numbers", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("parametric_json_noise_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
value = {"width": 799.9999, "ratio": 0.3333, "nested": [35.000001, True]}
normalized, changes = module.snap_near_integer(value, 0.01)
print(json.dumps({"normalized": normalized, "changes": changes}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const parsed = JSON.parse(result.stdout.toString());
  expect(parsed.normalized).toEqual({
    width: 800,
    ratio: 0.3333,
    nested: [35, true],
  });
  expect(parsed.changes).toHaveLength(2);
});
