import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/space_placement_recovery_candidate.py");

test("Space recovery accepts only integer translations with identity rotation", () => {
  const source = `
import importlib.util, json, pathlib, sys, numpy as np
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("space_placement_recovery_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
check = module.placement_evidence
valid = np.eye(4); valid[:3,3] = [-6600, -5300, 0]
decimal = valid.copy(); decimal[0,3] = -6599.8
rotated = valid.copy(); rotated[0,0] = 0.999
print(json.dumps({
  "valid": check(valid, 0.1),
  "decimal": check(decimal, 0.1),
  "rotated": check(rotated, 0.1),
}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const parsed = JSON.parse(result.stdout.toString());
  expect(parsed.valid[0]).toBe(true);
  expect(parsed.decimal[0]).toBe(false);
  expect(parsed.rotated[0]).toBe(false);
});
