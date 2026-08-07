import { expect, test } from "bun:test";
import { dirname, resolve } from "node:path";

const repositoryRoot = resolve(import.meta.dir, "../..");
const modulePath = resolve(repositoryRoot, "pipeline/scripts/c003_wall_candidate_compare.py");

function runPython(source: string): string {
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: repositoryRoot });
  expect(result.stderr.toString()).toBe("");
  expect(result.exitCode).toBe(0);
  return result.stdout.toString().trim();
}

function importPreamble(): string {
  return `
import importlib.util, json, sys
sys.path.insert(0, ${JSON.stringify(dirname(modulePath))})
spec = importlib.util.spec_from_file_location("c003_wall_candidate_compare", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
`;
}

test("C003 candidate snaps a near-cardinal proper rotation", () => {
  const output = runPython(`${importPreamble()}
import numpy as np
rotation = np.array([
  [-0.000204368432, 0.0, 0.999999979117],
  [0.999999979117, 0.0, 0.000204368432],
  [0.0, 1.0, 0.0],
])
print(json.dumps(module.nearest_proper_cardinal_rotation(rotation).tolist()))
`);
  expect(JSON.parse(output)).toEqual([
    [0, 0, 1],
    [1, 0, 0],
    [0, 1, 0],
  ]);
});

test("C003 representation rebase preserves the world transform", () => {
  const output = runPython(`${importPreamble()}
import numpy as np
source_product = np.eye(4)
source_product[:3, 3] = [-1899.8214006424, 3200.15978813171, 0.0]
source_solid = np.eye(4)
source_solid[:3, 3] = [40.0, -20.0, 5.0]
target_product = module.placement_target(source_product, cardinal=False)
target_solid = np.linalg.inv(target_product) @ source_product @ source_solid
print(np.max(np.abs(source_product @ source_solid - target_product @ target_solid)))
`);
  expect(Number(output)).toBeLessThan(1e-9);
});
