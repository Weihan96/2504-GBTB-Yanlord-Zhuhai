import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a105_floor_plan_candidate.py");

test("A105 top-plane fit reports slope and downhill direction", () => {
  const script = `
import importlib.util, json, numpy as np, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec=importlib.util.spec_from_file_location("a105", ${JSON.stringify(modulePath)})
module=importlib.util.module_from_spec(spec); sys.modules[spec.name]=module; spec.loader.exec_module(module)
vertices=np.array([[0,0,0],[1000,0,0],[1000,1000,10],[0,1000,10],[0,0,-20],[1000,0,-20],[1000,1000,-10],[0,1000,-10]],dtype=float)
faces=np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6]],dtype=int)
print(json.dumps(module.top_plane(vertices,faces)))
`;
  const result = Bun.spawnSync(["python3", "-c", script], { cwd: root });
  expect(result.exitCode).toBe(0);
  const plane = JSON.parse(result.stdout.toString());
  expect(plane.slope_percent).toBeCloseTo(1, 6);
  expect(plane.downhill_direction).toBe("S");
  expect(plane.maximum_fit_residual_mm).toBeLessThan(1e-8);
});

test("A105 source keeps every result read-only", async () => {
  const source = await Bun.file(modulePath).text();
  expect(source).toContain('"automatic_ifc_write_allowed": False');
  expect(source).toContain('"formal_ifc_write_allowed": "no"');
  expect(source).not.toContain("model.write(");
});
