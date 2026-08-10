import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a105_floor_plan_candidate.py");

function runPython(body: string) {
  return Bun.spawnSync(["python3", "-c", `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("a105_floor_plan_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
${body}
`], { cwd: root });
}

test("A105 uses the current IFC hash by default and rejects a mismatched caller freeze", async () => {
  const ifcPath = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
  const currentSha = createHash("sha256")
    .update(new Uint8Array(await Bun.file(ifcPath).arrayBuffer()))
    .digest("hex");
  const result = runPython(`
from pathlib import Path
values = [
    module.validate_source_sha(Path(${JSON.stringify(ifcPath)})),
    module.validate_source_sha(Path(${JSON.stringify(ifcPath)}), ${JSON.stringify(currentSha)}),
]
try:
    module.validate_source_sha(Path(${JSON.stringify(ifcPath)}), "0" * 64)
except RuntimeError as error:
    mismatch = str(error)
else:
    raise AssertionError("mismatched caller freeze was accepted")
print(json.dumps({"values": values, "mismatch": mismatch}))
`);
  expect(result.exitCode).toBe(0);
  const value = JSON.parse(result.stdout.toString());
  expect(value.values).toEqual([currentSha, currentSha]);
  expect(value.mismatch).toContain("formal IFC SHA-256 mismatch");
});

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

test("A105 slope arrow tip follows the mechanically fitted downhill vector", () => {
  const script = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec=importlib.util.spec_from_file_location("a105", ${JSON.stringify(modulePath)})
module=importlib.util.module_from_spec(spec); sys.modules[spec.name]=module; spec.loader.exec_module(module)
record={"bbox":{"centre_mm":[0,0,0],"dimensions_mm":[800,800,20]},"top_plane":{"a_dz_dx":0.0,"b_dz_dy":0.01}}
print(json.dumps(module.slope_arrow_geometry(record)))
`;
  const result = Bun.spawnSync(["python3", "-c", script], { cwd: root });
  expect(result.exitCode).toBe(0);
  const arrow = JSON.parse(result.stdout.toString());
  expect(arrow.end_mm[1]).toBeLessThan(arrow.start_mm[1]);
  expect(arrow.end_mm[0]).toBeCloseTo(arrow.start_mm[0], 9);
});

test("A105 source keeps every result read-only", async () => {
  const source = await Bun.file(modulePath).text();
  expect(source).toContain('"automatic_ifc_write_allowed": False');
  expect(source).toContain('"formal_ifc_write_allowed": "no"');
  expect(source).toContain('"wet_tile_slope_arrow_count"');
  expect(source).toContain('CONFIRMED_REFERENCE_MATERIALS');
  expect(source).toContain('A105-WET-SLOPE-DIRECTION-001');
  expect(source).toContain('"finish_reference_plane_confirmed"');
  expect(source).toContain('"confirmed_reference_count"');
  expect(source).toContain('Pset_A105FinishIntent');
  expect(source).not.toContain("model.write(");
});
