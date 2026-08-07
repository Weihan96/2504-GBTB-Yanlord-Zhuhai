import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/direction_noise_candidate.py");

function runPython(source: string): string {
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  return result.stdout.toString().trim();
}

test("direction cleanup distinguishes dimensionless direction noise from millimetres", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("direction_noise_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
near = module.nearest_cardinal_3d((2.0, 0.0000004, 0.0))
intentional = module.nearest_cardinal_3d((1.0, 0.01, 0.0))
outside = module.nearest_cardinal_3d((1.0, 0.0))
print(json.dumps({"near": near, "intentional": intentional, "outside": outside}))
`);
  const parsed = JSON.parse(output);
  expect(parsed.near[0]).toEqual([1, 0, 0]);
  expect(parsed.near[1]).toBeLessThan(1e-6);
  expect(parsed.intentional[1]).toBeGreaterThan(1e-6);
  expect(parsed.outside).toBeNull();
});

test("direction cleanup preserves explicit topology-sensitive exceptions", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
import ifcopenshell
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("direction_noise_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
model = ifcopenshell.file(schema="IFC4")
origin = model.create_entity("IfcCartesianPoint", Coordinates=(0.0, 0.0, 0.0))
direction = model.create_entity("IfcDirection", DirectionRatios=(1.0, 0.0000004, 0.0))
model.create_entity("IfcAxis2Placement3D", Location=origin, RefDirection=direction)
changes = module.apply_direction_cleanup(model, 1e-6, {direction.id()})
print(json.dumps({"changes": changes, "ratios": list(direction.DirectionRatios)}))
`);
  const parsed = JSON.parse(output);
  expect(parsed.changes).toEqual([]);
  expect(parsed.ratios).toEqual([1, 0.0000004, 0]);
});

test("direction cleanup writes only a uniquely owned Axis2Placement3D direction", () => {
  const output = runPython(`
import importlib.util, json, pathlib, sys
import ifcopenshell
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("direction_noise_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
model = ifcopenshell.file(schema="IFC4")
origin = model.create_entity("IfcCartesianPoint", Coordinates=(0.0, 0.0, 0.0))
direction = model.create_entity("IfcDirection", DirectionRatios=(1.0, 0.0000004, 0.0))
model.create_entity("IfcAxis2Placement3D", Location=origin, RefDirection=direction)
first = module.placement_usage(model, direction)
changes = module.apply_direction_cleanup(model, 1e-6)
after = list(direction.DirectionRatios)
second_origin = model.create_entity("IfcCartesianPoint", Coordinates=(1.0, 0.0, 0.0))
model.create_entity("IfcAxis2Placement3D", Location=second_origin, RefDirection=direction)
shared = module.placement_usage(model, direction)
print(json.dumps({"first": first, "changes": changes, "after": after, "shared": shared}))
`);
  const parsed = JSON.parse(output);
  expect(parsed.first).toMatchObject({ inverse_count: 1, eligible: true });
  expect(parsed.first.roles).toEqual([
    expect.objectContaining({ role: "RefDirection" }),
  ]);
  expect(parsed.changes).toHaveLength(1);
  expect(parsed.after).toEqual([1, 0, 0]);
  expect(parsed.shared).toMatchObject({ inverse_count: 2, eligible: false });
});

test("spatial hash accepts reordered vertices only within the millimetre gate", () => {
  const output = runPython(`
import importlib.util, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("direction_noise_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(module.vertex_hausdorff_within_tolerance_mm([(0, 0, 0), (1, 0, 0)], [(1.05, 0, 0), (0, 0, 0)], 0.1))
`);
  expect(Number(output)).toBeCloseTo(0.05, 9);
});

test("spatial hash rejects a vertex outside the millimetre gate", () => {
  const output = runPython(`
import importlib.util, math, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("direction_noise_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
value = module.vertex_hausdorff_within_tolerance_mm([(0, 0, 0)], [(0.11, 0, 0)], 0.1)
print("infinite" if math.isinf(value) else value)
`);
  expect(output).toBe("infinite");
});

test("entity comparison records a removed STEP entity", () => {
  const output = runPython(`
import importlib.util, pathlib, sys
import ifcopenshell
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("direction_noise_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
source = ifcopenshell.file(schema="IFC4")
point = source.create_entity("IfcCartesianPoint", Coordinates=(0.0, 0.0, 0.0))
candidate = ifcopenshell.file.from_string(source.to_string())
candidate.remove(candidate.by_id(point.id()))
print(module.entity_changes(source, candidate))
`);
  expect(output).toBe("[1]");
});
