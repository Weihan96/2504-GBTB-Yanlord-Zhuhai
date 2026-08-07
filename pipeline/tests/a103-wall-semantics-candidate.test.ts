import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a103_wall_semantics_candidate.py");

function runPython(body: string) {
  return Bun.spawnSync(["python3", "-c", `
import collections, importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("a103_wall_semantics_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
${body}
`], { cwd: root });
}

test("A103 confirmed wall boundary is 84 existing and 4 new", () => {
  const result = runPython(`
model = module.ifcopenshell.open("2504 GBTB Yanlord Zhuhai.ifc")
walls = module.final_built_walls(model)
records = [module.expected_wall_semantics(wall) for wall in walls]
print(json.dumps({
  "wall_count": len(walls),
  "demolition_count": len(model.by_type("IfcWall")) - len(walls),
  "phase": collections.Counter(record["status"] for record in records),
  "load": collections.Counter(record["load_bearing"] for record in records if record["status"] == "EXISTING"),
}))
`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    wall_count: 88,
    demolition_count: 13,
    phase: { EXISTING: 84, NEW: 4 },
    load: { "false": 64, "true": 20 },
  });
});

test("A103 kitchen pier normalization is idempotent at 154 mm", () => {
  const result = runPython(`
model = module.ifcopenshell.open("2504 GBTB Yanlord Zhuhai.ifc")
print(json.dumps(module.normalize_kitchen_pier_thickness(model)))
`);
  expect(result.exitCode).toBe(0);
  const output = JSON.parse(result.stdout.toString());
  expect(output.target_thickness_mm).toBe(154);
  expect(output.changed_coordinate_count).toBe(4);
  expect(output.maximum_intended_world_shift_mm).toBe(0);
});
