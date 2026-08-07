import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a103_wall_plan_candidate.py");

function runPython(body: string) {
  return Bun.spawnSync(["python3", "-c", `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("a103_wall_plan_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
${body}
`], { cwd: root });
}

test("A103 world-to-SVG mapping matches the 1:50 Wall Plan camera", () => {
  const result = runPython(`print(json.dumps([
    module.world_to_svg(-6600, 0),
    module.world_to_svg(3300, 0),
    module.world_to_svg(0, 4800),
  ]))`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual([
    [68, 200],
    [266, 200],
    [200, 104],
  ]);
});

test("A103 wall phase grouping requires confirmed IFC status and load-bearing semantics", () => {
  const result = runPython(`print(json.dumps([
    module.wall_status_candidate("NEW", None, ["Unknown"]),
    module.wall_status_candidate("EXISTING", False, ["Aircrete"]),
    module.wall_status_candidate("EXISTING", True, ["Concrete"]),
    module.wall_status_candidate(None, None, ["WhiteWall"]),
  ]))`);
  expect(result.exitCode).toBe(0);
  const values = JSON.parse(result.stdout.toString());
  expect(values.map((value: unknown[]) => value[0])).toEqual([
    "CONFIRMED_NEW",
    "CONFIRMED_EXISTING_NON_LOAD_BEARING",
    "CONFIRMED_EXISTING_LOAD_BEARING",
    "INVALID_SEMANTICS",
  ]);
  expect(values.slice(0, 3).every((value: unknown[]) => value[3] === "no")).toBe(true);
  expect(values[3][3]).toBe("yes");
});

test("A103 dimension chains close mechanically", () => {
  const result = runPython(`print(json.dumps(module.chain_closure([
    {"coordinate_mm": -6600, "segment_from_previous_mm": None},
    {"coordinate_mm": -4800, "segment_from_previous_mm": 1800},
    {"coordinate_mm": 3300, "segment_from_previous_mm": 8100},
    {"coordinate_mm": 6600, "segment_from_previous_mm": 3300},
  ])))`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    segment_sum_mm: 13200,
    overall_mm: 13200,
    residual_mm: 0,
    pass: true,
  });
});

test("A103 label placement refuses generated collisions", () => {
  const result = runPython(`
occupied = []
first = module.place_label(10, 10, "N01", occupied, offsets=[(0, 0)])
second = module.place_label(10, 10, "N02", occupied, offsets=[(8, 0)])
print(json.dumps({"count": len(occupied), "first": first[:2], "second": second[:2]}))
`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    count: 2,
    first: [10, 10],
    second: [18, 10],
  });
});

test("A103 colors every repeated SVG fragment exactly once", () => {
  const result = runPython(`
svg = '<g class="IfcWall WALL-A cut"></g><g class="WALL-A IfcWall cut"></g>'
print(module.add_wall_status_classes(svg, {"WALL-A": "CONFIRMED_NEW"}))
`);
  expect(result.exitCode).toBe(0);
  const output = result.stdout.toString();
  expect(output.match(/a103-confirmed-new/g)).toHaveLength(2);
  expect(output).not.toContain("a103-confirmed-new a103-confirmed-new");
});
