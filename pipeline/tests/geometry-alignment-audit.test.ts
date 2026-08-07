import { expect, test } from "bun:test";
import { resolve } from "node:path";

const repositoryRoot = resolve(import.meta.dir, "../..");
const modulePath = resolve(repositoryRoot, "pipeline/scripts/geometry_alignment_audit.py");

function runPython(source: string): string {
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: repositoryRoot });
  expect(result.exitCode).toBe(0);
  return result.stdout.toString().trim();
}

test("geometry audit keeps placement translation and rotation as separate units", () => {
  const output = runPython(`
import importlib.util, json, sys
spec = importlib.util.spec_from_file_location("geometry_alignment_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
first = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
second = [[1.0, 0.0, 0.0, 0.24], [0.0, 0.999, -0.001, 0.0], [0.0, 0.001, 0.999, 0.0], [0.0, 0.0, 0.0, 1.0]]
print(json.dumps(module.placement_deltas(first, second)))
`);
  const [translationMm, rotationDelta] = JSON.parse(output) as [number, number];
  expect(translationMm).toBeCloseTo(0.24, 9);
  expect(rotationDelta).toBeCloseTo(0.001, 9);
});

test("geometry audit measures a parallel wall-face gap in millimetres", () => {
  const output = runPython(`
import importlib.util, sys
spec = importlib.util.spec_from_file_location("geometry_alignment_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
first = module.Segment("A", None, (0.0, 0.0), (1000.0, 0.0), 1000.0)
second = module.Segment("B", None, (0.0, 0.056), (1000.0, 0.056), 1000.0)
print(module.parallel_line_gap_mm(first, second, module.segment_unit(first)))
`);
  expect(Number(output)).toBeCloseTo(0.056, 9);
});

test("geometry audit does not classify parallel edge extensions as junctions", () => {
  const output = runPython(`
import importlib.util, math, sys
spec = importlib.util.spec_from_file_location("geometry_alignment_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
tolerance = math.sin(math.radians(0.1))
print(module.segments_are_parallel((1.0, 0.0), (-1.0, 0.0), tolerance))
print(module.segments_are_parallel((1.0, 0.0), (0.0, 1.0), tolerance))
`);
  expect(output.split("\n")).toEqual(["True", "False"]);
});

test("geometry audit treats the minimum non-parallel wall-pair distance as the junction", () => {
  const output = runPython(`
import importlib.util, sys
spec = importlib.util.spec_from_file_location("geometry_alignment_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
records = {}
for gap in (1.0, 0.0, 0.05):
    key = ("A", "B")
    record = {"gap_mm": gap}
    module.retain_minimum_gap_record(records, key, record)
print(records[("A", "B")]["gap_mm"])
`);
  expect(Number(output)).toBe(0);
});

test("geometry audit excludes plan neighbours without vertical overlap", () => {
  const output = runPython(`
import importlib.util, sys
spec = importlib.util.spec_from_file_location("geometry_alignment_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
lower = module.Segment("LOW", None, (0.0, 0.0), (1000.0, 0.0), 1000.0, 0.0, 2800.0)
upper = module.Segment("UP", None, (0.0, 0.05), (1000.0, 0.05), 1000.0, 3000.0, 5800.0)
touching = module.Segment("TOUCH", None, (0.0, 0.05), (1000.0, 0.05), 1000.0, 2700.0, 5500.0)
print(module.vertical_overlap_mm(lower, upper), module.vertical_overlap_mm(lower, touching))
`);
  const [disjoint, overlapping] = output.split(" ").map(Number);
  expect(disjoint).toBe(0);
  expect(overlapping).toBe(100);
});

test("geometry audit clusters only unresolved wall-pair candidates", () => {
  const output = runPython(`
import importlib.util, json, sys
spec = importlib.util.spec_from_file_location("geometry_alignment_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
coplanar = [
    {"wall_a": "A", "wall_b": "B", "gap_mm": 0.3, "within_tolerance": False},
    {"wall_a": "B", "wall_b": "C", "gap_mm": 0.05, "within_tolerance": True},
    {"wall_a": "D", "wall_b": "E", "gap_mm": 0.4, "within_tolerance": False},
]
junctions = [
    {"wall_a": "B", "wall_b": "C", "gap_mm": 0.2, "within_tolerance": False},
]
print(json.dumps(module.build_review_clusters(coplanar, junctions)))
`);
  const clusters = JSON.parse(output) as Array<Record<string, unknown>>;
  expect(clusters).toHaveLength(2);
  expect(clusters[0]).toMatchObject({
    review_id: "C003-G1",
    display_order: 1,
    wall_count: 3,
    wall_global_ids: ["A", "B", "C"],
    pair_count: 2,
    coplanar_record_count: 1,
    junction_record_count: 1,
    minimum_gap_mm: 0.2,
    maximum_gap_mm: 0.3,
  });
  expect(String(clusters[0].stable_id)).toMatch(/^C003-S[0-9A-F]{8}$/);
  expect(clusters[1]).toMatchObject({
    review_id: "C003-G2",
    wall_count: 2,
    wall_global_ids: ["D", "E"],
    pair_count: 1,
  });
});
