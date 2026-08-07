import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

function runPython(expression: string) {
  return spawnSync("python3", ["-c", expression], {
    cwd: process.cwd(),
    encoding: "utf8",
  });
}

test("integer geometry anchor prefers the strongest project module", () => {
  const result = runPython(
    `import sys; sys.path.insert(0,"pipeline/scripts"); from integer_geometry_anchor_audit import preferred_module; print(preferred_module((100,-200,3000),0.01)); print(preferred_module((50,150,3000),0.01)); print(preferred_module((31,101,3000),0.01))`,
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim().split("\n")).toEqual(["100.0", "50.0", "1.0"]);
});

test("near-integer geometry vertex reports an exact lattice target", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from integer_geometry_anchor_audit import find_anchor; print(json.dumps(find_anchor([(100.04,199.96,300.03)],[],0.1)))`,
  );
  expect(result.status).toBe(0);
  const anchor = JSON.parse(result.stdout);
  expect(anchor.anchor_kind).toBe(
    "integer_point_within_tolerance_of_existing_vertex",
  );
  expect(anchor.point_mm).toEqual([100, 200, 300]);
  expect(anchor.source_geometry_point_mm[0]).toBeCloseTo(100.04);
  expect(anchor.source_geometry_point_mm[1]).toBeCloseTo(199.96);
  expect(anchor.source_geometry_point_mm[2]).toBeCloseTo(300.03);
  expect(anchor.snap_distance_mm).toBeCloseTo(0.04);
  expect(anchor.module_mm).toBe(100);
});

test("integer geometry anchor ignores a coplanar triangulation diagonal", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from integer_geometry_anchor_audit import physical_edges; v=[(0,0,0),(100,0,0),(100,100,0),(0,100,0)]; f=[(0,1,2),(0,2,3)]; print(json.dumps(physical_edges(v,f)))`,
  );
  expect(result.status).toBe(0);
  expect(JSON.parse(result.stdout)).toEqual([
    [0, 1],
    [0, 3],
    [1, 2],
    [2, 3],
  ]);
});

test("integer geometry anchor finds a module point on a physical axis edge", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from integer_geometry_anchor_audit import point_on_axis_aligned_edge; print(json.dumps(point_on_axis_aligned_edge((0,0,-55),(260,0,-55),0.01)))`,
  );
  expect(result.status).toBe(0);
  expect(JSON.parse(result.stdout)).toEqual([[130, 0, -55], 1]);
});

test("integer geometry anchor finds a module point on a diagonal physical edge", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from integer_geometry_anchor_audit import point_on_non_axis_physical_edge; print(json.dumps({"found": point_on_non_axis_physical_edge((0.2,0.2,0.2),(80.2,80.2,80.2),0.01), "missing": point_on_non_axis_physical_edge((0.2,0.4,0.6),(80.2,80.4,80.6),0.01)}))`,
  );
  expect(result.status).toBe(0);
  expect(JSON.parse(result.stdout)).toEqual({
    found: [[50, 50, 50], 50],
    missing: null,
  });
});

test("construction surface fallback finds an integer point inside a horizontal face", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from integer_geometry_anchor_audit import find_anchor; v=[(0.2,0.2,20),(200.2,0.2,20),(200.2,200.2,20),(0.2,200.2,20)]; f=[(0,1,2),(0,2,3)]; print(json.dumps(find_anchor(v,f,0.01,allow_axis_aligned_surface=True)))`,
  );
  expect(result.status).toBe(0);
  const anchor = JSON.parse(result.stdout);
  expect(anchor.anchor_kind).toBe("integer_point_on_axis_aligned_surface");
  expect(anchor.module_mm).toBe(20);
  expect(anchor.point_mm[2]).toBe(20);
  expect(anchor.point_mm.every((value: number) => value % 20 === 0)).toBe(true);
});

test("surface fallback stays disabled for semantically ambiguous meshes", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0,"pipeline/scripts"); from integer_geometry_anchor_audit import find_anchor; v=[(0.2,0.2,20),(200.2,0.2,20),(200.2,200.2,20),(0.2,200.2,20)]; f=[(0,1,2),(0,2,3)]; print(json.dumps(find_anchor(v,f,0.01)))`,
  );
  expect(result.status).toBe(0);
  expect(JSON.parse(result.stdout)).toBeNull();
});
