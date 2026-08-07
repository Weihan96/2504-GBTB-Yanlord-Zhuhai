import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

const script = "pipeline/scripts/beam_relationship_audit.py";

function runPython(expression: string) {
  return spawnSync("python3", ["-c", expression], {
    cwd: process.cwd(),
    encoding: "utf8",
  });
}

test("beam relationship audit classifies tolerance transitions", () => {
  const result = runPython(
    `import sys; sys.path.insert(0, "pipeline/scripts"); from beam_relationship_audit import transition; print(transition(0.05,0.2,0.1)); print(transition(0.2,0.05,0.1)); print(transition(0.05,0.08,0.1))`,
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim().split("\n")).toEqual([
    "regression",
    "improvement",
    "retained_within_tolerance",
  ]);
});

test("beam relationship audit distinguishes opposing and aligned faces", () => {
  const result = runPython(
    `import json,sys; sys.path.insert(0, "pipeline/scripts"); from beam_relationship_audit import face_metrics; print(json.dumps(face_metrics((0,0,0,200,1000,400),(200,0,0,400,1000,400),0,10)))`,
  );
  expect(result.status).toBe(0);
  const metric = JSON.parse(result.stdout);
  expect(metric.opposing_face).toBe("target_max_to_neighbor_min");
  expect(metric.opposing_residual_mm).toBe(0);
  expect(metric.aligned_residual_mm).toBe(200);
});

test("beam relationship audit rejects face comparisons without cross overlap", () => {
  const result = runPython(
    `import sys; sys.path.insert(0, "pipeline/scripts"); from beam_relationship_audit import face_metrics; print(face_metrics((0,0,0,200,100,100),(200,200,0,400,300,100),0,10))`,
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim()).toBe("None");
});
