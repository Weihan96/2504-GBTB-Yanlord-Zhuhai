import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("RA.LP evidence preserves its installation centre while rounding XY", () => {
  const source = `
import json,sys,numpy as np
sys.path.insert(0,"pipeline/scripts")
from light_fixture_anchor_candidate import light_fixture_evidence
r=light_fixture_evidence(
  "RA.LP","DIRECTIONSOURCE",
  np.array([5948.242188,1395.736814,2400.0]),
  np.array([5920.742188,1368.236814,2397.0]),
  np.array([5975.742188,1423.236814,2486.0]),
)
print(json.dumps(r))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.pass).toBe(true);
  expect(report.target_mm).toEqual([5948, 1396, 2400]);
  expect(report.translation_mm[0]).toBeCloseTo(-0.242188, 6);
  expect(report.translation_mm[1]).toBeCloseTo(0.263186, 6);
});

test("light evidence rejects a generic bbox-centre guess", () => {
  const source = `
import sys,numpy as np
sys.path.insert(0,"pipeline/scripts")
from light_fixture_anchor_candidate import light_fixture_evidence
r=light_fixture_evidence(
  "Other","DIRECTIONSOURCE",
  np.array([10.2,20.2,30.0]),
  np.array([0.0,0.0,0.0]),
  np.array([55.0,55.0,89.0]),
)
print(r["pass"])
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(result.stdout.toString().trim()).toBe("False");
});
