import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/a106_ceiling_device_candidate.py");
const output = resolve(root, "build/elec/a106-ceiling-device-candidate.test.json");
const svg = resolve(root, "build/elec/A106-ceiling-device-candidate.test.svg");

test("A-106 smoke and AP candidates pass known-geometry gates", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.candidate_count).toBe(5);
  expect(report.summary.smoke_alarm_candidates).toBe(3);
  expect(report.summary.bedroom_AP_candidates).toBe(2);
  expect(report.gates.known_geometry_pass).toBe(true);
  expect(report.gates.candidate_xy_coordinates_are_50mm_modular).toBe(true);
  expect(report.gates.bedroom_pair_separation_pass).toBe(true);
  expect(report.candidates.every((row: { nearest_high_level_obstacle_clearance_mm: number | null }) => row.nearest_high_level_obstacle_clearance_mm === null || row.nearest_high_level_obstacle_clearance_mm >= 500)).toBe(true);
  expect(report.gates.all_smoke_supply_air_clearances_verified).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(readFileSync(svg, "utf8")).toContain("a106-ceiling-device-candidate");
}, 60_000);

test("A-106 candidate has no formal IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("ifcopenshell.api.run");
  expect(source).not.toMatch(/\.write\s*\(/);
});
