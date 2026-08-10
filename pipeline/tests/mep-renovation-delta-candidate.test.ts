import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/mep_renovation_delta_candidate.py");
const output = resolve(root, "build/mep-positioning/mep-renovation-delta-candidate.test.json");

test("MEP delta distinguishes developer references from the current renovation design", () => {
  const run = spawnSync(["python3", script, "--output", output], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.developer_electrical_and_control_points).toBe(77);
  expect(report.summary.developer_plumbing_points).toBe(26);
  expect(report.summary.rooms).toBe(22);
  expect(report.summary.current_switch_instances).toBe(0);
  expect(report.summary.current_network_instances).toBe(0);
  expect(report.gates.all_22_spaces_programmed).toBe(true);
  expect(report.gates.developer_points_are_current_design).toBe(false);
  expect(report.gates.rough_in_connectors_inferred).toBe(false);
  expect(report.gates.whole_home_positioning_complete).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
}, 20_000);

test("MEP delta candidate has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});

test("MEP delta rejects a mismatched caller-frozen hash", () => {
  const run = spawnSync(["python3", script, "--expected-ifc-sha256", "0".repeat(64)], { cwd: root });
  expect(run.exitCode).not.toBe(0);
  expect(run.stderr.toString()).toContain("formal IFC hash changed");
});
