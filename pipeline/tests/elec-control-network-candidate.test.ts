import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_control_network_candidate.py");
const output = resolve(root, "build/elec/elec-control-network-candidate.test.json");
const svg = resolve(root, "build/elec/E302-E304-control-network-candidate.test.svg");

test("confirmed control and network roles compile as read-only coordination zones", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.control_coordination_zones).toBe(2);
  expect(report.summary.paired_two_way_control_groups).toBe(3);
  expect(report.summary.entrance_master_lighting_switches).toBe(1);
  expect(report.summary.bedroom_AP_zones).toBe(2);
  expect(report.summary.living_study_router_zones).toBe(1);
  expect(Object.values(report.gates).every((value) => value === true || value === false)).toBe(true);
  expect(report.gates.two_doorway_zones_present).toBe(true);
  expect(report.gates.three_two_way_groups_present).toBe(true);
  expect(report.gates.entrance_master_switch_present).toBe(true);
  expect(report.gates.two_bedroom_AP_zones_present).toBe(true);
  expect(report.gates.one_shared_router_no_AP_zone_present).toBe(true);
  expect(report.gates.all_positions_are_coordination_zones).toBe(true);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(readFileSync(svg, "utf8")).toContain("elec-control-network");
});

test("control/network candidate has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
