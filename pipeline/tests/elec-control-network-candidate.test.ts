import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_control_network_candidate.py");
const output = resolve(root, "build/elec/elec-control-network-candidate.test.json");
const svg = resolve(root, "build/elec/E302-E304-control-network-candidate.test.svg");

test("confirmed control, network, and safety roles compile as read-only coordination zones", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.control_coordination_zones).toBe(2);
  expect(report.summary.paired_two_way_control_groups).toBe(3);
  expect(report.summary.entrance_master_lighting_switches).toBe(1);
  expect(report.summary.bedroom_AP_candidates).toBe(2);
  expect(report.summary.living_study_router_zones).toBe(1);
  expect(report.summary.smoke_alarm_candidates).toBe(3);
  expect(report.summary.kitchen_fire_sensor_positions).toBe(1);
  expect(report.summary.kitchen_gas_alarm_room_zones).toBe(1);
  expect(Object.values(report.gates).every((value) => value === true || value === false)).toBe(true);
  expect(report.gates.two_doorway_zones_present).toBe(true);
  expect(report.gates.three_two_way_groups_present).toBe(true);
  expect(report.gates.entrance_master_switch_present).toBe(true);
  expect(report.gates.two_bedroom_AP_candidates_present).toBe(true);
  expect(report.gates.one_shared_router_no_AP_zone_present).toBe(true);
  expect(report.gates.three_smoke_candidates_present).toBe(true);
  expect(report.gates.smoke_positioning_constraints_present).toBe(true);
  expect(report.gates.one_formal_kitchen_fire_position_present).toBe(true);
  expect(report.gates.one_kitchen_gas_alarm_zone_present).toBe(true);
  expect(report.gates.gas_alarm_model_authority_pending).toBe(true);
  expect(report.gates.a106_exact_positions_imported).toBe(true);
  expect(report.gates.unconfirmed_positions_remain_coordination_zones).toBe(true);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  const renderedSvg = readFileSync(svg, "utf8");
  expect(renderedSvg).toContain("elec-control-network");
  expect(renderedSvg).toContain("A106-FIRE-R04");
  expect(renderedSvg).not.toContain("Wall Plan-underlay.png");
}, 30_000);

test("control/network candidate has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
