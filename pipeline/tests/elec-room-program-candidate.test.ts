import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_room_program_candidate.py");
const output = resolve(root, "build/elec/elec-room-program-candidate.test.json");
const svg = resolve(root, "build/elec/E302-E304-room-program-candidate.test.svg");

test("all rooms and modelled equipment receive read-only electrical program coverage", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.room_count).toBe(22);
  expect(report.summary.lighting_control_groups_min).toBe(23);
  expect(report.summary.general_power_groups_min).toBe(26);
  expect(report.summary.network_data_groups_min).toBe(8);
  expect(report.summary.manual_dedicated_power_min).toBe(9);
  expect(report.summary.cabinet_light_feed_zones_min).toBe(7);
  expect(report.summary.modelled_equipment_power_demands).toBe(16);
  expect(report.summary.formal_switch_instances).toBe(0);
  expect(report.summary.formal_network_instances).toBe(0);
  expect(report.design_rule_summary).toEqual({
    physical_wired_switch_only: true,
    panel_bottom_AFF_mm: {
      ordinary_socket: 300,
      television_point: 600,
      physical_switch: 1300,
    },
    bedroom_AP_count: 2,
    living_study_router_count: 2,
    paired_two_way_control_groups: 3,
    entrance_master_lighting_switches: 1,
    confirmed_island_dishwashers: 2,
    demolition_walls_hidden_in_general_reviews: true,
  });
  expect(report.rooms.map((row: { room_reference: string }) => row.room_reference)).toEqual(
    Array.from({ length: 22 }, (_, index) => `R${String(index + 1).padStart(2, "0")}`),
  );
  expect(report.gates.all_22_spaces_programmed).toBe(true);
  expect(report.gates.all_equipment_demands_have_evidence).toBe(true);
  expect(report.gates.developer_references_are_not_final_design).toBe(true);
  expect(report.gates.confirmed_design_rules_are_complete).toBe(true);
  expect(report.gates.equipment_centres_are_not_socket_positions).toBe(true);
  expect(report.gates.whole_home_switch_positioning_complete).toBe(false);
  expect(report.gates.whole_home_network_positioning_complete).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(readFileSync(svg, "utf8")).toContain("E-302/E-304 逐房间用电与弱电功能程序候选");
}, 30_000);

test("room electrical program compiler has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
