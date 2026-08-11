import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/a106_ceiling_device_candidate.py");
const output = resolve(root, "build/elec/a106-ceiling-device-candidate.test.json");
const svg = resolve(root, "build/elec/A106-ceiling-device-candidate.test.svg");

test("A-106 ceiling device candidates pass known-geometry gates", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.candidate_count).toBe(6);
  expect(report.summary.smoke_alarm_candidates).toBe(3);
  expect(report.summary.bedroom_AP_candidates).toBe(2);
  expect(report.summary.kitchen_fire_sensor_candidates).toBe(1);
  expect(report.gates.known_geometry_pass).toBe(true);
  expect(report.gates.candidate_xy_coordinates_are_integer_mm).toBe(true);
  expect(report.gates.candidate_xyz_coordinates_are_integer_mm).toBe(true);
  expect(report.gates.all_candidates_inside_space).toBe(true);
  expect(report.gates.all_candidates_on_ceiling).toBe(true);
  expect(report.gates.all_ap_light_edge_clearances_pass).toBe(true);
  expect(report.gates.all_ap_wall_beam_high_obstacle_clearances_pass).toBe(true);
  expect(report.gates.all_ap_same_room_smoke_clearances_pass).toBe(true);
  expect(report.gates.all_ap_mechanical_checks_pass).toBe(true);
  expect(report.gates.bedroom_pair_separation_pass).toBe(true);
  expect(report.gates.thirteen_demolition_walls_excluded).toBe(true);
  expect(report.gates.one_kitchen_fire_sensor_present).toBe(true);
  const smokeCandidates = report.candidates.filter((row: { device_role: string }) => row.device_role === "smoke_alarm");
  expect(smokeCandidates.every((row: { room_center_offset_mm: number }) => row.room_center_offset_mm <= 522.1)).toBe(true);
  expect(report.candidates.every((row: { nearest_high_level_obstacle_clearance_mm: number | null }) => row.nearest_high_level_obstacle_clearance_mm === null || row.nearest_high_level_obstacle_clearance_mm >= 500)).toBe(true);
  const primaryBedroomAp = report.candidates.find((row: { candidate_id: string }) => row.candidate_id === "A106-AP-R09");
  expect(primaryBedroomAp.position_mm).toEqual([-5840, -2844, 2720]);
  expect(primaryBedroomAp.ap_mechanical_assessment.candidate_reason).toContain("灯网空交点");
  expect(primaryBedroomAp.ap_mechanical_assessment.light_axis_alignment.matched_axis_count).toBe(2);
  expect(primaryBedroomAp.ap_mechanical_assessment.clearance_margins_mm.light_edge).toBeGreaterThan(500);
  expect(primaryBedroomAp.ap_mechanical_assessment.clearance_margins_mm.wall_boundary).toBeGreaterThan(250);
  expect(primaryBedroomAp.ap_mechanical_assessment.clearance_margins_mm.beam).toBeGreaterThan(1000);
  expect(primaryBedroomAp.ap_mechanical_assessment.clearance_margins_mm.same_room_smoke).toBeGreaterThan(700);
  expect(primaryBedroomAp.ap_mechanical_assessment.placement_checks).toEqual({
    inside_space: true,
    on_ceiling: true,
    integer_xyz_mm: true,
  });
  expect(report.gates.all_smoke_supply_air_clearances_verified).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  const renderedSvg = readFileSync(svg, "utf8");
  expect(renderedSvg).toContain("a106-ceiling-device-candidate");
  expect(renderedSvg).toContain("A106-FIRE-R04");
  expect(renderedSvg).toContain("Wall Plan-underlay.png");
  const wallGroupTags = renderedSvg.match(/<g\b[^>]*\bclass="[^"]*\bIfcWall\b[^"]*"[^>]*>/g) ?? [];
  for (const globalId of report.excluded_demolition_wall_global_ids) {
    const demolitionGroupTags = wallGroupTags.filter((tag) => tag.includes(globalId));
    expect(demolitionGroupTags.length).toBeGreaterThan(0);
    expect(demolitionGroupTags.every((tag) => tag.includes("a106-excluded-demolish"))).toBe(true);
  }
}, 300_000);

test("A-106 candidate has no formal IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("ifcopenshell.api.run");
  expect(source).not.toMatch(/\.write\s*\(/);
});
