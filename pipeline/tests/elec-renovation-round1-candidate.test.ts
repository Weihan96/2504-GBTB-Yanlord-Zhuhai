import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_renovation_round1_candidate.py");
const output = resolve(root, "build/elec/elec-renovation-round1-candidate.test.json");
const svg = resolve(root, "build/elec/E301-E303-renovation-round1-candidate.test.svg");

test("first-round renovation electrical demands stay read-only and complete", () => {
  const run = spawnSync([
    "python3",
    script,
    "--output",
    output,
    "--output-svg",
    svg,
  ], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.bedside_light_candidates).toBe(4);
  expect(report.summary.new_socket_candidates).toBe(2);
  expect(report.summary.cabinet_power_zones).toBe(7);
  expect(report.summary.kitchen_socket_rechecks).toBe(11);
  expect(report.summary.label_collision_count).toBe(0);
  expect(report.summary.developer_red_points_are_reference_only).toBe(true);
  expect(report.gates.four_bedside_lights_present).toBe(true);
  expect(report.gates.island_and_dining_bay_socket_present).toBe(true);
  expect(report.gates.illuminated_cabinet_power_is_grouped_not_fabricated).toBe(true);
  expect(report.gates.all_current_kitchen_sockets_reopened_for_review).toBe(true);
  expect(report.gates.label_collision_free).toBe(true);
  expect(report.gates.reverse_requirement_audit_planned).toBe(true);
  expect(report.gates.whole_home_electrical_positioning_complete).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(readFileSync(svg, "utf8")).toContain("elec-renovation-round1");
}, 30_000);

test("first-round renovation candidate has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
