import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/mep_positioning_source_audit.py");
const output = resolve(root, "build/mep-positioning/source-audit.test.json");

test("MEP source audit registers the latest kitchen drawing without authorizing IFC writes", () => {
  const run = spawnSync(["python3", script, "--output", output], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.source_drawings.older.pages).toBe(10);
  expect(report.source_drawings.latest.pages).toBe(10);
  expect(report.source_drawings.mechanical_page_difference).toEqual([6]);
  expect(report.source_requirements).toHaveLength(8);
  expect(report.current_ifc_inventory).toEqual({
    plum_registered_terminals: 27,
    plum_service_demand_candidates: 24,
    plum_non_service_components: 3,
    plum_existing_location_objects: 49,
    lights: 79,
    sockets: 11,
    typed_equipment: 8,
    unresolved_electrical_proxies: 9,
    switch_instances: 0,
    network_instances: 0,
  });
  expect(report.gates.latest_kitchen_source_registered).toBe(true);
  expect(report.gates.current_ifc_reports_match).toBe(true);
  expect(report.gates.whole_home_water_positioning_complete).toBe(false);
  expect(report.gates.whole_home_electrical_positioning_complete).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
}, 20_000);

test("MEP source audit contains no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});

test("MEP source audit rejects a mismatched caller-frozen hash", () => {
  const run = spawnSync(["python3", script, "--expected-ifc-sha256", "0".repeat(64)], { cwd: root });
  expect(run.exitCode).not.toBe(0);
  expect(run.stderr.toString()).toContain("formal IFC hash changed");
});
