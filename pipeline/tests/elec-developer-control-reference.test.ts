import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_developer_control_reference.py");
const output = resolve(root, "build/elec/elec-developer-control-reference.test.json");

test("developer HVAC and access controls remain visible references pending field confirmation", () => {
  const run = spawnSync(["python3", script, "--output", output], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary).toEqual({ hvac_control_panels: 5, video_intercoms: 1, doorbells: 1 });
  expect(report.hvac_control_panel_references).toHaveLength(5);
  expect(report.hvac_control_panel_references.every((row: any) => row.installation_height_mm === 1300)).toBe(true);
  expect(report.access_control_references.map((row: any) => row.device_role)).toEqual(["video_intercom", "doorbell"]);
  expect(report.access_control_references[0].position_mm).toEqual([3066.658536, -634.998805, 1400]);
  expect(report.access_control_references[1].position_mm).toEqual([6053.880075, -664.974809, 1300]);
  expect(report.gates.all_references_pending_field_confirmation).toBe(true);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
});
