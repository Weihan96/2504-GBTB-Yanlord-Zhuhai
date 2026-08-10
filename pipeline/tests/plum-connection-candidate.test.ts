import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/plum_connection_candidate.py");
const output = resolve(root, "build/plum/plum-connection-candidate.test.json");

test("PLUM connection candidate pairs current service objects only as read-only evidence", () => {
  const run = spawnSync(["python3", script, "--output", output], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.registered_sanitary_terminals).toBe(27);
  expect(report.summary.service_demand_candidates).toBe(24);
  expect(report.summary.non_service_components).toBe(3);
  expect(report.summary.derived_drainage_components).toBe(20);
  expect(report.summary.controlled_pvc110_branches).toBe(6);
  expect(report.records).toHaveLength(24);
  expect(report.records.every((row: any) => row.review_required)).toBe(true);
  expect(report.records.every((row: any) => !row.automatic_ifc_write_allowed)).toBe(true);
  expect(report.gates.all_service_candidates_have_nearest_branch_evidence).toBe(true);
  expect(report.gates.derived_points_are_formal_connectors).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
}, 120_000);

test("PLUM proximity script cannot write IFC", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});

test("PLUM proximity candidate rejects a mismatched caller-frozen hash", () => {
  const run = spawnSync([
    "python3", script, "--expected-ifc-sha256", "0".repeat(64), "--output", output,
  ], { cwd: root });
  expect(run.exitCode).not.toBe(0);
  expect(run.stderr.toString()).toContain("formal IFC hash changed");
});
