import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_positioning_candidate.py");
const output = resolve(root, "build/elec/elec-positioning-candidate.test.json");

test("ELEC positioning maps kitchen points to source-supported roles without writing IFC", () => {
  const run = spawnSync(["python3", script, "--output", output], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.existing_lights).toBe(79);
  expect(report.summary.existing_sockets).toBe(11);
  expect(report.summary.typed_equipment).toBe(8);
  expect(report.summary.proxy_handoffs).toBe(9);
  expect(report.summary.switch_instances).toBe(0);
  expect(report.summary.network_instances).toBe(0);
  expect(report.socket_candidates).toHaveLength(11);
  expect(report.proxy_identity_candidates).toHaveLength(9);
  const dishwashers = report.proxy_identity_candidates.filter(
    (row: { candidate_role: string }) => row.candidate_role === "island_dishwasher",
  );
  expect(dishwashers).toHaveLength(2);
  expect(dishwashers.every((row: { confidence: number }) => row.confidence === 1)).toBe(true);
  expect(report.gates.all_existing_kitchen_sockets_have_source_role_candidates).toBe(true);
  expect(report.gates.all_proxy_handoffs_have_identity_candidates).toBe(true);
  expect(report.gates.whole_home_socket_positioning_complete).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(report.gates.construction_release_ready).toBe(false);
});

test("ELEC positioning script contains no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});

test("ELEC positioning rejects a mismatched caller-frozen hash", () => {
  const run = spawnSync(["python3", script, "--expected-ifc-sha256", "0".repeat(64)], { cwd: root });
  expect(run.exitCode).not.toBe(0);
  expect(run.stderr.toString()).toContain("formal IFC hash changed");
});
