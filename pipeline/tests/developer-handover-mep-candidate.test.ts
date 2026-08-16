import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/developer_handover_mep_candidate.py");
const output = resolve(root, "build/mep-positioning/developer-handover-mep-candidate.test.json");
const elecSvg = resolve(root, "build/mep-positioning/E302-E304-developer-handover-reference.test.svg");
const plumSvg = resolve(root, "build/mep-positioning/P201-developer-handover-reference.test.svg");

test("developer handover MEP points map into the current IFC coordinate system without write authority", () => {
  const run = spawnSync([
    "python3",
    script,
    "--output",
    output,
    "--elec-svg",
    elecSvg,
    "--plum-svg",
    plumSvg,
  ], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.point_counts).toEqual({
    developer_electrical_point: 66,
    developer_floor_drain: 6,
    developer_switch_or_control: 11,
    developer_water_or_wc_drain: 20,
  });
  expect(report.records).toHaveLength(103);
  expect(report.gates.expected_point_counts_pass).toBe(true);
  expect(report.gates.source_matches_dwg_conversion).toBe(true);
  expect(report.gates.all_points_have_space_candidate).toBe(true);
  expect(report.gates.developer_reference_is_renovation_design).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  const renderedElec = readFileSync(elecSvg, "utf8");
  expect(renderedElec).toContain("developer-elec-reference");
  expect((renderedElec.match(/data-coordination-kind="door"/g) ?? []).length).toBe(8);
  expect((renderedElec.match(/data-coordination-kind="fixed_furniture"/g) ?? []).length).toBe(70);
  expect(readFileSync(plumSvg, "utf8")).toContain("developer-plum-reference");
}, 60_000);

test("developer handover candidate has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});

test("developer handover rejects a mismatched caller-frozen hash", () => {
  const run = spawnSync(["python3", script, "--expected-ifc-sha256", "0".repeat(64)], { cwd: root });
  expect(run.exitCode).not.toBe(0);
  expect(run.stderr.toString()).toContain("formal IFC hash changed");
});
