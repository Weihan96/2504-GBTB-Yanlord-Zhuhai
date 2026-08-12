import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/equipment_ssot.py");

test("equipment SSOT validates projections and covers every scoped IFC object", () => {
  const run = Bun.spawnSync(["python3", script, "all"], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(run.stdout.toString());
  expect(report.validate).toMatchObject({
    master_count: 146,
    requirement_count: 396,
    source_count: 89,
    schema_version: "1.0.0",
  });
  expect(report.projections).toMatchObject({
    appliance_rows: 18,
    furniture_rows: 9,
    elec_evidence_rows: 24,
    hvac_evidence_rows: 4,
    furniture_role_rows: 51,
    mode: "check",
  });
  expect(report.ifc_coverage).toMatchObject({
    scope_count: 162,
    covered_count: 162,
    missing_count: 0,
    duplicate_count: 0,
  });
  expect(report.ifc_coverage.class_counts).toEqual({
    IfcElectricAppliance: 19,
    IfcFurniture: 89,
    IfcSanitaryTerminal: 27,
    IfcWasteTerminal: 3,
    IfcSensor: 1,
    IfcDoor: 8,
    IfcWindow: 11,
    IfcElementAssembly: 3,
    IfcBuildingElementProxy: 1,
  });
}, 180_000);

test("equipment SSOT never writes the formal IFC", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifc.write(");
});

test("equipment SSOT can audit an IFC evidence hash refresh without writing", () => {
  const sourcePath = resolve(root, "pipeline/decisions/source-evidence-register.csv");
  const before = readFileSync(sourcePath);
  const run = Bun.spawnSync(["python3", script, "refresh-ifc-source-hashes", "--dry-run"], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const report = JSON.parse(run.stdout.toString());
  expect(report.scoped_ifc_object_count).toBe(162);
  expect(report.target_source_count).toBeGreaterThan(1);
  expect(report.dry_run).toBe(true);
  expect(readFileSync(sourcePath)).toEqual(before);
}, 30_000);
