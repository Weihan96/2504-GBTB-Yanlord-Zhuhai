import { expect, test } from "bun:test";
import { CryptoHasher } from "bun";
import {
  existsSync,
  mkdtempSync,
  readFileSync,
  readdirSync,
  writeFileSync,
} from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = join(root, "pipeline/scripts/a101_survey_input.py");
const template = join(root, "pipeline/templates/a101-survey-input.csv");
const formalIfc = join(root, "2504 GBTB Yanlord Zhuhai.ifc");

function sha256(path: string): string {
  const hasher = new CryptoHasher("sha256");
  hasher.update(readFileSync(path));
  return hasher.digest("hex");
}

function run(args: string[] = [], cwd = root) {
  return Bun.spawnSync(["python3", script, ...args], {
    cwd,
    stdout: "pipe",
    stderr: "pipe",
  });
}

function confirmedCsv(): string {
  const lines = readFileSync(template, "utf8").trimEnd().split("\n");
  const headers = lines[0].split(",");
  const index = Object.fromEntries(headers.map((header, column) => [header, column]));
  const rows = lines.slice(1).map((line) => {
    const values = line.split(",");
    const property = values[index.field_type] === "immovable_property_condition";
    values[index.value] = property ? "现场确认无新增不可移动物业条件" : "3000";
    values[index.status] = "confirmed";
    values[index.measurement_method] = property ? "visual_inspection" : "laser_distance_meter";
    values[index.observed_by] = "现场测量员";
    values[index.observed_at] = "2026-08-11";
    values[index.evidence_reference] = "evidence/A101-photo-set";
    values[index.confirmed_by] = "业主代表";
    values[index.confirmed_at] = "2026-08-11";
    return values.join(",");
  });
  return [lines[0], ...rows].join("\n") + "\n";
}

test("blank A-101 template is a valid pending list and remains construction-not-ready", () => {
  const temp = mkdtempSync(join(tmpdir(), "a101-pending-"));
  const beforeHash = sha256(formalIfc);
  const result = run([], temp);
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.schema_version).toBe("a101-survey-input/v1");
  expect(report.formal_ifc_write).toBe(false);
  expect(report.stdout_only_by_default).toBe(true);
  expect(report.source_ifc_sha256).toBe(beforeHash);
  expect(report.construction_ready).toBe(false);
  expect(report.summary.status_counts).toEqual({
    confirmed: 0,
    observed: 0,
    pending: 6,
  });
  expect(report.inputs.every((row: { effective_value: unknown }) => row.effective_value === null)).toBe(true);
  expect(readdirSync(temp)).toEqual([]);
  expect(sha256(formalIfc)).toBe(beforeHash);
});

test("confirmed survey rows compile normalized JSON and CSV only to explicit paths", () => {
  const temp = mkdtempSync(join(tmpdir(), "a101-confirmed-"));
  const input = join(temp, "survey.csv");
  const outputJson = join(temp, "normalized.json");
  const outputCsv = join(temp, "normalized.csv");
  const frozenHash = sha256(formalIfc);
  writeFileSync(input, confirmedCsv());
  const result = run([
    "--input", input,
    "--expected-ifc-sha256", frozenHash,
    "--output-json", outputJson,
    "--output-csv", outputCsv,
  ]);
  expect(result.exitCode).toBe(0);
  expect(existsSync(outputJson)).toBe(true);
  expect(existsSync(outputCsv)).toBe(true);
  const report = JSON.parse(readFileSync(outputJson, "utf8"));
  expect(report.caller_frozen_ifc_sha256).toBe(frozenHash);
  expect(report.construction_ready).toBe(true);
  expect(report.summary.blocking_input_ids).toEqual([]);
  expect(report.inputs.every((row: { effective_value: unknown }) => row.effective_value !== null)).toBe(true);
  expect(readFileSync(outputCsv, "utf8")).toContain("effective_value,effective_unit");
});

test("observed values remain recorded but never become construction-effective", () => {
  const temp = mkdtempSync(join(tmpdir(), "a101-observed-"));
  const lines = confirmedCsv().trimEnd().split("\n");
  const observed = lines.map((line, row) => row === 0
    ? line
    : line.replace(",confirmed,", ",observed,").replace(",业主代表,2026-08-11,", ",,,"));
  const input = join(temp, "observed.csv");
  writeFileSync(input, observed.join("\n") + "\n");
  const result = run(["--input", input]);
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.construction_ready).toBe(false);
  expect(report.inputs.every((row: { recorded_value: unknown }) => row.recorded_value !== null)).toBe(true);
  expect(report.inputs.every((row: { effective_value: unknown }) => row.effective_value === null)).toBe(true);
});

test("invalid unit is rejected before any report write", () => {
  const temp = mkdtempSync(join(tmpdir(), "a101-unit-"));
  const input = join(temp, "invalid.csv");
  const output = join(temp, "should-not-exist.json");
  writeFileSync(input, readFileSync(template, "utf8").replace(",mm,pending,", ",cm,pending,"));
  const result = run(["--input", input, "--output-json", output]);
  expect(result.exitCode).toBe(2);
  expect(result.stderr.toString()).toContain("unit must be 'mm'");
  expect(existsSync(output)).toBe(false);
});

test("unknown status is rejected", () => {
  const temp = mkdtempSync(join(tmpdir(), "a101-status-"));
  const input = join(temp, "invalid.csv");
  writeFileSync(input, readFileSync(template, "utf8").replace(",pending,", ",draft,"));
  const result = run(["--input", input]);
  expect(result.exitCode).toBe(2);
  expect(result.stderr.toString()).toContain("invalid status 'draft'");
});

test("missing schema field is rejected", () => {
  const temp = mkdtempSync(join(tmpdir(), "a101-header-"));
  const input = join(temp, "invalid.csv");
  const rows = readFileSync(template, "utf8").trimEnd().split("\n");
  writeFileSync(input, rows.map((row) => row.split(",").slice(0, -1).join(",")).join("\n") + "\n");
  const result = run(["--input", input]);
  expect(result.exitCode).toBe(2);
  expect(result.stderr.toString()).toContain("input headers changed; expected exact schema");
});

test("caller-frozen IFC mismatch is a hard failure and IFC remains unchanged", () => {
  const beforeHash = sha256(formalIfc);
  const result = run(["--expected-ifc-sha256", "0".repeat(64)]);
  expect(result.exitCode).toBe(2);
  expect(result.stderr.toString()).toContain("IFC SHA-256 mismatch");
  expect(sha256(formalIfc)).toBe(beforeHash);
  expect(readFileSync(script, "utf8")).not.toContain("ifcopenshell");
});
