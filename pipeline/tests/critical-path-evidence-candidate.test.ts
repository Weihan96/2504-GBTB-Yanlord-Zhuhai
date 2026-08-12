import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/critical_path_evidence_candidate.py");

test("ELEC/RCP1 critical-path evidence validates current reports and publish hashes", () => {
  const temporary = mkdtempSync(join(tmpdir(), "critical-path-evidence-"));
  const output = join(temporary, "report.json");
  const result = Bun.spawnSync(["python3", script, "--root", root, "--output", output]);
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.scope_sheet_count).toBe(6);
  expect(report.summary.critical_report_count).toBe(9);
  expect(report.summary.current_critical_report_count).toBe(9);
  expect(report.summary.publish_bundle_count).toBe(4);
  expect(report.summary.mechanical_publish_bundle_pass_count).toBe(4);
  expect(report.summary.excluded_preview_count).toBe(4);
  expect(report.summary.stale_excluded_preview_count).toBeGreaterThan(0);
  expect(report.gates.reviewed_candidate_evidence_pass).toBe(true);
  expect(report.gates.construction_release_ready).toBe(false);
  expect(report.gates.formal_ifc_write_allowed).toBe(false);
  expect(report.publish_bundles.every((bundle: any) =>
    Object.values(bundle.files).every((file: any) => file.inside_project)
  )).toBe(true);
  expect(report.excluded_previews.every((row: any) =>
    !row.registered_for_release && !row.used_as_publish_target && !row.eligible_for_release
  )).toBe(true);
});

test("critical-path evidence fails closed on a publish hash mismatch", () => {
  const temporary = mkdtempSync(join(tmpdir(), "critical-path-stale-"));
  const render = JSON.parse(readFileSync(
    resolve(root, "build/elec/E-301-E-303-renovation-round1-render.json"),
    "utf8",
  ));
  render.output_pdf_sha256 = "0".repeat(64);
  const badRender = join(temporary, "bad-render.json");
  const output = join(temporary, "report.json");
  writeFileSync(badRender, JSON.stringify(render));
  const result = Bun.spawnSync([
    "python3",
    script,
    "--root",
    root,
    "--e301-e303-render",
    badRender,
    "--output",
    output,
  ]);
  expect(result.exitCode).toBe(1);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.gates.all_publish_bundles_mechanically_match).toBe(false);
  expect(report.gates.reviewed_candidate_evidence_pass).toBe(false);
});

test("critical-path evidence has no Blender or IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("import bpy");
  expect(source).not.toContain("ifcopenshell");
  expect(source).not.toContain("model.write(");
  expect(source).toContain('"formal_ifc_write_allowed": False');
});
