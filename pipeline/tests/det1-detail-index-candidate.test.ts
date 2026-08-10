import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { mkdtempSync, readFileSync, statSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/det1_detail_index_candidate.py");
const ifc = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");

function sha256(path: string): string {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("DET1 candidate records current hashes, unresolved variables, and PNG proof", async () => {
  const temporary = mkdtempSync(join(tmpdir(), "det1-index-"));
  const reportPath = join(temporary, "detail-index-candidate.json");
  const reviewCsv = join(temporary, "det1-detail-review.csv");
  const svg = join(temporary, "D601-D602-detail-index-candidate.svg");
  const png = join(temporary, "D601-D602-detail-index-candidate.png");
  const issues = join(temporary, "open-issues.md");
  const process = Bun.spawn([
    "python3", script,
    "--root", root,
    "--review-csv", reviewCsv,
    "--output-svg", svg,
    "--report", reportPath,
    "--open-issues", issues,
    "--proof-png", png,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  const [exitCode, stderr] = await Promise.all([
    process.exited,
    new Response(process.stderr).text(),
  ]);
  expect(exitCode).toBe(0);
  expect(stderr).not.toContain("Traceback");

  const report = await Bun.file(reportPath).json();
  expect(report.source.ifc_sha256).toBe(sha256(ifc));
  for (const input of Object.values(report.source.inputs) as Array<{ path: string; sha256: string }>) {
    expect(input.sha256).toBe(sha256(input.path));
  }
  for (const output of Object.values(report.outputs) as Array<{ path: string; sha256: string }>) {
    expect(output.sha256).toBe(sha256(output.path));
  }
  expect(report.summary).toMatchObject({
    node_count: 12,
    d601_node_count: 6,
    d602_node_count: 6,
    unresolved_field_count: 36,
    invented_dimension_count: 0,
  });
  expect(report.records.every((row: any) =>
    row.variable_parameters &&
    row.unresolved_material_or_product &&
    row.unresolved_manufacturer &&
    row.unresolved_construction &&
    row.review_status === "candidate_pending_review" &&
    row.automatic_ifc_write_allowed === "false" &&
    row.construction_release_ready === "false"
  )).toBe(true);
  expect(report.gates).toMatchObject({
    input_hashes_current: true,
    unresolved_items_explicit: true,
    svg_nonempty: true,
    png_nonempty: true,
    automatic_ifc_write_allowed: false,
    construction_release_ready: false,
  });
  expect(statSync(svg).size).toBeGreaterThan(1_000);
  expect(statSync(png).size).toBeGreaterThan(10_000);
  expect(readFileSync(png).subarray(0, 8).toString("hex")).toBe("89504e470d0a1a0a");
  expect(readFileSync(svg, "utf8")).toContain("不写 IFC");
  expect(readFileSync(reviewCsv, "utf8")).not.toContain("\r\n");
}, 30_000);

test("DET1 candidate rejects a stale caller-frozen IFC hash", () => {
  const temporary = mkdtempSync(join(tmpdir(), "det1-stale-"));
  const run = Bun.spawnSync([
    "python3", script,
    "--root", root,
    "--expected-ifc-sha256", "0".repeat(64),
    "--review-csv", join(temporary, "review.csv"),
    "--output-svg", join(temporary, "candidate.svg"),
    "--report", join(temporary, "report.json"),
    "--open-issues", join(temporary, "issues.md"),
    "--proof-png", join(temporary, "proof.png"),
  ], { cwd: root });
  expect(run.exitCode).not.toBe(0);
  expect(run.stderr.toString()).toContain("formal IFC SHA mismatch");
});

test("DET1 generator has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
  expect(source).not.toContain("save_ifc_file");
  expect(source).not.toContain("bpy");
});
