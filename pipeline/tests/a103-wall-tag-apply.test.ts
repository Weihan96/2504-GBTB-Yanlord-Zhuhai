import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/a103_wall_tag_apply.py");
const sourceIfc = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");

function buildTaglessFixture(temporary: string) {
  const baseline = join(temporary, "tagless.ifc");
  const candidate = join(temporary, "candidate-report.json");
  const result = Bun.spawnSync([
    "python3", "-c", `
import hashlib, ifcopenshell, json, pathlib, sys
source, original_report, baseline, candidate = map(pathlib.Path, sys.argv[1:])
data = json.loads(original_report.read_text(encoding="utf-8"))
model = ifcopenshell.open(source)
for record in data["records"]:
    model.by_guid(record["global_id"]).Tag = None
model.write(baseline)
digest = hashlib.sha256(baseline.read_bytes()).hexdigest()
data["source_ifc_sha256"] = digest
data["pass"] = True
candidate.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
print(digest)
`, sourceIfc, resolve(root, "build/a103/a103-wall-tag-candidate.json"), baseline, candidate,
  ], { cwd: root });
  expect(result.exitCode).toBe(0);
  return { baseline, candidate, sourceHash: result.stdout.toString().trim() };
}

test("A-103 formal apply rejects an absent human approval", () => {
  const temporary = mkdtempSync(join(tmpdir(), "a103-tag-reject-"));
  try {
    const result = Bun.spawnSync([
      "python3", script,
      "--input", sourceIfc,
      "--output", join(temporary, "candidate.ifc"),
      "--report", join(temporary, "report.json"),
      "--expected-ifc-sha256", "0".repeat(64),
      "--approval-token", "NOT-APPROVED",
    ], { cwd: root });
    expect(result.exitCode).not.toBe(0);
    expect(result.stderr.toString()).toContain("exact explicit approval token");
  } finally {
    rmSync(temporary, { recursive: true, force: true });
  }
});

test("A-103 approved apply validates tags, IDS, and zero geometry change on a separate output", () => {
  const temporary = mkdtempSync(join(tmpdir(), "a103-tag-apply-"));
  try {
    const fixture = buildTaglessFixture(temporary);
    const output = join(temporary, "candidate.ifc");
    const reportPath = join(temporary, "report.json");
    const result = Bun.spawnSync([
      "python3", script,
      "--input", fixture.baseline,
      "--output", output,
      "--candidate-report", fixture.candidate,
      "--report", reportPath,
      "--expected-ifc-sha256", fixture.sourceHash,
      "--approval-token", "APPROVE-A103-EW-NW",
    ], { cwd: root });
    expect(result.exitCode).toBe(0);
    const report = JSON.parse(readFileSync(reportPath, "utf8"));
    expect(report.mode).toBe("validated_output");
    expect(report.formal_ifc_write_performed).toBe(false);
    expect(report.gates.target_count).toBe(88);
    expect(report.gates.tags_correct_after_reload).toBe(true);
    expect(report.gates.entity_counts_unchanged).toBe(true);
    expect(report.gates.all_root_attributes_except_target_tags_unchanged).toBe(true);
    expect(report.gates.all_product_physical_graphs_unchanged).toBe(true);
    expect(report.gates.maximum_world_geometry_delta_mm).toBe(0);
    expect(report.gates.maximum_placement_matrix_delta).toBe(0);
    expect(report.gates.fills_voids_relationships_unchanged).toBe(true);
    expect(report.gates.source_ids_pass).toBe(499);
    expect(report.gates.postwrite_ids_pass).toBe(587);
    expect(report.gates.postwrite_ids_fail).toBe(6);
    expect(report.pass).toBe(true);
  } finally {
    rmSync(temporary, { recursive: true, force: true });
  }
}, 120_000);

test("A-103 apply refuses an in-place path without the replace-source guard", () => {
  const result = Bun.spawnSync([
    "python3", script,
    "--input", sourceIfc,
    "--output", sourceIfc,
    "--report", join(tmpdir(), "a103-should-not-write.json"),
    "--expected-ifc-sha256", "0".repeat(64),
    "--approval-token", "APPROVE-A103-EW-NW",
  ], { cwd: root });
  expect(result.exitCode).not.toBe(0);
  expect(result.stderr.toString()).toContain("--replace-source must be supplied");
});

test("A-103 apply source exposes no unapproved in-place write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).toContain('APPROVAL_TOKEN = "APPROVE-A103-EW-NW"');
  expect(source).toContain('parser.add_argument("--expected-ifc-sha256", required=True)');
  expect(source).toContain('parser.add_argument("--approval-token", required=True)');
  expect(source).toContain("replacing_source != args.replace_source");
  expect(source).toContain("os.replace(temporary, output_path)");
});
