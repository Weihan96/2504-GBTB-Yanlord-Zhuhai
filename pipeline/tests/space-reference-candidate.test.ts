import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/space_reference_candidate.py");

test("Space Reference generator stays read-only and fixes the agreed ordering invariant", async () => {
  const source = await Bun.file(script).text();
  expect(source).toContain('foyers = [record for record in records if record["long_name"] == "玄关"]');
  expect(source).toContain('record["candidate_reference"] = f"R{index:02d}"');
  expect(source).toContain('"automatic_ifc_write_allowed": False');
  expect(source).not.toContain("model.write(");
});

test("Space Reference candidate produces 22 unique codes with foyer first", async () => {
  const temp = resolve(root, "build/test-space-reference");
  const result = Bun.spawnSync([
    "python3", script,
    "--input", resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
    "--register", resolve(temp, "register.csv"),
    "--report", resolve(temp, "report.json"),
  ], { cwd: root });
  expect(result.exitCode).toBe(0);
  const report = await Bun.file(resolve(temp, "report.json")).json();
  expect(report.space_count).toBe(22);
  expect(report.unique_reference_count).toBe(22);
  expect(report.r01_long_name).toBe("玄关");
  expect([0, 22]).toContain(report.formal_ifc_reference_count);
  expect(report.existing_references_match_candidate).toBe(true);
  expect(report.qa.formal_ifc_unchanged).toBe(true);
});
