import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("owner input sync is dry-run by default and never writes IFC", () => {
  const temp = mkdtempSync(join(tmpdir(), "owner-inputs-"));
  const report = join(temp, "report.json");
  const openItems = join(temp, "open.md");
  const decisions = join(root, "pipeline/decisions/owner-input-register.csv");
  const before = readFileSync(decisions, "utf8");
  const result = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py",
    "--input", decisions,
    "--report", report,
    "--open-items", openItems,
  ], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(readFileSync(decisions, "utf8")).toBe(before);
  const parsed = JSON.parse(readFileSync(report, "utf8"));
  expect(parsed.mode).toBe("dry-run");
  expect(parsed.formal_ifc_write).toBe(false);
  expect(parsed.summary.release_blocking_open_count).toBeGreaterThan(0);
});

test("owner input sync rejects a custom confirmation without a value", () => {
  const temp = mkdtempSync(join(tmpdir(), "owner-inputs-invalid-"));
  const source = readFileSync(join(root, "pipeline/decisions/owner-input-register.csv"), "utf8");
  const invalid = source.replace(",待填写,,门口控制应位于真实墙面", ",自定义确认,,门口控制应位于真实墙面");
  const input = join(temp, "invalid.csv");
  writeFileSync(input, invalid);
  const result = Bun.spawnSync([
    "python3", "pipeline/scripts/sync_owner_inputs.py",
    "--input", input,
    "--report", join(temp, "report.json"),
    "--open-items", join(temp, "open.md"),
  ], { cwd: root });
  expect(result.exitCode).toBe(2);
  expect(result.stderr.toString()).toContain("custom confirmation requires user_value");
});
