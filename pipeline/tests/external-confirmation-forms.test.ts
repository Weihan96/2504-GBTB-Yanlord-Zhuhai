import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("current external forms consume SSOT context and retain explicit ownership", () => {
  const result = Bun.spawnSync([
    "python3",
    "pipeline/scripts/validate_external_confirmation_forms.py",
  ], { cwd: root });
  expect(result.stderr.toString()).toBe("");
  expect(result.exitCode).toBe(0);
  expect(result.stdout.toString()).toContain("11 current forms");
});
