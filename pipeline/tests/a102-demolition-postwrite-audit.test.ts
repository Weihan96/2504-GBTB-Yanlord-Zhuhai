import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

test("A102 postwrite audit protects the approved scope", async () => {
  const source = await Bun.file(
    "pipeline/scripts/a102_demolition_postwrite_audit.py",
  ).text();
  expect(source).toContain('"baseline_product_count"');
  expect(source).toContain('"demolition_records_pass"');
  expect(source).toContain('"protected_products_max_world_geometry_change_mm"');
  expect(source).toContain('"POSITION_DIRECTION_LENGTH_APPROXIMATE_NOT_SURVEY_GRADE"');
});

test("A102 postwrite audit keeps the geometry tolerance configurable", () => {
  const result = spawnSync(
    "python3",
    ["pipeline/scripts/a102_demolition_postwrite_audit.py", "--help"],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout).toContain("--tolerance-mm");
  expect(result.stdout).toContain("--expected-ifc-sha256");
});

test("A102 postwrite audit rejects a mismatched caller-frozen IFC hash", () => {
  const result = spawnSync(
    "python3",
    [
      "pipeline/scripts/a102_demolition_postwrite_audit.py",
      "--formal", "2504 GBTB Yanlord Zhuhai.ifc",
      "--candidate", "build/candidates/2504-GBTB-a102-demolition.ifc",
      "--register", "pipeline/decisions/a102-demolition-review.csv",
      "--report", "build/a102/should-not-write.json",
      "--expected-ifc-sha256", "0".repeat(64),
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).not.toBe(0);
  expect(result.stderr).toContain("differs from caller-frozen hash");
});
