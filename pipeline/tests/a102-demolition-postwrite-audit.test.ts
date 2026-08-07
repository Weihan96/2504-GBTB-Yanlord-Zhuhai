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
});
