import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

test("C003 postwrite audit checks exact placement and PVC110 geometry", async () => {
  const source = await Bun.file(
    "pipeline/scripts/c003_pending_origin_batch_postwrite_audit.py",
  ).text();
  expect(source).toContain('"target_placement_mismatches"');
  expect(source).toContain('"all_product_placement_mismatches"');
  expect(source).toContain(
    '"pvc110_max_world_geometry_change_from_prewrite_baseline_mm"',
  );
  expect(source).toContain('record["formal_body_item_count"] == 3');
});

test("C003 postwrite audit keeps the mechanical tolerance configurable", () => {
  const result = spawnSync(
    "python3",
    ["pipeline/scripts/c003_pending_origin_batch_postwrite_audit.py", "--help"],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout).toContain("--tolerance-mm");
});
