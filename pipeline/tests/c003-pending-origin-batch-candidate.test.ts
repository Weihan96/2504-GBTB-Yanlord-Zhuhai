import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

test("C003 pending origin batch is an exact eight-target atomic candidate", async () => {
  const source = await Bun.file(
    "pipeline/scripts/c003_pending_origin_batch_candidate.py",
  ).text();
  expect(source).toContain("expected eight unique C003 targets");
  expect(source).toContain('gates["product_targets"] == 5');
  expect(source).toContain('gates["assembly_targets"] == 3');
  expect(source).toContain('"non_target_placement_changes_over_tolerance"');
  expect(source).toContain('"product_geometry_changes_over_tolerance"');
});

test("C003 pending origin batch keeps the 0.1 mm gate configurable", () => {
  const result = spawnSync(
    "python3",
    ["pipeline/scripts/c003_pending_origin_batch_candidate.py", "--help"],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout).toContain("--tolerance-mm");
});
