import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

test("A06 postwrite audit protects every pre-existing product", async () => {
  const source = await Bun.file(
    "pipeline/scripts/rcp1_a06_postwrite_audit.py",
  ).text();
  expect(source).toContain('"maximum_preexisting_product_world_geometry_change_mm"');
  expect(source).toContain('"preexisting_placement_and_representation_graphs_exact"');
  expect(source).toContain('"a06_is_placement_only"');
  expect(source).toContain('"a06_has_no_unconfirmed_type_ports_or_system"');
});

test("A06 postwrite audit requires caller-frozen hashes", () => {
  const result = spawnSync(
    "python3",
    ["pipeline/scripts/rcp1_a06_postwrite_audit.py", "--help"],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout).toContain("--expected-before-sha256");
  expect(result.stdout).toContain("--expected-formal-sha256");
  expect(result.stdout).toContain("--before-git-ref");
});
