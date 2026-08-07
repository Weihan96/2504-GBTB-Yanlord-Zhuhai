import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

const targets =
  "pipeline/decisions/c003-bay-window-assembly-origin-reset-targets.csv";
const cabinetTargets =
  "pipeline/decisions/c003-cabinet-assembly-origin-reset-targets.csv";

test("bay-window assembly targets are two exact integer construction anchors", () => {
  const result = spawnSync(
    "python3",
    [
      "-c",
      `import sys; sys.path.insert(0,"pipeline/scripts"); from pathlib import Path; from element_assembly_origin_reset_candidate import read_targets; rows=read_targets(Path(sys.argv[1])); print([(row["global_id"],row["target_mm"].tolist()) for row in rows])`,
      targets,
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim()).toBe(
    "[('14da0BQgfDmfHZBuWs7SnK', [-6500.0, -4600.0, 3000.0]), ('2$hL31sEL9xAyH$TmcJ1pN', [-6500.0, -4600.0, 0.0])]",
  );
});

test("assembly placement QA compares reopened source and candidate files", () => {
  const result = spawnSync(
    "python3",
    [
      "-c",
      `import sys; sys.path.insert(0,"pipeline/scripts"); import element_assembly_origin_reset_candidate as module; print(module.maximum_descendant_world_placement_residual.__name__)`,
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim()).toBe(
    "maximum_descendant_world_placement_residual",
  );
});

test("cabinet assembly reuses the confirmed SB03 installation anchor", () => {
  const result = spawnSync(
    "python3",
    [
      "-c",
      `import sys; sys.path.insert(0,"pipeline/scripts"); from pathlib import Path; from element_assembly_origin_reset_candidate import read_targets; row=read_targets(Path(sys.argv[1]))[0]; print(row["global_id"]); print(row["target_mm"].tolist()); print(row["anchor_kind"])`,
      cabinetTargets,
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim().split("\n")).toEqual([
    "3NIZVJBir6rQlVwhKs84Ca",
    "[2900.0, -1400.0, 0.0]",
    "confirmed_descendant_installation_anchor",
  ]);
});
