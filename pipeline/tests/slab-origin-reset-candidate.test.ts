import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";

test("slab origin reset reuses the no-world-movement product implementation", () => {
  const result = spawnSync(
    "python3",
    [
      "-c",
      `import sys; sys.path.insert(0,"pipeline/scripts"); import slab_origin_reset_candidate as slab; from beam_origin_reset_candidate import reset_product_origin; print(slab.reset_product_origin is reset_product_origin)`,
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim()).toBe("True");
});
