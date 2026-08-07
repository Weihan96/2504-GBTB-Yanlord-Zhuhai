import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("stacked appliance reset uses exact bottom-support surface targets", () => {
  const source = `
import json,sys,pathlib
sys.path.insert(0,"pipeline/scripts")
from structural_origin_reset_candidate import read_targets
rows=read_targets(pathlib.Path("pipeline/decisions/c003-stacked-appliance-origin-reset-targets.csv"))
print(json.dumps([{**r,"target_mm":r["target_mm"].tolist()} for r in rows]))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const rows = JSON.parse(result.stdout.toString());
  expect(rows).toHaveLength(2);
  expect(rows).toEqual(
    expect.arrayContaining([
      expect.objectContaining({
        global_id: "0UOnmuAdP1MPy6p3olwiEU",
        target_mm: [2863, -2376, 970],
        anchor_kind: "existing_bottom_support_surface_point",
      }),
      expect.objectContaining({
        global_id: "3PQOXKxgj6IftqWXFXMQXG",
        target_mm: [2859, -2376, 370],
        anchor_kind: "existing_bottom_support_surface_point",
      }),
    ]),
  );
});
