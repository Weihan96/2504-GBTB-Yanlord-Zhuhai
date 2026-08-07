import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("HD01 origin reset is one exact approved geometry-surface target", () => {
  const source = `
import json,sys,pathlib
sys.path.insert(0,"pipeline/scripts")
from structural_origin_reset_candidate import read_targets
rows=read_targets(pathlib.Path("pipeline/decisions/c003-electric-appliance-origin-reset-targets.csv"))
print(json.dumps([{**r,"target_mm":r["target_mm"].tolist()} for r in rows]))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const rows = JSON.parse(result.stdout.toString());
  expect(rows).toHaveLength(1);
  expect(rows[0]).toMatchObject({
    expected_class: "IfcElectricAppliance",
    global_id: "288GLY62v8kPPydA1lAK8W",
    target_mm: [3180, -4468, 1437],
    anchor_kind: "existing_bottom_surface_point",
    confidence: 1,
  });
});
