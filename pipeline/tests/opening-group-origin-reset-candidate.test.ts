import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/opening_group_origin_reset_candidate.py");
const targets = resolve(root, "pipeline/decisions/c003-opening-group-origin-reset-targets.csv");

test("shared Opening target table defines two exact groups", () => {
  const code = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(script)}).parent))
spec = importlib.util.spec_from_file_location("opening_group", ${JSON.stringify(script)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
rows = module.read_targets(pathlib.Path(${JSON.stringify(targets)}))
print(json.dumps([{**row, "member_global_ids": sorted(row["member_global_ids"]), "target_mm": row["target_mm"].tolist()} for row in rows]))
`;
  const result = Bun.spawnSync(["python3", "-c", code], { cwd: root });
  expect(result.exitCode).toBe(0);
  const rows = JSON.parse(result.stdout.toString());
  expect(rows).toHaveLength(2);
  expect(rows.map((row: any) => row.member_global_ids.length)).toEqual([3, 2]);
  expect(rows.map((row: any) => row.target_mm)).toEqual([
    [-1815, 728, 0],
    [6000, -575, 0],
  ]);
});
