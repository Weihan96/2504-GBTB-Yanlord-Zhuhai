import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/structural_origin_reset_candidate.py");
const targets = resolve(root, "pipeline/decisions/c003-surface-edge-origin-reset-targets.csv");

test("surface edge targets are three unique integer boundary anchors", () => {
  const code = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(script)}).parent))
spec = importlib.util.spec_from_file_location("structural", ${JSON.stringify(script)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
rows = module.read_targets(pathlib.Path(${JSON.stringify(targets)}))
print(json.dumps([{"class": row["expected_class"], "global_id": row["global_id"], "target": row["target_mm"].tolist()} for row in rows]))
`;
  const result = Bun.spawnSync(["python3", "-c", code], { cwd: root });
  expect(result.exitCode).toBe(0);
  const rows = JSON.parse(result.stdout.toString());
  expect(rows).toHaveLength(3);
  expect(new Set(rows.map((row: any) => row.global_id)).size).toBe(3);
  expect(rows.map((row: any) => row.class).sort()).toEqual([
    "IfcCovering",
    "IfcCovering",
    "IfcSlab",
  ]);
  expect(rows.every((row: any) => row.target.every(Number.isInteger))).toBe(true);
});
