import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/beam_anchor_candidate.py");
const targetPath = resolve(
  root,
  "pipeline/decisions/c003-beam-anchor-targets.csv",
);

test("beam anchor targets are explicit integer-millimetre decisions", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("beam_anchor_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(json.dumps(module.read_targets(pathlib.Path(${JSON.stringify(targetPath)}))))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const rows = JSON.parse(result.stdout.toString());
  expect(rows).toHaveLength(6);
  for (const row of rows) {
    expect(row.target_mm.every(Number.isInteger)).toBe(true);
    expect(row.confidence).toBeGreaterThanOrEqual(0.95);
    expect(row.basis.length).toBeGreaterThan(0);
  }
  const safeRows = rows.filter((row: { review_required: boolean }) =>
    !row.review_required
  );
  expect(safeRows).toHaveLength(3);
});
