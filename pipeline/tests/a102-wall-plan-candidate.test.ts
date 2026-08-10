import { expect, test } from "bun:test";
import { spawnSync } from "node:child_process";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/a102_wall_plan_candidate.py");

function runPython(body: string) {
  return Bun.spawnSync(["python3", "-c", `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("a102_wall_plan_candidate", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
${body}
`], { cwd: root });
}

test("A102 bbox overlays use the Wall Plan 1:50 mapping", () => {
  const result = runPython(`print(json.dumps(module.rect_for_bbox({
    "min_mm": [-6500, 0, 0],
    "max_mm": [-5800, 100, 3000],
    "dimensions_mm": [700, 100, 3000],
})))`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual([70, 198, 14, 2]);
});

test("A102 status boundary distinguishes both demolition source states", () => {
  const result = runPython(`print(json.dumps(module.EXPECTED_STATUS_COUNTS))`);
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    EXISTING: 84,
    NEW: 4,
    DEMOLISH_ALREADY_REMOVED: 11,
    DEMOLISH_PLANNED: 2,
  });
});

test("A102 side panel preserves the approximate-not-survey boundary", async () => {
  const source = await Bun.file(modulePath).text();
  expect(source).toContain("不是现场测量或施工放线依据");
  expect(source).toContain("CONFIRMED_APPROXIMATE");
  expect(source).toContain('"protected_products_max_world_geometry_change_mm"');
});

test("A102 plan rejects a mismatched caller-frozen IFC hash", () => {
  const result = spawnSync(
    "python3",
    [
      "pipeline/scripts/a102_wall_plan_candidate.py",
      "--input", "2504 GBTB Yanlord Zhuhai.ifc",
      "--source-svg", "drawings/Wall Plan.svg",
      "--review-register", "pipeline/decisions/a102-demolition-review.csv",
      "--postwrite-report", "build/a102/a102-demolition-postwrite.json",
      "--output-svg", "build/a102/should-not-write.svg",
      "--report", "build/a102/should-not-write.json",
      "--expected-ifc-sha256", "0".repeat(64),
    ],
    { cwd: root, encoding: "utf8" },
  );
  expect(result.status).not.toBe(0);
  expect(result.stderr).toContain("differs from caller-frozen hash");
});
