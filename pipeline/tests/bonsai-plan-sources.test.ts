import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const statusScript = resolve(root, "pipeline/scripts/bonsai_plan_sources.py");
const refreshScript = resolve(root, "pipeline/scripts/refresh_bonsai_plan_sources.py");
const drawingScript = resolve(root, "pipeline/scripts/drawing_pipeline.py");

test("Bonsai plan source register has five unique native Drawings", () => {
  const rows = readFileSync(
    resolve(root, "pipeline/decisions/bonsai-plan-source-register.csv"),
    "utf8",
  ).trim().split("\n");
  expect(rows).toHaveLength(6);
  const names = rows.slice(1).map((row) => row.split(",")[0]);
  expect(names).toEqual([
    "Wall Plan",
    "Furniture Plan",
    "Sanitary Plan",
    "FFL PLAN",
    "Wall Finish Plan",
  ]);
  expect(new Set(names).size).toBe(names.length);
});

test("underlay status is a read-only view over source manifests", () => {
  const result = spawnSync(["python3", statusScript, "status"], { cwd: root });
  expect(result.exitCode, result.stderr.toString()).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.mode).toBe("read_only_bonsai_plan_source_status");
  expect(report.drawing_count).toBe(5);
  expect(report.current_count + report.stale_count).toBe(5);
  expect(report.drawings.map((item: { drawing_name: string }) => item.drawing_name)).toEqual([
    "Wall Plan",
    "Furniture Plan",
    "Sanitary Plan",
    "FFL PLAN",
    "Wall Finish Plan",
  ]);
});

test("native refresh loads Bonsai and never writes the formal IFC", () => {
  const source = readFileSync(refreshScript, "utf8");
  expect(source).toContain("bpy.ops.bim.load_project");
  expect(source).toContain("bpy.ops.bim.activate_drawing");
  expect(source).toContain("bpy.ops.bim.create_drawing");
  expect(source).toContain("should_use_underlay_cache = False");
  expect(source).toContain("formal IFC changed");
  expect(source).toContain("bpy.app.timers.register");
  expect(source).toContain("bpy.context.temp_override");
  expect(source).toContain("apply_official_elevation_index");
  expect(source).toContain("os._exit(0)");
  expect(source).not.toMatch(/tool\.Ifc\.get\(\)\.write/);
  expect(source).not.toMatch(/model\.write\s*\(/);
  const wrapper = readFileSync(statusScript, "utf8");
  expect(wrapper).not.toContain('"--background"');
  expect(wrapper).toContain("capture_output=True");
});

test("single-sheet dry run resolves its native underlays and grouped build", () => {
  const result = spawnSync(
    ["python3", drawingScript, "--dry-run", "--sheet", "I-501"],
    { cwd: root },
  );
  expect(result.exitCode, result.stderr.toString()).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.requested_sheets).toEqual(["I-501"]);
  expect(report.selected_sheets).toEqual(["I-501", "I-502", "I-503", "I-504"]);
  expect(report.underlays).toEqual(["Furniture Plan", "Sanitary Plan"]);
  expect(report.commands).toEqual([["bun", "run", "pipeline:int1-candidate"]]);
  expect(report.dry_run).toBe(true);
});
