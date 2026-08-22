import { expect, test } from "bun:test";
import { copyFileSync, mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/plan_elevation_index.py");
const source = resolve(root, "drawings/Wall Plan.svg");
const register = resolve(root, "pipeline/decisions/int1-elevation-view-register.csv");

test("plan elevation index exposes twelve official anchors and 36 directions", () => {
  const temporary = mkdtempSync(join(tmpdir(), "plan-elevation-index-"));
  const candidate = join(temporary, "Wall Plan.svg");
  copyFileSync(source, candidate);
  const result = Bun.spawnSync([
    "python3", script, candidate, "--register", register,
  ], { cwd: root });
  expect(result.exitCode, result.stderr.toString()).toBe(0);
  const svg = readFileSync(candidate, "utf8");
  expect(svg.match(/class="official-elevation-anchor"/g)).toHaveLength(12);
  expect(svg.match(/class="official-elevation-direction"/g)).toHaveLength(36);
  expect(svg).not.toContain('xlink:href="#elevation-arrow"');
  expect(svg).not.toContain('class="ELEVATION"');
  expect(svg).not.toContain("EL-P01");
  expect(svg).not.toContain("EL-P02");
  const ids = [...svg.matchAll(/data-view-id="(\d{2})"/g)].map((match) => match[1]);
  expect(new Set(ids)).toEqual(new Set(Array.from({ length: 36 }, (_, index) => String(index + 1).padStart(2, "0"))));

  const check = Bun.spawnSync([
    "python3", script, candidate, "--register", register, "--check",
  ], { cwd: root });
  expect(check.exitCode, check.stderr.toString()).toBe(0);
});

test("native plan refresh always applies the managed index", () => {
  const source = readFileSync(
    resolve(root, "pipeline/scripts/refresh_bonsai_plan_sources.py"),
    "utf8",
  );
  expect(source).toContain("apply_official_elevation_index");
  expect(source).toContain("ELEVATION_INDEX_REGISTER");
});
