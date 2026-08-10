import { expect, test } from "bun:test";
import {
  copyFileSync,
  mkdtempSync,
  readFileSync,
  rmSync,
  symlinkSync,
  writeFileSync,
} from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const helper = resolve(root, "pipeline/scripts/svg_audit_underlay.py");
const ifc = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");

function validate(sourceSvg: string) {
  return spawnSync(
    [
      "python3",
      "-c",
      [
        "from pathlib import Path",
        "from svg_audit_underlay import validate_wall_plan_source",
        `source_svg = Path(${JSON.stringify(sourceSvg)})`,
        `validate_wall_plan_source(source_svg.read_text(encoding='utf-8'), source_svg, Path(${JSON.stringify(ifc)}))`,
      ].join(";"),
    ],
    {
      cwd: root,
      env: { ...process.env, PYTHONPATH: resolve(root, "pipeline/scripts") },
      stdout: "pipe",
      stderr: "pipe",
    },
  );
}

test("Wall Plan SVG and raster underlay require a current camera-state manifest", () => {
  const source = resolve(root, "drawings/Wall Plan.svg");
  expect(validate(source).exitCode).toBe(0);

  const temp = mkdtempSync(join(tmpdir(), "wall-plan-source-test-"));
  try {
    const testSvg = join(temp, "Wall Plan.svg");
    const testPng = join(temp, "Wall Plan-underlay.png");
    const testManifest = join(temp, "Wall Plan-source.json");
    copyFileSync(source, testSvg);
    symlinkSync(resolve(root, "drawings/Wall Plan-underlay.png"), testPng);
    const manifest = JSON.parse(
      readFileSync(resolve(root, "drawings/Wall Plan-source.json"), "utf8"),
    );
    manifest.viewport_perspectives = ["PERSP"];
    writeFileSync(testManifest, JSON.stringify(manifest));
    const stale = validate(testSvg);
    expect(stale.exitCode).not.toBe(0);
    expect(stale.stderr.toString()).toContain("invalid Wall Plan camera/texture state");
  } finally {
    rmSync(temp, { recursive: true, force: true });
  }
});
