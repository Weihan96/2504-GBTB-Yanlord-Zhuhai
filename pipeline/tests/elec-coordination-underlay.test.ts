import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_coordination_underlay.py");
const svg = resolve(root, "drawings/Electrical Coordination Plan.svg");
const manifest = resolve(root, "drawings/Electrical Coordination Plan-source.json");

test("electrical coordination underlay carries current walls, doors, and fixed joinery", () => {
  const run = spawnSync(["python3", script], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const report = JSON.parse(readFileSync(manifest, "utf8"));
  expect(report.counts).toMatchObject({ walls: 101, doors: 8, fixed_furniture: 70 });
  expect(report.wall_plan_svg_path).toBe("drawings/Wall Plan.svg");
  expect(report.int1_report_path).toBe("build/int1/int1-existing-report.json");
  expect(report.coordination_svg_path).toBe("drawings/Electrical Coordination Plan.svg");
  expect(JSON.stringify(report)).not.toContain("/Users/");
  expect(report.gates).toEqual({
    wall_plan_current: true,
    doors_from_current_ifc: true,
    fixed_furniture_from_current_ifc_and_int1_role_register: true,
    stale_furniture_plan_not_used: true,
    automatic_ifc_write_allowed: false,
  });
  const rendered = readFileSync(svg, "utf8");
  expect(rendered.match(/class="official-elevation-anchor"/g)).toHaveLength(12);
  expect(rendered.match(/class="official-elevation-direction"/g)).toHaveLength(36);
  expect((rendered.match(/data-coordination-kind="door"/g) ?? []).length).toBe(8);
  expect((rendered.match(/data-coordination-kind="fixed_furniture"/g) ?? []).length).toBe(70);
  expect(rendered).toContain('data-electrical-coordination="current-ifc"');
  expect(rendered).toContain("Wall Plan-underlay.png");
}, 40_000);
