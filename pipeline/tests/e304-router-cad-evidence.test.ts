import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/e304_router_cad_evidence.py");
const output = resolve(root, "build/elec/e304-router-cad-evidence.test.json");
const svg = resolve(root, "build/elec/E304-router-entry-cabinet-evidence.test.svg");

test("official CAD pins the router to the entry weak-current cabinet without inventing Z", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.evidence_register.rows.map((row: { evidence_id: string }) => row.evidence_id)).toEqual([
    "E304-CAD-001",
    "E304-CAD-002",
    "E304-USER-001",
  ]);
  expect(report.text_evidence.map((row: { handle: string }) => row.handle)).toEqual(["2598C0", "224270", "224295"]);
  expect(report.leader_evidence.handle).toBe("224271");
  expect(report.source.coordinate_transform.viewport_handle).toBe("224238");
  expect(report.weak_current_box.ifc_plan_position_mm[0]).toBeCloseTo(4600.016493, 6);
  expect(report.weak_current_box.ifc_plan_position_mm[1]).toBeCloseTo(-735.368749, 6);
  expect(report.weak_current_box.weak_current_box_bottom_aff_mm).toBe(350);
  expect(report.router_decision.installation_z_mm).toBeNull();
  expect(report.gates.router_z_not_inferred_from_weak_box_datum).toBe(true);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(readFileSync(svg, "utf8")).toContain("H+350 是弱电箱底边，不是路由器安装高度");
}, 20_000);

test("CAD evidence extractor has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
