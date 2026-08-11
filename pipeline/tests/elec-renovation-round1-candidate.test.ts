import { afterAll, expect, test } from "bun:test";
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_renovation_round1_candidate.py");
const currentIfcHash = "6c2fd8da9e9ad7ddbc2b63415a27f1c979e8995b880d8fce210a2dda2ef2aab6";
const temp = mkdtempSync(join(tmpdir(), "elec-renovation-round1-test-"));
const output = join(temp, "candidate.json");
const svg = join(temp, "candidate.svg");

afterAll(() => rmSync(temp, { recursive: true, force: true }));

test("first-round renovation electrical demands stay read-only and complete", () => {
  const existing = JSON.parse(readFileSync(resolve(root, "build/elec/elec-existing-candidate.json"), "utf8"));
  const positioning = JSON.parse(readFileSync(resolve(root, "build/elec/elec-positioning-candidate.json"), "utf8"));
  const int1 = JSON.parse(readFileSync(resolve(root, "build/int1/int1-existing-report.json"), "utf8"));
  existing.source.sha256 = currentIfcHash;
  positioning.source_ifc_sha256 = currentIfcHash;
  int1.source.ifc_sha256 = currentIfcHash;
  const existingPath = join(temp, "existing.json");
  const positioningPath = join(temp, "positioning.json");
  const int1Path = join(temp, "int1.json");
  writeFileSync(existingPath, JSON.stringify(existing));
  writeFileSync(positioningPath, JSON.stringify(positioning));
  writeFileSync(int1Path, JSON.stringify(int1));
  const run = spawnSync([
    "python3",
    script,
    "--elec-existing",
    existingPath,
    "--elec-positioning",
    positioningPath,
    "--int1",
    int1Path,
    "--output",
    output,
    "--output-svg",
    svg,
  ], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.bedside_light_candidates).toBe(4);
  expect(report.summary.new_socket_candidates).toBe(2);
  expect(report.summary.new_socket_known_connected_load_w).toEqual({ "NS-01": 0, "NS-02": 0 });
  expect(report.summary.new_socket_planning_envelopes["NS-01"]).toMatchObject({
    listed_device_connected_minimum_w: 2600,
    listed_device_connected_maximum_w: 4900,
    simultaneous_design_minimum_w: 1800,
    simultaneous_design_maximum_w: 3700,
    planning_voltage_v: 220,
    single_socket_class_current_a: 16,
    minimum_independent_circuit_count_candidate: 2,
    minimum_connection_positions_candidate: 3,
    calculation_status: "listed_and_simultaneous_planning_envelopes_not_confirmed_product_load",
  });
  expect(report.summary.new_socket_planning_envelopes["NS-02"]).toMatchObject({
    listed_device_connected_minimum_w: 2350,
    listed_device_connected_maximum_w: 4000,
    simultaneous_design_minimum_w: 1000,
    simultaneous_design_maximum_w: 2200,
    minimum_independent_circuit_count_candidate: 1,
    minimum_connection_positions_candidate: 3,
  });
  expect(report.summary.cabinet_power_zones).toBe(7);
  expect(report.summary.kitchen_socket_rechecks).toBe(11);
  expect(report.summary.label_collision_count).toBe(0);
  expect(report.summary.developer_red_points_are_reference_only).toBe(true);
  expect(report.gates.four_bedside_lights_present).toBe(true);
  expect(report.gates.island_and_dining_bay_socket_present).toBe(true);
  expect(report.gates.socket_use_lists_compiled_without_fabricated_load).toBe(true);
  expect(report.gates.circuit_planning_candidates_match_confirmed_use).toBe(true);
  expect(report.gates.illuminated_cabinet_power_is_grouped_not_fabricated).toBe(true);
  expect(report.gates.all_current_kitchen_sockets_reopened_for_review).toBe(true);
  expect(report.gates.label_collision_free).toBe(true);
  expect(report.gates.reverse_requirement_audit_planned).toBe(true);
  expect(report.gates.whole_home_electrical_positioning_complete).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  const renderedSvg = readFileSync(svg, "utf8");
  expect(renderedSvg).toContain("elec-renovation-round1");
  expect(renderedSvg).toContain("Wall Plan-underlay.png");
  expect(renderedSvg).toContain("NS-01 全连接 2.6–4.9｜同时 1.8–3.7 kW");
  expect(renderedSvg).toContain("NS-02 同时使用：1.00–2.2 kW");
  expect(renderedSvg).toContain("三面各 1 个隐藏盖板位｜2 回路候选");
  expect(renderedSvg).toContain("线性轨道或自制翻盖｜1 回路候选");
  expect(renderedSvg).toContain("全连接/同时包络 ≠ 未购设备铭牌功率");
  const ns01 = report.new_socket_candidates.find((row: { candidate_id: string }) => row.candidate_id === "NS-01");
  const ns02 = report.new_socket_candidates.find((row: { candidate_id: string }) => row.candidate_id === "NS-02");
  expect(ns01.appliance_context.items.map((row: { appliance_name: string }) => row.appliance_name)).toEqual([
    "火锅电器", "搅拌机", "Sous-vide 棒",
  ]);
  expect(ns02.appliance_context.items.map((row: { appliance_name: string }) => row.appliance_name)).toEqual([
    "咖啡机", "手冲电热水壶", "磨豆机",
  ]);
  expect(ns01.appliance_context.socket_form_and_circuit_sizing_ready).toBe(false);
  expect(ns01.socket_form_candidate).toContain("three concealed covered");
  expect(ns02.circuit_strategy_candidate).toContain("not used together");
}, 30_000);

test("first-round renovation candidate has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
