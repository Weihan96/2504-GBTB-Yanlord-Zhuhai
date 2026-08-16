import { afterAll, expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_renovation_round1_candidate.py");
const formalIfc = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
const currentIfcHash = createHash("sha256").update(readFileSync(formalIfc)).digest("hex");
const temp = mkdtempSync(join(tmpdir(), "elec-renovation-round1-test-"));
const output = join(temp, "candidate.json");
const svg = join(temp, "candidate.svg");
const ownerDecisions = resolve(root, "pipeline/decisions/owner-input-register.csv");
const loadScenarios = resolve(root, "pipeline/decisions/e303-load-scenarios.csv");
const circuitDesign = resolve(root, "pipeline/decisions/e303-circuit-design.csv");

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

function runCandidate(ownerDecisionPath = ownerDecisions, outputPath = output, svgPath = svg, circuitDesignPath = circuitDesign) {
  return spawnSync([
    "python3",
    script,
    "--elec-existing",
    existingPath,
    "--elec-positioning",
    positioningPath,
    "--int1",
    int1Path,
    "--owner-decisions",
    ownerDecisionPath,
    "--load-scenarios",
    loadScenarios,
    "--circuit-design",
    circuitDesignPath,
    "--output",
    outputPath,
    "--output-svg",
    svgPath,
  ], { cwd: root });
}

function sha256(path: string) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

afterAll(() => rmSync(temp, { recursive: true, force: true }));

test("first-round renovation electrical demands stay read-only and complete", () => {
  const run = runCandidate();
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.owner_inputs.sha256).toBe(sha256(ownerDecisions));
  expect(report.load_scenarios.sha256).toBe(sha256(loadScenarios));
  expect(report.owner_inputs.e303_circuit_decisions).toEqual({
    "E303-NS01-CIRCUIT": { status: "自定义确认" },
    "E303-NS02-CIRCUIT": { status: "自定义确认" },
  });
  expect(report.summary.bedside_light_candidates).toBe(4);
  expect(report.summary.new_socket_candidates).toBe(2);
  expect(report.summary.new_socket_known_connected_load_w).toEqual({ "NS-01": 0, "NS-02": 1200 });
  expect(report.summary.new_socket_planning_envelopes["NS-01"]).toMatchObject({
    design_load_w: 4400,
    design_current_a: 20,
    planning_voltage_v: 220,
    minimum_independent_circuit_count_candidate: 2,
    minimum_connection_positions_candidate: 3,
    breaker_curve_rating: "C16",
    conductor_cross_section_mm2: 2.5,
    rcbo_required: true,
    calculation_status: "research_conclusion_not_product_nameplate_load",
  });
  expect(report.summary.new_socket_planning_envelopes["NS-02"]).toMatchObject({
    design_load_w: 2600,
    design_current_a: 11.82,
    minimum_independent_circuit_count_candidate: 2,
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
  expect(report.gates.use_confirmed).toBe(true);
  expect(report.gates.planning_envelope_compiled).toBe(true);
  expect(report.gates.product_and_circuit_fixed).toBe(false);
  expect(report.gates.illuminated_cabinet_power_is_grouped_not_fabricated).toBe(true);
  expect(report.gates.all_current_kitchen_sockets_reopened_for_review).toBe(true);
  expect(report.gates.label_collision_free).toBe(true);
  expect(report.gates.reverse_requirement_audit_planned).toBe(true);
  expect(report.gates.whole_home_electrical_positioning_complete).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  const renderedSvg = readFileSync(svg, "utf8");
  expect(renderedSvg).toContain("elec-renovation-round1");
  expect(renderedSvg).toContain("Wall Plan-underlay.png");
  expect((renderedSvg.match(/data-coordination-kind="door"/g) ?? []).length).toBe(8);
  expect((renderedSvg.match(/data-coordination-kind="fixed_furniture"/g) ?? []).length).toBe(70);
  expect(renderedSvg).toContain("NS-01 两个 10A 端口同时：20A / 4.4kW");
  expect(renderedSvg).toContain("NS-02 咖啡候选上限：11.82A / 2.6kW");
  expect(renderedSvg).toContain("E-303 容量研究结论｜220V / C16 RCBO");
  expect(renderedSvg).toContain("三面各 1 个隐藏盖板位｜2 回路已闭合");
  expect(renderedSvg).toContain("线性轨道或自制翻盖｜2 回路（1咖啡专用+1辅助）");
  expect(renderedSvg).toContain("设计容量 ≠ 未购设备铭牌功率");
  const ns01 = report.new_socket_candidates.find((row: { candidate_id: string }) => row.candidate_id === "NS-01");
  const ns02 = report.new_socket_candidates.find((row: { candidate_id: string }) => row.candidate_id === "NS-02");
  expect(ns01.appliance_context.items.map((row: { appliance_name: string }) => row.appliance_name)).toEqual([
    "火锅电器", "搅拌机", "Sous-vide 棒",
  ]);
  expect(ns02.appliance_context.items.map((row: { appliance_name: string }) => row.appliance_name)).toEqual([
    "咖啡机", "手冲电热水壶", "磨豆机",
  ]);
  expect(ns01.appliance_context.explicit_load_scenarios.map((row: { scenario_label: string }) => row.scenario_label)).toEqual([
    "火锅+搅拌机", "火锅+sous-vide",
  ]);
  expect(ns02.appliance_context.explicit_load_scenarios.map((row: { scenario_label: string }) => row.scenario_label)).toEqual([
    "咖啡机+磨豆机", "手冲壶单独",
  ]);
  expect(report.e303_circuit_semantic_gates).toEqual({
    "NS-01": {
      use_confirmed: true,
      planning_envelope_compiled: true,
      product_and_circuit_fixed: false,
    },
    "NS-02": {
      use_confirmed: true,
      planning_envelope_compiled: true,
      product_and_circuit_fixed: false,
    },
  });
  expect(ns01.appliance_context.socket_form_and_circuit_sizing_ready).toBe(true);
  expect(ns01.socket_form_candidate).toContain("three concealed covered");
  expect(ns01.circuit_strategy_candidate).toContain("two independent C16 RCBO");
  expect(ns02.circuit_strategy_candidate).toContain("dedicated C16 RCBO coffee circuit");
}, 30_000);

test("E-303 rejects a circuit row that is not an adopted research conclusion", () => {
  const source = readFileSync(circuitDesign, "utf8");
  const invalid = source.replace(",research_conclusion,ELEC-RESEARCH-20260815-001,", ",candidate,ELEC-RESEARCH-20260815-001,");
  const invalidPath = join(temp, "circuit-design-invalid.csv");
  writeFileSync(invalidPath, invalid);
  const run = runCandidate(ownerDecisions, join(temp, "invalid.json"), join(temp, "invalid.svg"), invalidPath);
  expect(run.exitCode).not.toBe(0);
  expect(run.stderr.toString()).toContain("circuit design is not adopted");
}, 30_000);

test("first-round renovation candidate has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
