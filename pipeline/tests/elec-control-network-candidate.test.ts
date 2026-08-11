import { expect, test } from "bun:test";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/elec_control_network_candidate.py");
const output = resolve(root, "build/elec/elec-control-network-candidate.test.json");
const svg = resolve(root, "build/elec/E302-E304-control-network-candidate.test.svg");
const productGateNames = [
  "exact_sku_certificate_match",
  "control_role_closed",
  "physical_wired_two_way_verified",
  "neutral_and_wiring_diagram_verified",
  "rated_load_schedule_verified",
  "box_and_joinery_interface_verified",
  "matter_infrastructure_verified",
  "ecosystem_behavior_acceptance",
];

test("confirmed control, network, and safety roles compile as read-only coordination zones", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.summary.control_coordination_zones).toBe(2);
  expect(report.summary.control_wall_side_options).toBe(4);
  expect(report.summary.control_panel_candidates).toBe(3);
  expect(report.summary.paired_two_way_control_groups).toBe(3);
  expect(report.summary.entrance_master_lighting_switches).toBe(1);
  expect(report.summary.bedroom_AP_candidates).toBe(2);
  expect(report.summary.living_study_router_zones).toBe(1);
  expect(report.summary.smoke_alarm_candidates).toBe(3);
  expect(report.summary.kitchen_fire_sensor_positions).toBe(1);
  expect(report.summary.kitchen_gas_alarm_room_zones).toBe(1);
  expect(report.summary.network_requirement_links).toBe(6);
  expect(report.summary.network_downstream_endpoints).toBe(5);
  expect(report.summary.network_AP_power_method_pending_endpoints).toBe(2);
  expect(report.summary.switch_product_review_panels).toBe(3);
  expect(report.summary.switch_product_release_gates).toBe(24);
  expect(report.summary.switch_product_release_blockers).toBe(24);
  expect(Object.values(report.gates).every((value) => value === true || value === false)).toBe(true);
  expect(report.gates.two_doorway_zones_present).toBe(true);
  expect(report.gates.four_wall_side_options_present).toBe(true);
  expect(report.gates.master_a_and_b_distinct_roles_pending_a104).toBe(true);
  expect(report.gates.entry_wall_side_not_auto_closed).toBe(true);
  expect(report.gates.entry_candidate_not_mislabeled_selected).toBe(true);
  expect(report.gates.controlled_fixture_group_mapping_complete).toBe(false);
  expect(report.gates.three_two_way_groups_present).toBe(true);
  expect(report.gates.master_a_internal_lighting_separate).toBe(true);
  expect(report.gates.entrance_master_switch_present).toBe(true);
  expect(report.gates.two_bedroom_AP_candidates_present).toBe(true);
  expect(report.gates.one_shared_router_no_AP_zone_present).toBe(true);
  expect(report.gates.router_entry_cabinet_plan_position_verified).toBe(true);
  expect(report.gates.six_network_requirement_links_present).toBe(true);
  expect(report.gates.five_downstream_endpoints_present).toBe(true);
  expect(report.gates.two_AP_power_methods_pending).toBe(true);
  expect(report.gates.physical_ports_remain_unassigned).toBe(true);
  expect(report.gates.main_bedroom_AP_clearance_reserve_pass).toBe(true);
  expect(report.gates.guest_bedroom_AP_clearance_reserve_pass).toBe(true);
  const mainBedroomAp = report.network_coordination_zones.find(
    (candidate: any) => candidate.candidate_id === "A106-AP-R09",
  );
  expect(mainBedroomAp.position_mm).toEqual([-5840, -2844, 2720]);
  expect(mainBedroomAp.clearance_margin_mm).toBeCloseTo(578.5, 3);
  expect(report.gates.gateway_identity_complete).toBe(false);
  expect(report.gates.cabinet_dimensions_complete).toBe(false);
  expect(report.gates.thermal_test_complete).toBe(false);
  expect(report.gates.cable_continuity_complete).toBe(false);
  expect(report.gates.ap_power_method_complete).toBe(false);
  expect(report.gates.three_smoke_candidates_present).toBe(true);
  expect(report.gates.smoke_positioning_constraints_present).toBe(true);
  expect(report.gates.one_formal_kitchen_fire_position_present).toBe(true);
  expect(report.gates.one_kitchen_gas_alarm_zone_present).toBe(true);
  expect(report.gates.gas_alarm_model_authority_pending).toBe(true);
  expect(report.gates.a106_exact_positions_imported).toBe(true);
  expect(report.gates.unconfirmed_positions_remain_coordination_zones).toBe(true);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(report.gates.switch_product_release_ready).toBe(false);
  expect(report.gates.construction_release_ready).toBe(false);
  const renderedSvg = readFileSync(svg, "utf8");
  expect(renderedSvg).toContain("elec-control-network");
  expect(renderedSvg).toContain("A106-FIRE-R04");
  expect(renderedSvg).toContain("CTRL-MASTER-A USER/A104");
  expect(renderedSvg).toContain("CTRL-MASTER-B USER/A104");
  expect(renderedSvg).toContain("CTRL-ENTRY-A USER/PENDING");
  expect(renderedSvg).toContain("CTRL-ENTRY-B OPTION");
  expect(renderedSvg).toContain("主卧 AP 灯具净距余量 578.5mm｜机械候选");
  expect(renderedSvg).toContain("玄关高柜路由器平面柜位");
  expect(renderedSvg).toContain("3 面板 × 8 门＝24 项，全部 BLOCK");
  expect(renderedSvg).toContain("Matter 网络类型按准确 SKU；Matter ≠ KNX");
  expect(renderedSvg).toContain("Wall Plan-underlay.png");
}, 30_000);

test("Entry A, Master A, and Master B each fail every named product release gate with a reason", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.switch_product_review.gate_names).toEqual(productGateNames);
  expect(report.switch_product_review.release_ready).toBe(false);
  expect(report.switch_product_review.protocol_position).toEqual({
    candidate_protocols: {
      "CTRL-ENTRY-A": "Matter over Thread candidate",
      "CTRL-MASTER-A": "Matter over Wi-Fi candidate",
      "CTRL-MASTER-B": "unassigned",
    },
    not_equivalent_to: "KNX",
    infrastructure_verified: false,
  });

  const panels = new Map(report.switch_product_review.panels.map((panel: any) => [panel.panel_name, panel]));
  expect([...panels.keys()]).toEqual(["Entry A", "Master A", "Master B"]);
  for (const panelName of ["Entry A", "Master A", "Master B"]) {
    const panel: any = panels.get(panelName);
    expect(Object.keys(panel.gates)).toEqual(productGateNames);
    expect(Object.values(panel.gates)).toEqual(productGateNames.map(() => false));
    expect(panel.release_ready).toBe(false);
    expect(panel.release_blockers).toHaveLength(8);
    for (const gate of productGateNames) {
      const blocker = panel.release_blockers.find((item: any) => item.gate === gate);
      expect(blocker.panel_name).toBe(panelName);
      expect(blocker.failure_reason.length).toBeGreaterThan(0);
      expect(blocker.required_evidence.length).toBeGreaterThan(0);
    }
  }
  expect(report.release_blockers).toHaveLength(24);
}, 30_000);

test("product evidence limits and panel-specific closeout reasons remain explicit", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  const byName = new Map(report.switch_product_review.panels.map((panel: any) => [panel.panel_name, panel]));
  const entry: any = byName.get("Entry A");
  const masterA: any = byName.get("Master A");
  const masterB: any = byName.get("Master B");

  expect(entry.product_evidence_scope).toBe("csa_product_family_only_exact_taobao_variant_unmatched");
  expect(entry.release_blockers.find((item: any) => item.gate === "exact_sku_certificate_match").failure_reason)
    .toContain("淘宝 EGG 具体变体");
  expect(masterA.product_evidence_scope).toBe("manufacturer_parameters_present_official_csa_association_missing");
  expect(masterA.release_blockers.find((item: any) => item.gate === "exact_sku_certificate_match").failure_reason)
    .toContain("2.5D Neo 厂家参数存在");
  const joineryReason = masterA.release_blockers
    .find((item: any) => item.gate === "box_and_joinery_interface_verified").failure_reason;
  for (const condition of ["阻燃背盒", "固定基层", "散热", "可检修"]) {
    expect(joineryReason).toContain(condition);
  }
  expect(masterB.product_evidence_scope).toBe("no_exact_product_assigned_master_a_2_5d_neo_evidence_not_inherited");
  expect(masterB.release_blockers.find((item: any) => item.gate === "exact_sku_certificate_match").failure_reason)
    .toContain("不得自动继承");
  for (const panel of [entry, masterA, masterB]) {
    expect(panel.release_blockers.find((item: any) => item.gate === "matter_infrastructure_verified").failure_reason)
      .toContain("Matter 不等同于 KNX");
  }
}, 30_000);

test("control/network candidate has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
