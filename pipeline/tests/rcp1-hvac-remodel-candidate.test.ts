import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/rcp1_hvac_remodel_candidate.py");
const output = resolve(root, "build/rcp1/hvac-remodel-candidate.test.json");

test("RCP1 remodel HVAC candidate keeps proximity separate from connection", async () => {
  const source = await Bun.file(script).text();
  expect(source).toContain('"formal_ifc_write_allowed": False');
  expect(source).toContain('"proximity_is_connection": False');
  expect(source).not.toContain("model.write(");

  const process = Bun.spawn(
    [
      "python3",
      script,
      "--input",
      resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
      "--rcp1-report",
      resolve(root, "build/rcp1/rcp1-existing-candidate.json"),
      "--coordination-report",
      resolve(root, "build/rcp1/coordination-report.json"),
      "--m401-report",
      resolve(root, "build/rcp1/m401-existing-report.json"),
      "--legacy-report",
      resolve(root, "build/rcp1/legacy-base-audit.json"),
      "--review",
      resolve(root, "pipeline/decisions/rcp1-hvac-remodel-review.csv"),
      "--output",
      output,
      "--tolerance-mm",
      "0.1",
    ],
    { cwd: root, stdout: "pipe", stderr: "pipe" },
  );
  const [exitCode, stderr] = await Promise.all([
    process.exited,
    new Response(process.stderr).text(),
  ]);
  expect(stderr).toBe("");
  expect(exitCode).toBe(0);

  const report = await Bun.file(output).json();
  expect(report.summary.formal_ac_candidates).toBe(5);
  expect(report.summary.legacy_east_ac_candidates_without_global_id).toBe(1);
  expect(report.summary.developer_openings).toBe(7);
  expect(report.summary.legacy_pipe_products).toBe(5);
  expect(report.summary.legacy_pipe_independent_components).toBe(12);
  expect(report.summary.confirmed_outlet_identity_candidates).toBe(2);
  expect(report.formal_equipment_pairing).toHaveLength(5);
  expect(
    report.formal_equipment_pairing.every(
      (item: any) => item.service_space_probe.direct_space_candidate,
    ),
  ).toBe(true);
  expect(report.developer_opening_pairing).toHaveLength(7);
  expect(report.airside_pairing).toHaveLength(2);
  expect(report.human_review_bundle).toHaveLength(3);
  expect(report.gates.candidate_ready_for_blender_review).toBe(true);
  expect(report.gates.service_space_probe_ready_for_review).toBe(true);
  expect(
    report.service_space_summary.living_room_R20_has_no_direct_equipment_or_confirmed_outlet_candidate,
  ).toBe(true);
  expect(report.gates.living_room_service_resolved).toBe(false);
  expect(report.gates.fixed_equipment_positions_confirmed).toBe(true);
  expect(report.gates.fixed_equipment_airside_paths_mechanically_diagnosed).toBe(true);
  expect(report.fixed_equipment_airside_diagnostics).toHaveLength(2);
  const [a05Diagnostic, a06Diagnostic] = report.fixed_equipment_airside_diagnostics;
  expect(a05Diagnostic.diagnostic_id).toBe("RCP1-AIR-FIXED-A05-TO-R20");
  expect(a05Diagnostic.equipment_position_status).toBe("confirmed_fixed");
  expect(a05Diagnostic.formal_ifc_identity_required).toBe(false);
  expect(a05Diagnostic.route.turn_from_current_airside_degrees).toBeLessThan(1e-6);
  expect(
    a05Diagnostic.route.demolition_wall_crossings.some(
      (wall: any) => wall.global_id === "12lp8aIu9LTeewHdYAmHs2",
    ),
  ).toBe(true);
  expect(a05Diagnostic.route.permanent_wall_crossings).toHaveLength(0);
  expect(a06Diagnostic.diagnostic_id).toBe("RCP1-AIR-FIXED-A06-TO-R20");
  expect(a06Diagnostic.equipment_position_status).toBe("confirmed_fixed");
  expect(a06Diagnostic.formal_ifc_identity_required).toBe(true);
  expect(a06Diagnostic.route.turn_from_current_airside_degrees).toBeGreaterThan(90);
  expect(report.gates.formal_ifc_write_allowed).toBe(false);
  expect(report.gates.hvac_design_ready).toBe(false);
}, 30_000);
