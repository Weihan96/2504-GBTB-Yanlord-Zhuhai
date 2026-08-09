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
  expect(report.gates.formal_ifc_write_allowed).toBe(false);
  expect(report.gates.hvac_design_ready).toBe(false);
}, 30_000);
