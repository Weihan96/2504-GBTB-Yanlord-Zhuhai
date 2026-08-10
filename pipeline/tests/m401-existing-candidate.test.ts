import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/m401_existing_candidate.py");
const ifcPath = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
const currentIfcHash = createHash("sha256").update(readFileSync(ifcPath)).digest("hex");

test("M-401 inventory separates instances, types, context and missing inputs", () => {
  const temporary = mkdtempSync(join(tmpdir(), "m401-existing-"));
  const result = Bun.spawnSync(
    [
      "python3",
      script,
      "--input",
      ifcPath,
      "--decision-csv",
      join(temporary, "m401-existing-review.csv"),
      "--output-dir",
      join(temporary, "build"),
    ],
    { cwd: root },
  );
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(
    readFileSync(join(temporary, "build/m401-existing-report.json"), "utf8"),
  );
  expect(report.mode).toBe("read_only_m401_existing_candidate");
  expect(report.source_ifc_sha256).toBe(currentIfcHash);
  expect(report.source.ifc_sha256).toBe(currentIfcHash);
  expect(report.legacy_evidence.status).toBe("stale_not_refreshed");
  expect(report.legacy_evidence.current_formal_ifc).toBe(false);
  expect(report.summary.actual_instances).toBe(28);
  expect(report.summary.instance_role_counts).toEqual({
    assigned_ac_equipment_instance: 5,
    legacy_base_condensate_geometry: 1,
    legacy_base_refrigerant_gas_geometry: 2,
    legacy_base_refrigerant_liquid_geometry: 2,
    high_level_ceiling_or_led_coordination_context: 15,
    named_embedded_ac_diffuser_proxy: 2,
    named_flue_check_valve_proxy: 1,
  });
  expect(report.summary.type_definitions).toBe(3);
  expect(report.summary.missing_input_blocks).toBe(5);
  expect(report.summary.human_review_queue).toBe(18);
  expect(report.summary.ifc_air_terminal_instances).toBe(0);
  expect(report.summary.ifc_fan_instances).toBe(0);
  expect(report.summary.ifc_sensor_instances).toBe(1);
  expect(report.summary.ifc_alarm_instances).toBe(0);
  expect(report.summary.ifc_distribution_ports).toBe(0);
  expect(report.summary.ifc_systems).toBe(0);
  expect(report.summary.ifc_port_connections).toBe(0);
  expect(report.gates.inventory_pass).toBe(true);
  expect(report.gates.formal_ifc_write_allowed).toBe(false);
  expect(report.gates.m401_design_ready).toBe(false);
  expect(report.gates.rcp1_completion_pass).toBe(false);

  const instances = report.records.filter(
    (item: { record_kind: string }) => item.record_kind === "actual_instance",
  );
  expect(new Set(instances.map((item: { global_id: string }) => item.global_id)).size).toBe(28);
  const valve = instances.find(
    (item: { global_id: string }) => item.global_id === "1faflkXXH6M9cnYPE9Liir",
  );
  expect(valve.observable_role).toBe("named_flue_check_valve_proxy");
  expect(valve.review_status).toBe("HUMAN_REVIEW_REQUIRED");
  for (const globalId of ["16Ey9Flj9BK9VRun$ozzjH", "3Bv_Kl3jDC5RvaUMdyge1U"]) {
    const diffuser = instances.find(
      (item: { global_id: string }) => item.global_id === globalId,
    );
    expect(diffuser.observable_role).toBe("named_embedded_ac_diffuser_proxy");
  }
  const confirmedCornerDiffuser = instances.find(
    (item: { global_id: string }) => item.global_id === "16Ey9Flj9BK9VRun$ozzjH",
  );
  expect(confirmedCornerDiffuser.review_status).toBe(
    "IDENTITY_CONFIRMED_REMODEL_DESIGN_PENDING",
  );
  expect(confirmedCornerDiffuser.confidence).toBe(1);
  const confirmedLegacyDiffuser = instances.find(
    (item: { global_id: string }) => item.global_id === "3Bv_Kl3jDC5RvaUMdyge1U",
  );
  expect(confirmedLegacyDiffuser.review_status).toBe(
    "LEGACY_SCHEME_CONFIRMED_REDESIGN_REQUIRED",
  );
  expect(confirmedLegacyDiffuser.stop_condition).toContain("legacy design base");
  const confirmedFlowRoles = instances.filter(
    (item: { observable_role: string }) => item.observable_role.startsWith("legacy_base_"),
  );
  expect(confirmedFlowRoles).toHaveLength(5);
  for (const flow of confirmedFlowRoles) {
    expect(flow.review_status).toBe("LEGACY_BASE_CONFIRMED_REMODEL_DESIGN_PENDING");
    expect(flow.human_review_required).toBe(true);
    expect(flow.stop_condition).toContain("legacy design base only");
  }
  expect(report.records.filter(
    (item: { review_status: string }) => item.review_status === "BLOCK",
  )).toHaveLength(5);
}, 20_000);

test("M-401 candidate contains no write or inferred-design path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run(");
  for (const forbidden of ["airflow_cfm", "duct_diameter", "circuit_number", "automatic_connection"]){
    expect(source).not.toContain(forbidden);
  }
});
