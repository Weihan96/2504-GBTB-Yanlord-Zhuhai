import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/m401_existing_candidate.py");

test("M-401 inventory separates instances, types, context and missing inputs", () => {
  const temporary = mkdtempSync(join(tmpdir(), "m401-existing-"));
  const result = Bun.spawnSync(
    [
      "python3",
      script,
      "--input",
      resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
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
  expect(report.summary.actual_instances).toBe(28);
  expect(report.summary.instance_role_counts).toEqual({
    assigned_ac_equipment_instance: 5,
    high_level_ceiling_or_led_coordination_context: 15,
    named_embedded_ac_diffuser_proxy: 2,
    named_flue_check_valve_proxy: 1,
    untyped_hvac_service_mesh: 5,
  });
  expect(report.summary.type_definitions).toBe(3);
  expect(report.summary.missing_input_blocks).toBe(5);
  expect(report.summary.human_review_queue).toBe(18);
  expect(report.summary.ifc_air_terminal_instances).toBe(0);
  expect(report.summary.ifc_fan_instances).toBe(0);
  expect(report.summary.ifc_sensor_instances).toBe(0);
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
