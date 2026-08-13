import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync, statSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/m401_coordination_sheet_candidate.py");
const ifc = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");

function sha256(path: string): string {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("M-401 coordination sheet preserves evidence boundaries and emits PNG proof", async () => {
  const temporary = mkdtempSync(join(tmpdir(), "m401-sheet-"));
  const svg = join(temporary, "M-401-candidate.svg");
  const png = join(temporary, "M-401-candidate.png");
  const reportPath = join(temporary, "M-401-candidate.json");
  const before = sha256(ifc);
  const process = Bun.spawn([
    "python3", script,
    "--root", root,
    "--expected-ifc-sha256", before,
    "--output-svg", svg,
    "--proof-png", png,
    "--report", reportPath,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  const [exitCode, stderr] = await Promise.all([
    process.exited,
    new Response(process.stderr).text(),
  ]);
  expect(exitCode).toBe(0);
  expect(stderr).not.toContain("Traceback");
  expect(sha256(ifc)).toBe(before);

  const report = await Bun.file(reportPath).json();
  expect(report.mode).toBe("read_only_m401_hvac_safety_coordination_sheet_candidate");
  expect(report.source_ifc_sha256).toBe(before);
  expect(report.summary).toMatchObject({
    actual_instance_count: 29,
    ac_type_count: 3,
    assigned_ac_instance_count: 5,
    placement_only_ac_instance_count: 1,
    missing_input_block_count: 5,
    confirmed_route_constraint_count: 6,
    confirmed_waypoint_count: 7,
    ifc_sensor_instances: 1,
    ifc_alarm_instances: 0,
    legacy_routes_declared_final_count: 0,
  });
  expect(report.summary.instance_role_counts).toEqual({
    assigned_ac_equipment_instance: 5,
    placement_only_ac_equipment_instance: 1,
    high_level_ceiling_or_led_coordination_context: 15,
    legacy_base_condensate_geometry: 1,
    legacy_base_refrigerant_gas_geometry: 2,
    legacy_base_refrigerant_liquid_geometry: 2,
    named_embedded_ac_diffuser_proxy: 2,
    named_flue_check_valve_proxy: 1,
  });
  expect(report.equipment_types.map((item: any) => [item.type_name, item.type_occurrence_count]).sort()).toEqual([
    ["AC1180", 1], ["AC700", 2], ["AC700F", 2],
  ]);
  expect(report.route_constraints.every((item: any) =>
    item.status === "user_confirmed" &&
    item.formal_ifc_write_allowed === false &&
    item.final_remodel_route === false
  )).toBe(true);
  expect(report.route_constraints.find((item: any) => item.route_id === "RCP1-SERVICE-A03").waypoint_order).toEqual(["A03", "H04", "H02"]);
  expect(report.release_blockers).toHaveLength(5);
  expect(report.release_blockers.every((item: any) => item.review_status === "BLOCK")).toBe(true);
  expect(report.safety_context).toMatchObject({
    ifc_sensor_instances: 1,
    ifc_alarm_instances: 0,
    sensor_name: "A106-FIRE-R04",
    sensor_predefined_type: "FIRESENSOR",
    review_status: "position_confirmed_type_pending",
    final_type_product_power_communication_pending: true,
    final_release_pass: false,
  });
  expect(report.gates).toEqual({
    caller_frozen_ifc_hash_match: true,
    source_hashes_current: true,
    inventory_counts_match: true,
    route_constraints_confirmed: true,
    legacy_routes_marked_nonfinal: true,
    one_ifc_sensor_preserved: true,
    kitchen_fire_position_only_closed: true,
    five_release_blockers_explicit: true,
    automatic_ifc_write_allowed: false,
    construction_release_ready: false,
  });
  for (const source of Object.values(report.source.inputs) as Array<{ path: string; sha256: string }>) {
    expect(source.sha256).toBe(sha256(source.path));
  }
  expect(statSync(svg).size).toBeGreaterThan(5_000);
  expect(statSync(png).size).toBeGreaterThan(20_000);
  expect(readFileSync(png).subarray(0, 8).toString("hex")).toBe("89504e470d0a1a0a");
  const svgText = readFileSync(svg, "utf8");
  expect(svgText).toContain("旧紫色管线不是装修后最终路线");
  expect(svgText).toContain("A03 → H04 → H02");
  expect(svgText).toContain("5 BLOCK");
  expect(svgText).toContain("construction_release_ready=false");
}, 30_000);

test("M-401 coordination sheet rejects a stale caller-frozen IFC hash before outputs", () => {
  const temporary = mkdtempSync(join(tmpdir(), "m401-sheet-stale-"));
  const svg = join(temporary, "candidate.svg");
  const png = join(temporary, "candidate.png");
  const report = join(temporary, "candidate.json");
  const result = Bun.spawnSync([
    "python3", script,
    "--root", root,
    "--expected-ifc-sha256", "0".repeat(64),
    "--output-svg", svg,
    "--proof-png", png,
    "--report", report,
  ], { cwd: root });
  expect(result.exitCode).not.toBe(0);
  expect(result.stderr.toString()).toContain("formal IFC SHA mismatch");
  expect(existsSync(svg)).toBe(false);
  expect(existsSync(png)).toBe(false);
  expect(existsSync(report)).toBe(false);
});

test("M-401 reports an omitted caller freeze as not evaluated", async () => {
  const temporary = mkdtempSync(join(tmpdir(), "m401-sheet-unfrozen-"));
  const reportPath = join(temporary, "candidate.json");
  const process = Bun.spawn([
    "python3", script,
    "--root", root,
    "--output-svg", join(temporary, "candidate.svg"),
    "--proof-png", join(temporary, "candidate.png"),
    "--report", reportPath,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  const [exitCode, stderr] = await Promise.all([
    process.exited,
    new Response(process.stderr).text(),
  ]);
  expect(exitCode, stderr).toBe(0);
  const report = await Bun.file(reportPath).json();
  expect(report.gates.caller_frozen_ifc_hash_match).toBeNull();
}, 30_000);

test("M-401 coordination generator has no IFC or Blender write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
  expect(source).not.toContain("save_ifc_file");
  expect(source).not.toContain("bpy");
});
