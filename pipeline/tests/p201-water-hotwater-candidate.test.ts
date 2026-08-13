import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { mkdtempSync, readFileSync, statSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/p201_water_hotwater_candidate.py");
const ifc = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
const endpoints = resolve(root, "build/plum/p201-demand-endpoints.json");

function sha256(path: string): string {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("P-201 renders all demand endpoints without inventing water topology", async () => {
  const temporary = mkdtempSync(join(tmpdir(), "p201-water-demand-"));
  const svg = join(temporary, "P-201-water-hotwater-demand-candidate.svg");
  const png = join(temporary, "P-201-water-hotwater-demand-candidate.png");
  const reportPath = join(temporary, "p201-water-hotwater-demand-candidate.json");
  const frozenHash = sha256(ifc);
  const process = Bun.spawn([
    "python3", script,
    "--root", root,
    "--expected-ifc-sha256", frozenHash,
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

  const report = await Bun.file(reportPath).json();
  expect(report.source_ifc_sha256).toBe(frozenHash);
  expect(report.source.ifc.sha256).toBe(frozenHash);
  expect(report.summary).toMatchObject({
    registered_object_count: 27,
    service_demand_candidate_count: 24,
    non_service_component_count: 3,
    ifc_distribution_port_count: 0,
    ifc_system_count: 0,
    ifc_distribution_system_count: 0,
    ifc_rel_connects_ports_count: 0,
    invented_pipe_route_count: 0,
    invented_pipe_size_count: 0,
    invented_pressure_value_count: 0,
    invented_equipment_interface_count: 0,
  });
  expect(report.summary.connection_requirement_counts).toEqual({ unknown: 27 });
  expect(report.summary.service_media_status_counts).toEqual({
    confirmed_exact_model: 3,
    unknown: 21,
    not_applicable: 3,
  });
  expect(report.gates).toMatchObject({
    source_hashes_current: true,
    caller_frozen_hash_checked: true,
    all_registered_objects_drawn: true,
    service_and_non_service_split_explicit: true,
    unknown_connection_coordinates_preserved: true,
    confirmed_service_media_has_exact_model_evidence: true,
    formal_distribution_topology_present: false,
    ifc_unchanged_during_generation: true,
    automatic_ifc_write_allowed: false,
    construction_release_ready: false,
  });
  expect(report.records).toHaveLength(27);
  expect(report.records.filter((row: any) => row.service_demand_candidate)).toHaveLength(24);
  expect(report.records.filter((row: any) => !row.service_demand_candidate)).toHaveLength(3);
  expect(report.records.every((row: any) => row.connection_requirement === "unknown")).toBe(true);

  const sourceRecords = JSON.parse(readFileSync(endpoints, "utf8")).demand_endpoints;
  const rendered = readFileSync(svg, "utf8");
  for (const row of sourceRecords) expect(rendered).toContain(`data-global-id="${row.global_id}"`);
  expect((rendered.match(/data-service-demand="true"/g) ?? [])).toHaveLength(24);
  expect((rendered.match(/data-service-demand="false"/g) ?? [])).toHaveLength(3);
  expect(rendered).toContain("冷水 是｜热水 否｜排水 是");
  expect(rendered).toContain("冷水 否｜热水 否｜排水 是");
  expect(rendered).toContain("冷热排需求：unknown");
  expect(rendered).toContain("无 IfcDistributionPort / System / 正式管线拓扑");
  expect(rendered).not.toContain("<polyline");
  expect(rendered).not.toContain("<path");
  expect(statSync(svg).size).toBeGreaterThan(10_000);
  expect(statSync(png).size).toBeGreaterThan(10_000);
  expect(readFileSync(png).subarray(0, 8).toString("hex")).toBe("89504e470d0a1a0a");
  expect(sha256(ifc)).toBe(frozenHash);
}, 30_000);

test("P-201 rejects a stale caller-frozen IFC hash", () => {
  const temporary = mkdtempSync(join(tmpdir(), "p201-water-stale-"));
  const run = Bun.spawnSync([
    "python3", script,
    "--root", root,
    "--expected-ifc-sha256", "0".repeat(64),
    "--output-svg", join(temporary, "candidate.svg"),
    "--proof-png", join(temporary, "proof.png"),
    "--report", join(temporary, "report.json"),
  ], { cwd: root });
  expect(run.exitCode).not.toBe(0);
  expect(run.stderr.toString()).toContain("caller-frozen hash");
});

test("P-201 generator has no IFC or Blender write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).toContain("exact_media_evidence_complete(records)");
  expect(source).not.toContain("confirmed_media_count == 3");
  expect(source).not.toContain("ifcopenshell");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
  expect(source).not.toContain("save_ifc_file");
  expect(source).not.toContain("bpy");
  expect(source).not.toContain("Blender");
});
