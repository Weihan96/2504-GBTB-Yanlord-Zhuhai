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
  expect(report.summary.cabinet_power_zones).toBe(7);
  expect(report.summary.kitchen_socket_rechecks).toBe(11);
  expect(report.summary.label_collision_count).toBe(0);
  expect(report.summary.developer_red_points_are_reference_only).toBe(true);
  expect(report.gates.four_bedside_lights_present).toBe(true);
  expect(report.gates.island_and_dining_bay_socket_present).toBe(true);
  expect(report.gates.illuminated_cabinet_power_is_grouped_not_fabricated).toBe(true);
  expect(report.gates.all_current_kitchen_sockets_reopened_for_review).toBe(true);
  expect(report.gates.label_collision_free).toBe(true);
  expect(report.gates.reverse_requirement_audit_planned).toBe(true);
  expect(report.gates.whole_home_electrical_positioning_complete).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  const renderedSvg = readFileSync(svg, "utf8");
  expect(renderedSvg).toContain("elec-renovation-round1");
  expect(renderedSvg).not.toContain("Wall Plan-underlay.png");
}, 30_000);

test("first-round renovation candidate has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
