import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/rcp1_route_readiness_candidate.py");
const output = resolve(root, "build/rcp1/route-readiness-candidate.test.json");

test("RCP1 route readiness keeps fixed positions and exposes missing endpoints", async () => {
  const source = await Bun.file(script).text();
  expect(source).toContain('"equipment_positions_may_move": False');
  expect(source).toContain('"demolition_walls_are_permanent_obstacles": False');
  expect(source).not.toContain("model.write(");

  const process = Bun.spawn(
    [
      "python3",
      script,
      "--input",
      resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
      "--hvac-report",
      resolve(root, "build/rcp1/hvac-remodel-candidate.json"),
      "--output",
      output,
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
  expect(report.equipment).toHaveLength(6);
  expect(report.equipment.every((item: any) => item.position_status === "confirmed_fixed")).toBe(true);
  expect(report.openings).toHaveLength(7);
  expect(report.high_confidence_equipment_opening_pairings.map((item: any) => [item.equipment_id, item.opening_id])).toEqual([
    ["A03", "H04"],
    ["A04", "H06"],
    ["A05", "H07"],
  ]);
  expect(report.global_equipment_opening_assignment.pairs.map((item: any) => [item.equipment_id, item.opening_id])).toEqual([
    ["A01", "H05"],
    ["A02", "H02"],
    ["A03", "H04"],
    ["A04", "H06"],
    ["A05", "H07"],
    ["A06", "H03"],
  ]);
  expect(report.global_equipment_opening_assignment.unused_opening_ids).toEqual(["H01"]);
  expect(report.global_equipment_opening_assignment.best_to_second_margin_mm).toBeGreaterThan(800);
  expect(report.shared_or_route_opening_evidence[0].opening_id).toBe("H05");
  expect(report.refrigerant_and_condensate_readiness.distribution_port_count).toBe(0);
  expect(report.refrigerant_and_condensate_readiness.outdoor_condenser_or_compressor_count).toBe(0);
  expect(report.minimum_required_inputs).toHaveLength(3);
  expect(report.gates.global_six_equipment_assignment_ready).toBe(true);
  expect(report.gates.demolition_walls_excluded_from_permanent_obstacles).toBe(true);
  expect(report.gates.formal_ifc_write_allowed).toBe(false);
}, 30_000);
