import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/rcp1_route_readiness_candidate.py");
const output = resolve(root, "build/rcp1/route-readiness-candidate.test.json");

test("RCP1 route graph accepts shared openings and ordered multi-hop paths", async () => {
  const source = await Bun.file(script).text();
  expect(source).toContain('"blender_preview_is_source_of_truth": False');
  expect(source).toContain('"equipment_positions_may_move": False');
  expect(source).not.toContain("itertools.permutations");
  expect(source).not.toContain("model.write(");

  const process = Bun.spawn(
    [
      "python3",
      script,
      "--input",
      resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
      "--hvac-report",
      resolve(root, "build/rcp1/hvac-remodel-candidate.json"),
      "--legacy-report",
      resolve(root, "build/rcp1/legacy-base-audit.json"),
      "--delivery-dwg",
      resolve(root, "../图纸/矩阵纵横/D户型交付竣工图.dwg"),
      "--routes",
      resolve(root, "pipeline/decisions/rcp1-hvac-route-register.csv"),
      "--waypoints",
      resolve(root, "pipeline/decisions/rcp1-hvac-route-waypoints.csv"),
      "--equipment-register",
      resolve(root, "pipeline/decisions/equipment-register.csv"),
      "--requirements",
      resolve(root, "pipeline/decisions/equipment-installation-requirements.csv"),
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
  expect(report.openings).toHaveLength(7);
  expect(report.airside_readiness.a02_and_a06_both_serve_dining).toBe(true);
  const routes = Object.fromEntries(report.confirmed_route_graph.map((item: any) => [item.route_id, item]));
  expect(routes["RCP1-SERVICE-A02"].waypoints.map((item: any) => item.anchor_id)).toEqual(["A02", "H03"]);
  expect(routes["RCP1-SERVICE-A03"].waypoints.map((item: any) => item.anchor_id)).toEqual(["A03", "H04", "H02"]);
  expect(report.interface_roles.H01.role).toBe("outdoor_unit_interface_at_existing_opening");
  expect(report.interface_roles.H07.role).toBe("condensate_discharge_interface_at_existing_opening");
  expect(report.interface_roles.H01.ifc_entity_remains).toBe("IfcOpeningElement");
  expect(report.interface_roles.H07.ifc_entity_remains).toBe("IfcOpeningElement");
  expect(report.legacy_reference.component_topology_matches).toBe(true);
  expect(report.legacy_reference.maximum_component_dimension_difference_mm).toBeLessThan(0.016);
  expect(report.gates.confirmed_shared_and_multihop_graph_ready).toBe(true);
  expect(report.gates.manufacturer_constraints_consumed_from_ssot).toBe(true);
  expect(report.gates.a01_a04_refrigerant_and_condensate_nominals_available).toBe(true);
  const a02Requirements = Object.fromEntries(report.manufacturer_interface_inputs.by_equipment.A02.requirements
    .map((item: any) => [item.requirement_name, item.value]));
  expect(a02Requirements.gas_pipe_od).toBe(12.7);
  expect(a02Requirements.liquid_pipe_od).toBe(6.35);
  expect(a02Requirements.drain_pipe_od).toBe(32);
  const a05Requirements = Object.fromEntries(report.manufacturer_interface_inputs.by_equipment.A05.requirements
    .map((item: any) => [item.requirement_name, item.value]));
  expect(a05Requirements.service_access_width_min).toBe(450);
  expect(a05Requirements.supply_duct_length_min).toBe(1000);
  expect(report.manufacturer_interface_inputs.by_equipment.A06.release_blocking_requirements)
    .toContain("nominal_body_width");
  expect(report.gates.formal_ifc_write_allowed).toBe(false);
}, 30_000);
