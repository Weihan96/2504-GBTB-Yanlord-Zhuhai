import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/rcp1_hvac_preview_audit.py");
const routeReportPath = resolve(root, "build/rcp1/route-readiness-candidate.json");
const previewPath = resolve(root, "build/rcp1/hvac-route-preview-audit.test.input.json");
const outputPath = resolve(root, "build/rcp1/hvac-route-preview-audit.test.json");

function buildOrthogonalPoints(anchors: any[]) {
  const points = [anchors[0].world_mm];
  for (const target of anchors.slice(1)) {
    const current = points.at(-1)!;
    const xPoint = [target.world_mm[0], current[1], current[2]];
    const yPoint = [target.world_mm[0], target.world_mm[1], current[2]];
    const zPoint = [...target.world_mm];
    for (const point of [xPoint, yPoint, zPoint]) {
      if (point.some((value, index) => value !== points.at(-1)![index])) {
        points.push(point);
      }
    }
  }
  return points;
}

test("HVAC preview audit preserves IFC anchors and allows reviewed bends", async () => {
  const source = await Bun.file(script).text();
  expect(source).not.toContain("model.write(");
  const report = await Bun.file(routeReportPath).json();
  const routes = report.confirmed_route_graph
    .filter((route: any) => route.waypoints.length >= 2)
    .map((route: any) => {
      const anchors = route.waypoints.map((row: any, index: number) => ({
        anchor_id: row.anchor_id,
        anchor_kind: row.anchor_kind,
        global_id: row.anchor_global_id || null,
        world_mm: index === 0 && row.anchor_kind === "equipment"
          ? [row.centre_mm[0] + 400, row.centre_mm[1], row.centre_mm[2]]
          : row.centre_mm,
        placement_basis: row.anchor_kind === "equipment"
          ? "equipment_local_positive_x_service_face"
          : "ifc_object_bbox_centre",
        reference_world_mm: row.centre_mm,
        service_local_mm: row.anchor_kind === "equipment" ? [800, 200, 0] : null,
        service_axis_world: row.anchor_kind === "equipment" ? [1, 0, 0] : null,
      }));
      return {
        route_id: route.route_id,
        status: "blender_preview_candidate_human_review_required",
        anchors,
        orthogonal_points_world_mm: buildOrthogonalPoints(anchors),
      };
    });
  routes[0].anchors.splice(1, 0, {
    anchor_id: "TEST_BEND_01",
    anchor_kind: "blender_bend",
    global_id: null,
    world_mm: [-1500, -500, 2700],
  });
  routes[0].orthogonal_points_world_mm = buildOrthogonalPoints(routes[0].anchors);
  await Bun.write(previewPath, JSON.stringify({
    mode: "blender_hvac_route_preview_candidate",
    source_ifc_sha256: report.source.ifc_sha256,
    routes,
    preview_parameters: { fillet_radius_m: 0.12 },
    formal_ifc_write_allowed: false,
  }));

  const process = Bun.spawn([
    "python3",
    script,
    "--input",
    resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
    "--route-report",
    routeReportPath,
    "--preview",
    previewPath,
    "--output",
    outputPath,
    "--tolerance-mm",
    "0.1",
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  const [exitCode, stderr] = await Promise.all([
    process.exited,
    new Response(process.stderr).text(),
  ]);
  expect(stderr).toBe("");
  expect(exitCode).toBe(0);
  const audit = await Bun.file(outputPath).json();
  expect(audit.summary.route_count).toBe(2);
  expect(audit.summary.fixed_anchor_check_count).toBe(5);
  expect(audit.summary.temporary_bend_count).toBe(1);
  expect(audit.summary.maximum_fixed_anchor_deviation_mm).toBe(0);
  expect(audit.summary.maximum_equipment_service_offset_from_reference_mm).toBe(400);
  expect(audit.summary.all_routes_orthogonal).toBe(true);
  expect(audit.summary.all_route_anchors_on_path_in_order).toBe(true);
  expect(audit.gates.fixed_anchor_order_matches).toBe(true);
  expect(audit.gates.fixed_anchors_within_tolerance).toBe(true);
  expect(audit.gates.equipment_service_axes_cardinal).toBe(true);
  expect(audit.gates.all_preview_segments_axis_aligned).toBe(true);
  expect(audit.gates.all_route_anchors_on_path_in_order).toBe(true);
  expect(audit.gates.geometry_nodes_fillet_present).toBe(true);
  expect(audit.gates.temporary_bends_require_human_review).toBe(true);
  expect(audit.gates.formal_ifc_write_allowed).toBe(false);
}, 30_000);
