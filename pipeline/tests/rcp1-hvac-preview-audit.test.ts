import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/rcp1_hvac_preview_audit.py");
const routeReportPath = resolve(root, "build/rcp1/route-readiness-candidate.json");
const previewPath = resolve(root, "build/rcp1/hvac-route-preview-audit.test.input.json");
const outputPath = resolve(root, "build/rcp1/hvac-route-preview-audit.test.json");

test("HVAC preview audit preserves IFC anchors and allows reviewed bends", async () => {
  const source = await Bun.file(script).text();
  expect(source).not.toContain("model.write(");
  const report = await Bun.file(routeReportPath).json();
  const routes = report.confirmed_route_graph
    .filter((route: any) => route.waypoints.length >= 2)
    .map((route: any) => ({
      route_id: route.route_id,
      status: "blender_preview_candidate_human_review_required",
      anchors: route.waypoints.map((row: any) => ({
        anchor_id: row.anchor_id,
        anchor_kind: row.anchor_kind,
        global_id: row.anchor_global_id || null,
        world_mm: row.centre_mm,
      })),
    }));
  routes[0].anchors.splice(1, 0, {
    anchor_id: "TEST_BEND_01",
    anchor_kind: "blender_bend",
    global_id: null,
    world_mm: [-1500, -500, 2700],
  });
  await Bun.write(previewPath, JSON.stringify({
    mode: "blender_hvac_route_preview_candidate",
    source_ifc_sha256: report.source.ifc_sha256,
    routes,
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
  expect(audit.gates.fixed_anchor_order_matches).toBe(true);
  expect(audit.gates.fixed_anchors_within_tolerance).toBe(true);
  expect(audit.gates.temporary_bends_require_human_review).toBe(true);
  expect(audit.gates.formal_ifc_write_allowed).toBe(false);
}, 30_000);
