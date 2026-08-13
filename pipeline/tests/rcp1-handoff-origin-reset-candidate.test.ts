import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const targetPath = resolve(
  root,
  "pipeline/decisions/rcp1-handoff-origin-reset-targets.csv",
);
const reportPath = resolve(
  root,
  "build/coordinate-normalization/rcp1-handoff-origin-reset-candidate.test.json",
);
const candidatePath = resolve(
  root,
  "build/candidates/2504-GBTB-rcp1-handoff-origin-reset.test.ifc",
);

test(
  "RCP1 confirmed handoffs keep world geometry while resetting origins",
  async () => {
    const result = Bun.spawnSync(
      [
        "python3",
        "pipeline/scripts/structural_origin_reset_candidate.py",
        "--input",
        "2504 GBTB Yanlord Zhuhai.ifc",
        "--targets",
        targetPath,
        "--output",
        candidatePath,
        "--report",
        reportPath,
        "--tolerance-mm",
        "0.1",
      ],
      { cwd: root, stdout: "pipe", stderr: "pipe" },
    );
    expect(result.exitCode).toBe(0);

    const report = await Bun.file(reportPath).json();
    expect(report.gates.pass).toBe(true);
    expect(report.gates.target_products).toBe(2);
    expect(report.gates.targets_by_class).toEqual({
      IfcBuildingElementProxy: 2,
    });
    expect(report.gates.target_placement_mismatches).toEqual([]);
    expect(report.gates.anchors_over_tolerance).toBe(0);
    expect(report.gates.product_geometry_over_tolerance).toBe(0);
    expect(report.gates.root_global_ids_equal).toBe(true);
    expect(report.gates.schema_equal).toBe(true);
    expect(report.gates.entity_count_delta_matches_expected).toBe(true);

    const byId = new Map(
      report.results.map((item: { global_id: string }) => [item.global_id, item]),
    );
    expect(byId.get("0zWtSQZzjFQg_PORjlssbe").target_mm).toEqual([
      5800, 0, 2400,
    ]);
    expect(byId.get("16Ey9Flj9BK9VRun$ozzjH").target_mm).toEqual([
      5800, 2970, 2500,
    ]);
    expect(
      Math.max(...report.all_product_geometry.records.map(
        (item: { world_corresponding_vertex_max_delta_mm: number | null }) =>
          item.world_corresponding_vertex_max_delta_mm ?? 0,
      )),
    ).toBeLessThanOrEqual(0.1);
  },
  600_000,
);
