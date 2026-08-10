import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const scriptPath = resolve(root, "pipeline/scripts/rcp1_existing_candidate.py");
const reviewPath = resolve(root, "pipeline/decisions/rcp1-existing-review.csv");
const reportPath = resolve(root, "build/rcp1/rcp1-existing-candidate.test.json");
const coordinationPath = resolve(root, "build/rcp1/coordination-report.test.json");
const ifcPath = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
const currentIfcHash = createHash("sha256").update(readFileSync(ifcPath)).digest("hex");

test("RCP1 source protects handoffs and does not write or invent systems", async () => {
  const source = await Bun.file(scriptPath).text();
  const review = await Bun.file(reviewPath).text();
  for (const globalId of ["16Ey9Flj9BK9VRun$ozzjH", "0zWtSQZzjFQg_PORjlssbe"]) {
    expect(source).toContain(globalId);
    expect(review).toContain(globalId);
  }
  expect(source).toContain('"automatic_ifc_write_allowed": False');
  expect(source).toContain('gates["construction_release_ready"] = False');
  expect(source).toContain('"IfcAirTerminal"');
  expect(source).toContain('"IfcSensor"');
  expect(source).toContain('"IfcAlarm"');
  expect(source).toContain("aabb_distance_mm");
  expect(source).toContain("clash_collision_many");
  expect(source).toContain("clash_intersection_many");
  expect(source).toContain("clash_clearance_many");
  expect(source).toContain('default=0.1');
  expect(source).toContain("expected_light_ceiling_embedding_candidate");
  expect(source).toContain("expected_opening_void_host");
  expect(source).toContain("LEGACY_BASE_PAIR_IDS");
  expect(review).toContain("RCP1-LEGACY-BASE-001");
  expect(review).toContain("2504_lowpoly.blend");
  expect(source).not.toContain("model.write(");
});

test(
  "RCP1 candidate inventories current high-level geometry and explicit gaps",
  async () => {
    const result = Bun.spawnSync(
      [
        "python3",
        scriptPath,
        "--input",
        resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
        "--review",
        reviewPath,
        "--output",
        reportPath,
        "--coordination-output",
        coordinationPath,
      ],
      { cwd: root, stdout: "pipe", stderr: "pipe" },
    );
    expect(result.exitCode, result.stderr.toString()).toBe(0);
    const report = JSON.parse(await Bun.file(reportPath).text());
    expect(report.source.sha256).toBe(currentIfcHash);
    expect(report.source_ifc_sha256).toBe(currentIfcHash);
    expect(report.legacy_evidence.status).toBe("current_external_audit");
    expect(report.legacy_evidence.current_formal_ifc).toBe(true);
    expect(report.gates.candidate_pass).toBe(true);
    expect(report.gates.construction_release_ready).toBe(false);
    expect(report.gates.space_count).toBe(22);
    expect(report.gates.spaces_with_reference).toBe(22);
    expect(report.gates.ceiling_covering_count).toBe(2);
    expect(report.gates.light_fixture_count).toBe(79);
    expect(report.gates.dcl_proxy_count).toBe(18);
    expect(report.gates.name_only_air_outlet_proxy_count).toBe(2);
    expect(report.gates.typed_high_equipment_count).toBe(5);
    expect(report.gates.named_high_opening_count).toBe(7);
    expect(report.gates.high_flow_segment_count).toBe(5);
    expect(report.gates.other_high_proxy_count).toBe(1);
    expect(report.gates.duplicate_inventory_ids).toEqual([]);
    expect(report.gates.protected_handoffs_over_0_1_mm).toBe(0);
    expect(report.gates.missing_instances.IfcSensor).toBe(1);
    expect(Object.entries(report.gates.missing_instances).every(([key, value]) => key === "IfcSensor" || value === 0)).toBe(true);
    expect(Object.values(report.gates.topology).every((value) => value === 0)).toBe(true);
    for (const group of Object.values(report.inventory) as any[][]) {
      for (const item of group) {
        expect(item.bbox.dimensions_mm.length).toBe(3);
        expect(item.elevation_mm.bottom).toBeNumber();
        expect(item.basis.length).toBeGreaterThan(0);
        expect(item.confidence).toBeNumber();
        expect(["yes", "no"]).toContain(item.review_required);
        expect(item.formal_ifc_write_allowed).toBe("no");
      }
    }
    const coordination = JSON.parse(await Bun.file(coordinationPath).text());
    expect(coordination.source.sha256).toBe(report.source.sha256);
    expect(coordination.source_ifc_sha256).toBe(currentIfcHash);
    expect(coordination.method.tolerance_mm).toBe(0.1);
    expect(coordination.summary.object_count).toBe(117);
    expect(coordination.summary.total_pair_count).toBe(6786);
    expect(coordination.summary.aabb_candidate_pair_count).toBeGreaterThan(0);
    expect(coordination.summary.aabb_proven_far_pair_count).toBeGreaterThan(0);
    expect(coordination.summary.geometry_state_counts).toEqual({
      contacting: 12,
      intersecting: 106,
      separated: 6668,
    });
    expect(coordination.summary.human_review_pair_count).toBe(0);
    expect(coordination.summary.human_review_state_counts).toEqual({});
    expect(coordination.summary.unresolved_conflict_candidate_count).toBe(0);
    expect(coordination.summary.legacy_base_pair_count).toBe(10);
    expect(coordination.gates.all_pairs_classified).toBe(true);
    expect(coordination.gates.review_pairs_have_exact_ids_and_basis).toBe(true);
    expect(coordination.gates.review_pairs_exact_mesh_tested).toBe(true);
    expect(coordination.gates.expected_light_or_host_relations_not_conflicts).toBe(true);
    expect(coordination.gates.legacy_base_pairs_preserve_exact_intersection_scope).toBe(true);
    expect(coordination.gates.automatic_ifc_write_allowed).toBe(false);
    expect(coordination.gates.candidate_pass).toBe(true);
    expect(coordination.pairs.length).toBe(6786);
    for (const pair of coordination.human_review_pairs) {
      expect(pair.pair.length).toBe(2);
      expect(pair.review_required).toBe("yes");
      expect(pair.basis).toContain("formal IFC world Body meshes");
      expect(pair.conflict_candidate).toBeBoolean();
      expect(pair.minimum_clearance_candidate_mm).toBeNumber();
      expect(pair.aabb_candidate).toBe(true);
      expect(pair.minimum_clearance_is_lower_bound).toBe(false);
      expect(pair.geometry_state).toBe("intersecting");
      expect(pair.expected_relation_rule).toBe("none");
    }
    expect(coordination.pairs.filter(
      (pair: { geometry_state: string; review_required: string }) =>
        pair.geometry_state === "separated" && pair.review_required === "yes",
    )).toHaveLength(0);
    const legacyBasePairs = coordination.pairs.filter(
      (pair: { expected_relation_rule: string }) =>
        pair.expected_relation_rule === "legacy_hvac_base_intersection_pending_redesign",
    );
    expect(legacyBasePairs).toHaveLength(10);
    for (const pair of legacyBasePairs) {
      expect(pair.geometry_state).toBe("intersecting");
      expect(pair.review_required).toBe("no");
      expect(pair.conflict_candidate).toBe(false);
      expect(pair.basis).toContain("rcp1_legacy_base_audit.py");
    }
  },
  90_000,
);

test("RCP1 rejects a caller-frozen IFC hash mismatch before auditing geometry", () => {
  const result = Bun.spawnSync(
    ["python3", scriptPath, "--input", ifcPath, "--expected-ifc-sha256", "0".repeat(64)],
    { cwd: root, stdout: "pipe", stderr: "pipe" },
  );
  expect(result.exitCode).not.toBe(0);
  expect(result.stderr.toString()).toContain("formal IFC SHA drift");
});

test("RCP1 legacy audit uses a caller freeze instead of a historical IFC constant", async () => {
  const source = await Bun.file(resolve(root, "pipeline/scripts/rcp1_legacy_base_audit.py")).text();
  expect(source).toContain('"--expected-ifc-sha256"');
  expect(source).toContain('"caller_frozen_ifc_hash_match"');
  expect(source).not.toContain("EXPECTED_IFC_SHA256");
});
