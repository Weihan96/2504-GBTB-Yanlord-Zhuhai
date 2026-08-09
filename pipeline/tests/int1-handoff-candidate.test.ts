import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("INT1 handoff classifies the exact 16 proxies without authorizing writes", () => {
  const temporary = mkdtempSync(join(tmpdir(), "int1-handoff-"));
  const output = join(temporary, "report.json");
  const csv = join(temporary, "review.csv");
  const result = Bun.spawnSync([
    "python3",
    "pipeline/scripts/int1_handoff_candidate.py",
    "--input",
    "2504 GBTB Yanlord Zhuhai.ifc",
    "--p0-review",
    "pipeline/decisions/p0-review.csv",
    "--origin-review",
    "build/coordinate-normalization/remaining-origin-review.json",
    "--anchor-audit",
    "build/coordinate-normalization/integer-geometry-anchor-audit.json",
    "--decision-csv",
    csv,
    "--output",
    output,
  ], { cwd: root });
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.mode).toBe("read_only_int1_c003_handoff_candidate");
  expect(report.summary).toEqual({
    handoff_object_count: 16,
    records_by_sheet: { "I-501": 11, "I-502": 1, "I-504": 4 },
    named_worktop_or_countertop_candidates: 3,
    named_bathroom_pipe_wall_candidates: 1,
    legacy_cad_references: 4,
    unresolved_joinery_proxies: 8,
    shape_error_count: 3,
    integer_geometry_anchor_count: 2,
  });
  expect(report.gates.source_reports_match_formal_ifc).toBe(true);
  expect(report.gates.handoff_set_matches_exactly).toBe(true);
  expect(report.gates.all_objects_classified).toBe(true);
  expect(report.gates.all_objects_have_fabrication_identity).toBe(false);
  expect(report.gates.all_geometry_readable_as_solid_mesh).toBe(false);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(report.gates.fabrication_dimensions_ready).toBe(false);
  expect(report.records.every((row: any) => row.review_required === true)).toBe(true);
  expect(report.records.every((row: any) => row.automatic_ifc_write_allowed === false)).toBe(true);
  expect(readFileSync(csv, "utf8").split("\n").filter(Boolean)).toHaveLength(17);
});

test("INT1 handoff candidate contains no IFC write path", () => {
  const source = readFileSync(resolve(root, "pipeline/scripts/int1_handoff_candidate.py"), "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run(");
});
