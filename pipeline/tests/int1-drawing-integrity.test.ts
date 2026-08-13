import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/int1_drawing_integrity.py");

test("INT1 locks all 49 formal IFC Drawings and their final provenance", () => {
  const temporary = mkdtempSync(join(tmpdir(), "int1-drawing-integrity-"));
  const reportPath = join(temporary, "report.json");
  const run = Bun.spawnSync(
    ["python3", script, "--root", root, "--output", reportPath],
    { cwd: root, stdout: "pipe", stderr: "pipe" },
  );
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const report = JSON.parse(readFileSync(reportPath, "utf8"));
  expect(report.mode).toBe("read_only_int1_drawing_integrity");
  expect(report.status).toBe(true);
  expect(report.summary).toMatchObject({
    drawing_count: 49,
    plan_count: 5,
    official_elevation_count: 36,
    project_compiled_elevation_count: 8,
    drawing_document_relation_count: 49,
    final_svg_count: 49,
    source_json_count: 44,
    final_provenance_is_aggregate_report: true,
  });
  expect(report.output_records).toHaveLength(93);
  expect(report.drawings).toHaveLength(49);
  expect(report.drawings.filter((row: any) => row.source_json)).toHaveLength(44);
  expect(report.drawings.every((row: any) => row.svg_sha256.length === 64)).toBe(true);
  expect(report.gates).toMatchObject({
    drawing_integrity_pass: true,
    all_document_relations_resolve: true,
    all_final_svg_artifacts_hashed: true,
    all_elevation_source_json_artifacts_hashed: true,
    construction_release_ready: false,
    fabrication_dimension_ready: false,
  });
});

test("INT1 Drawing integrity generator is read-only with respect to IFC and SVG", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("write_text(svg");
  expect(source).toContain("read_only_int1_drawing_integrity");
});
