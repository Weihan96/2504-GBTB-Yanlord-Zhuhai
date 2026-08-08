import { expect, test } from "bun:test";
import { mkdtempSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/p202_existing_candidate.py");
const source = await Bun.file(script).text();

test("P-202 preserves the observable object inventory", () => {
  expect(source).toContain('(len(sanitary), len(waste), len(assemblies), len(drainage)) != (27, 3, 3, 16)');
  expect(source).toContain('qa["sanitary_marker_count"] == 27');
  expect(source).toContain('qa["waste_marker_count"] == 3');
  expect(source).toContain('qa["drainage_marker_count"] == 16');
  expect(source).toContain('qa["assembly_marker_count"] == 3');
});

test("P-202 labels placement markers as non-connectors", () => {
  expect(source).toContain('existing-object-placement-not-connector');
  expect(source).toContain('not a connector or rough-in point');
  expect(source).toContain('"connectivity_qa_passed": False');
  expect(source).toContain('qa["construction_release_pass"] = False');
  expect(source).toContain('return [float(value) for value in matrix[:3, 3]]');
});

test("P-202 protects PVC110 products and six branches", () => {
  expect(source).toContain('qa["pvc110_product_count"] == 2');
  expect(source).toContain('qa["pvc110_branch_count"] == 6');
  expect(source).toContain('qa["pvc110_world_geometry_unchanged"]');
});

test("P-202 candidate renders and passes its mechanical gate", async () => {
  const temporary = mkdtempSync(join(tmpdir(), "p202-candidate-"));
  const reportPath = join(temporary, "report.json");
  const process = Bun.spawn([
    "python3", script,
    "--root", root,
    "--output-svg", join(temporary, "candidate.svg"),
    "--output-pdf", join(temporary, "candidate.pdf"),
    "--report", reportPath,
    "--render-report", join(temporary, "render-report.json"),
    "--proof-png", join(temporary, "proof.png"),
  ], {cwd: root, stdout: "pipe", stderr: "pipe"});
  const [exitCode, stdout, stderr] = await Promise.all([
    process.exited,
    new Response(process.stdout).text(),
    new Response(process.stderr).text(),
  ]);
  expect(exitCode).toBe(0);
  expect(stderr).not.toContain("Traceback");
  expect(stdout).toContain('"candidate_mechanical_pass": true');
  const report = await Bun.file(reportPath).json();
  expect(report.qa.candidate_mechanical_pass).toBe(true);
  expect(report.qa.construction_release_pass).toBe(false);
  expect(report.qa.connectivity_status).toBe("data_missing");
  expect(report.qa.pdf_page_pass).toBe(true);
  expect(report.qa.off_plan_waste_terminal_ids).toEqual(["3IVqCnhGr51hY4LrOq_5G_"]);
}, 60_000);
