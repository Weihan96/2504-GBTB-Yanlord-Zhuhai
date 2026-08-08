import { expect, test } from "bun:test";

const script = "pipeline/scripts/plum_existing_candidate.py";
const source = await Bun.file(script).text();

test("PLUM candidate protects observable counts and PVC110 identity", () => {
  expect(source).toContain('len(sanitary) != 27');
  expect(source).toContain('len(waste) != 3');
  expect(source).toContain('len(assemblies) != 3');
  expect(source).toContain('len(drainage) != 16');
  expect(source).toContain('"178mqyyzzFowLcbXcH6prO"');
  expect(source).toContain('"0bfVg4Ys1CevZs$qxhkXTo"');
  expect(source).toContain('qa["pvc110_branch_count"] == 6');
  expect(source).toContain('qa["pvc110_world_geometry_unchanged"]');
});

test("P-201 endpoints never become inferred connectors", () => {
  expect(source).toContain('"connection_requirement": "unknown"');
  expect(source).toContain('"is_ifc_distribution_port": False');
  expect(source).toContain('"candidate_is_write_authority": False');
  expect(source).toContain('"cold-water connection inference"');
  expect(source).toContain('"hot-water connection inference"');
  expect(source).toContain('"world ObjectPlacement for existing-object registration only; not a connector or rough-in point"');
  expect(source).toContain('return [float(value) for value in matrix[:3, 3]]');
  expect(source).not.toContain('value * 1000.0');
});

test("missing system data is disclosed instead of passed", () => {
  expect(source).toContain('"status": "data_missing"');
  expect(source).toContain('"connectivity_qa_passed": False');
  expect(source).toContain('qa["construction_release_pass"] = False');
  expect(source).toContain('"candidate_registry_pass"');
});

test("PLUM candidate runs against the frozen formal IFC", async () => {
  const process = Bun.spawn(["python3", script], { stdout: "pipe", stderr: "pipe" });
  const [exitCode, stdout, stderr] = await Promise.all([
    process.exited,
    new Response(process.stdout).text(),
    new Response(process.stderr).text(),
  ]);
  expect(stderr).toBe("");
  expect(exitCode).toBe(0);
  expect(stdout).toContain('"p201_demand_endpoint_count": 27');
  expect(stdout).toContain('"p202_existing_object_count": 49');
  const report = await Bun.file("build/plum/plum-report.json").json();
  expect(report.source.ifc_sha256).toBe(
    "c7295688003f3f36775a25f6adc2c9878e203c52c980f8e46faed66a3537c4a8",
  );
  expect(report.qa.candidate_registry_pass).toBe(true);
  expect(report.qa.construction_release_pass).toBe(false);
  expect(report.qa.distribution_data.status).toBe("data_missing");
  expect(report.qa.pvc110_world_geometry_unchanged).toBe(true);
  const p202 = await Bun.file("build/plum/p202-existing-location-register.json").json();
  expect(
    Math.max(...p202.objects.flatMap((row: any) => row.object_origin_mm.map(Math.abs))),
  ).toBeLessThan(20_000);
}, 30_000);
