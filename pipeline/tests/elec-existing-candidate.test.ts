import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const scriptPath = resolve(root, "pipeline/scripts/elec_existing_candidate.py");
const reviewPath = resolve(root, "pipeline/decisions/elec-existing-review.csv");
const reportPath = resolve(root, "build/elec/elec-existing-candidate.test.json");

test("ELEC candidate protects exact socket exceptions and proxy handoffs", async () => {
  const source = await Bun.file(scriptPath).text();
  const review = await Bun.file(reviewPath).text();
  for (const globalId of [
    "0laejMoxn8Lu_X3FZaCXmi",
    "27MTenki57DQsfMryX_1U0",
    "2OOjqQDMHDjRcXQCniWXnp",
    "3KXtmVvejA78j_iAS$kydj",
  ]) {
    expect(source).toContain(globalId);
    expect(review).toContain(globalId);
  }
  expect((review.match(/,proxy_handoff,/g) ?? []).length).toBe(9);
  expect(source).toContain('"automatic_ifc_write_allowed": False');
  expect(source).toContain('"construction_release_ready"] = False');
  expect(source).not.toContain("model.write(");
});

test("ELEC candidate separates unused type definitions from instances", async () => {
  const source = await Bun.file(scriptPath).text();
  expect(source).toContain('"used": used_types, "unused": unused_types');
  expect(source).toContain('"available_uninstantiated_type_definitions": control_types');
  expect(source).toContain('"available_uninstantiated_type_definitions": network_types');
  expect(source).toContain('gates["e302_switch_instances"] == 0');
  expect(source).toContain('gates["e304_network_instances"] == 0');
});

test(
  "ELEC candidate inventories the formal IFC without writing it",
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
      ],
      { cwd: root, stdout: "pipe", stderr: "pipe" },
    );
    expect(result.exitCode, result.stderr.toString()).toBe(0);
    const report = JSON.parse(await Bun.file(reportPath).text());
    expect(report.source.sha256).toBe(
      "c7295688003f3f36775a25f6adc2c9878e203c52c980f8e46faed66a3537c4a8",
    );
    expect(report.gates.candidate_pass).toBe(true);
    expect(report.gates.construction_release_ready).toBe(false);
    expect(report.gates.light_count).toBe(79);
    expect(report.gates.socket_count).toBe(11);
    expect(report.gates.typed_equipment_count).toBe(8);
    expect(report.gates.proxy_handoff_count).toBe(9);
    expect(report.gates.used_appliance_type_count).toBe(9);
    expect(report.gates.unused_appliance_type_count).toBe(33);
    expect(report.gates.e302_switch_instances).toBe(0);
    expect(report.gates.e304_network_instances).toBe(0);
    expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  },
  30_000,
);
