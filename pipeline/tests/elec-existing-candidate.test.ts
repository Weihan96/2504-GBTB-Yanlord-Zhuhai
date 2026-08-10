import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const scriptPath = resolve(root, "pipeline/scripts/elec_existing_candidate.py");
const reviewPath = resolve(root, "pipeline/decisions/elec-existing-review.csv");
const reportPath = resolve(root, "build/elec/elec-existing-candidate.test.json");
const e301TestSvg = resolve(root, "build/elec/E301-lighting-location-candidate.test.svg");
const e303TestSvg = resolve(root, "build/elec/E303-socket-equipment-location-candidate.test.svg");

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
  expect(source).toContain("不含灯组、回路、功率或控制推断");
  expect(source).toContain("不含回路、功率、防水或接口推断");
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
        "--source-svg",
        resolve(root, "drawings/Wall Plan.svg"),
        "--e301-svg",
        e301TestSvg,
        "--e303-svg",
        e303TestSvg,
      ],
      { cwd: root, stdout: "pipe", stderr: "pipe" },
    );
    expect(result.exitCode, result.stderr.toString()).toBe(0);
    const report = JSON.parse(await Bun.file(reportPath).text());
    expect(report.source.sha256).toBe(
      "6c2fd8da9e9ad7ddbc2b63415a27f1c979e8995b880d8fce210a2dda2ef2aab6",
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
    expect(report.gates.drawing_candidate_pass).toBe(true);
    expect(report.drawing_gates["E-301"].marker_count).toBe(79);
    expect(report.drawing_gates["E-301"].label_collisions).toBe(0);
    expect(report.drawing_gates["E-303"].socket_markers).toBe(11);
    expect(report.drawing_gates["E-303"].typed_equipment_markers).toBe(8);
    expect(report.drawing_gates["E-303"].proxy_markers).toBe(9);
    expect(report.drawing_gates["E-303"].label_collisions).toBe(0);
    const e301 = await Bun.file(e301TestSvg).text();
    const e303 = await Bun.file(e303TestSvg).text();
    expect((e301.match(/data-elec-kind="light"/g) ?? []).length).toBe(79);
    expect((e303.match(/data-elec-kind="socket(?:-exception)?"/g) ?? []).length).toBe(11);
    expect((e303.match(/data-elec-kind="equipment"/g) ?? []).length).toBe(8);
    expect((e303.match(/data-elec-kind="proxy"/g) ?? []).length).toBe(9);
    expect(e301).toContain('width="500mm"');
    expect(e303).toContain('viewBox="0 0 500 400"');
    expect(e301).toContain("Wall Plan-underlay.png");
    expect(e303).toContain("Wall Plan-underlay.png");
  },
  30_000,
);
