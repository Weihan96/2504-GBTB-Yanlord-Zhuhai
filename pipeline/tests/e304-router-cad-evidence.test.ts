import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { spawnSync } from "bun";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/e304_router_cad_evidence.py");
const output = resolve(root, "build/elec/e304-router-cad-evidence.test.json");
const svg = resolve(root, "build/elec/E304-router-entry-cabinet-evidence.test.svg");

test("official CAD pins the router to the entry weak-current cabinet without inventing Z", () => {
  const run = spawnSync(["python3", script, "--output", output, "--output-svg", svg], { cwd: root });
  expect(run.exitCode).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  const evidenceIds = report.evidence_register.rows.map((row: { evidence_id: string }) => row.evidence_id);
  expect(evidenceIds).toEqual([
    "E304-CAD-001",
    "E304-CAD-002",
    "E304-USER-001",
    "E304-PHOTO-001",
    "E304-PHOTO-002",
    "E304-PHOTO-003",
    "E304-PHOTO-004",
    "E304-PHOTO-005",
    "E304-USER-002",
  ]);
  expect(report.evidence_register.rows.find((row: { evidence_id: string }) => row.evidence_id === "E304-PHOTO-003").source_sha256).toBe(
    "3d01b1292d7620389e3920cfbd3420a3f0bb1fd77d754a0f59c48ec56dc4d384",
  );
  expect(report.evidence_register.rows.find((row: { evidence_id: string }) => row.evidence_id === "E304-PHOTO-005").source_sha256).toBe(
    "24ef14dd6871b469a74ee167302d894214b0df30cc31eb4029df2a18f4d8a41b",
  );
  const photo003 = report.evidence_register.rows.find(
    (row: { evidence_id: string }) => row.evidence_id === "E304-PHOTO-003",
  );
  expect(photo003.proves).toContain("设备和线缆在弱电箱内存在");
  expect(photo003.does_not_prove).toContain("准确型号角色");
  expect(photo003.status).toBe("verified_user_photo");
  expect(photo003.confidence).toBe("1.00");
  expect(photo003.review_required).toBe("yes");
  expect(photo003.formal_ifc_write_allowed).toBe("no");
  expect(report.text_evidence.map((row: { handle: string }) => row.handle)).toEqual(["2598C0", "224270", "224295"]);
  expect(report.leader_evidence.handle).toBe("224271");
  expect(report.source.coordinate_transform.viewport_handle).toBe("224238");
  expect(report.weak_current_box.ifc_plan_position_mm[0]).toBeCloseTo(4600.016493, 6);
  expect(report.weak_current_box.ifc_plan_position_mm[1]).toBeCloseTo(-735.368749, 6);
  expect(report.weak_current_box.weak_current_box_bottom_aff_mm).toBe(350);
  expect(report.router_decision.installation_z_mm).toBeNull();
  expect(report.gates.router_z_not_inferred_from_weak_box_datum).toBe(true);
  expect(report.gates.all_site_photo_hashes_and_dimensions_match).toBe(true);
  const photos = report.evidence_register.rows.filter((row: { source_kind: string }) => row.source_kind === "user_site_photo");
  expect(photos).toHaveLength(5);
  expect(photos.every((row: { verified_file?: { sha256: string; pixel_dimensions: number[] } }) =>
    Boolean(row.verified_file?.sha256) && row.verified_file!.pixel_dimensions.length === 2)).toBe(true);
  expect(report.gates.automatic_ifc_write_allowed).toBe(false);
  expect(readFileSync(svg, "utf8")).toContain("H+350 是弱电箱底边，不是路由器安装高度");
}, 20_000);

test("site-photo evidence fails closed when a registered hash drifts", () => {
  const temp = mkdtempSync(join(tmpdir(), "e304-photo-gate-"));
  try {
    const register = resolve(root, "pipeline/decisions/elec-source-evidence.csv");
    const drifted = readFileSync(register, "utf8").replace(
      "5f428bf20247e8a131ac4315e1f0c3958ab0e0e2823966f96cd7d80f4e1df938",
      "0".repeat(64),
    );
    const driftedRegister = join(temp, "evidence.csv");
    writeFileSync(driftedRegister, drifted);
    const run = spawnSync([
      "python3", script,
      "--evidence-register", driftedRegister,
      "--output", join(temp, "out.json"),
      "--output-svg", join(temp, "out.svg"),
    ], { cwd: root });
    expect(run.exitCode).not.toBe(0);
  } finally {
    rmSync(temp, { recursive: true, force: true });
  }
}, 20_000);

test("CAD evidence extractor has no IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
