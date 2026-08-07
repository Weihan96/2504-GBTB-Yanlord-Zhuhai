import { expect, test } from "bun:test";
import { mkdtempSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(
  root,
  "pipeline/scripts/furniture_group_origin_candidate.py",
);

test("furniture anchor targets require an integer approved destination", () => {
  const dir = mkdtempSync(join(tmpdir(), "furniture-anchor-"));
  const csv = join(dir, "targets.csv");
  const header =
    "expected_class,global_id,source_anchor_x_mm,source_anchor_y_mm,source_anchor_z_mm,target_x_mm,target_y_mm,target_z_mm,anchor_kind,contact_group,basis,confidence,human_review_required,status\n";
  writeFileSync(
    csv,
    `${header}IfcFurniture,A,0.2,-0.4,0,0,0,0,existing_vertex,G1,test,0.99,no,approved\n`,
  );
  const ok = Bun.spawnSync(
    [
      "python3",
      "-c",
      `import sys; from pathlib import Path; sys.path.insert(0,"pipeline/scripts"); from furniture_group_origin_candidate import read_targets; row=read_targets(Path(sys.argv[1]))[0]; print(row["translation_mm"].tolist())`,
      csv,
    ],
    { cwd: root },
  );
  expect(ok.exitCode).toBe(0);
  expect(ok.stdout.toString().trim()).toBe("[-0.2, 0.4, 0.0]");
  writeFileSync(
    csv,
    `${header}IfcFurniture,A,10,-20,30,10,-20,30,existing_surface_point,G1,test,1.0,no,approved\n`,
  );
  const surfacePoint = Bun.spawnSync(
    [
      "python3",
      "-c",
      `import sys; from pathlib import Path; sys.path.insert(0,"pipeline/scripts"); from furniture_group_origin_candidate import read_targets; row=read_targets(Path(sys.argv[1]))[0]; print(row["translation_mm"].tolist())`,
      csv,
    ],
    { cwd: root },
  );
  expect(surfacePoint.exitCode).toBe(0);
  expect(surfacePoint.stdout.toString().trim()).toBe("[0.0, 0.0, 0.0]");
  writeFileSync(
    csv,
    `${header}IfcFurniture,A,0.2,-0.4,0,0.1,0,0,existing_vertex,G1,test,0.99,no,approved\n`,
  );
  const invalid = Bun.spawnSync(
    [
      "python3",
      "-c",
      `import sys; from pathlib import Path; sys.path.insert(0,"pipeline/scripts"); from furniture_group_origin_candidate import read_targets; read_targets(Path(sys.argv[1]))`,
      csv,
    ],
    { cwd: root },
  );
  expect(invalid.exitCode).not.toBe(0);
  rmSync(dir, { recursive: true, force: true });
});

test("contact changes distinguish allowed closure from regressions", () => {
  const source = `
import json,sys
sys.path.insert(0,"pipeline/scripts")
from furniture_group_origin_candidate import classify_contact_changes
A=("A","B"); C=("C","D"); E=("E","F")
print(json.dumps(classify_contact_changes({A:0.02,C:0.4},{A:0.03,C:0.0,E:0.0},0.1,{C})))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.regressions).toHaveLength(0);
  expect(report.new_contacts).toHaveLength(2);
  expect(report.unexpected_new_contacts).toHaveLength(1);
  expect(report.missing_allowed_contacts).toHaveLength(0);
});
