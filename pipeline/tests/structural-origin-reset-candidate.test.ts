import { expect, test } from "bun:test";
import { mkdtempSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { spawnSync } from "node:child_process";

const claddingTargets =
  "pipeline/decisions/c003-cladding-origin-reset-targets.csv";

test("structural origin reset parses mixed classes and rejects duplicates", () => {
  const dir = mkdtempSync(join(tmpdir(), "structural-origin-"));
  const csv = join(dir, "targets.csv");
  writeFileSync(
    csv,
    "expected_class,global_id,target_x_mm,target_y_mm,target_z_mm,anchor_kind,basis,confidence\nIfcBeam,A,0,0,0,vertex,test,1\nIfcSlab,B,100,200,0,edge,test,1\nIfcCovering,C,-100,0,0,edge,test,1\nIfcFurniture,D,200,-300,0,edge,test,1\n",
  );
  const ok = spawnSync(
    "python3",
    [
      "-c",
      `import sys; sys.path.insert(0,"pipeline/scripts"); from pathlib import Path; from structural_origin_reset_candidate import read_targets; print([(r["expected_class"],r["global_id"]) for r in read_targets(Path(sys.argv[1]))])`,
      csv,
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(ok.status).toBe(0);
  expect(ok.stdout.trim()).toBe(
    "[('IfcBeam', 'A'), ('IfcSlab', 'B'), ('IfcCovering', 'C'), ('IfcFurniture', 'D')]",
  );
  writeFileSync(
    csv,
    "expected_class,global_id,target_x_mm,target_y_mm,target_z_mm,anchor_kind,basis,confidence\nIfcBeam,A,0,0,0,vertex,test,1\nIfcSlab,A,100,200,0,edge,test,1\n",
  );
  const duplicate = spawnSync(
    "python3",
    [
      "-c",
      `import sys; sys.path.insert(0,"pipeline/scripts"); from pathlib import Path; from structural_origin_reset_candidate import read_targets; read_targets(Path(sys.argv[1]))`,
      csv,
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(duplicate.status).not.toBe(0);
  rmSync(dir, { recursive: true, force: true });
});

test("cladding origin targets are five unique IfcCovering anchors", () => {
  const result = spawnSync(
    "python3",
    [
      "-c",
      `import sys; sys.path.insert(0,"pipeline/scripts"); from pathlib import Path; from structural_origin_reset_candidate import read_targets; rows=read_targets(Path(sys.argv[1])); print(len(rows)); print(sorted({row["expected_class"] for row in rows})); print(len({row["global_id"] for row in rows}))`,
      claddingTargets,
    ],
    { cwd: process.cwd(), encoding: "utf8" },
  );
  expect(result.status).toBe(0);
  expect(result.stdout.trim().split("\n")).toEqual([
    "5",
    "['IfcCovering']",
    "5",
  ]);
});
