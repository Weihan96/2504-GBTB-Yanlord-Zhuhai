import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/wfin_semantics_candidate.py");

test("WFIN semantics remain a read-only project intent boundary", async () => {
  const source = await Bun.file(script).text();
  expect(source).toContain('PSET_NAME = "Pset_WallFinishIntent"');
  expect(source).toContain("MULTI_FINISH_SPLIT_REQUIRED");
  expect(source).toContain("CandidateSegmentSummary");
  expect(source).toContain('"standard_status": "project-specific design-intent Pset; not a buildingSMART standard Pset"');
  expect(source).toContain('"ExistingMaterialAssociationPreserved": True');
  expect(source).toContain('"FormalIfcWriteAllowed": False');
  expect(source).toContain('"automatic_formal_ifc_write_allowed": False');
  expect(source).not.toContain("save_and_load_ifc(");
});

test("WFIN open issues preserve the exact C003 handoff boundary", async () => {
  const source = await Bun.file(script).text();
  expect(source).toContain('"0WdsBKMo545Ryicomy6Mqg"');
  expect(source).toContain('"0moBrr1cf5WgtBzqXWKtqn"');
  expect(source).toContain('"2ntxn4aYnB0xQOSuraf1r2"');
  expect(source).toContain('"WFIN-R01"');
  expect(source).toContain('"WFIN-R04"');
  expect(source).toContain('"finish-boundary-and-missing-face-review"');
  expect(source).toContain("06wFwLoDD6ie5iCTnc_yad");
  expect(source).toContain("0e0XOb$L18ZBVYJiQJrQ1p");
  expect(source).toContain("3bMoS7bIT8wBmjqRvdPunu");
});
