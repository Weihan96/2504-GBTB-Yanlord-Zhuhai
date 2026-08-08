import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/p0_semantics_batch_candidate.py");

test("P0 semantic batch is one read-only A104/A105 candidate", async () => {
  const source = await Bun.file(modulePath).text();
  expect(source).toContain("a104.apply_semantics(candidate, a104_rows)");
  expect(source).toContain("a105.apply_semantics(candidate, a105_rows)");
  expect(source).toContain("candidate output must not overwrite the formal IFC");
  expect(source).toContain('"new_root_count"] == 54');
  expect(source).toContain('"formal_ifc_write_allowed": False');
});
