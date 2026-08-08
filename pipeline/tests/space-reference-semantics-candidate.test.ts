import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const candidateScript = resolve(root, "pipeline/scripts/space_reference_semantics_candidate.py");
const postwriteScript = resolve(root, "pipeline/scripts/space_reference_postwrite_audit.py");

test("S003 candidate uses a deterministic exact write boundary", async () => {
  const source = await Bun.file(candidateScript).text();
  expect(source).toContain('GUID_NAMESPACE = uuid.UUID(');
  expect(source).toContain('properties={"Reference": row["candidate_reference"]}');
  expect(source).toContain('"automatic_formal_ifc_write_allowed": False');
  expect(source).not.toContain("save_and_load_ifc(");
});

test("S003 postwrite audit checks semantic identity and geometry", async () => {
  const source = await Bun.file(postwriteScript).text();
  expect(source).toContain('formal_roots != candidate_roots');
  expect(source).toContain('current_reference(formal_space)');
  expect(source).toContain('maximum_space_world_vertex_change_mm');
  expect(source).not.toContain("model.write(");
});
