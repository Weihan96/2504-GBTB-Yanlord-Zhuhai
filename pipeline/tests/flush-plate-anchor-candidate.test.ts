import { describe, expect, test } from "bun:test";
import { readFileSync } from "node:fs";

const script = readFileSync(
  "pipeline/scripts/flush_plate_anchor_candidate.py",
  "utf8",
);

describe("Geberit 115.770 installation anchor", () => {
  test("scope is the exact two thin wall plates", () => {
    expect(script).toContain('"2gFgcOYEXEaQWAzcKulTFt"');
    expect(script).toContain('"2lDPsdQevFSfeOThtjSlPG"');
    expect(script).toContain('EXPECTED_TYPE_NAME = "Geberit 115.770"');
    expect(script).toContain("EXPECTED_SIZE_SORTED_MM");
  });

  test("candidate uses the installation face and does not rewrite semantics", () => {
    expect(script).toContain('"integer_installation_face_centre"');
    expect(script).toContain('"semantic_write_performed": False');
    expect(script).toContain('"contact_regressions"');
    expect(script).toContain('"non_target_geometry_changes"');
  });
});
