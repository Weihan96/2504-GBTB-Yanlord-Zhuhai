import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const output = resolve(root, "build/rcp1/hvac-interface-candidate.test.json");

test("HVAC interface candidate validates types and preserves unresolved gates", async () => {
  const process = Bun.spawn([
    "python3",
    "pipeline/scripts/rcp1_hvac_interface_candidate.py",
    "--input",
    "2504 GBTB Yanlord Zhuhai.ifc",
    "--m401-review",
    "pipeline/decisions/m401-existing-review.csv",
    "--evidence",
    "pipeline/decisions/rcp1-hvac-equipment-interface-evidence.csv",
    "--manual-dir",
    "tmp/pdfs",
    "--output",
    output,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  const exitCode = await process.exited;
  expect(await new Response(process.stderr).text()).toBe("");
  expect(exitCode).toBe(0);

  const report = await Bun.file(output).json();
  expect(report.mode).toBe("read_only_rcp1_hvac_manufacturer_interface_candidate");
  expect(report.summary).toEqual({
    equipment_count: 6,
    formal_identity_present_count: 5,
    formal_identity_missing_count: 1,
    official_model_family_match_count: 4,
    model_text_mismatch_count: 1,
    precise_port_coordinate_ready_count: 0,
    formal_ifc_write_allowed_count: 0,
  });
  expect(report.gates.source_hash_matches_m401_review).toBe(true);
  expect(report.gates.formal_ifc_type_assignments_match).toBe(true);
  expect(report.gates.available_local_manual_hashes_match).toBe(true);
  expect(report.gates.all_equipment_has_formal_identity).toBe(false);
  expect(report.gates.all_model_texts_have_exact_official_match).toBe(false);
  expect(report.gates.precise_port_coordinates_ready).toBe(false);
  expect(report.gates.formal_ifc_write_allowed).toBe(false);
  expect(report.occurrences.every((row: any) => row.formal_ifc_write_allowed !== true)).toBe(true);
});
