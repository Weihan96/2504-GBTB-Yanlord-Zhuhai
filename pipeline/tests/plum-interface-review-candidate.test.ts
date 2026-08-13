import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { mkdtempSync, readFileSync, statSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/plum_interface_review_candidate.py");
const ifc = resolve(root, "2504 GBTB Yanlord Zhuhai.ifc");
const requirementsCsv = resolve(root, "pipeline/decisions/equipment-installation-requirements.csv");

function sha256(path: string): string {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

function parseCsv(text: string): string[][] {
  const rows: string[][] = [];
  let row: string[] = [];
  let field = "";
  let quoted = false;
  const source = text.replace(/^\uFEFF/, "");
  for (let index = 0; index < source.length; index += 1) {
    const character = source[index];
    if (quoted) {
      if (character === '"' && source[index + 1] === '"') {
        field += '"';
        index += 1;
      } else if (character === '"') {
        quoted = false;
      } else {
        field += character;
      }
    } else if (character === '"') {
      quoted = true;
    } else if (character === ",") {
      row.push(field);
      field = "";
    } else if (character === "\n") {
      row.push(field.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      field = "";
    } else {
      field += character;
    }
  }
  if (field || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  return rows;
}

function csvRecords(path: string): Record<string, string>[] {
  const [header, ...rows] = parseCsv(readFileSync(path, "utf8"));
  return rows
    .filter((row) => row.some(Boolean))
    .map((row) => Object.fromEntries(header.map((name, index) => [name, row[index] ?? ""])));
}

test("PLUM interface review projects every PLUM requirement without touching IFC", async () => {
  const temporary = mkdtempSync(join(tmpdir(), "plum-interface-review-"));
  const outputJson = join(temporary, "review.json");
  const outputCsv = join(temporary, "review.csv");
  const outputSvg = join(temporary, "review.svg");
  const outputPng = join(temporary, "review.png");
  const frozenIfcHash = sha256(ifc);
  const sourceRequirements = csvRecords(requirementsCsv).filter((row) => row.discipline.includes("PLUM"));
  const expectedEquipmentIds = new Set(sourceRequirements.map((row) => row.equipment_id));
  const expectedBlocking = sourceRequirements.filter((row) => row.blocks_release === "yes").length;

  const process = Bun.spawn([
    "python3", script,
    "--root", root,
    "--expected-ifc-sha256", frozenIfcHash,
    "--output-json", outputJson,
    "--output-csv", outputCsv,
    "--output-svg", outputSvg,
    "--output-png", outputPng,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  const [exitCode, stdout, stderr] = await Promise.all([
    process.exited,
    new Response(process.stdout).text(),
    new Response(process.stderr).text(),
  ]);
  expect(exitCode).toBe(0);
  expect(stderr).not.toContain("Traceback");
  expect(stdout).toContain('"automatic_ifc_write_allowed": false');

  const report = await Bun.file(outputJson).json();
  expect(report.summary.equipment_count).toBe(expectedEquipmentIds.size);
  expect(report.summary.requirement_count).toBe(sourceRequirements.length);
  expect(report.summary.blocking_requirement_count).toBe(expectedBlocking);
  expect(report.equipment.map((row: any) => row.equipment_id).sort()).toEqual([...expectedEquipmentIds].sort());
  expect(report.equipment.flatMap((row: any) => row.requirements)).toHaveLength(sourceRequirements.length);
  expect(new Set(report.equipment.flatMap((row: any) => row.requirements.map((item: any) => item.requirement_id))).size)
    .toBe(sourceRequirements.length);

  expect(Object.keys(report.summary.ifc_mapping_counts).sort()).toEqual([
    "linked_single_occurrence",
    "multiple_occurrences_review_required",
    "no_formal_ifc_occurrence",
  ]);
  expect(report.summary.requirement_classification_counts.exact_model_confirmed).toBeGreaterThan(0);
  expect(report.summary.requirement_classification_counts.family_or_project_candidate).toBeGreaterThan(0);
  expect(report.summary.requirement_classification_counts.pending_blocking).toBeGreaterThan(0);

  expect(report.source_validation.all_referenced_source_ids_exist).toBe(true);
  expect(report.source_validation.missing_source_ids).toEqual([]);
  expect(report.source_validation.local_evidence_hash_check_count).toBeGreaterThan(0);
  expect(report.source_validation.local_evidence_hashes_match).toBe(true);
  expect(report.source_validation.local_evidence_hash_mismatches).toEqual([]);

  expect(report.component_relationship_candidates).toHaveLength(3);
  expect(report.component_relationship_candidates.every((row: any) =>
    row.candidate_only === true && row.automatic_ifc_write_allowed === false,
  )).toBe(true);
  const relationships = JSON.stringify(report.component_relationship_candidates).toLowerCase();
  expect(relationships).not.toContain("coordinate");
  expect(relationships).not.toContain("direction");
  expect(relationships).not.toContain("placement");

  expect(report.gates).toMatchObject({
    all_plum_equipment_projected: true,
    all_plum_requirements_projected: true,
    all_referenced_source_ids_exist: true,
    local_evidence_hashes_match: true,
    component_relationships_are_candidate_only: true,
    contains_connector_coordinates: false,
    contains_route_directions: false,
    automatic_ifc_write_allowed: false,
    construction_release_ready: false,
    ifc_unchanged_during_generation: true,
  });
  expect(report.source_ifc_sha256).toBe(frozenIfcHash);
  expect(report.source.ifc.sha256).toBe(frozenIfcHash);
  expect(report.source.ifc.ending_sha256).toBe(frozenIfcHash);
  expect(sha256(ifc)).toBe(frozenIfcHash);

  const csvRows = csvRecords(outputCsv);
  expect(csvRows).toHaveLength(sourceRequirements.length);
  expect(csvRows.every((row) =>
    row.automatic_ifc_write_allowed === "false" && row.construction_release_ready === "false",
  )).toBe(true);
  const svg = readFileSync(outputSvg, "utf8");
  expect(svg).toContain("PLUM 接口审核候选矩阵");
  expect(svg).not.toContain("<image");
  expect(statSync(outputSvg).size).toBeGreaterThan(10_000);
  expect(statSync(outputPng).size).toBeGreaterThan(50_000);
  expect(readFileSync(outputPng).subarray(0, 8).toString("hex")).toBe("89504e470d0a1a0a");
});

test("PLUM interface review rejects a stale frozen IFC hash", () => {
  const temporary = mkdtempSync(join(tmpdir(), "plum-interface-stale-"));
  const process = Bun.spawnSync([
    "python3", script,
    "--root", root,
    "--expected-ifc-sha256", "0".repeat(64),
    "--output-json", join(temporary, "review.json"),
    "--output-csv", join(temporary, "review.csv"),
    "--output-svg", join(temporary, "review.svg"),
    "--output-png", join(temporary, "review.png"),
  ], { cwd: root });
  expect(process.exitCode).not.toBe(0);
  expect(process.stderr.toString()).toContain("formal IFC hash changed");
});

test("PLUM interface review has no Blender or IFC write path", () => {
  const source = readFileSync(script, "utf8");
  expect(source).not.toContain("ifcopenshell");
  expect(source).not.toContain("bpy");
  expect(source).not.toContain("save_ifc_file");
  expect(source).not.toContain("model.write(");
  expect(source).not.toContain("ifcopenshell.api.run");
});
