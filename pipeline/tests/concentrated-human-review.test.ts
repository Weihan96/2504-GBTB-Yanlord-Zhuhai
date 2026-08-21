import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import {
  appendFileSync,
  copyFileSync,
  existsSync,
  mkdirSync,
  mkdtempSync,
  readFileSync,
  writeFileSync,
} from "node:fs";
import { tmpdir } from "node:os";
import { basename, join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/concentrated_human_review.py");
const sources = {
  decisions: resolve(root, "pipeline/decisions/owner-input-register.csv"),
  rules: resolve(root, "pipeline/decisions/owner-input-closeout-rules.csv"),
  equipment: resolve(root, "pipeline/decisions/equipment-register.csv"),
  requirements: resolve(root, "pipeline/decisions/equipment-installation-requirements.csv"),
  drawings: resolve(root, "pipeline/decisions/drawing-register.csv"),
  reviewed: resolve(root, "build/release/reviewed-candidate-current.json"),
  construction: resolve(root, "build/release/construction-candidate-current.json"),
  ifc: resolve(root, "2504 GBTB Yanlord Zhuhai.ifc"),
};

function sha256(path: string): string {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

function parseCsv(path: string): Record<string, string>[] {
  const text = readFileSync(path, "utf8").replace(/^\uFEFF/, "");
  const lines: string[][] = [];
  let row: string[] = [];
  let field = "";
  let quoted = false;
  for (let index = 0; index < text.length; index += 1) {
    const char = text[index];
    if (char === '"') {
      if (quoted && text[index + 1] === '"') {
        field += '"';
        index += 1;
      } else quoted = !quoted;
    } else if (char === "," && !quoted) {
      row.push(field);
      field = "";
    } else if ((char === "\n" || char === "\r") && !quoted) {
      if (char === "\r" && text[index + 1] === "\n") index += 1;
      row.push(field);
      if (row.some((value) => value.length > 0)) lines.push(row);
      row = [];
      field = "";
    } else field += char;
  }
  if (field || row.length) {
    row.push(field);
    lines.push(row);
  }
  const [header, ...records] = lines;
  return records.map((values) => Object.fromEntries(header.map((key, index) => [key, values[index] ?? ""])));
}

function runCompiler(extra: string[] = []) {
  const outputRoot = mkdtempSync(join(tmpdir(), "concentrated-review-output-"));
  const jsonOutput = join(outputRoot, "review.json");
  const markdownOutput = join(outputRoot, "review.md");
  const process = Bun.spawnSync(
    [
      "python3",
      script,
      "--root",
      root,
      "--json-output",
      jsonOutput,
      "--markdown-output",
      markdownOutput,
      ...extra,
    ],
    { cwd: root, stdout: "pipe", stderr: "pipe" },
  );
  return { process, jsonOutput, markdownOutput };
}

test("compiler preserves decision status and evidence-gated closeout status", () => {
  const sourceHashes = Object.fromEntries(Object.entries(sources).map(([key, path]) => [key, sha256(path)]));
  const { process, jsonOutput, markdownOutput } = runCompiler();
  expect(process.exitCode, process.stderr.toString()).toBe(0);
  const report = JSON.parse(readFileSync(jsonOutput, "utf8"));
  const decisions = parseCsv(sources.decisions);
  const rules = new Map(parseCsv(sources.rules).map((row) => [row.input_id, row]));
  const reportItems = new Map(report.decision_review_items.map((row: any) => [row.input_id, row]));
  expect(report.decision_review_items).toHaveLength(decisions.length);
  for (const decision of decisions) {
    const output: any = reportItems.get(decision.input_id);
    expect(output.decision_status).toBe(decision.status);
    const rule = rules.get(decision.input_id)!;
    if (decision.status === "采用候选" && rule.automatic_close_allowed === "no") {
      expect(output.closeout_status).toBe("decision_confirmed_evidence_pending");
      expect(output.closeout_status).not.toBe("verified");
    }
    if (decision.status === "不适用") {
      expect(output.closeout_status).toBe("not_applicable");
    }
  }
  expect(report.formal_ifc.sha256).toBe(sha256(sources.ifc));
  expect(report.source_ifc_sha256).toBe(report.formal_ifc.sha256);
  expect(report.formal_ifc.write_prohibited).toBe(true);
  expect(report.gates.all_release_blocking_requirements_mapped_once).toBe(true);
  expect(report.gates.all_blocking_requirements_mapped).toBe(true);
  expect(report.gates.construction_release_ready).toBe(false);
  expect(report.summary.unmapped_blocker_count).toBe(0);
  expect(report.summary.review_item_count).toBe(report.review_items.length);
  expect(report.summary.open_root_review_item_count).toBe(report.review_items.length);
  expect(report.summary.recovery_stream_counts).toEqual({
    authority_review: 2,
    owner_preference: 4,
    project_internal: 0,
    site_evidence: 15,
    vendor_review: 25,
  });
  expect(report.gates.all_open_owner_inputs_routed).toBe(true);
  expect(report.gates.all_requirement_packages_have_delivery_target).toBe(true);
  expect(report.gates.all_root_items_have_delivery_target).toBe(true);
  for (const item of report.review_items) {
    expect(item.delivery_target.length).toBeGreaterThan(0);
    if (item.review_item_kind === "owner_input") expect(item.recovery_stream.length).toBeGreaterThan(0);
    if (item.review_item_kind === "requirement_package") {
      expect(item.delivery_target.startsWith("sheet:")).toBe(false);
    }
  }
  expect(readFileSync(markdownOutput, "utf8")).toContain(report.formal_ifc.sha256);
  for (const [key, path] of Object.entries(sources)) expect(sha256(path)).toBe(sourceHashes[key]);
});

test("all current INT1 blockers form twenty-two stable and lossless review packages", () => {
  const { process, jsonOutput } = runCompiler();
  expect(process.exitCode, process.stderr.toString()).toBe(0);
  const report = JSON.parse(readFileSync(jsonOutput, "utf8"));
  const currentInt1 = parseCsv(sources.requirements)
    .filter((row) => row.blocks_release === "yes" && row.discipline.toUpperCase().split("/").includes("INT1"))
    .map((row) => row.requirement_id)
    .sort();
  const currentInt1Set = new Set(currentInt1);
  const packages = report.requirement_review_packages.filter((row: any) => row.root_cause_id.startsWith("INT1-"));
  const projected = packages
    .flatMap((row: any) => row.requirement_ids)
    .filter((requirementId: string) => currentInt1Set.has(requirementId))
    .sort();
  expect(packages).toHaveLength(22);
  expect(projected).toEqual(currentInt1);
  expect(new Set(projected).size).toBe(projected.length);
  for (const reviewPackage of packages) {
    expect(reviewPackage.equipment_ids.length).toBeGreaterThan(0);
    expect(reviewPackage.affected_sheets.length).toBeGreaterThan(0);
    expect(reviewPackage.responsible_party.length).toBeGreaterThan(0);
    expect(reviewPackage.required_evidence.length).toBeGreaterThan(0);
    expect(reviewPackage.delivery_target.length).toBeGreaterThan(0);
  }
});

test("an unmapped release-blocking root cause fails closed before writing views", () => {
  const fixture = mkdtempSync(join(tmpdir(), "concentrated-review-fixture-"));
  const fixtureDecisions = join(fixture, basename(sources.decisions));
  const fixtureRules = join(fixture, basename(sources.rules));
  const fixtureEquipment = join(fixture, basename(sources.equipment));
  const fixtureRequirements = join(fixture, basename(sources.requirements));
  const fixtureDrawings = join(fixture, basename(sources.drawings));
  for (const [source, target] of [
    [sources.decisions, fixtureDecisions],
    [sources.rules, fixtureRules],
    [sources.equipment, fixtureEquipment],
    [sources.requirements, fixtureRequirements],
    [sources.drawings, fixtureDrawings],
  ]) copyFileSync(source, target);
  appendFileSync(
    fixtureRequirements,
    "REQ-FAIL-CLOSED,APP-001,UNKNOWN,new_unmapped_root_cause,,,,,pending,pending,,,yes,fixture\n",
  );
  const fixtureIfc = join(fixture, "formal.ifc");
  writeFileSync(fixtureIfc, "read-only IFC hash fixture\n");
  const hash = sha256(fixtureIfc);
  const reviewed = join(fixture, "reviewed.json");
  const construction = join(fixture, "construction.json");
  writeFileSync(reviewed, JSON.stringify({ source: { formal_ifc_sha256: hash }, result: {} }));
  writeFileSync(construction, JSON.stringify({ source: { formal_ifc_sha256: hash }, result: {} }));
  const outputDir = join(fixture, "output");
  mkdirSync(outputDir);
  const jsonOutput = join(outputDir, "review.json");
  const markdownOutput = join(outputDir, "review.md");
  const process = Bun.spawnSync(
    [
      "python3",
      script,
      "--root",
      fixture,
      "--ifc",
      fixtureIfc,
      "--owner-input-register",
      fixtureDecisions,
      "--owner-input-closeout-rules",
      fixtureRules,
      "--equipment-register",
      fixtureEquipment,
      "--requirements",
      fixtureRequirements,
      "--drawing-register",
      fixtureDrawings,
      "--reviewed-release",
      reviewed,
      "--construction-release",
      construction,
      "--json-output",
      jsonOutput,
      "--markdown-output",
      markdownOutput,
    ],
    { cwd: root, stdout: "pipe", stderr: "pipe" },
  );
  expect(process.exitCode).toBe(1);
  expect(process.stderr.toString()).toContain("unmapped release-blocking root causes");
  expect(existsSync(jsonOutput)).toBe(false);
  expect(existsSync(markdownOutput)).toBe(false);
  expect(sha256(fixtureIfc)).toBe(hash);
});

test("missing derived release snapshots do not prevent canonical review compilation", () => {
  const missingRoot = mkdtempSync(join(tmpdir(), "concentrated-review-missing-release-"));
  const { process, jsonOutput } = runCompiler([
    "--reviewed-release", join(missingRoot, "reviewed.json"),
    "--construction-release", join(missingRoot, "construction.json"),
  ]);
  expect(process.exitCode, process.stderr.toString()).toBe(0);
  const report = JSON.parse(readFileSync(jsonOutput, "utf8"));
  expect(report.gates.current_release_reports_match_formal_ifc).toBe(false);
  expect(report.gates.all_release_blocking_requirements_mapped_once).toBe(true);
  expect(report.release_context.reviewed_candidate.available).toBe(false);
  expect(report.release_context.construction_release_candidate.available).toBe(false);
});
