import { expect, test } from "bun:test";
import {
  mkdirSync,
  mkdtempSync,
  readFileSync,
  readdirSync,
  writeFileSync,
} from "node:fs";
import { createHash } from "node:crypto";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/release_candidate_gate.py");

type Fixture = {
  root: string;
  ifc: string;
  register: string;
  report: string;
  reportRegister: string;
};

function sha256(path: string): string {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

function makeFixture(options: { stale?: boolean; duplicate?: boolean; malformedIfc?: boolean } = {}): Fixture {
  const fixtureRoot = mkdtempSync(join(tmpdir(), "release-gate-"));
  const ifc = join(fixtureRoot, "formal.ifc");
  const register = join(fixtureRoot, "drawing-register.csv");
  const report = join(fixtureRoot, "professional-report.json");
  const reportRegister = join(fixtureRoot, "release-report-register.csv");
  mkdirSync(join(fixtureRoot, "output"));
  if (options.malformedIfc) {
    writeFileSync(ifc, "ISO-10303-21;\nEND-ISO-10303-21;\n");
  } else {
    const createIfc = Bun.spawnSync([
      "python3",
      "-c",
      "import ifcopenshell,sys; f=ifcopenshell.file(schema='IFC4'); f.create_entity('IfcProject',GlobalId=ifcopenshell.guid.new(),Name='Test'); f.write(sys.argv[1])",
      ifc,
    ]);
    if (createIfc.exitCode !== 0) throw new Error(createIfc.stderr.toString());
  }
  writeFileSync(join(fixtureRoot, "output", "A-001.pdf"), "candidate");
  writeFileSync(join(fixtureRoot, "output", "report-proof.png"), "proof");
  const duplicate = options.duplicate
    ? "A-001,Duplicate,candidate,output/A-001.pdf,duplicate\n"
    : "";
  writeFileSync(
    register,
    [
      "sheet_number,title,status,publish_target,notes",
      "A-001,Index,candidate,output/A-001.pdf,ready",
      duplicate.trimEnd(),
      "D-601,Details,planned,,requires review",
    ]
      .filter(Boolean)
      .join("\n") + "\n",
  );
  writeFileSync(
    report,
    JSON.stringify({
      source_ifc_sha256: options.stale ? "0".repeat(64) : sha256(ifc),
      outputs: { proof: "output/report-proof.png" },
    }),
  );
  writeFileSync(
    reportRegister,
    "report_id,workstream,report_path,required_stage\nTEST,TEST,professional-report.json,reviewed-candidate\n",
  );
  return { root: fixtureRoot, ifc, register, report, reportRegister };
}

function snapshot(path: string): string[] {
  const visit = (current: string, prefix = ""): string[] =>
    readdirSync(current, { withFileTypes: true }).flatMap((entry) => {
      const absolute = join(current, entry.name);
      const relative = join(prefix, entry.name);
      return entry.isDirectory()
        ? visit(absolute, relative)
        : [`${relative}:${sha256(absolute)}`];
    });
  return visit(path).sort();
}

function runGate(fixture: Fixture, stage: string) {
  const process = Bun.spawnSync(
    [
      "python3",
      script,
      "--stage",
      stage,
      "--root",
      fixture.root,
      "--ifc",
      fixture.ifc,
      "--drawing-register",
      fixture.register,
      "--report-register",
      fixture.reportRegister,
    ],
    { cwd: root, stdout: "pipe", stderr: "pipe" },
  );
  return {
    exitCode: process.exitCode,
    stdout: process.stdout.toString(),
    stderr: process.stderr.toString(),
  };
}

test("reviewed candidate discloses planned drawings and writes no view file", () => {
  const fixture = makeFixture();
  const before = snapshot(fixture.root);
  const run = runGate(fixture, "reviewed-candidate");
  const after = snapshot(fixture.root);
  expect(run.exitCode, run.stderr).toBe(0);
  expect(run.stderr).toBe("");
  const output = JSON.parse(run.stdout);
  expect(output.read_only).toBe(true);
  expect(output.stage).toBe("reviewed-candidate");
  expect(output.result.reviewed_candidate_ready).toBe(true);
  expect(output.source.professional_reports).toHaveLength(1);
  expect(output.result.construction_release_candidate_ready).toBe(false);
  expect(output.summary.planned_count).toBe(1);
  expect(output.unresolved[0].sheet_number).toBe("D-601");
  expect(before).toEqual(after);
});

test("construction release rejects planned drawings and missing targets", () => {
  const fixture = makeFixture();
  const run = runGate(fixture, "construction-release-candidate");
  expect(run.exitCode).toBe(1);
  const output = JSON.parse(run.stdout);
  expect(output.result.pass).toBe(false);
  expect(output.result.construction_release_candidate_ready).toBe(false);
  expect(output.checks.find((check: any) => check.id === "PLANNED-DRAWINGS").status).toBe("fail");
  expect(output.checks.find((check: any) => check.id === "DRAWING-OUTPUTS").status).toBe("fail");
});

test("every stage rejects a professional report with an old IFC hash", () => {
  const fixture = makeFixture({ stale: true });
  const run = runGate(fixture, "reviewed-candidate");
  expect(run.exitCode).toBe(1);
  const output = JSON.parse(run.stdout);
  const check = output.checks.find((item: any) => item.id === "PROFESSIONAL-REPORT-HASHES");
  expect(check.status).toBe("fail");
  expect(check.details.reports[0].current_ifc_hash).toBe(false);
  expect(check.details.reports[0].ifc_hashes).toEqual(["0".repeat(64)]);
});

test("release gate atomically refreshes an explicit JSON snapshot", () => {
  const fixture = makeFixture();
  const output = join(fixture.root, "build", "release", "reviewed-current.json");
  const result = Bun.spawnSync([
    "python3",
    script,
    "--stage",
    "reviewed-candidate",
    "--root",
    fixture.root,
    "--ifc",
    fixture.ifc,
    "--drawing-register",
    fixture.register,
    "--report-register",
    fixture.reportRegister,
    "--output",
    output,
  ]);
  expect(result.exitCode).toBe(0);
  const stdout = JSON.parse(result.stdout.toString());
  const saved = JSON.parse(readFileSync(output, "utf8"));
  expect(saved).toEqual(stdout);
  expect(saved.source.formal_ifc_sha256).toBe(sha256(fixture.ifc));
});

test("duplicate drawing numbers fail the gate", () => {
  const fixture = makeFixture({ duplicate: true });
  const run = runGate(fixture, "reviewed-candidate");
  expect(run.exitCode).toBe(1);
  const output = JSON.parse(run.stdout);
  const check = output.checks.find((item: any) => item.id === "DRAWING-NUMBERS");
  expect(check.status).toBe("fail");
  expect(check.details.duplicates).toEqual(["A-001"]);
});

test("every stage rejects a malformed IFC even when its file hash is current", () => {
  const fixture = makeFixture({ malformedIfc: true });
  const run = runGate(fixture, "reviewed-candidate");
  expect(run.exitCode).toBe(1);
  const output = JSON.parse(run.stdout);
  const check = output.checks.find((item: any) => item.id === "FORMAL-IFC");
  expect(check.status).toBe("fail");
  expect(check.message).toContain("cannot be parsed");
});

test("canonical professional reports cannot be omitted by the caller", () => {
  const fixture = makeFixture({ stale: true });
  const run = runGate(fixture, "reviewed-candidate");
  expect(run.exitCode).toBe(1);
  const output = JSON.parse(run.stdout);
  expect(output.summary.report_count).toBe(1);
  expect(output.checks.find((item: any) => item.id === "PROFESSIONAL-REPORT-INVENTORY").status).toBe("pass");
  expect(output.checks.find((item: any) => item.id === "PROFESSIONAL-REPORT-HASHES").status).toBe("fail");
});

test("construction release rejects an intrinsically failing professional report", () => {
  const fixture = makeFixture();
  writeFileSync(
    fixture.register,
    "sheet_number,title,status,publish_target,notes\nA-001,Index,candidate,output/A-001.pdf,ready\n",
  );
  writeFileSync(
    fixture.report,
    JSON.stringify({ source_ifc_sha256: sha256(fixture.ifc), status: false }),
  );
  const run = runGate(fixture, "construction-release-candidate");
  expect(run.exitCode).toBe(1);
  const output = JSON.parse(run.stdout);
  const check = output.checks.find((item: any) => item.id === "PROFESSIONAL-REPORT-READINESS");
  expect(check.status).toBe("fail");
  expect(check.details.blocked_reports[0].blockers).toEqual(["status=false"]);
});

test("reviewed candidates disclose explicit construction readiness blockers", () => {
  const fixture = makeFixture();
  writeFileSync(
    fixture.register,
    "sheet_number,title,status,publish_target,notes\nA-001,Index,candidate,output/A-001.pdf,ready\n",
  );
  writeFileSync(
    fixture.report,
    JSON.stringify({
      source_ifc_sha256: sha256(fixture.ifc),
      gates: { construction_release_ready: false },
    }),
  );
  const run = runGate(fixture, "reviewed-candidate");
  expect(run.exitCode).toBe(0);
  const output = JSON.parse(run.stdout);
  const check = output.checks.find((item: any) => item.id === "PROFESSIONAL-REPORT-READINESS");
  expect(check.status).toBe("disclosed");
  expect(check.details.blocked_reports[0].blockers).toEqual([
    "gates.construction_release_ready=false",
  ]);
  expect(output.result.construction_release_candidate_ready).toBe(false);
});

test("construction release becomes ready only in the construction stage with no blockers", () => {
  const fixture = makeFixture();
  writeFileSync(
    fixture.register,
    "sheet_number,title,status,publish_target,notes\nA-001,Index,candidate,output/A-001.pdf,ready\n",
  );
  const run = runGate(fixture, "construction-release-candidate");
  expect(run.exitCode, run.stderr).toBe(0);
  const output = JSON.parse(run.stdout);
  expect(output.result.reviewed_candidate_ready).toBe(true);
  expect(output.result.construction_release_candidate_ready).toBe(true);
});

test("reviewed candidates reject a report artifact whose bytes drift from its declared hash", () => {
  const fixture = makeFixture();
  const artifact = join(fixture.root, "output", "A-001.pdf");
  writeFileSync(
    fixture.report,
    JSON.stringify({
      source_ifc_sha256: sha256(fixture.ifc),
      output_records: [
        {
          path: "output/A-001.pdf",
          passes: true,
          sha256: sha256(artifact),
        },
      ],
    }),
  );
  writeFileSync(artifact, "stale candidate bytes");
  const run = runGate(fixture, "reviewed-candidate");
  expect(run.exitCode).toBe(1);
  const output = JSON.parse(run.stdout);
  const check = output.checks.find((item: any) => item.id === "REPORT-ARTIFACT-INTEGRITY");
  expect(check.status).toBe("fail");
  expect(check.details.records[0].exists).toBe(true);
  expect(check.details.records[0].hash_matches).toBe(false);
  expect(check.details.errors[0]).toContain("artifact SHA-256 does not match");
});

test("reviewed candidates reject a report whose declared source dependency drifts", () => {
  const fixture = makeFixture();
  const dependency = join(fixture.root, "source.csv");
  writeFileSync(dependency, "current source\n");
  writeFileSync(
    fixture.report,
    JSON.stringify({
      source_ifc_sha256: sha256(fixture.ifc),
      source_dependencies: [
        {
          path: "source.csv",
          sha256: sha256(dependency),
        },
      ],
    }),
  );
  writeFileSync(dependency, "drifted source\n");
  const run = runGate(fixture, "reviewed-candidate");
  expect(run.exitCode).toBe(1);
  const output = JSON.parse(run.stdout);
  const check = output.checks.find(
    (item: any) => item.id === "REPORT-SOURCE-DEPENDENCY-INTEGRITY",
  );
  expect(check.status).toBe("fail");
  expect(check.details.records[0].exists).toBe(true);
  expect(check.details.records[0].hash_matches).toBe(false);
  expect(check.details.errors[0]).toContain("source dependency SHA-256 does not match");
});

test("construction release recognizes every declared blocker shape used by project reports", () => {
  const fixture = makeFixture();
  writeFileSync(
    fixture.register,
    "sheet_number,title,status,publish_target,notes\nA-001,Index,candidate,output/A-001.pdf,ready\n",
  );
  writeFileSync(
    fixture.report,
    JSON.stringify({
      source_ifc_sha256: sha256(fixture.ifc),
      summary: { releasable: false, construction_release_pass: false },
      gates: [
        { id: "QA-01", status: "block" },
        { fabrication_dimension_ready: false },
        { whole_home_switch_positioning_complete: false },
      ],
      candidates: [{ final_release_pass: false, release_blocker: "manufacturer review pending" }],
      blockers: [{ issue_id: "INT1-BLOCK-001" }],
      release_blockers: [{ issue_id: "M401-MISS-001" }],
      release_blocks: [{ issue_id: "PLUM-MISS-001" }],
      open_release_items: ["confirm service interface"],
      blocking_input_ids: ["A101-MAIN-BAY-01"],
    }),
  );
  const run = runGate(fixture, "construction-release-candidate");
  expect(run.exitCode).toBe(1);
  const output = JSON.parse(run.stdout);
  const blockers = output.checks.find(
    (item: any) => item.id === "PROFESSIONAL-REPORT-READINESS",
  ).details.blocked_reports[0].blockers;
  expect(blockers).toContain("summary.releasable=false");
  expect(blockers).toContain("summary.construction_release_pass=false");
  expect(blockers).toContain("gates[0].status=block");
  expect(blockers).toContain("gates[1].fabrication_dimension_ready=false");
  expect(blockers).toContain("gates[2].whole_home_switch_positioning_complete=false");
  expect(blockers).toContain("candidates[0].final_release_pass=false");
  expect(blockers).toContain("candidates[0].release_blocker=nonempty");
  expect(blockers).toContain("blockers=nonempty");
  expect(blockers).toContain("release_blockers=nonempty");
  expect(blockers).toContain("release_blocks=nonempty");
  expect(blockers).toContain("open_release_items=nonempty");
  expect(blockers).toContain("blocking_input_ids=nonempty");
});
