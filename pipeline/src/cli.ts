import { mkdir } from "node:fs/promises";
import { dirname, resolve } from "node:path";
import { snapshotIfc } from "./ifc-step";
import { buildQaReport, renderQaMarkdown, type DecisionSnapshot, type DrawingSnapshot } from "./qa";

interface DrawingConfig {
  id: string;
  path: string;
  role: string;
  expectedScale: string;
}

interface ProjectConfig {
  projectId: string;
  projectName: string;
  ifcPath: string;
  ifcSchema: string;
  baselineCounts: Record<string, number>;
  drawings: DrawingConfig[];
  decisionFiles: string[];
  output: {
    snapshot: string;
    qaJson: string;
    qaMarkdown: string;
    candidateDirectory: string;
    releaseDirectory: string;
  };
}

const repositoryRoot = resolve(import.meta.dir, "../..");
const configPath = resolve(repositoryRoot, "pipeline/project.json");

async function writeJson(path: string, value: unknown): Promise<void> {
  await mkdir(dirname(path), { recursive: true });
  await Bun.write(path, `${JSON.stringify(value, null, 2)}\n`);
}

async function inspectDrawing(config: DrawingConfig): Promise<DrawingSnapshot> {
  const absolutePath = resolve(repositoryRoot, config.path);
  const file = Bun.file(absolutePath);
  if (!(await file.exists())) {
    return { id: config.id, path: config.path, role: config.role, exists: false, byteSize: 0, scale: null, dimensionLineCount: 0, textCount: 0 };
  }
  const text = await file.text();
  return {
    id: config.id,
    path: config.path,
    role: config.role,
    exists: true,
    byteSize: file.size,
    scale: /data-scale="([^"]+)"/.exec(text)?.[1] ?? null,
    dimensionLineCount: (text.match(/<line\b[^>]*class="[^"]*PredefinedType-DIMENSION[^"]*"/g) ?? []).length,
    textCount: (text.match(/<text\b/g) ?? []).length,
  };
}

export function splitCsvLine(line: string): string[] {
  const fields: string[] = [];
  let value = "";
  let quoted = false;
  for (let index = 0; index < line.length; index += 1) {
    const character = line[index];
    if (character === '"') {
      if (quoted && line[index + 1] === '"') {
        value += '"';
        index += 1;
      } else {
        quoted = !quoted;
      }
      continue;
    }
    if (character === "," && !quoted) {
      fields.push(value);
      value = "";
      continue;
    }
    value += character;
  }
  fields.push(value);
  return fields;
}

async function inspectDecisionFile(path: string): Promise<DecisionSnapshot> {
  const absolutePath = resolve(repositoryRoot, path);
  const text = await Bun.file(absolutePath).text();
  const lines = text.split(/\r?\n/).filter((line) => line.trim().length > 0);
  const headers = splitCsvLine(lines[0] ?? "");
  const records = lines.slice(1).map((line) => Object.fromEntries(headers.map((header, index) => [header, splitCsvLine(line)[index] ?? ""])));
  const unresolved: DecisionSnapshot["unresolved"] = [];
  const invalid: DecisionSnapshot["invalid"] = [];
  for (const record of records) {
    const decisionId = record.decision_id ?? "<missing>";
    const confidence = Number(record.confidence);
    if (!record.basis) invalid.push({ decisionId, reason: "missing basis" });
    if (!Number.isFinite(confidence) || confidence < 0 || confidence > 1) invalid.push({ decisionId, reason: "confidence must be between 0 and 1" });
    if (!new Set(["pending", "confirmed", "rejected", "implemented"]).has(record.status)) invalid.push({ decisionId, reason: `unsupported status: ${record.status}` });
    if (record.review_required === "yes" && record.status === "pending") {
      unresolved.push({ decisionId, scope: record.scope, objectGuid: record.object_guid, status: record.status });
    }
  }
  return {
    path,
    total: records.length,
    reviewRequired: records.filter((record) => record.review_required === "yes").length,
    unresolved,
    invalid,
  };
}

async function main(): Promise<void> {
  const command = process.argv[2] ?? "help";
  if (!new Set(["snapshot", "check"]).has(command)) {
    console.log("Usage: bun run pipeline/src/cli.ts <snapshot|check> [--report-only]");
    process.exitCode = command === "help" ? 0 : 64;
    return;
  }

  const config = (await Bun.file(configPath).json()) as ProjectConfig;
  const ifcPath = resolve(repositoryRoot, config.ifcPath);
  const snapshot = await snapshotIfc(ifcPath);
  const drawingSnapshots = await Promise.all(config.drawings.map(inspectDrawing));
  const decisionSnapshots = await Promise.all(config.decisionFiles.map(inspectDecisionFile));
  const snapshotArtifact = { projectId: config.projectId, projectName: config.projectName, ifc: snapshot, drawings: drawingSnapshots, decisions: decisionSnapshots };
  const snapshotPath = resolve(repositoryRoot, config.output.snapshot);
  await writeJson(snapshotPath, snapshotArtifact);
  console.log(`snapshot: ${config.output.snapshot}`);
  console.log(`ifc sha256: ${snapshot.source.sha256}`);
  console.log(`entities: ${snapshot.integrity.stepEntityCount}`);

  if (command === "snapshot") return;

  const report = buildQaReport(snapshot, drawingSnapshots, config, decisionSnapshots);
  const qaJsonPath = resolve(repositoryRoot, config.output.qaJson);
  const qaMarkdownPath = resolve(repositoryRoot, config.output.qaMarkdown);
  await writeJson(qaJsonPath, report);
  await mkdir(dirname(qaMarkdownPath), { recursive: true });
  await Bun.write(qaMarkdownPath, renderQaMarkdown(report));
  console.log(`qa json: ${config.output.qaJson}`);
  console.log(`qa markdown: ${config.output.qaMarkdown}`);
  console.log(`pass/warn/block: ${report.summary.pass}/${report.summary.warn}/${report.summary.block}`);

  if (!report.summary.releasable && !process.argv.includes("--report-only")) process.exitCode = 2;
}

if (import.meta.main) await main();
