import type { IfcSnapshot } from "./ifc-step";

export type GateStatus = "pass" | "warn" | "block";

export interface DrawingSnapshot {
  id: string;
  path: string;
  role: string;
  exists: boolean;
  byteSize: number;
  scale: string | null;
  dimensionLineCount: number;
  textCount: number;
}

export interface DecisionSnapshot {
  path: string;
  total: number;
  reviewRequired: number;
  unresolved: Array<{ decisionId: string; scope: string; objectGuid: string; status: string }>;
  invalid: Array<{ decisionId: string; reason: string }>;
}

export interface GateResult {
  id: string;
  status: GateStatus;
  summary: string;
  evidence: Record<string, unknown>;
}

export interface QaReport {
  generatedAt: string;
  sourceSha256: string;
  summary: { pass: number; warn: number; block: number; releasable: boolean };
  gates: GateResult[];
}

interface ProjectConfig {
  ifcSchema: string;
  baselineCounts: Record<string, number>;
}

function gate(id: string, status: GateStatus, summary: string, evidence: Record<string, unknown> = {}): GateResult {
  return { id, status, summary, evidence };
}

export function buildQaReport(snapshot: IfcSnapshot, drawings: DrawingSnapshot[], config: ProjectConfig, decisions: DecisionSnapshot[] = []): QaReport {
  const gates: GateResult[] = [];
  gates.push(
    gate(
      "SOURCE-SCHEMA",
      snapshot.schema === config.ifcSchema ? "pass" : "block",
      snapshot.schema === config.ifcSchema ? `IFC schema is ${snapshot.schema}.` : `Expected ${config.ifcSchema}, found ${snapshot.schema ?? "unknown"}.`,
      { expected: config.ifcSchema, actual: snapshot.schema },
    ),
  );

  const integrityOk =
    snapshot.integrity.duplicateStepIds.length === 0 &&
    snapshot.integrity.duplicateGlobalIds.length === 0 &&
    snapshot.integrity.hasValidFooter;
  gates.push(
    gate("IFC-INTEGRITY", integrityOk ? "pass" : "block", integrityOk ? "STEP ids, GlobalIds, and footer passed mechanical checks." : "IFC mechanical integrity failed.", snapshot.integrity),
  );

  const countDrift = Object.entries(config.baselineCounts)
    .filter(([type, expected]) => (snapshot.entityCounts[type] ?? 0) !== expected)
    .map(([type, expected]) => ({ type, expected, actual: snapshot.entityCounts[type] ?? 0 }));
  gates.push(gate("BASELINE-COUNTS", countDrift.length === 0 ? "pass" : "warn", countDrift.length === 0 ? "Critical entity counts match the recorded baseline." : "Critical entity counts drifted from the recorded baseline.", { drift: countDrift }));

  const spaceIdentityOk = snapshot.spaces.missingName === 0 && snapshot.spaces.missingLongName === 0;
  gates.push(gate("SPACE-IDENTITY", spaceIdentityOk ? "pass" : "block", spaceIdentityOk ? "All spaces have Name and LongName." : "Some spaces are missing Name or LongName.", snapshot.spaces));
  gates.push(gate("SPACE-REVIEW", snapshot.spaces.provisional === 0 && snapshot.spaces.pendingDelete === 0 ? "pass" : "block", snapshot.spaces.provisional === 0 && snapshot.spaces.pendingDelete === 0 ? "All space semantics are release-ready." : "Provisional or pending-delete spaces remain.", { provisional: snapshot.spaces.provisional, pendingDelete: snapshot.spaces.pendingDelete }));
  gates.push(gate("SPACE-STOREY-RELATION", snapshot.spaces.aggregatedUnderStorey === snapshot.spaces.total ? "pass" : "block", snapshot.spaces.aggregatedUnderStorey === snapshot.spaces.total ? "Every space is aggregated under a storey." : "Spaces are not represented through the expected storey aggregation relationship.", { total: snapshot.spaces.total, contained: snapshot.spaces.containedInStorey, aggregated: snapshot.spaces.aggregatedUnderStorey }));

  gates.push(gate("DOOR-LOCATION-DATA", snapshot.doors.missingTag === 0 && snapshot.doors.missingOverallSize === 0 && snapshot.doors.missingOpeningFill === 0 ? "pass" : "block", "Door numbering, overall size, and opening-fill readiness for A-104.", { total: snapshot.doors.total, missingTag: snapshot.doors.missingTag, missingOverallSize: snapshot.doors.missingOverallSize, missingOpeningFill: snapshot.doors.missingOpeningFill }));
  gates.push(gate("WINDOW-LOCATION-DATA", snapshot.windows.missingTag === 0 && snapshot.windows.missingOverallSize === 0 && snapshot.windows.missingOpeningFill === 0 ? "pass" : "block", "Window numbering, overall size, and opening-fill readiness for A-104.", { total: snapshot.windows.total, missingTag: snapshot.windows.missingTag, missingOverallSize: snapshot.windows.missingOverallSize, missingOpeningFill: snapshot.windows.missingOpeningFill }));

  const formalDimensions = snapshot.annotations.byObjectType.DIMENSION ?? 0;
  gates.push(gate("FORMAL-DIMENSIONS", formalDimensions > 0 ? "pass" : "block", formalDimensions > 0 ? `${formalDimensions} formal IFC dimension annotations found.` : "No formal IFC dimension annotations found.", { formalDimensions }));

  const missingDrawings = drawings.filter((drawing) => !drawing.exists);
  const wrongScale = drawings.filter((drawing) => drawing.exists && drawing.scale !== "1:50");
  gates.push(gate("DRAWING-SOURCES", missingDrawings.length === 0 && wrongScale.length === 0 ? "pass" : "block", missingDrawings.length === 0 && wrongScale.length === 0 ? "Configured drawing sources exist at 1:50." : "Drawing source files or scales are invalid.", { missing: missingDrawings.map((drawing) => drawing.path), wrongScale: wrongScale.map((drawing) => ({ path: drawing.path, scale: drawing.scale })) }));

  const candidate = drawings.find((drawing) => drawing.id === "A-103-CANDIDATE");
  gates.push(gate("A103-CANDIDATE", candidate && candidate.dimensionLineCount > 0 ? "warn" : "block", candidate && candidate.dimensionLineCount > 0 ? `A-103 contains ${candidate.dimensionLineCount} SVG candidate dimension lines; they remain non-IFC and require site verification.` : "A-103 candidate dimensions are missing.", { candidate }));

  const unresolvedDecisions = decisions.flatMap((decision) => decision.unresolved);
  const invalidDecisions = decisions.flatMap((decision) => decision.invalid);
  const decisionStatus: GateStatus = invalidDecisions.length > 0 || unresolvedDecisions.length > 0 ? "block" : "pass";
  gates.push(
    gate(
      "DECISION-QUEUE",
      decisionStatus,
      decisionStatus === "pass" ? "All required human reviews are resolved and decision records are valid." : `${unresolvedDecisions.length} required reviews remain unresolved; ${invalidDecisions.length} records are invalid.`,
      { files: decisions, unresolved: unresolvedDecisions, invalid: invalidDecisions },
    ),
  );

  const summary = {
    pass: gates.filter((item) => item.status === "pass").length,
    warn: gates.filter((item) => item.status === "warn").length,
    block: gates.filter((item) => item.status === "block").length,
    releasable: !gates.some((item) => item.status === "block"),
  };
  return { generatedAt: new Date().toISOString(), sourceSha256: snapshot.source.sha256, summary, gates };
}

export function renderQaMarkdown(report: QaReport): string {
  const lines = [
    "# IFC → 施工图机械 QA",
    "",
    `- 生成时间：${report.generatedAt}`,
    `- IFC SHA-256：\`${report.sourceSha256}\``,
    `- 结果：${report.summary.releasable ? "可发布" : "不可发布"}`,
    `- PASS/WARN/BLOCK：${report.summary.pass}/${report.summary.warn}/${report.summary.block}`,
    "",
    "| Gate | 状态 | 结论 |",
    "| --- | --- | --- |",
  ];
  for (const item of report.gates) {
    lines.push(`| ${item.id} | ${item.status.toUpperCase()} | ${item.summary.replaceAll("|", "\\|")} |`);
  }
  lines.push("", "## 证据", "");
  for (const item of report.gates) {
    lines.push(`### ${item.id}`, "", "```json", JSON.stringify(item.evidence, null, 2), "```", "");
  }
  return `${lines.join("\n")}\n`;
}
