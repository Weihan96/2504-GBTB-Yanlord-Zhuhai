import {
  accessSync,
  constants,
  existsSync,
  lstatSync,
  readdirSync,
  readFileSync,
  realpathSync,
} from "node:fs";
import { isAbsolute, relative, resolve, sep } from "node:path";

const SOURCE_REGISTER = "pipeline/decisions/source-evidence-register.csv";
const DECISIONS_DIR = "pipeline/decisions";
const EVIDENCE_DIR = "drawings/evidence";

const SOURCE_FIELDS = new Set([
  "local_path",
  "source_document",
  "source_url",
  "evidence",
  "evidence_reference",
  "source_locator",
  "source_file",
  "source_basis",
  "source_references",
  "official_reference",
  "current_evidence",
  "report_path",
  "source_view",
]);

const SHA256_PATTERN = /^[0-9a-f]{64}$/;
const HTTP_URL_PATTERN = /^https?:\/\//i;
const FILE_URL_PATTERN = /file:\/\//i;
const EMBEDDED_ABSOLUTE_PATH_PATTERN =
  /(?:^|[\s"'(=;,\[])(?:~\/|\/(?:Users|home|Volumes|private|var|tmp|opt|mnt|etc|Applications|Library|System)(?:\/|$)|[A-Za-z]:[\\/])/;
const TRANSIENT_PATH_PATTERN =
  /(?:^|[\\/])(?:Downloads|Desktop|TemporaryItems|tmp)(?:[\\/]|$)/i;
const SECONDARY_OUTPUT_PATTERN = /(?:^|[\\/])outputs(?:[\\/]|$)/i;

export type PortabilityIssue = {
  code: string;
  file: string;
  row_number?: number;
  row_id?: string;
  field?: string;
  value?: string;
  message: string;
};

export type PortabilityReport = {
  mode: "read_only_evidence_portability_scan";
  root: string;
  ok: boolean;
  canonical_data: string;
  summary: {
    csv_files_scanned: number;
    source_records: number;
    local_evidence_checked: number;
    evidence_files_found: number;
    registered_evidence_files: number;
    allowlisted_unregistered_files: number;
    error_count: number;
  };
  allowlisted_unregistered: Array<{ path: string; reason: string }>;
  errors: PortabilityIssue[];
};

type CsvRecord = {
  rowNumber: number;
  values: Record<string, string>;
};

function normalizeRepoPath(path: string): string {
  return path.replaceAll("\\", "/").replace(/^\.\//, "");
}

function isInsideRoot(root: string, candidate: string): boolean {
  const rel = relative(root, candidate);
  return rel === "" || (!rel.startsWith(`..${sep}`) && rel !== ".." && !isAbsolute(rel));
}

function fileSha256(path: string): string {
  const hasher = new Bun.CryptoHasher("sha256");
  hasher.update(readFileSync(path));
  return hasher.digest("hex");
}

function parseCsv(text: string): CsvRecord[] {
  const input = text.replace(/^\uFEFF/, "");
  const rows: Array<{ line: number; cells: string[] }> = [];
  let cells: string[] = [];
  let cell = "";
  let quoted = false;
  let line = 1;
  let rowStart = 1;

  const pushRow = () => {
    cells.push(cell);
    if (cells.some((value) => value.length > 0)) rows.push({ line: rowStart, cells });
    cells = [];
    cell = "";
    rowStart = line + 1;
  };

  for (let index = 0; index < input.length; index += 1) {
    const char = input[index];
    if (quoted) {
      if (char === '"') {
        if (input[index + 1] === '"') {
          cell += '"';
          index += 1;
        } else {
          quoted = false;
        }
      } else {
        cell += char;
        if (char === "\n") line += 1;
      }
      continue;
    }
    if (char === '"') {
      quoted = true;
    } else if (char === ",") {
      cells.push(cell);
      cell = "";
    } else if (char === "\n") {
      pushRow();
      line += 1;
    } else if (char !== "\r") {
      cell += char;
    }
  }
  if (cell.length > 0 || cells.length > 0) pushRow();
  if (quoted) throw new Error("unterminated quoted CSV field");
  if (rows.length === 0) return [];

  const headers = rows[0].cells;
  for (const row of rows.slice(1)) {
    if (row.cells.length !== headers.length) {
      throw new Error(
        `row ${row.line} has ${row.cells.length} fields; expected ${headers.length}`,
      );
    }
  }
  return rows.slice(1).map(({ line: rowNumber, cells: rowCells }) => ({
    rowNumber,
    values: Object.fromEntries(headers.map((header, index) => [header, rowCells[index] ?? ""])),
  }));
}

function listFilesRecursively(directory: string): string[] {
  if (!existsSync(directory)) return [];
  const result: string[] = [];
  for (const entry of readdirSync(directory, { withFileTypes: true })) {
    const path = resolve(directory, entry.name);
    if (entry.isDirectory()) result.push(...listFilesRecursively(path));
    else if (entry.isFile() || entry.isSymbolicLink()) result.push(path);
  }
  return result.sort();
}

function pathPortabilityCodes(value: string): string[] {
  const trimmed = value.trim();
  if (!trimmed) return [];
  const withoutHttpUrls = trimmed.replace(/https?:\/\/\S+/gi, "");
  const codes: string[] = [];
  if (FILE_URL_PATTERN.test(trimmed)) codes.push("file_uri_reference");
  if (isAbsolute(trimmed) || EMBEDDED_ABSOLUTE_PATH_PATTERN.test(trimmed)) {
    codes.push("absolute_external_reference");
  }
  if (TRANSIENT_PATH_PATTERN.test(withoutHttpUrls)) codes.push("transient_local_reference");
  if (SECONDARY_OUTPUT_PATTERN.test(withoutHttpUrls)) codes.push("secondary_output_reference");
  return [...new Set(codes)];
}

function isSourceField(field: string): boolean {
  return SOURCE_FIELDS.has(field)
    || /(?:^|_)(?:source|evidence|reference|document|locator|path)(?:_|$)/i.test(field);
}

function issueMessage(code: string): string {
  const messages: Record<string, string> = {
    file_uri_reference: "file:// reference is machine-local and not portable",
    absolute_external_reference: "absolute local path is outside the portable repository contract",
    transient_local_reference: "Downloads/Desktop/TemporaryItems/tmp reference is transient and must be archived",
    secondary_output_reference: "outputs/ is not an allowed source or output root; use output/",
  };
  return messages[code] ?? code;
}

function addIssue(errors: PortabilityIssue[], issue: PortabilityIssue): void {
  const key = [issue.code, issue.file, issue.row_number, issue.row_id, issue.field, issue.value].join("\u0000");
  if (!errors.some((item) =>
    [item.code, item.file, item.row_number, item.row_id, item.field, item.value].join("\u0000") === key
  )) errors.push(issue);
}

export function scanEvidencePortability(options: {
  root: string;
  allowUnregistered?: string[];
}): PortabilityReport {
  const root = realpathSync(options.root);
  const errors: PortabilityIssue[] = [];
  const decisionsPath = resolve(root, DECISIONS_DIR);
  const csvFiles = existsSync(decisionsPath)
    ? readdirSync(decisionsPath)
      .filter((name) => name.endsWith(".csv"))
      .map((name) => resolve(decisionsPath, name))
      .sort()
    : [];

  for (const csvPath of csvFiles) {
    const csvRelative = normalizeRepoPath(relative(root, csvPath));
    let records: CsvRecord[];
    try {
      records = parseCsv(readFileSync(csvPath, "utf8"));
    } catch (error) {
      addIssue(errors, {
        code: "csv_parse_error",
        file: csvRelative,
        message: error instanceof Error ? error.message : String(error),
      });
      continue;
    }
    for (const record of records) {
      const rowId = record.values.source_id
        || record.values.evidence_id
        || record.values.input_id
        || record.values.appliance_id
        || record.values.equipment_id
        || "";
      for (const [field, value] of Object.entries(record.values)) {
        if (!isSourceField(field) || !value.trim()) continue;
        for (const code of pathPortabilityCodes(value)) {
          addIssue(errors, {
            code,
            file: csvRelative,
            row_number: record.rowNumber,
            row_id: rowId,
            field,
            value,
            message: issueMessage(code),
          });
        }
      }
    }
  }

  const registerPath = resolve(root, SOURCE_REGISTER);
  let sourceRecords: CsvRecord[] = [];
  if (!existsSync(registerPath)) {
    addIssue(errors, {
      code: "missing_source_register",
      file: SOURCE_REGISTER,
      message: "canonical source-evidence register is missing",
    });
  } else {
    try {
      sourceRecords = parseCsv(readFileSync(registerPath, "utf8"));
    } catch (error) {
      addIssue(errors, {
        code: "csv_parse_error",
        file: SOURCE_REGISTER,
        message: error instanceof Error ? error.message : String(error),
      });
    }
  }

  const registeredEvidence = new Set<string>();
  let localEvidenceChecked = 0;
  for (const record of sourceRecords) {
    const row = record.values;
    const sourceId = row.source_id || "";
    const sourceKind = (row.source_kind || "").trim();
    const localPath = (row.local_path || "").trim();
    const document = (row.source_document || "").trim();
    const sourceUrl = (row.source_url || "").trim();
    const expectedHash = (row.sha256 || "").trim();
    const remoteOnly = !localPath && (HTTP_URL_PATTERN.test(sourceUrl) || HTTP_URL_PATTERN.test(document));
    const conversationProvenance = !localPath
      && sourceKind === "user_confirmation"
      && expectedHash === "not_applicable_user_confirmation"
      && /^Codex task [0-9a-f-]{36} \(\d{4}-\d{2}-\d{2}\)$/.test(document);

    if (!SHA256_PATTERN.test(expectedHash) && !expectedHash.startsWith("not_applicable_")) {
      addIssue(errors, {
        code: "invalid_sha256",
        file: SOURCE_REGISTER,
        row_number: record.rowNumber,
        row_id: sourceId,
        field: "sha256",
        value: expectedHash,
        message: "sha256 must be 64 lowercase hexadecimal characters or an explicit not_applicable_* marker",
      });
    }

    if (!localPath) {
      if (!remoteOnly && !conversationProvenance) {
        addIssue(errors, {
          code: "missing_project_copy",
          file: SOURCE_REGISTER,
          row_number: record.rowNumber,
          row_id: sourceId,
          field: "local_path",
          value: document,
          message: "local evidence has no repository-relative archived copy",
        });
      }
      continue;
    }

    if (HTTP_URL_PATTERN.test(localPath) || FILE_URL_PATTERN.test(localPath)) {
      addIssue(errors, {
        code: "invalid_local_path_scheme",
        file: SOURCE_REGISTER,
        row_number: record.rowNumber,
        row_id: sourceId,
        field: "local_path",
        value: localPath,
        message: "local_path must be a repository-relative filesystem path",
      });
      continue;
    }
    if (isAbsolute(localPath)) {
      addIssue(errors, {
        code: "nonportable_local_path",
        file: SOURCE_REGISTER,
        row_number: record.rowNumber,
        row_id: sourceId,
        field: "local_path",
        value: localPath,
        message: "local_path must not be absolute",
      });
      continue;
    }

    const normalized = normalizeRepoPath(localPath);
    const diskPath = resolve(root, normalized);
    if (!isInsideRoot(root, diskPath)) {
      addIssue(errors, {
        code: "local_path_escapes_repo",
        file: SOURCE_REGISTER,
        row_number: record.rowNumber,
        row_id: sourceId,
        field: "local_path",
        value: localPath,
        message: "local_path resolves outside the repository",
      });
      continue;
    }
    if (!existsSync(diskPath) || !lstatSync(diskPath).isFile() && !lstatSync(diskPath).isSymbolicLink()) {
      addIssue(errors, {
        code: "missing_local_evidence",
        file: SOURCE_REGISTER,
        row_number: record.rowNumber,
        row_id: sourceId,
        field: "local_path",
        value: localPath,
        message: "registered local evidence file is missing",
      });
      continue;
    }
    const realDiskPath = realpathSync(diskPath);
    if (!isInsideRoot(root, realDiskPath)) {
      addIssue(errors, {
        code: "local_path_symlink_escapes_repo",
        file: SOURCE_REGISTER,
        row_number: record.rowNumber,
        row_id: sourceId,
        field: "local_path",
        value: localPath,
        message: "registered evidence symlink resolves outside the repository",
      });
      continue;
    }
    try {
      accessSync(realDiskPath, constants.R_OK);
    } catch {
      addIssue(errors, {
        code: "unreadable_local_evidence",
        file: SOURCE_REGISTER,
        row_number: record.rowNumber,
        row_id: sourceId,
        field: "local_path",
        value: localPath,
        message: "registered local evidence is not readable",
      });
      continue;
    }
    localEvidenceChecked += 1;
    if (normalized === EVIDENCE_DIR || normalized.startsWith(`${EVIDENCE_DIR}/`)) {
      registeredEvidence.add(normalized);
    }
    if (!SHA256_PATTERN.test(expectedHash)) {
      addIssue(errors, {
        code: "local_evidence_requires_sha256",
        file: SOURCE_REGISTER,
        row_number: record.rowNumber,
        row_id: sourceId,
        field: "sha256",
        value: expectedHash,
        message: "repository-local evidence requires an exact lowercase SHA-256",
      });
    } else {
      const actualHash = fileSha256(realDiskPath);
      if (actualHash !== expectedHash) {
        addIssue(errors, {
          code: "hash_mismatch",
          file: SOURCE_REGISTER,
          row_number: record.rowNumber,
          row_id: sourceId,
          field: "sha256",
          value: expectedHash,
          message: `registered SHA-256 does not match disk (${actualHash})`,
        });
      }
    }
  }

  const evidenceRoot = resolve(root, EVIDENCE_DIR);
  const evidenceFiles = listFilesRecursively(evidenceRoot)
    .map((path) => normalizeRepoPath(relative(root, path)));
  const allowlisted = new Map<string, string>();
  for (const entry of options.allowUnregistered || []) {
    const separator = entry.indexOf("::");
    const allowlistedPath = normalizeRepoPath(separator < 0 ? entry : entry.slice(0, separator));
    const reason = separator < 0 ? "" : entry.slice(separator + 2).trim();
    if (!allowlistedPath.startsWith(`${EVIDENCE_DIR}/`) || !reason) {
      addIssue(errors, {
        code: "invalid_unregistered_allowlist",
        file: allowlistedPath,
        message: `--allow-unregistered must use ${EVIDENCE_DIR}/<path>::<reason>`,
      });
    } else {
      allowlisted.set(allowlistedPath, reason);
    }
  }
  for (const evidencePath of evidenceFiles) {
    if (registeredEvidence.has(evidencePath) || allowlisted.has(evidencePath)) continue;
    addIssue(errors, {
      code: "unregistered_evidence",
      file: evidencePath,
      message: "evidence file exists in the repository but is absent from the canonical source register",
    });
  }

  const pluralOutput = resolve(root, "outputs");
  if (existsSync(pluralOutput)) {
    addIssue(errors, {
      code: "plural_output_root",
      file: "outputs/",
      message: "repository root outputs/ is forbidden; output/ is the only formal output root",
    });
  }

  errors.sort((a, b) =>
    a.file.localeCompare(b.file)
    || (a.row_number || 0) - (b.row_number || 0)
    || a.code.localeCompare(b.code)
  );
  const allowlistedEntries = evidenceFiles
    .filter((path) => allowlisted.has(path))
    .map((path) => ({ path, reason: allowlisted.get(path)! }));
  return {
    mode: "read_only_evidence_portability_scan",
    root,
    ok: errors.length === 0,
    canonical_data: SOURCE_REGISTER,
    summary: {
      csv_files_scanned: csvFiles.length,
      source_records: sourceRecords.length,
      local_evidence_checked: localEvidenceChecked,
      evidence_files_found: evidenceFiles.length,
      registered_evidence_files: registeredEvidence.size,
      allowlisted_unregistered_files: allowlistedEntries.length,
      error_count: errors.length,
    },
    allowlisted_unregistered: allowlistedEntries,
    errors,
  };
}

export function renderPortabilityReport(report: PortabilityReport): string {
  const lines = [
    `Evidence portability: ${report.ok ? "PASS" : "FAIL"}`,
    `Canonical source register: ${report.canonical_data}`,
    `Scanned ${report.summary.csv_files_scanned} CSV files, ${report.summary.source_records} source records, and ${report.summary.evidence_files_found} evidence files.`,
    `Checked ${report.summary.local_evidence_checked} local evidence hashes; errors=${report.summary.error_count}.`,
  ];
  for (const issue of report.errors) {
    const location = [issue.file, issue.row_number ? `row ${issue.row_number}` : "", issue.row_id || "", issue.field || ""]
      .filter(Boolean)
      .join(":");
    lines.push(`ERROR [${issue.code}] ${location} — ${issue.message}`);
  }
  return `${lines.join("\n")}\n`;
}

function parseArgs(args: string[]): { root: string; json: boolean; allowUnregistered: string[] } {
  let root = process.cwd();
  let json = false;
  const allowUnregistered: string[] = [];
  for (let index = 0; index < args.length; index += 1) {
    const arg = args[index];
    if (arg === "--json") json = true;
    else if (arg === "--root") root = args[++index] || "";
    else if (arg.startsWith("--root=")) root = arg.slice("--root=".length);
    else if (arg === "--allow-unregistered") allowUnregistered.push(args[++index] || "");
    else if (arg.startsWith("--allow-unregistered=")) allowUnregistered.push(arg.slice("--allow-unregistered=".length));
    else throw new Error(`unknown argument: ${arg}`);
  }
  if (!root) throw new Error("--root requires a path");
  return { root, json, allowUnregistered };
}

if (import.meta.main) {
  try {
    const options = parseArgs(process.argv.slice(2));
    const report = scanEvidencePortability(options);
    process.stdout.write(options.json ? `${JSON.stringify(report, null, 2)}\n` : renderPortabilityReport(report));
    process.exitCode = report.ok ? 0 : 1;
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    process.stderr.write(`Evidence portability scan could not run: ${message}\n`);
    process.exitCode = 2;
  }
}
