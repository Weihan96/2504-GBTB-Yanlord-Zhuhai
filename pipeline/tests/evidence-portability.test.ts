import { afterEach, expect, test } from "bun:test";
import {
  mkdirSync,
  mkdtempSync,
  readdirSync,
  readFileSync,
  rmSync,
  writeFileSync,
} from "node:fs";
import { tmpdir } from "node:os";
import { dirname, resolve } from "node:path";
import { scanEvidencePortability } from "../src/evidence-portability";

const script = resolve(import.meta.dir, "../src/evidence-portability.ts");
const temporaryRoots: string[] = [];

afterEach(() => {
  while (temporaryRoots.length) rmSync(temporaryRoots.pop()!, { recursive: true, force: true });
});

function sha256(value: string | Buffer): string {
  const hasher = new Bun.CryptoHasher("sha256");
  hasher.update(value);
  return hasher.digest("hex");
}

function csvCell(value: string): string {
  return /[",\n]/.test(value) ? `"${value.replaceAll('"', '""')}"` : value;
}

function write(path: string, value: string): void {
  mkdirSync(dirname(path), { recursive: true });
  writeFileSync(path, value);
}

function makeFixture(options: {
  localPath?: string;
  expectedHash?: string;
  extraEvidence?: boolean;
  outputs?: boolean;
  sourceDocument?: string;
  sourceKind?: string;
  sourceUrl?: string;
} = {}): string {
  const root = mkdtempSync(resolve(tmpdir(), "evidence-portability-"));
  temporaryRoots.push(root);
  const localPath = options.localPath ?? "drawings/evidence/manual.pdf";
  const content = "official evidence\n";
  if (localPath && !localPath.startsWith("/") && !localPath.startsWith("file://")) {
    write(resolve(root, localPath), content);
  }
  if (options.extraEvidence) write(resolve(root, "drawings/evidence/unregistered.jpg"), "photo\n");
  if (options.outputs) write(resolve(root, "outputs/task-id/old.xlsx"), "old export\n");

  const header = [
    "source_id", "source_kind", "source_document", "source_url", "local_path", "sha256", "locator",
  ];
  const row = [
    "SRC-001",
    options.sourceKind ?? "official_product_pdf",
    options.sourceDocument ?? "manual.pdf",
    options.sourceUrl ?? "https://manufacturer.example/manual",
    localPath,
    options.expectedHash ?? sha256(content),
    "official manual page 1",
  ];
  const remote = [
    "SRC-URL-001",
    "official_product_web",
    "https://manufacturer.example/product",
    "https://manufacturer.example/product",
    "",
    "not_applicable_live_official_web",
    "official product page",
  ];
  write(
    resolve(root, "pipeline/decisions/source-evidence-register.csv"),
    `${header.join(",")}\n${row.map(csvCell).join(",")}\n${remote.map(csvCell).join(",")}\n`,
  );
  write(
    resolve(root, "pipeline/decisions/owner-input-register.csv"),
    "input_id,evidence_reference,source_basis\nIN-001,https://example.com/item,owner input\n",
  );
  return root;
}

function snapshotTree(root: string): Record<string, string> {
  const result: Record<string, string> = {};
  const walk = (directory: string) => {
    for (const entry of readdirSync(directory, { withFileTypes: true })) {
      const path = resolve(directory, entry.name);
      const relativePath = path.slice(root.length + 1);
      if (entry.isDirectory()) {
        result[`${relativePath}/`] = "directory";
        walk(path);
      } else {
        result[relativePath] = sha256(readFileSync(path));
      }
    }
  };
  walk(root);
  return result;
}

test("healthy portable evidence repository passes without writing files", () => {
  const root = makeFixture();
  const before = snapshotTree(root);
  const run = Bun.spawnSync(["bun", "run", script, "--root", root, "--json"]);
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const report = JSON.parse(run.stdout.toString());
  expect(report).toMatchObject({
    mode: "read_only_evidence_portability_scan",
    ok: true,
    summary: { source_records: 2, local_evidence_checked: 1, error_count: 0 },
  });
  expect(snapshotTree(root)).toEqual(before);
});

test("absolute, file URI, and transient source references fail portability", () => {
  const root = makeFixture({
    localPath: "/Users/example/Downloads/manual.pdf",
    sourceDocument: "file:///private/var/folders/TemporaryItems/manual.pdf",
  });
  const report = scanEvidencePortability({ root });
  expect(report.ok).toBe(false);
  const codes = new Set(report.errors.map((error) => error.code));
  expect(codes).toContain("nonportable_local_path");
  expect(codes).toContain("absolute_external_reference");
  expect(codes).toContain("file_uri_reference");
  expect(codes).toContain("transient_local_reference");
});

test("registered hash drift fails", () => {
  const root = makeFixture({ expectedHash: "0".repeat(64) });
  const report = scanEvidencePortability({ root });
  expect(report.ok).toBe(false);
  expect(report.errors.map((error) => error.code)).toContain("hash_mismatch");
});

test("unregistered evidence fails by default", () => {
  const root = makeFixture({ extraEvidence: true });
  const report = scanEvidencePortability({ root });
  expect(report.ok).toBe(false);
  expect(report.errors).toContainEqual(expect.objectContaining({
    code: "unregistered_evidence",
    file: "drawings/evidence/unregistered.jpg",
  }));
});

test("an explicit unregistered allowlist requires and reports a reason", () => {
  const root = makeFixture({ extraEvidence: true });
  const invalid = scanEvidencePortability({
    root,
    allowUnregistered: ["drawings/evidence/unregistered.jpg"],
  });
  expect(invalid.ok).toBe(false);
  expect(invalid.errors.map((error) => error.code)).toContain("invalid_unregistered_allowlist");

  const allowed = scanEvidencePortability({
    root,
    allowUnregistered: ["drawings/evidence/unregistered.jpg::fixture-only derived preview"],
  });
  expect(allowed.ok).toBe(true);
  expect(allowed.allowlisted_unregistered).toEqual([{
    path: "drawings/evidence/unregistered.jpg",
    reason: "fixture-only derived preview",
  }]);
});

test("plural outputs root fails", () => {
  const root = makeFixture({ outputs: true });
  const report = scanEvidencePortability({ root });
  expect(report.ok).toBe(false);
  expect(report.errors).toContainEqual(expect.objectContaining({
    code: "plural_output_root",
    file: "outputs/",
  }));
});

test("thread-addressed owner confirmation is valid non-file provenance", () => {
  const root = makeFixture({
    localPath: "",
    expectedHash: "not_applicable_user_confirmation",
    sourceKind: "user_confirmation",
    sourceDocument: "Codex task 019fea35-19f8-7a92-916c-575d2cd8aaf6 (2026-08-13)",
    sourceUrl: "",
  });
  const report = scanEvidencePortability({ root });
  expect(report.ok, JSON.stringify(report.errors)).toBe(true);
  expect(report.summary.local_evidence_checked).toBe(0);
});

test("CSV rows with extra or missing fields fail instead of being truncated", () => {
  const root = makeFixture();
  const register = resolve(root, "pipeline/decisions/source-evidence-register.csv");
  writeFileSync(register, `${readFileSync(register, "utf8")}BROKEN,ROW,WITH,TOO,MANY,FIELDS,EXTRA,CELL\n`);
  const report = scanEvidencePortability({ root });
  expect(report.ok).toBe(false);
  expect(report.errors).toContainEqual(expect.objectContaining({
    code: "csv_parse_error",
    file: "pipeline/decisions/source-evidence-register.csv",
    message: expect.stringContaining("fields; expected"),
  }));
});
