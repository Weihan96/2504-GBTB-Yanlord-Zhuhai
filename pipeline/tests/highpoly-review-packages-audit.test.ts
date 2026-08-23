import { expect, test } from "bun:test";
import { createHash } from "node:crypto";
import { existsSync, mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const formalIfc = join(root, "2504 GBTB Yanlord Zhuhai.ifc");
const formalHash = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c";

function sha256(path: string): string {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

test("all high-poly review packages obey approval, context and commit gates", () => {
  const temporary = mkdtempSync(join(tmpdir(), "highpoly-package-audit-"));
  const output = join(temporary, "queue.json");
  const before = sha256(formalIfc);
  const run = Bun.spawnSync([
    "python3",
    "pipeline/scripts/audit_highpoly_review_packages.py",
    "--output",
    output,
  ], { cwd: root });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const report = JSON.parse(readFileSync(output, "utf8"));
  expect(report.formal_ifc_sha256).toBe(formalHash);
  expect(report.formal_ifc_bytes_unchanged).toBe(true);
  expect(report.stash_mutation_performed).toBe(false);
  expect(report.stashes_observed_read_only).toEqual(expect.arrayContaining([
    expect.stringContaining("Libelle background studies from other task"),
    expect.stringContaining("superseded highpoly review candidates"),
  ]));
  expect(report.summary).toMatchObject({
    product_count: 40,
    completed_count: 1,
    fully_packaged_pending_count: 38,
    pending_human_approval_count: 38,
    excluded_other_task_count: 1,
    derived_ifc_product_count: 1,
    independently_committed_product_count: 39,
    package_mechanical_pass: true,
    stash_recovery_pass: true,
    goal_complete: false,
  });
  expect(report.stash_recovery).toMatchObject({
    source_stash_commit: "8620d0cd0ee91f0f6694b98d5099fb518aa65aaa",
    recovered_file_count: 14,
    candidate_product_slugs: ["falper-sorgente", "geberit-146-140"],
    pass: true,
    protected_other_task_stash: {
      content_inspected: false,
      mutated: false,
    },
  });

  const pending = report.products.filter(
    (product: any) => product.status === "review_ready_pending_approval",
  );
  expect(pending).toHaveLength(38);
  for (const product of pending) {
    expect(product.mechanical_pass, product.slug).toBe(true);
    expect(product.package.missing_files).toEqual([]);
    expect(product.package.contact_sheet).toBeString();
    expect(product.bonsai).toMatchObject({
      actual_ifc_body_camera_render: true,
      saved_camera_count: 4,
      one_product_only: true,
      whole_model_render: false,
      pass: true,
    });
    expect(product.project_context).toMatchObject({
      project_context_retained: true,
      walls_and_surrounding_project_elements_retained: true,
      pass: true,
    });
    expect(product.approval.hash_matches).toBe(true);
    expect(["pending", "visual_review_pending"]).toContain(product.approval.status);
    expect(product.approval.write_allowed).toBe(false);
    expect(product.derived_ifcs).toEqual([]);
    expect(product.commits.length, product.slug).toBeGreaterThanOrEqual(1);
    expect(product.commit_gate_pass, product.slug).toBe(true);
    expect(product.drawing_source_gate.pass, product.slug).toBe(true);
  }

  const officialCad = pending.filter((product: any) => product.official_cad_used);
  expect(officialCad.map((product: any) => product.slug).sort()).toEqual([
    "geberit-115-770",
    "geberit-146-140",
    "geberit-154-446-ks-1",
    "geberit-duofix-sigma-224-212",
    "gessi316-54038",
    "gessi316-54093",
    "gessi316-54145",
    "gessi316-54146",
    "gessi316-54294",
    "marilyn-01",
    "marilyn-02",
  ]);
  expect(officialCad.every((product: any) =>
    product.drawing_source_gate.mode === "official_native_cad"
  )).toBe(true);

  const fallback = pending.filter((product: any) => !product.official_cad_used);
  expect(fallback).toHaveLength(27);
  for (const product of fallback) {
    expect(product.source_kind, product.slug).toBe("geometry_derived_simplified_proxy");
    expect(product.drawing_source_gate).toMatchObject({
      mode: "geometry_derived_fallback",
      source_label_zh: "基于原始高模几何生成的简化图纸表达",
      official_cad_used: false,
      pass: true,
    });
    expect(Object.values(product.drawing_source_gate.checks).every(Boolean), product.slug).toBe(true);
  }

  const falper = report.products.find((product: any) => product.slug === "falper-sorgente");
  expect(falper).toMatchObject({
    status: "completed",
    source_kind: "native_dwg",
    official_cad_used: true,
    drawing_source_gate: {
      mode: "approved_official_native_cad",
      pass: true,
    },
    mechanical_pass: true,
  });
  expect(falper.project_context).toMatchObject({
    view_count: 2,
    walls_and_surrounding_project_elements_retained: true,
    pass: true,
  });
  expect(falper.approval).toMatchObject({ status: "approved", hash_matches: true });
  expect(falper.derived_ifcs.length).toBeGreaterThanOrEqual(1);
  expect(falper.commits).toContain(
    "567dfb8e61f390fec5ea49691b8f4dee110f59b0",
  );

  const libelle = report.products.find((product: any) => product.slug === "libelle");
  expect(libelle).toMatchObject({
    status: "excluded_other_task_do_not_touch",
    excluded: true,
    do_not_touch: true,
    mechanical_pass: true,
  });
  expect(existsSync(output.replace(/\.json$/, ".html"))).toBe(true);
  expect(sha256(formalIfc)).toBe(before);
}, 30_000);
