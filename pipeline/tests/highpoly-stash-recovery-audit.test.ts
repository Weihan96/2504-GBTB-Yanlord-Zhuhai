import { expect, test } from "bun:test";
import { mkdtempSync, readFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("stash high-poly candidates are recovered without applying Libelle changes", () => {
  const temporary = mkdtempSync(join(tmpdir(), "highpoly-stash-audit-"));
  const output = join(temporary, "audit.json");
  const before = Bun.spawnSync([
    "git", "stash", "list", "--format=%H%x09%gd%x09%gs",
  ], { cwd: root }).stdout.toString();
  const run = Bun.spawnSync([
    "python3",
    "pipeline/scripts/audit_highpoly_stash_recovery.py",
    "--output",
    output,
  ], { cwd: root, stdout: "pipe", stderr: "pipe" });
  expect(run.exitCode, run.stderr.toString()).toBe(0);
  const audit = JSON.parse(readFileSync(output, "utf8"));
  expect(audit.pass).toBe(true);
  expect(audit.recovery_strategy).toBe(
    "read_only_extract_untracked_parent_then_supersede_with_audited_packages",
  );
  expect(audit.source_stash.commit).toBe("8620d0cd0ee91f0f6694b98d5099fb518aa65aaa");
  expect(audit.candidate_product_slugs).toEqual(["falper-sorgente", "geberit-146-140"]);
  expect(audit.untracked_snapshot).toHaveLength(14);
  expect(audit.untracked_snapshot.every((record: any) => record.current_exists)).toBe(true);
  expect(audit.untracked_snapshot.find((record: any) =>
    record.path === "pipeline/tests/int1-highpoly-type-review.test.ts"
  )).toMatchObject({
    current_path: "pipeline/tests/falper-sorgente-highpoly-review.test.ts",
    successor_mapping_used: true,
    recovery_state: "recovered_and_superseded_by_audited_package",
  });
  expect(audit.tracked_delta_not_restored.map((record: any) => record.path).sort()).toEqual([
    "output/images/libelle-background-studies/libelle-preview-warm-silver-white.png",
    "output/images/libelle-background-studies/libelle-scene-transparent-background.png",
  ]);
  expect(audit.protected_other_task_stash).toMatchObject({
    ref_observed_from_stash_list_only: "stash@{0}",
    content_inspected: false,
    mutated: false,
  });
  expect(audit.stash_apply_performed).toBe(false);
  expect(audit.stash_pop_performed).toBe(false);
  expect(audit.stash_drop_performed).toBe(false);
  expect(audit.stash_mutation_performed).toBe(false);
  expect(audit.inventory_product_count_after_recovery).toBe(40);
  const after = Bun.spawnSync([
    "git", "stash", "list", "--format=%H%x09%gd%x09%gs",
  ], { cwd: root }).stdout.toString();
  expect(after).toBe(before);
});
