import { expect, test } from "bun:test";
import { mkdtempSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const script = resolve(root, "pipeline/scripts/pm_progress_check.py");

function fixture(summaryDone: number) {
  const directory = mkdtempSync(join(tmpdir(), "pm-progress-"));
  const path = join(directory, "pm.md");
  writeFileSync(path, `
| PM 进度（派生） | 2 项中 ${summaryDone} 项完成（50%） |

\`\`\`mermaid
gantt
    已完成 :done, A001, 2026-08-04, 1d
    进行中 :active, A002, 2026-08-05, 1d
\`\`\`

## 5. PM 任务台账

- [x] **A001｜已完成**
- [ ] **A002｜进行中**

### 任务台账操作规则
`);
  return path;
}

function run(path: string) {
  return Bun.spawnSync(["python3", script, "--pm", path], {
    cwd: root,
    stdout: "pipe",
    stderr: "pipe",
  });
}

test("PM summary and Gantt match the authoritative ledger", () => {
  const result = run(fixture(1));
  expect(result.exitCode, result.stderr.toString()).toBe(0);
  expect(JSON.parse(result.stdout.toString()).summary).toEqual({
    total: 2,
    done: 1,
    open: 1,
    percent: 50,
  });
});

test("stale PM summary fails closed", () => {
  const result = run(fixture(0));
  expect(result.exitCode).toBe(1);
  expect(JSON.parse(result.stdout.toString()).errors[0]).toContain("does not match ledger");
});
