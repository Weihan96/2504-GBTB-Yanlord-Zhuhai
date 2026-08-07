import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");

test("SOC04 evidence identifies the existing mounting-face centre", () => {
  const source = `
import json,sys,numpy as np
sys.path.insert(0,"pipeline/scripts")
from socket_anchor_candidate import socket_evidence
r=socket_evidence(
  "SOC04",
  np.array([3129.773378,-2289.012909,1691.785216]),
  np.array([3119.773378,-2339.012909,1641.785216]),
  np.array([3129.773378,-2239.012909,1741.785216]),
)
print(json.dumps(r))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.pass).toBe(true);
  expect(report.target_mm).toEqual([3130, -2289, 1692]);
});

test("socket evidence rejects an unrelated appliance type", () => {
  const source = `
import sys,numpy as np
sys.path.insert(0,"pipeline/scripts")
from socket_anchor_candidate import socket_evidence
r=socket_evidence("OV01",np.zeros(3),np.zeros(3),np.array([10,100,100]))
print(r["pass"])
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(result.stdout.toString().trim()).toBe("False");
});

test("socket contact exceptions and allowed closures are exact sets", () => {
  const source = `
import json,sys
sys.path.insert(0,"pipeline/scripts")
from socket_anchor_candidate import CONTROLLED_EXCEPTION_IDS, ALLOWED_NEW_CONTEXT_CONTACTS
print(json.dumps({"exceptions": sorted(CONTROLLED_EXCEPTION_IDS), "contacts": sorted(ALLOWED_NEW_CONTEXT_CONTACTS)}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const report = JSON.parse(result.stdout.toString());
  expect(report.exceptions).toEqual([
    "0laejMoxn8Lu_X3FZaCXmi",
    "27MTenki57DQsfMryX_1U0",
    "2OOjqQDMHDjRcXQCniWXnp",
    "3KXtmVvejA78j_iAS$kydj",
  ]);
  expect(report.contacts).toHaveLength(3);
});
