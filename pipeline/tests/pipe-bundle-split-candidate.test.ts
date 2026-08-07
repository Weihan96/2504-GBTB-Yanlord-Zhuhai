import { describe, expect, test } from "bun:test";
import { readFileSync } from "node:fs";

const script = readFileSync(
  "pipeline/scripts/pipe_bundle_split_candidate.py",
  "utf8",
);

describe("PVC110 bundle split candidate", () => {
  test("scope is the exact two disjoint three-pipe bundles", () => {
    expect(script).toContain('"178mqyyzzFowLcbXcH6prO"');
    expect(script).toContain('"0bfVg4Ys1CevZs$qxhkXTo"');
    expect(script).toContain("EXPECTED_BUNDLE_SIZE = 3");
    expect(script).toContain("minimum_pair_bbox_distance_mm");
  });

  test("candidate preserves union geometry and adds four roots", () => {
    expect(script).toContain("union_vertex_hausdorff_mm");
    expect(script).toContain("root_delta_matches_expected");
    expect(script).toContain('"candidate_pipe_segments"] == 6');
    expect(script).toContain('"formal_write_allowed": False');
    expect(script).toContain('"decision_status": "rejected_for_formal_write"');
  });
});
