#!/usr/bin/env python3
"""Audit safe recovery of high-poly candidates from the stash untracked tree."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
INVENTORY = ROOT / "pipeline/decisions/highpoly-product-inventory.json"
DEFAULT_OUTPUT = ROOT / "pipeline/decisions/highpoly-stash-recovery-audit.json"
TARGET_MESSAGE = "superseded highpoly review candidates"
PROTECTED_MESSAGE = "Libelle background studies from other task"
EXPECTED_TRACKED_OTHER_TASK_PATHS = {
    "output/images/libelle-background-studies/libelle-preview-warm-silver-white.png",
    "output/images/libelle-background-studies/libelle-scene-transparent-background.png",
}
SUCCESSOR_PATHS = {
    "pipeline/tests/int1-highpoly-type-review.test.ts":
        "pipeline/tests/falper-sorgente-highpoly-review.test.ts",
}


def git(*arguments: str) -> str:
    result = subprocess.run(
        ["git", *arguments],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stash_records() -> list[dict[str, str]]:
    records = []
    for line in git("stash", "list", "--format=%H%x09%gd%x09%gs").splitlines():
        commit, ref, subject = line.split("\t", 2)
        records.append({"commit": commit, "ref": ref, "subject": subject})
    return records


def unique_stash(records: list[dict[str, str]], message: str) -> dict[str, str]:
    matches = [record for record in records if message in record["subject"]]
    if len(matches) != 1:
        raise RuntimeError(f"expected one stash matching {message!r}, found {len(matches)}")
    return matches[0]


def ls_tree(tree: str) -> list[dict[str, str]]:
    records = []
    for line in git("ls-tree", "-r", tree).splitlines():
        metadata, path = line.split("\t", 1)
        mode, object_type, object_id = metadata.split()
        records.append(
            {
                "mode": mode,
                "object_type": object_type,
                "stash_blob_oid": object_id,
                "path": path,
            }
        )
    return records


def tracked_delta(stash_commit: str) -> list[dict[str, str]]:
    records = []
    for line in git("diff", "--name-status", f"{stash_commit}^1", stash_commit).splitlines():
        status, path = line.split("\t", 1)
        records.append({"status": status, "path": path})
    return records


def current_record(record: dict[str, str]) -> dict[str, Any]:
    current_path = SUCCESSOR_PATHS.get(record["path"], record["path"])
    path = ROOT / current_path
    exists = path.is_file()
    current_oid = git("hash-object", "--", current_path) if exists else None
    exact = (
        exists
        and current_path == record["path"]
        and current_oid == record["stash_blob_oid"]
    )
    return {
        **record,
        "current_path": current_path,
        "successor_mapping_used": current_path != record["path"],
        "current_exists": exists,
        "current_blob_oid": current_oid,
        "current_sha256": sha256(path) if exists else None,
        "recovery_state": (
            "recovered_exact"
            if exact
            else "recovered_and_superseded_by_audited_package"
            if exists
            else "missing"
        ),
    }


def product_slug(path: str) -> str | None:
    prefix = "output/review/highpoly-types/"
    if not path.startswith(prefix):
        return None
    remainder = path[len(prefix) :]
    return remainder.split("/", 1)[0] if "/" in remainder else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    before = stash_records()
    target = unique_stash(before, TARGET_MESSAGE)
    protected = unique_stash(before, PROTECTED_MESSAGE)
    parents = git("rev-list", "--parents", "-n", "1", target["commit"]).split()
    if len(parents) != 4:
        raise RuntimeError("target stash must contain base, index, and untracked parents")
    untracked_tree = f"{target['commit']}^3"
    recovered = [current_record(record) for record in ls_tree(untracked_tree)]
    tracked = tracked_delta(target["commit"])
    inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
    inventory_slugs = {Path(product["folder"]).name for product in inventory["products"]}
    candidate_slugs = sorted(
        slug for slug in {product_slug(record["path"]) for record in recovered} if slug
    )
    after = stash_records()
    checks = {
        "target_stash_found_by_subject": target["commit"] == "8620d0cd0ee91f0f6694b98d5099fb518aa65aaa",
        "target_has_untracked_parent": len(parents) == 4,
        "untracked_snapshot_has_expected_file_count": len(recovered) == 14,
        "all_untracked_snapshot_files_recovered": all(record["current_exists"] for record in recovered),
        "candidate_products_are_in_inventory": set(candidate_slugs) <= inventory_slugs,
        "candidate_products_recovered": candidate_slugs == ["falper-sorgente", "geberit-146-140"],
        "tracked_delta_is_other_task_only": {record["path"] for record in tracked} == EXPECTED_TRACKED_OTHER_TASK_PATHS,
        "protected_stash_is_distinct": protected["commit"] != target["commit"],
        "stash_list_unchanged": before == after,
    }
    payload = {
        "schema_version": 1,
        "generator": "pipeline/scripts/audit_highpoly_stash_recovery.py",
        "recovery_strategy": "read_only_extract_untracked_parent_then_supersede_with_audited_packages",
        "top_level_stash_command_performed": None,
        "stash_apply_performed": False,
        "stash_pop_performed": False,
        "stash_drop_performed": False,
        "stash_mutation_performed": False,
        "protected_other_task_stash": {
            "ref_observed_from_stash_list_only": protected["ref"],
            "commit": protected["commit"],
            "subject": protected["subject"],
            "content_inspected": False,
            "mutated": False,
        },
        "source_stash": target,
        "source_untracked_tree": untracked_tree,
        "tracked_delta_not_restored": tracked,
        "tracked_delta_reason": "Libelle background images belong to another task; applying the stash would violate task isolation.",
        "untracked_snapshot": recovered,
        "candidate_product_slugs": candidate_slugs,
        "inventory_product_count_after_recovery": len(inventory["products"]),
        "checks": checks,
        "pass": all(checks.values()),
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    try:
        output_label = str(output.relative_to(ROOT))
    except ValueError:
        output_label = str(output)
    print(json.dumps({"output": output_label, "pass": payload["pass"]}, ensure_ascii=False))
    if not payload["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
