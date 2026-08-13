#!/usr/bin/env python3
"""Prove the placement-only A06 IFC write changed no pre-existing product geometry."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from collections import deque
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.util.placement


A06_GLOBAL_ID = "0YzEUom7522RIg1TonOQXn"
A06_CONTAINER_GLOBAL_ID = "33eWFtzoP21ffkl7wZ6sk3"
A06_WORLD_TRANSLATION_MM = (-1068.494529, -905.287176, 2620.516300)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    before_source = parser.add_mutually_exclusive_group(required=True)
    before_source.add_argument("--before", type=Path)
    before_source.add_argument(
        "--before-git-ref",
        help="Git ref whose formal IFC is the frozen prewrite baseline",
    )
    parser.add_argument("--formal", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--expected-before-sha256", required=True)
    parser.add_argument("--expected-formal-sha256", required=True)
    return parser.parse_args()


def materialize_before_git_ref(
    git_ref: str, formal_path: Path, destination: Path
) -> dict[str, str]:
    repository = Path(
        subprocess.check_output(
            ["git", "-C", str(formal_path.parent), "rev-parse", "--show-toplevel"],
            text=True,
        ).strip()
    ).resolve()
    try:
        relative_formal = formal_path.relative_to(repository).as_posix()
    except ValueError as error:
        raise RuntimeError(
            f"formal IFC is outside its Git repository: {formal_path}"
        ) from error
    with destination.open("wb") as handle:
        result = subprocess.run(
            ["git", "-C", str(repository), "show", f"{git_ref}:{relative_formal}"],
            stdout=handle,
            stderr=subprocess.PIPE,
            check=False,
        )
    if result.returncode != 0:
        message = result.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(
            f"cannot read prewrite IFC from {git_ref}:{relative_formal}: {message}"
        )
    return {
        "source": "git_ref",
        "repository": str(repository),
        "git_ref": git_ref,
        "path": relative_formal,
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def forward_entities(roots: tuple[Any, ...]) -> dict[int, str]:
    queue = deque(root for root in roots if root is not None)
    entities: dict[int, str] = {}
    while queue:
        entity = queue.popleft()
        if not hasattr(entity, "id") or entity.id() in entities:
            continue
        entities[entity.id()] = str(entity)
        for value in entity:
            if hasattr(value, "id"):
                queue.append(value)
            elif isinstance(value, (tuple, list)):
                queue.extend(item for item in value if hasattr(item, "id"))
    return entities


def physical_graph(product: Any) -> dict[int, str]:
    representations = (
        tuple(product.Representation.Representations) if product.Representation else ()
    )
    return forward_entities((product.ObjectPlacement, *representations))


def matrix(model_product: Any) -> list[list[float]]:
    return [
        [float(value) for value in row]
        for row in ifcopenshell.util.placement.get_local_placement(
            model_product.ObjectPlacement
        )
    ]


def matrix_max_delta(first: list[list[float]], second: list[list[float]]) -> float:
    return max(
        abs(first[row][column] - second[row][column])
        for row in range(4)
        for column in range(4)
    )


def products_by_global_id(model: ifcopenshell.file) -> dict[str, Any]:
    return {
        product.GlobalId: product
        for product in model.by_type("IfcProduct")
        if getattr(product, "GlobalId", None)
    }


def main() -> None:
    args = parse_args()
    formal_path = args.formal.resolve()
    with tempfile.TemporaryDirectory(prefix="yanlord-a06-audit-") as temp_dir:
        if args.before:
            before_path = args.before.resolve()
            before_provenance = {"source": "filesystem", "path": str(before_path)}
        else:
            before_path = Path(temp_dir) / "prewrite.ifc"
            before_provenance = materialize_before_git_ref(
                args.before_git_ref, formal_path, before_path
            )
        run_audit(args, before_path, before_provenance, formal_path)


def run_audit(
    args: argparse.Namespace,
    before_path: Path,
    before_provenance: dict[str, str],
    formal_path: Path,
) -> None:
    before_sha = sha256(before_path)
    formal_sha = sha256(formal_path)
    if before_sha != args.expected_before_sha256:
        raise RuntimeError(
            f"before IFC hash drift: expected {args.expected_before_sha256}, found {before_sha}"
        )
    if formal_sha != args.expected_formal_sha256:
        raise RuntimeError(
            f"formal IFC hash drift: expected {args.expected_formal_sha256}, found {formal_sha}"
        )

    before = ifcopenshell.open(before_path)
    formal = ifcopenshell.open(formal_path)
    before_products = products_by_global_id(before)
    formal_products = products_by_global_id(formal)
    added_global_ids = sorted(set(formal_products) - set(before_products))
    removed_global_ids = sorted(set(before_products) - set(formal_products))

    changed_records = []
    maximum_placement_delta_mm = 0.0
    for global_id, before_product in before_products.items():
        formal_product = formal_products.get(global_id)
        if formal_product is None:
            continue
        before_matrix = matrix(before_product)
        formal_matrix = matrix(formal_product)
        placement_delta = matrix_max_delta(before_matrix, formal_matrix)
        maximum_placement_delta_mm = max(maximum_placement_delta_mm, placement_delta)
        graph_equal = physical_graph(before_product) == physical_graph(formal_product)
        class_equal = before_product.is_a() == formal_product.is_a()
        if placement_delta != 0.0 or not graph_equal or not class_equal:
            changed_records.append(
                {
                    "global_id": global_id,
                    "before_class": before_product.is_a(),
                    "formal_class": formal_product.is_a(),
                    "placement_delta_mm": placement_delta,
                    "placement_and_representation_graph_equal": graph_equal,
                }
            )

    a06 = formal.by_guid(A06_GLOBAL_ID)
    a06_matrix = matrix(a06)
    a06_translation = tuple(a06_matrix[index][3] for index in range(3))
    a06_translation_delta = max(
        abs(actual - expected)
        for actual, expected in zip(a06_translation, A06_WORLD_TRANSLATION_MM)
    )
    containers = [
        relation.RelatingStructure.GlobalId
        for relation in getattr(a06, "ContainedInStructure", ())
    ]
    typed_by = list(getattr(a06, "IsTypedBy", ()))
    ports = list(getattr(a06, "HasPorts", ()))
    systems = [
        relation
        for relation in formal.by_type("IfcRelAssignsToGroup")
        if a06 in relation.RelatedObjects
    ]
    a06_record = {
        "ifc_id": a06.id(),
        "global_id": a06.GlobalId,
        "ifc_class": a06.is_a(),
        "predefined_type": a06.PredefinedType,
        "name": a06.Name,
        "tag": a06.Tag,
        "description": a06.Description,
        "world_translation_mm": list(a06_translation),
        "maximum_translation_delta_mm": a06_translation_delta,
        "container_global_ids": containers,
        "has_representation": a06.Representation is not None,
        "type_assignment_count": len(typed_by),
        "port_count": len(ports),
        "system_assignment_count": len(systems),
    }

    entity_delta = len(list(formal)) - len(list(before))
    gates = {
        "before_product_count": len(before_products),
        "formal_product_count": len(formal_products),
        "entity_count_delta": entity_delta,
        "added_global_ids": added_global_ids,
        "removed_global_ids": removed_global_ids,
        "preexisting_product_changes": len(changed_records),
        "maximum_preexisting_product_world_geometry_change_mm": (
            0.0 if not changed_records else maximum_placement_delta_mm
        ),
        "preexisting_placement_and_representation_graphs_exact": not changed_records,
        "a06_semantics_exact": (
            a06.is_a() == "IfcUnitaryEquipment"
            and a06.PredefinedType == "AIRCONDITIONINGUNIT"
            and a06.Name == "A06"
            and a06.Tag == "A06"
        ),
        "a06_placement_exact": a06_translation_delta == 0.0,
        "a06_container_exact": containers == [A06_CONTAINER_GLOBAL_ID],
        "a06_is_placement_only": a06.Representation is None,
        "a06_has_no_unconfirmed_type_ports_or_system": (
            not typed_by and not ports and not systems
        ),
    }
    gates["pass"] = (
        gates["before_product_count"] == 924
        and gates["formal_product_count"] == 925
        and gates["entity_count_delta"] == 7
        and added_global_ids == [A06_GLOBAL_ID]
        and not removed_global_ids
        and not changed_records
        and gates["maximum_preexisting_product_world_geometry_change_mm"] == 0.0
        and gates["a06_semantics_exact"]
        and gates["a06_placement_exact"]
        and gates["a06_container_exact"]
        and gates["a06_is_placement_only"]
        and gates["a06_has_no_unconfirmed_type_ports_or_system"]
    )
    report = {
        "mode": "read-only-rcp1-a06-postwrite-audit",
        "before": {**before_provenance, "sha256": before_sha},
        "formal": {"path": str(formal_path), "sha256": formal_sha},
        "a06": a06_record,
        "changed_preexisting_products": changed_records,
        "gates": gates,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"report": str(args.report), **gates}, ensure_ascii=False))
    if not gates["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
