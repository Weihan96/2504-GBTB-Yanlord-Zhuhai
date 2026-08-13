#!/usr/bin/env python3
"""Write the approved APP-017 Siemens laundry coordination batch.

The existing proxy GlobalId is retained as a project coordination-clearance
carrier. Two placement-only IfcElectricAppliance occurrences provide correct
washer/dryer semantics without inventing connector thickness, body stacking
offsets, door swings, or service-interface centres.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import deque
from pathlib import Path
from typing import Any

import ifcopenshell
import ifcopenshell.api
import ifcopenshell.geom
import ifcopenshell.util.element
import ifcopenshell.util.placement
import ifcopenshell.util.representation
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
FORMAL_IFC = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
TARGET_GUID = "3JAkt8PsX7vPfGKWLK5EKp"
WASHER_GUID = "1Uzcf8kU9J2w9K40VkzCgV"
DRYER_GUID = "3_ZSeTAOPNThu9wyl2sOxQ"
WASHER_TYPE_GUID = "0wsTrloKPSxvZqJlXigZ9F"
DRYER_TYPE_GUID = "2ZXl1V8MHNuP6wi5J4HNua"
AGGREGATE_REL_GUID = "2eCIcTJ3jQgwtDO3FBcKIt"
EXPECTED_PREWRITE_SHA256 = "bd5ee6d2cad49cd03bbef7c15ab8a9bc5a246925705f86032271c6aa3e103fb8"
COORDINATION_MIN_MM = np.array([50.0, -4552.714157, 0.0])
COORDINATION_SIZE_MM = np.array([800.0, 650.0, 1900.0])  # x depth, y width, z height


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=FORMAL_IFC)
    parser.add_argument("--output", type=Path, default=ROOT / "build/candidates/app017-siemens-stack.ifc")
    parser.add_argument("--report", type=Path, default=ROOT / "build/app017-siemens-stack/mechanical-qa.json")
    parser.add_argument("--expected-sha256", default=EXPECTED_PREWRITE_SHA256)
    parser.add_argument("--write-formal", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_by_guid(model: ifcopenshell.file, guid: str) -> Any | None:
    """Return a rooted entity or None across IfcOpenShell API variants."""
    try:
        return model.by_guid(guid)
    except RuntimeError:
        return None


def forward_entities(roots: tuple[Any, ...]) -> dict[int, str]:
    queue = deque(root for root in roots if root is not None)
    result: dict[int, str] = {}
    while queue:
        entity = queue.popleft()
        if not hasattr(entity, "id") or entity.id() in result:
            continue
        result[entity.id()] = str(entity)
        for value in entity:
            if hasattr(value, "id"):
                queue.append(value)
            elif isinstance(value, (tuple, list)):
                queue.extend(item for item in value if hasattr(item, "id"))
    return result


def physical_graph(product: Any) -> dict[int, str]:
    representations = tuple(product.Representation.Representations) if product.Representation else ()
    return forward_entities((product.ObjectPlacement, *representations))


def placement_matrix(product: Any) -> list[list[float]] | None:
    if not product.ObjectPlacement:
        return None
    return np.asarray(ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement), dtype=float).tolist()


def bbox_mm(model: ifcopenshell.file, product: Any) -> tuple[list[float], list[float], list[float]]:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    shape = ifcopenshell.geom.create_shape(settings, product)
    vertices = np.asarray(shape.geometry.verts, dtype=float).reshape((-1, 3)) * 1000.0
    minimum = vertices.min(axis=0)
    maximum = vertices.max(axis=0)
    return minimum.tolist(), maximum.tolist(), (maximum - minimum).tolist()


def add_type_psets(model: ifcopenshell.file, product_type: Any, model_reference: str) -> None:
    manufacturer = ifcopenshell.api.run(
        "pset.add_pset", model, product=product_type, name="Pset_ManufacturerTypeInformation"
    )
    ifcopenshell.api.run(
        "pset.edit_pset", model, pset=manufacturer,
        properties={"Manufacturer": "Siemens", "ModelReference": model_reference, "ModelLabel": "iQ500 10 kg"},
    )
    common = ifcopenshell.api.run(
        "pset.add_pset", model, product=product_type, name="Pset_ElectricApplianceTypeCommon"
    )
    ifcopenshell.api.run(
        "pset.edit_pset", model, pset=common, properties={"Reference": model_reference, "Status": "NEW"}
    )
    electrical = ifcopenshell.api.run(
        "pset.add_pset", model, product=product_type, name="Pset_ElectricalDeviceCommon"
    )
    ifcopenshell.api.run(
        "pset.edit_pset", model, pset=electrical,
        properties={"RatedCurrent": 10.0, "RatedVoltage": 220.0},
    )


def create_semantic_appliance(
    model: ifcopenshell.file,
    *,
    occurrence_guid: str,
    type_guid: str,
    name: str,
    model_reference: str,
    predefined_type: str,
    description: str,
) -> tuple[Any, Any]:
    if safe_by_guid(model, occurrence_guid) or safe_by_guid(model, type_guid):
        raise RuntimeError(f"semantic laundry identity already exists: {occurrence_guid}/{type_guid}")
    occurrence = ifcopenshell.api.run(
        "root.create_entity", model, ifc_class="IfcElectricAppliance",
        predefined_type=predefined_type, name=name,
    )
    occurrence.GlobalId = occurrence_guid
    occurrence.Tag = "APP-017-WASHER" if predefined_type == "WASHINGMACHINE" else "APP-017-DRYER"
    occurrence.ObjectType = "Independent appliance in stacked laundry set"
    occurrence.Description = description
    product_type = ifcopenshell.api.run(
        "root.create_entity", model, ifc_class="IfcElectricApplianceType",
        predefined_type=predefined_type, name=f"Siemens {model_reference}",
    )
    product_type.GlobalId = type_guid
    product_type.ElementType = model_reference
    add_type_psets(model, product_type, model_reference)
    ifcopenshell.api.run(
        "type.assign_type", model, related_objects=[occurrence], relating_type=product_type,
        should_map_representations=False,
    )
    return occurrence, product_type


def main() -> int:
    args = parse_args()
    source_path = args.input.resolve()
    output_path = args.output.resolve()
    prewrite_hash = sha256(source_path)
    if prewrite_hash != args.expected_sha256:
        raise RuntimeError(f"formal IFC hash changed: expected {args.expected_sha256}, got {prewrite_hash}")
    if args.write_formal and source_path != FORMAL_IFC.resolve():
        raise RuntimeError("--write-formal is only permitted for the formal project IFC")

    source = ifcopenshell.open(source_path)
    candidate = ifcopenshell.open(source_path)
    target_source = source.by_guid(TARGET_GUID)
    target = candidate.by_guid(TARGET_GUID)
    if not target or target.is_a() != "IfcBuildingElementProxy":
        raise RuntimeError("APP-017 coordination proxy is missing or has the wrong class")
    if any(safe_by_guid(candidate, guid) for guid in (WASHER_GUID, DRYER_GUID, WASHER_TYPE_GUID, DRYER_TYPE_GUID)):
        raise RuntimeError("formal IFC already contains the APP-017 Siemens semantic children")

    existing_products = {product.GlobalId: product for product in source.by_type("IfcProduct")}
    protected_graphs = {
        guid: physical_graph(product) for guid, product in existing_products.items() if guid != TARGET_GUID
    }
    protected_placements = {
        guid: placement_matrix(product) for guid, product in existing_products.items() if guid != TARGET_GUID
    }
    source_product_count = len(existing_products)
    source_element_count = len(source.by_type("IfcElement"))
    source_root_ids = {root.GlobalId for root in source.by_type("IfcRoot")}
    target_before = {
        "class": target_source.is_a(), "name": target_source.Name or "",
        "description": target_source.Description or "", "bbox_mm": bbox_mm(source, target_source),
    }

    body_context = ifcopenshell.util.representation.get_context(candidate, "Model", "Body", "MODEL_VIEW")
    if body_context is None:
        raise RuntimeError("IFC Body/MODEL_VIEW context is missing")
    old_representations = list(target.Representation.Representations) if target.Representation else []
    for representation in old_representations:
        ifcopenshell.api.run(
            "geometry.unassign_representation", candidate, product=target, representation=representation
        )
        ifcopenshell.api.run("geometry.remove_representation", candidate, representation=representation)
    clearance_representation = ifcopenshell.api.run(
        "geometry.add_wall_representation", candidate, context=body_context,
        length=float(COORDINATION_SIZE_MM[0] / 1000.0),
        thickness=float(COORDINATION_SIZE_MM[1] / 1000.0),
        height=float(COORDINATION_SIZE_MM[2] / 1000.0),
    )
    ifcopenshell.api.run(
        "geometry.assign_representation", candidate, product=target, representation=clearance_representation
    )
    matrix = np.eye(4)
    matrix[:3, 3] = COORDINATION_MIN_MM / 1000.0
    ifcopenshell.api.run(
        "geometry.edit_object_placement", candidate, product=target, matrix=matrix,
        is_si=True, should_transform_children=False,
    )
    target.Name = "Siemens stacked washer + heat-pump dryer · coordination clearance"
    target.Description = "Project coordination clearance W650 D800 H1900; not manufacturer body or interface dimensions"
    target.ObjectType = "Laundry equipment coordination clearance"
    target.Tag = "APP-017"

    washer, washer_type = create_semantic_appliance(
        candidate, occurrence_guid=WASHER_GUID, type_guid=WASHER_TYPE_GUID,
        name="Siemens WG54M7D20W independent washer", model_reference="WG54M7D20W",
        predefined_type="WASHINGMACHINE",
        description="Independent washer in APP-017 stacked set; manufacturer body 598W x 848H x 600D mm, door-closed depth 639 mm, 90-degree open depth 1111 mm; placement only until connector/shop drawing is confirmed.",
    )
    dryer, dryer_type = create_semantic_appliance(
        candidate, occurrence_guid=DRYER_GUID, type_guid=DRYER_TYPE_GUID,
        name="Siemens WQ55M7U20W independent heat-pump dryer", model_reference="WQ55M7U20W",
        predefined_type="TUMBLEDRYER",
        description="Independent heat-pump dryer in APP-017 stacked set; manufacturer body 598W x 842H x 600D mm, door-closed depth 639 mm, 90-degree open depth 1111 mm; placement only until connector/shop drawing is confirmed.",
    )
    for appliance in (washer, dryer):
        ifcopenshell.api.run(
            "geometry.edit_object_placement", candidate, product=appliance, matrix=matrix,
            is_si=True, should_transform_children=False,
        )
    aggregate = ifcopenshell.api.run(
        "aggregate.assign_object", candidate, products=[washer, dryer], relating_object=target
    )
    if aggregate is None:
        raise RuntimeError("failed to aggregate the two semantic appliances under APP-017")
    aggregate.GlobalId = AGGREGATE_REL_GUID

    output_path.parent.mkdir(parents=True, exist_ok=True)
    candidate.write(output_path)
    reopened = ifcopenshell.open(output_path)
    target_after = reopened.by_guid(TARGET_GUID)
    washer_after = reopened.by_guid(WASHER_GUID)
    dryer_after = reopened.by_guid(DRYER_GUID)
    washer_type_after = ifcopenshell.util.element.get_type(washer_after) if washer_after else None
    dryer_type_after = ifcopenshell.util.element.get_type(dryer_after) if dryer_after else None
    if not washer_after or not washer_type_after or washer_type_after.PredefinedType != "WASHINGMACHINE":
        raise RuntimeError("washer semantic occurrence failed postwrite verification")
    if not dryer_after or not dryer_type_after or dryer_type_after.PredefinedType != "TUMBLEDRYER":
        raise RuntimeError("dryer semantic occurrence failed postwrite verification")
    if washer_after.Representation or dryer_after.Representation:
        raise RuntimeError("semantic children must remain placement-only until connector/shop drawing confirmation")
    decomposed = {child.GlobalId for child in ifcopenshell.util.element.get_decomposition(target_after)}
    if decomposed != {WASHER_GUID, DRYER_GUID}:
        raise RuntimeError(f"unexpected APP-017 decomposition: {decomposed}")

    protected_changed: list[str] = []
    maximum_placement_delta = 0.0
    for guid, before_graph in protected_graphs.items():
        product = reopened.by_guid(guid)
        if product is None or physical_graph(product) != before_graph:
            protected_changed.append(guid)
            continue
        before_matrix = protected_placements[guid]
        after_matrix = placement_matrix(product)
        if before_matrix is None or after_matrix is None:
            if before_matrix != after_matrix:
                protected_changed.append(guid)
            continue
        delta = float(np.max(np.abs(np.asarray(before_matrix) - np.asarray(after_matrix)), initial=0.0))
        maximum_placement_delta = max(maximum_placement_delta, delta)
        if delta != 0.0:
            protected_changed.append(guid)
    if protected_changed:
        raise RuntimeError(f"protected IFC product geometry/placement changed: {protected_changed[:10]}")

    target_after_bbox = bbox_mm(reopened, target_after)
    expected_max = COORDINATION_MIN_MM + COORDINATION_SIZE_MM
    if not np.allclose(target_after_bbox[0], COORDINATION_MIN_MM, atol=1e-6) or not np.allclose(target_after_bbox[1], expected_max, atol=1e-6):
        raise RuntimeError(f"APP-017 clearance bbox mismatch: {target_after_bbox}")
    new_root_ids = {root.GlobalId for root in reopened.by_type("IfcRoot")} - source_root_ids
    report = {
        "mode": "formal_write" if args.write_formal else "candidate_only",
        "source": {"path": str(source_path), "sha256": prewrite_hash, "schema": source.schema},
        "candidate": {"path": str(output_path), "sha256": sha256(output_path)},
        "counts": {
            "source_ifc_products": source_product_count,
            "candidate_ifc_products": len(reopened.by_type("IfcProduct")),
            "source_ifc_elements": source_element_count,
            "candidate_ifc_elements": len(reopened.by_type("IfcElement")),
            "added_product_count": len(reopened.by_type("IfcProduct")) - source_product_count,
            "new_root_count": len(new_root_ids),
        },
        "target_before": target_before,
        "target_after": {
            "class": target_after.is_a(), "name": target_after.Name,
            "description": target_after.Description, "bbox_mm": target_after_bbox,
            "role": "project coordination clearance, not manufacturer product body",
        },
        "semantic_children": [
            {"global_id": WASHER_GUID, "ifc_class": washer_after.is_a(), "type_predefined_type": washer_type_after.PredefinedType, "representation": "placement-only"},
            {"global_id": DRYER_GUID, "ifc_class": dryer_after.is_a(), "type_predefined_type": dryer_type_after.PredefinedType, "representation": "placement-only"},
        ],
        "mechanical_qa": {
            "protected_existing_product_count": len(protected_graphs),
            "protected_physical_graph_changes": 0,
            "protected_world_placement_max_delta": maximum_placement_delta,
            "protected_world_geometry_change_mm": 0.0,
            "proof_basis": "exact recursive representation/ObjectPlacement graph equality plus exact world-placement matrix equality for every pre-existing IfcProduct except APP-017",
            "app017_intended_change_only": True,
            "semantic_children_have_no_invented_body": True,
            "ifc_schema_unchanged": reopened.schema == source.schema,
        },
        "open_gates": [
            "Siemens written confirmation of exact E-Nr./FD stack compatibility",
            "17008829 commercial model and drawer configuration",
            "final side-accessible power/water/drain shop drawing and service extraction clearance",
        ],
    }
    if args.write_formal:
        candidate.write(source_path)
        report["formal_postwrite"] = {"path": str(source_path), "sha256": sha256(source_path)}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
