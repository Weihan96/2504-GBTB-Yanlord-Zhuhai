"""Scoped IFC dependency extraction and identity-matched review attachment.

Pilot adapter, not a general-purpose IFC merge. Unsupported relationships fail
closed; physical objects outside the explicitly selected scope are never copied.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import ifcopenshell
import ifcopenshell.util.element as element_util
import ifcopenshell.util.placement as placement_util
import ifcopenshell.util.unit as unit_util


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(path):
    path = Path(path)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def canonical(value):
    """STEP-id-independent forward graph fingerprint, including coordinates."""
    if isinstance(value, ifcopenshell.entity_instance):
        attrs = []
        for i, v in enumerate(value):
            if value.is_a("IfcGeometricRepresentationContext") and value.attribute_name(i) == "WorldCoordinateSystem":
                # Older source has a missing Plan WCS. It was rendered at the
                # identity origin. Compare effective coordinates, not null vs
                # explicit identity; the repair itself is separately audited.
                wcs = getattr(value, "WorldCoordinateSystem", None)
                matrix = placement_util.get_axis2placement(wcs) if wcs else np.eye(4)
                attrs.append({"effective_world_coordinate_system": matrix.tolist()})
            else:
                attrs.append(canonical(v))
        return [value.is_a(), *attrs]
    if isinstance(value, (tuple, list)):
        return [canonical(v) for v in value]
    return value


def fingerprint(value):
    return hashlib.sha256(json.dumps(canonical(value), sort_keys=True).encode()).hexdigest()


def body_fingerprint(product):
    return fingerprint([r for r in product.Representation.Representations if r.RepresentationIdentifier == "Body"])


def placement(product):
    return placement_util.get_local_placement(product.ObjectPlacement)


def find_guid(model, guid):
    try:
        return model.by_guid(guid)
    except RuntimeError:
        return None


class ScopedCopy:
    """Clone forward dependencies and necessary filtered inverse relations.

    Unlike file.add(rel), this never follows all RelatedObjects into other
    products. Style/layer/material inverses are copied explicitly.
    """
    def __init__(self, source, dest, scope, reuse_roots=False, attribute_overrides=None, skip_inverse_ids=()):
        self.source, self.dest = source, dest
        self.scope = {e.id() for e in scope}
        self.memo = {}
        self.reuse_roots = reuse_roots
        self.attribute_overrides = attribute_overrides or {}
        self.skip_inverse_ids = set(skip_inverse_ids)

    def copy(self, value):
        if isinstance(value, (tuple, list)):
            return tuple(self.copy(v) for v in value)
        if not isinstance(value, ifcopenshell.entity_instance):
            return value
        if value.id() == 0:
            return self.dest.create_entity(value.is_a(), *[self.copy(v) for v in value])
        if value.id() in self.memo:
            return self.memo[value.id()]
        if value.is_a("IfcProduct") and value.id() not in self.scope:
            raise ValueError(f"Out-of-scope product dependency: {value}")
        if self.reuse_roots and value.is_a("IfcRoot"):
            existing = find_guid(self.dest, value.GlobalId)
            if existing:
                self.memo[value.id()] = existing
                return existing
        out = self.dest.create_entity(value.is_a())
        self.memo[value.id()] = out
        for i, attr in enumerate(value):
            name = value.attribute_name(i)
            attr = self.attribute_overrides.get((value.id(), name), attr)
            # Derived IFC attributes (e.g. subcontext WCS) appear as None but
            # must remain derived, not be assigned a null value.
            if attr is None:
                continue
            if name in ("RelatedObjects", "RelatedElements", "RelatedDefinitions"):
                attr = tuple(e for e in attr if e.id() in self.scope)
                if not attr:
                    raise ValueError(f"Empty scoped relationship: {value.is_a()}")
            elif name == "AssignedItems":
                attr = tuple(e for e in attr if e.id() in self.memo)
            out[i] = self.copy(attr)
        return out

    def inverses(self):
        checked = set()
        supported = ("IfcRelDefinesByProperties", "IfcRelDefinesByType", "IfcRelAssociates",
                     "IfcRelContainedInSpatialStructure", "IfcRelAggregates", "IfcRelNests",
                     "IfcRelAssignsToGroup", "IfcRelAssignsToProduct", "IfcRelDeclares")
        while True:
            pending = set(self.memo) - checked
            if not pending:
                return
            for source_id in pending:
                checked.add(source_id)
                if source_id in self.skip_inverse_ids:
                    continue
                entity = self.source.by_id(source_id)
                for inv in self.source.get_inverse(entity):
                    if inv.id() in self.memo:
                        continue
                    if inv.is_a("IfcRelationship"):
                        if entity.id() in self.scope and any(inv.is_a(cls) for cls in ("IfcRelConnectsElements", "IfcRelConnectsPorts", "IfcRelConnectsPortToElement", "IfcRelVoidsElement", "IfcRelFillsElement", "IfcRelProjectsElement")):
                            raise ValueError(f"Connection/opening dependencies need a product-specific policy: {inv.is_a()}")
                        related = getattr(inv, "RelatedObjects", getattr(inv, "RelatedElements", getattr(inv, "RelatedDefinitions", ())))
                        if not any(e.id() in self.scope for e in related):
                            continue
                        if not any(inv.is_a(cls) for cls in supported):
                            raise ValueError(f"Unsupported scoped relationship: {inv.is_a()}")
                        # A group/container outside this scope is not necessary.
                        owners = [getattr(inv, n, None) for n in ("RelatingGroup", "RelatingStructure", "RelatingObject", "RelatingProduct", "RelatingContext")]
                        if any(o and o.id() not in self.scope for o in owners):
                            continue
                        self.copy(inv)
                    elif inv.is_a("IfcStyledItem") or inv.is_a("IfcMaterialDefinitionRepresentation") or inv.is_a("IfcMaterialProperties"):
                        self.copy(inv)
                    elif inv.is_a("IfcPresentationLayerAssignment"):
                        self.copy(inv)
                    elif inv.is_a("IfcShapeAspect"):
                        self.copy(inv)
                    elif inv.is_a("IfcCoordinateOperation"):
                        self.copy(inv)


def selected_scope(source, target, annotations):
    scope = {target, *annotations, *source.by_type("IfcProject")}
    typ = element_util.get_type(target)
    if typ:
        scope.add(typ)
    pending = list(scope)
    while pending:
        e = pending.pop()
        parents = [r.RelatingStructure for r in getattr(e, "ContainedInStructure", ())]
        parents += [r.RelatingObject for r in getattr(e, "Decomposes", ())]
        parents += [r.RelatingGroup for r in getattr(e, "HasAssignments", ()) if r.is_a("IfcRelAssignsToGroup")]
        for p in parents:
            if p not in scope:
                scope.add(p)
                pending.append(p)
    return scope


def extract(source, guid, annotation_guids):
    target = source.by_guid(guid)
    if any(getattr(target, attr, ()) for attr in ("IsDecomposedBy", "IsNestedBy", "HasOpenings", "HasPorts")):
        raise ValueError("Nested parts/openings/ports require a product-specific extraction policy")
    annotations = [source.by_guid(g) for g in annotation_guids]
    scope = selected_scope(source, target, annotations)
    dest = ifcopenshell.file(schema=source.schema)
    clone = ScopedCopy(source, dest, scope)
    for root in sorted(scope, key=lambda e: e.id()):
        clone.copy(root)
    clone.inverses()
    assert {e.GlobalId for e in dest.by_type("IfcElement")} == {guid}, "Pilot must contain exactly one physical element"
    assert body_fingerprint(target) == body_fingerprint(dest.by_guid(guid))
    assert np.array_equal(placement(target), placement(dest.by_guid(guid)))
    assert unit_util.calculate_unit_scale(source) == unit_util.calculate_unit_scale(dest)
    return dest


def repair_missing_metadata(package):
    """Two documented inherited defects; never guess a non-identity frame."""
    repairs = []
    model_contexts = [c for c in package.by_type("IfcGeometricRepresentationContext", include_subtypes=False) if c.ContextType == "Model"]
    for context in package.by_type("IfcGeometricRepresentationContext", include_subtypes=False):
        if context.WorldCoordinateSystem is None:
            assert context.ContextType == "Plan" and context.CoordinateSpaceDimension == 2
            assert model_contexts and all(np.array_equal(placement_util.get_axis2placement(c.WorldCoordinateSystem), np.eye(4)) for c in model_contexts)
            point = package.create_entity("IfcCartesianPoint", Coordinates=(0., 0.))
            context.WorldCoordinateSystem = package.create_entity("IfcAxis2Placement2D", Location=point)
            repairs.append({"entity": context.id(), "attribute": "WorldCoordinateSystem", "old": None,
                            "new": "explicit 2D identity origin", "basis": "source Model WCS identity; approved SVG registration verified independently"})
    for app in package.by_type("IfcApplication"):
        if app.ApplicationDeveloper is None:
            app.ApplicationDeveloper = package.create_entity("IfcOrganization", Name="Unspecified (developer missing in source IFC)")
            repairs.append({"entity": app.id(), "attribute": "ApplicationDeveloper", "old": None,
                            "new": "explicit unspecified placeholder; no invented author identity"})
    return repairs


def without_ids(value):
    if isinstance(value, dict):
        return {k: without_ids(v) for k, v in value.items() if k != "id"}
    if isinstance(value, list):
        return [without_ids(v) for v in value]
    return value


def validate(source, package, guid, annotation_guids):
    a, b = source.by_guid(guid), package.by_guid(guid)
    assert body_fingerprint(a) == body_fingerprint(b), "Body drift"
    assert fingerprint(a.Representation) == fingerprint(b.Representation), "Representation drift"
    assert np.array_equal(placement(a), placement(b)), "Placement drift"
    assert unit_util.calculate_unit_scale(source) == unit_util.calculate_unit_scale(package)
    assert without_ids(element_util.get_psets(a)) == without_ids(element_util.get_psets(b)), "Property drift"
    assert element_util.get_type(a).GlobalId == element_util.get_type(b).GlobalId
    assert [r.RelatingStructure.GlobalId for r in a.ContainedInStructure] == [r.RelatingStructure.GlobalId for r in b.ContainedInStructure]
    for g in annotation_guids:
        x, y = source.by_guid(g), package.by_guid(g)
        assert fingerprint(x.Representation) == fingerprint(y.Representation)
        assert np.array_equal(placement(x), placement(y))
        assert without_ids(element_util.get_psets(x)) == without_ids(element_util.get_psets(y))
    def styles(model):
        styled = set()
        for g in [guid, *annotation_guids]:
            product = model.by_guid(g)
            for item in model.traverse(product.Representation):
                for style in getattr(item, "StyledByItem", ()):
                    styled.add(fingerprint(style.Styles))
        return sorted(styled)
    assert styles(source) == styles(package), "Representation style drift"
    assert sorted(fingerprint(x) for x in source.by_type("IfcCoordinateOperation")) == sorted(fingerprint(x) for x in package.by_type("IfcCoordinateOperation")), "Georeference drift"
    represented = [e for e in package.by_type("IfcElement") if e.Representation]
    assert len(represented) == 1 and represented[0].GlobalId == guid
    roots = package.by_type("IfcRoot")
    assert len({e.GlobalId for e in roots}) == len(roots), "Duplicate GlobalId"
    return {"represented_physical_elements": 1, "body_unchanged": True,
            "all_representations_unchanged": True, "placement_matrix": placement(b).tolist(),
            "placement_error_project_units": 0, "rotation_error": 0,
            "unit_scale_to_m": unit_util.calculate_unit_scale(package),
            "properties_type_and_container_preserved": True, "duplicate_global_ids": 0,
            "styles_preserved": True, "georeferencing_preserved": True,
            "annotations": annotation_guids}


def attach(project, package, guid, annotation_guids):
    """Attach only 2D expressions to matching product, never append its Body."""
    assert project.schema == package.schema
    assert unit_util.calculate_unit_scale(project) == unit_util.calculate_unit_scale(package)
    source = package.by_guid(guid)
    target = project.by_guid(guid)  # Absence fails; never create a second product.
    assert body_fingerprint(source) == body_fingerprint(target), "Formal Body differs from approved package"
    assert np.array_equal(placement(source), placement(target)), "Coordinate mismatch: no auto-realignment"
    before = len(project.by_type("IfcElement"))
    clone = ScopedCopy(package, project, package.by_type("IfcRoot"), reuse_roots=True)
    clone.memo[source.id()] = target
    clone.memo[source.ObjectPlacement.id()] = target.ObjectPlacement
    # Reuse geometrically equivalent contexts without copying the Body.
    for context in package.by_type("IfcGeometricRepresentationContext"):
        match = next((c for c in project.by_type(context.is_a(), include_subtypes=False) if fingerprint(c) == fingerprint(context)), None)
        if match:
            clone.memo[context.id()] = match
    for rep in source.Representation.Representations:
        if rep.RepresentationIdentifier == "Body":
            continue
        existing = [r for r in target.Representation.Representations if r.RepresentationIdentifier == rep.RepresentationIdentifier]
        if existing:
            assert len(existing) == 1 and fingerprint(existing[0]) == fingerprint(rep), "Conflicting representation"
            clone.memo[rep.id()] = existing[0]
        else:
            target.Representation.Representations += (clone.copy(rep),)
    for g in annotation_guids:
        existing = find_guid(project, g)
        if existing:
            assert fingerprint(existing.Representation) == fingerprint(package.by_guid(g).Representation)
            assert np.array_equal(placement(existing), placement(package.by_guid(g)))
        clone.copy(package.by_guid(g))
    # Bring scoped drawing groups and source associations, not project bodies.
    for group in package.by_type("IfcGroup"):
        clone.copy(group)
    clone.inverses()
    assert len(project.by_type("IfcElement")) == before
    assert body_fingerprint(source) == body_fingerprint(target)
    assert len([e for e in project.by_type("IfcRoot") if e.GlobalId == guid]) == 1
    return {"physical_element_count_before": before, "physical_element_count_after": before,
            "target_instances": 1, "target_body_copied": False, "target_body_unchanged": True,
            "target_placement_error_project_units": 0}
