"""Pure 3D+three-view IFC with an external JSON scene recipe.

Migration input may be a legacy project IFC; runtime scene reconstruction uses
only the pure IFC, this JSON recipe, and a disposable formal-project copy.
"""
import numpy as np
import ifcopenshell
import ifcopenshell.util.element as eu
import review_product_package as pkg


def spatial_scope(source, target):
    scope = {target, *source.by_type("IfcProject")}
    pending = [target]
    while pending:
        item = pending.pop()
        parents = [r.RelatingStructure for r in getattr(item, "ContainedInStructure", ())]
        parents += [r.RelatingObject for r in getattr(item, "Decomposes", ())]
        for parent in parents:
            if parent not in scope:
                scope.add(parent)
                pending.append(parent)
    return scope


def no_spatial_geometry(scope):
    return {(e.id(), "Representation"): None for e in scope if e.is_a("IfcSpatialElement")}


def representation_content(rep):
    # A representation may be renamed from Annotation to ApprovedPlan without
    # changing its local geometry, type, context or inverse presentation style.
    return pkg.fingerprint((rep.ContextOfItems, rep.RepresentationType, rep.Items))


def graph_to_json(model):
    def value(x):
        if isinstance(x, ifcopenshell.entity_instance):
            if x.id():
                return {"ref": x.id()}
            return {"inline": x.is_a(), "values": [value(v) for v in x]}
        if isinstance(x, (tuple, list)):
            return [value(v) for v in x]
        return x
    return {"schema": model.schema, "entities": [{"id": e.id(), "type": e.is_a(), "attributes": [value(v) for v in e]} for e in model]}


def graph_from_json(graph):
    model = ifcopenshell.file(schema=graph["schema"])
    mapping = {r["id"]: model.create_entity(r["type"]) for r in graph["entities"]}
    def value(x):
        if isinstance(x, dict):
            if "ref" in x:
                return mapping[x["ref"]]
            return model.create_entity(x["inline"], *[value(v) for v in x["values"]])
        if isinstance(x, list):
            return tuple(value(v) for v in x)
        return x
    for row in graph["entities"]:
        for i, x in enumerate(row["attributes"]):
            if x is not None:
                mapping[row["id"]][i] = value(x)
    return model


def build_pure_package(source, target_guid, views):
    """views maps plan/front/side to annotation_guid, drawing_guid and optional
    representation_identifier. Existing approved product reps are reused when
    explicitly selected; otherwise the corresponding annotation rep is copied.
    All original Body reps are retained. No annotation entity is stored in IFC.
    """
    assert set(views) == {"plan", "front", "side"}
    target = source.by_guid(target_guid)
    if any(getattr(target, attr, ()) for attr in ("IsDecomposedBy", "IsNestedBy", "HasOpenings", "HasPorts")):
        raise ValueError("Nested parts/openings/ports require a product-specific policy")
    scope = spatial_scope(source, target)
    typ = eu.get_type(target)
    if typ:
        scope.add(typ)
    body = [r for r in target.Representation.Representations if r.RepresentationIdentifier == "Body"]
    assert body
    selected = list(body)
    selected_view_reps = {}
    for view, spec in views.items():
        identifier = spec.get("representation_identifier")
        if identifier:
            reps = [r for r in target.Representation.Representations if r.RepresentationIdentifier == identifier]
            assert len(reps) == 1
            selected.append(reps[0])
            selected_view_reps[view] = reps[0]
    overrides = no_spatial_geometry(scope)
    overrides[(target.Representation.id(), "Representations")] = tuple(selected)
    pure = ifcopenshell.file(schema=source.schema)
    clone = pkg.ScopedCopy(source, pure, scope, attribute_overrides=overrides)
    for e in sorted(scope, key=lambda e: e.id()):
        clone.copy(e)
    product = pure.by_guid(target_guid)
    details = {}
    for view, spec in views.items():
        if view in selected_view_reps:
            original = selected_view_reps[view]
            rep = clone.memo[original.id()]
        else:
            annotation = source.by_guid(spec["annotation_guid"])
            assert np.allclose(pkg.placement(annotation), pkg.placement(target), atol=1e-7), "Annotation placement differs; explicit coordinate conversion required"
            reps = annotation.Representation.Representations
            assert len(reps) == 1, "Multiple annotation representations need explicit selection"
            original = reps[0]
            rep = clone.copy(original)
            if rep in product.Representation.Representations:
                # IFC may legally share one geometric representation across
                # views. Keep independent view identifiers without duplicating
                # the underlying geometry or renaming a previously added view.
                rep = pure.create_entity(rep.is_a(), **{
                    rep.attribute_name(i): value for i, value in enumerate(rep)})
            rep.RepresentationIdentifier = f"Approved{view.title()}"
            product.Representation.Representations += (rep,)
        details[view] = {"representation_identifier": rep.RepresentationIdentifier,
                         "content_fingerprint": representation_content(rep),
                         "source_annotation_guid": spec.get("annotation_guid")}
        assert representation_content(original) == representation_content(rep)
    clone.inverses()
    repairs = pkg.repair_missing_metadata(pure)
    assert pkg.body_fingerprint(target) == pkg.body_fingerprint(product)
    assert np.array_equal(pkg.placement(target), pkg.placement(product))
    assert not pure.by_type("IfcAnnotation") and not pure.by_type("IfcGroup")
    assert all(not e.Representation for e in pure.by_type("IfcSpatialElement"))
    assert len(pure.by_type("IfcElement")) == 1
    # No scene camera, Include filter or drawing pset is allowed in the package.
    assert not [p for p in pure.by_type("IfcPropertySet") if p.Name == "EPset_Drawing"]
    assert not [p for p in pure.by_type("IfcPropertySingleValue") if p.Name in ("Include", "Exclude")]

    annotations = [source.by_guid(spec[key]) for spec in views.values() for key in ("annotation_guid", "drawing_guid")]
    recipe_scope = spatial_scope(source, target) | set(annotations)
    for e in annotations:
        recipe_scope.update(r.RelatingGroup for r in getattr(e, "HasAssignments", ()) if r.is_a("IfcRelAssignsToGroup"))
    recipe_file = ifcopenshell.file(schema=source.schema)
    recipe_overrides = no_spatial_geometry(recipe_scope)
    recipe_overrides[(target.id(), "Representation")] = None
    for spec in views.values():
        recipe_overrides[(source.by_guid(spec["annotation_guid"]).id(), "Representation")] = None
    copier = pkg.ScopedCopy(source, recipe_file, recipe_scope, attribute_overrides=recipe_overrides, skip_inverse_ids=[target.id()])
    for e in sorted(recipe_scope, key=lambda e: e.id()):
        copier.copy(e)
    for view, spec in views.items():
        context = source.by_guid(spec["annotation_guid"]).Representation.Representations[0].ContextOfItems
        copied_context = copier.copy(context)
        details[view]["scene_context_fingerprint"] = pkg.fingerprint(copied_context)
    copier.inverses()
    assert all(not e.Representation for e in recipe_file.by_type("IfcElement")), "External recipe must not hide another copy of 3D Body"
    recipe = {"schema_version": 1, "target_global_id": target_guid, "views": views,
              "format": "scoped-ifc-scene-dependency-graph", "body_geometry_included": False,
              "linework_geometry_included": False, "pure_view_representations": details,
              "runtime_inputs": ["pure_product_ifc", "scene_recipe_json", "formal_project_temporary_copy"],
              "graph": graph_to_json(recipe_file)}
    audit = {"target_global_id": target_guid, "body_unchanged": True, "body_representation_count": len(body),
             "placement_matrix_project_units": pkg.placement(product).tolist(), "views": details,
             "ifc_annotation_count": 0, "ifc_group_count": 0, "represented_products": 1,
             "spatial_geometry_count": 0, "scene_settings_in_ifc": False, "inherited_schema_repairs": repairs}
    return pure, recipe, audit


def attach_recipe(project, pure, recipe):
    """Mutates only the caller-supplied disposable in-memory project copy."""
    guid = recipe["target_global_id"]
    report = pkg.attach(project, pure, guid, [])
    data = graph_from_json(recipe["graph"])
    # The pure IFC is the sole geometry source: recipe holds only annotation
    # identities, placements and settings, never a second copy of approved 2D.
    importer = pkg.ScopedCopy(pure, data, pure.by_type("IfcRoot"), reuse_roots=True)
    for view, spec in recipe["views"].items():
        identifier = recipe["pure_view_representations"][view]["representation_identifier"]
        reps = [r for r in pure.by_guid(guid).Representation.Representations if r.RepresentationIdentifier == identifier]
        assert len(reps) == 1
        representation = importer.copy(reps[0])
        representation.RepresentationIdentifier = "Annotation"
        context_fingerprint = recipe["pure_view_representations"][view].get("scene_context_fingerprint")
        if context_fingerprint:
            contexts = [c for c in data.by_type("IfcGeometricRepresentationContext") if pkg.fingerprint(c) == context_fingerprint]
            assert contexts, "External recipe lost its scene annotation context"
            representation.ContextOfItems = contexts[0]
        data.by_guid(spec["annotation_guid"]).Representation = data.create_entity("IfcProductDefinitionShape", Representations=[representation])
    importer.inverses()
    roots = data.by_type("IfcRoot")
    clone = pkg.ScopedCopy(data, project, roots, reuse_roots=True, skip_inverse_ids=[data.by_guid(guid).id()])
    clone.memo[data.by_guid(guid).id()] = project.by_guid(guid)
    clone.memo[data.by_guid(guid).ObjectPlacement.id()] = project.by_guid(guid).ObjectPlacement
    for context in data.by_type("IfcGeometricRepresentationContext"):
        existing = next((c for c in project.by_type(context.is_a(), include_subtypes=False) if pkg.fingerprint(c) == pkg.fingerprint(context)), None)
        if existing:
            clone.memo[context.id()] = existing
    for e in [*data.by_type("IfcAnnotation"), *data.by_type("IfcGroup")]:
        clone.copy(e)
    clone.inverses()
    report["linework_geometry_source"] = "pure_product_ifc_representations"
    assert len(project.by_type("IfcElement")) == report["physical_element_count_before"]
    for spec in recipe["views"].values():
        for key in ("annotation_guid", "drawing_guid"):
            expected = data.by_guid(spec[key])
            actual = project.by_guid(spec[key])
            assert np.array_equal(pkg.placement(expected), pkg.placement(actual))
            assert representation_content(expected.Representation.Representations[0]) == representation_content(actual.Representation.Representations[0])
    return report
