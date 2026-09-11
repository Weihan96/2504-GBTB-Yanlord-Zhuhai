"""Pure IFC/recipe boundary tests independent of Blender rendering."""
import json
import unittest
import ifcopenshell
import ifcopenshell.api.root
from test_review_product_package import fixture
import pure_product_package as pure
import review_product_package as pkg


def three_views():
    src, formal, guid, annotation = fixture()
    target = src.by_guid(guid)
    views = {}
    for view in ("plan", "front", "side"):
        ann = ifcopenshell.api.root.create_entity(src, ifc_class="IfcAnnotation", name=view)
        ann.ObjectPlacement = target.ObjectPlacement
        ann.Representation = src.by_guid(annotation).Representation
        drawing = ifcopenshell.api.root.create_entity(src, ifc_class="IfcAnnotation", name=view + " camera")
        drawing.ObjectPlacement = target.ObjectPlacement
        drawing.Representation = src.by_guid(annotation).Representation
        views[view] = {"annotation_guid": ann.GlobalId, "drawing_guid": drawing.GlobalId}
    return src, formal, guid, views


class PurePackageTests(unittest.TestCase):
    def test_pure_boundary_and_json_roundtrip(self):
        src, _, guid, views = three_views()
        model, recipe, audit = pure.build_pure_package(src, guid, views)
        model = ifcopenshell.file.from_string(model.to_string())
        self.assertEqual(len(model.by_type("IfcElement")), 1)
        self.assertFalse(model.by_type("IfcAnnotation"))
        self.assertFalse(model.by_type("IfcGroup"))
        self.assertEqual(len(model.by_guid(guid).Representation.Representations), 4)
        graph = pure.graph_from_json(json.loads(json.dumps(recipe))["graph"])
        self.assertIsNone(graph.by_guid(guid).Representation)
        for spec in views.values():
            self.assertIsNone(graph.by_guid(spec["annotation_guid"]).Representation)
        self.assertFalse(recipe["body_geometry_included"])
        self.assertFalse(recipe["linework_geometry_included"])
        self.assertEqual(pkg.body_fingerprint(src.by_guid(guid)), pkg.body_fingerprint(model.by_guid(guid)))

    def test_attach_uses_pure_geometry_and_is_idempotent(self):
        src, formal, guid, views = three_views()
        model, recipe, _ = pure.build_pure_package(src, guid, views)
        result = pure.attach_recipe(formal, model, recipe)
        count = len(list(formal))
        pure.attach_recipe(formal, model, recipe)
        self.assertEqual(count, len(list(formal)))
        self.assertEqual(len(formal.by_type("IfcElement")), 2)
        self.assertEqual(result["linework_geometry_source"], "pure_product_ifc_representations")
        self.assertEqual(pkg.body_fingerprint(src.by_guid(guid)), pkg.body_fingerprint(formal.by_guid(guid)))

    def test_annotation_coordinate_drift_fails_closed(self):
        src, _, guid, views = three_views()
        ann = src.by_guid(views["side"]["annotation_guid"])
        ann.ObjectPlacement = src.create_entity("IfcLocalPlacement", RelativePlacement=src.create_entity(
            "IfcAxis2Placement3D", Location=src.create_entity("IfcCartesianPoint", Coordinates=(0., 0., 0.))))
        with self.assertRaisesRegex(AssertionError, "Annotation placement differs"):
            pure.build_pure_package(src, guid, views)


if __name__ == "__main__":
    unittest.main()
