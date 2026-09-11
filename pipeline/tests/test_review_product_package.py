"""Run: python3 -m unittest discover -s pipeline/tests -p test_review_product_package.py"""
import sys
from pathlib import Path
import unittest
import numpy as np
import ifcopenshell
import ifcopenshell.guid
import ifcopenshell.api.root
import ifcopenshell.api.context
import ifcopenshell.api.unit
import ifcopenshell.api.geometry
import ifcopenshell.api.aggregate
import ifcopenshell.api.spatial
import ifcopenshell.api.type
import ifcopenshell.api.pset

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import review_product_package as pkg


def fixture():
    m = ifcopenshell.file(schema="IFC4")
    root = ifcopenshell.api.root.create_entity(m, ifc_class="IfcProject", name="Pilot")
    ifcopenshell.api.unit.assign_unit(m)
    context = ifcopenshell.api.context.add_context(m, context_type="Model")
    floor = ifcopenshell.api.root.create_entity(m, ifc_class="IfcBuildingStorey", name="FFL")
    ifcopenshell.api.aggregate.assign_object(m, products=[floor], relating_object=root)
    target = ifcopenshell.api.root.create_entity(m, ifc_class="IfcFurniture", name="Target")
    other = ifcopenshell.api.root.create_entity(m, ifc_class="IfcFurniture", name="Do not extract")
    ifcopenshell.api.spatial.assign_container(m, products=[target, other], relating_structure=floor)
    typ = ifcopenshell.api.root.create_entity(m, ifc_class="IfcFurnitureType", name="Type")
    ifcopenshell.api.type.assign_type(m, related_objects=[target, other], relating_type=typ, should_map_representations=False)
    matrix = np.array([[0, 1, 0, 5.897], [-1, 0, 0, .896], [0, 0, 1, .048], [0, 0, 0, 1.]])
    ifcopenshell.api.geometry.edit_object_placement(m, product=target, matrix=matrix)
    def rep(identifier):
        pts = [m.create_entity("IfcCartesianPoint", Coordinates=p) for p in [(0., 0., 0.), (100., 0., 0.)]]
        line = m.create_entity("IfcPolyline", Points=pts)
        color = m.create_entity("IfcColourRgb", Red=0., Green=0., Blue=0.)
        style = m.create_entity("IfcCurveStyle", CurveColour=color)
        m.create_entity("IfcStyledItem", Item=line, Styles=[style])
        return m.create_entity("IfcShapeRepresentation", ContextOfItems=context,
                               RepresentationIdentifier=identifier, RepresentationType="Curve3D", Items=[line])
    target.Representation = m.create_entity("IfcProductDefinitionShape", Representations=[rep("Body")])
    formal = ifcopenshell.file.from_string(m.to_string())
    target.Representation.Representations += (rep("ApprovedPlan"),)
    prop = ifcopenshell.api.pset.add_pset(m, product=target, name="ReviewSource")
    ifcopenshell.api.pset.edit_pset(m, pset=prop, properties={"Approved": True, "Source": "geometry_derived"})
    ann = ifcopenshell.api.root.create_entity(m, ifc_class="IfcAnnotation", name="ApprovedPlan")
    ann.ObjectPlacement = target.ObjectPlacement
    ann.Representation = m.create_entity("IfcProductDefinitionShape", Representations=[rep("Annotation")])
    return m, formal, target.GlobalId, ann.GlobalId


class ProductPackageTests(unittest.TestCase):
    def test_roundtrip_keeps_identity_coordinates_and_dependencies(self):
        src, _, guid, ann = fixture()
        package = pkg.extract(src, guid, [ann])
        package = ifcopenshell.file.from_string(package.to_string())
        pkg.validate(src, package, guid, [ann])
        self.assertEqual(len(package.by_type("IfcElement")), 1)
        self.assertEqual(len(package.by_type("IfcStyledItem")), 3)
        self.assertEqual(pkg.body_fingerprint(src.by_guid(guid)), pkg.body_fingerprint(package.by_guid(guid)))

    def test_attach_idempotent_without_duplicate_product(self):
        src, formal, guid, ann = fixture()
        package = pkg.extract(src, guid, [ann])
        count = len(formal.by_type("IfcElement"))
        pkg.attach(formal, package, guid, [ann])
        entities = len(list(formal))
        pkg.attach(formal, package, guid, [ann])
        self.assertEqual(len(list(formal)), entities)
        self.assertEqual(len(formal.by_type("IfcElement")), count)
        self.assertEqual(len(formal.by_type("IfcAnnotation")), 1)

    def test_coordinate_mismatch_fails_closed(self):
        src, formal, guid, ann = fixture()
        package = pkg.extract(src, guid, [ann])
        formal.by_guid(guid).ObjectPlacement.RelativePlacement.Location.Coordinates = (0., 0., 0.)
        with self.assertRaisesRegex(AssertionError, "Coordinate mismatch"):
            pkg.attach(formal, package, guid, [ann])

    def test_body_mismatch_fails_closed(self):
        src, formal, guid, ann = fixture()
        package = pkg.extract(src, guid, [ann])
        formal.by_guid(guid).Representation.Representations[0].Items[0].Points[1].Coordinates = (120., 0., 0.)
        with self.assertRaisesRegex(AssertionError, "Body differs"):
            pkg.attach(formal, package, guid, [ann])

    def test_missing_target_never_appended(self):
        src, formal, guid, ann = fixture()
        package = pkg.extract(src, guid, [ann])
        formal.by_guid(guid).GlobalId = ifcopenshell.guid.new()
        with self.assertRaises(RuntimeError):
            pkg.attach(formal, package, guid, [ann])

    def test_units_mismatch_fails_closed(self):
        src, formal, guid, ann = fixture()
        package = pkg.extract(src, guid, [ann])
        for u in formal.by_type("IfcSIUnit"):
            if u.UnitType == "LENGTHUNIT":
                u.Prefix = None
        with self.assertRaises(AssertionError):
            pkg.attach(formal, package, guid, [ann])

    def test_nested_components_require_explicit_policy(self):
        src, _, guid, ann = fixture()
        child = ifcopenshell.api.root.create_entity(src, ifc_class="IfcFurniture", name="Nested child")
        ifcopenshell.api.aggregate.assign_object(src, products=[child], relating_object=src.by_guid(guid))
        with self.assertRaisesRegex(ValueError, "Nested parts"):
            pkg.extract(src, guid, [ann])

    def test_missing_plan_context_repair_preserves_effective_frame(self):
        src, _, guid, ann = fixture()
        plan = src.create_entity("IfcGeometricRepresentationContext", ContextType="Plan", CoordinateSpaceDimension=2, Precision=1e-5)
        src.by_guid(ann).Representation.Representations[0].ContextOfItems = plan
        project = src.by_type("IfcProject")[0]
        project.RepresentationContexts += (plan,)
        package = pkg.extract(src, guid, [ann])
        repairs = pkg.repair_missing_metadata(package)
        self.assertTrue(any(r['attribute'] == 'WorldCoordinateSystem' for r in repairs))
        pkg.validate(src, package, guid, [ann])


if __name__ == "__main__":
    unittest.main()
