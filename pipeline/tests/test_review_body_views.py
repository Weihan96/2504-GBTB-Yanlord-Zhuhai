import importlib.util
from pathlib import Path
import unittest
import numpy as np
import ifcopenshell
import ifcopenshell.api.project
import ifcopenshell.api.root
import ifcopenshell.api.unit

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("body_views", ROOT / "pipeline/addons/highpoly_review_library/body_views.py")
views = importlib.util.module_from_spec(spec)
spec.loader.exec_module(views)


class DirectionTests(unittest.TestCase):
    def setUp(self):
        self.roles = {"Plan": (None, np.array([0,0,1.])), "Front": (None, np.array([0,-1.,0])), "Side": (None, np.array([-1.,0,0]))}

    def test_front_side_not_overlaid(self):
        self.assertEqual(views.choose_role(self.roles, [0,-1,0], "ELEVATION_VIEW"), "Front")
        self.assertEqual(views.choose_role(self.roles, [-1,0,0], "ELEVATION_VIEW"), "Side")

    def test_opposite_uses_same_plane(self):
        self.assertEqual(views.choose_role(self.roles, [0,1,0], "ELEVATION_VIEW"), "Front")

    def test_plan_does_not_select_front(self):
        self.assertEqual(views.choose_role(self.roles, [0,0,1], "PLAN_VIEW"), "Plan")
        self.assertIsNone(views.choose_role(self.roles, [0,-1,0], "PLAN_VIEW"))

    def test_oblique_and_section_keep_model_body(self):
        self.assertIsNone(views.choose_role(self.roles, [1,1,0], "ELEVATION_VIEW"))
        self.assertIsNone(views.choose_role(self.roles, [0,-1,0], "SECTION_VIEW"))

    def test_rotated_instance_direction(self):
        rotation = np.array([[0,-1,0],[1,0,0],[0,0,1.]])
        self.assertEqual(views.choose_role(self.roles, np.linalg.solve(rotation,[1,0,0]), "ELEVATION_VIEW"), "Front")

    def test_invalid_camera_is_rejected(self):
        with self.assertRaises(ValueError):
            views.choose_role(self.roles,[0,0,0],"PLAN_VIEW")

    def test_shape_aspects_survive_append_context_merge_and_reload(self):
        path = ROOT / "output/review/approved-product-library/runtime/ifc/gessi316-54294.ifc"
        source = ifcopenshell.open(str(path))
        product = source.by_type("IfcElement")[0]
        before = views.validate(product)
        target = ifcopenshell.file(schema="IFC4")
        ifcopenshell.api.root.create_entity(target, ifc_class="IfcProject")
        ifcopenshell.api.unit.assign_unit(target)
        appended = ifcopenshell.api.project.append_asset(target, library=source, element=product)
        # Keep source alive while serializing (IfcOpenShell 0.8.4 borrowed data).
        result = views.validate(appended)
        for role in views.ROLES:
            self.assertTrue(np.array_equal(result[role][1],before[role][1]))
            self.assertEqual(len(result[role][0].Items),len(before[role][0].Items))
        reopened = ifcopenshell.file.from_string(target.to_string())
        views.validate(reopened.by_guid(appended.GlobalId))


if __name__ == "__main__":
    unittest.main()
