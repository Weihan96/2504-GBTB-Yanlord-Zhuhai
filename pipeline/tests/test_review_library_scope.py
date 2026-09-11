"""Inclusion in a verified library must never silently approve a scene."""
import copy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import build_approved_product_library as library


class ApprovalScopeTests(unittest.TestCase):
    def setUp(self):
        self.scope = json.loads((library.OUT / "acceptance-scope.json").read_text())
        self.authority = json.loads((library.OUT / "migration-authorization.json").read_text())

    def test_explicit_scope_is_fifteen_with_two_tiers(self):
        products = library.authorized_scope(self.scope, self.authority)
        self.assertEqual(len(products), 15)
        self.assertEqual(len(library.ACCEPTED), 8)
        self.assertEqual(len(library.SINGLE_ACCEPTED), 7)

    def test_each_tier_has_unambiguous_label(self):
        for key in library.ACCEPTED:
            self.assertEqual(library.approval_fields(key)["scene_approval_status"], "approved")
        for key in library.SINGLE_ACCEPTED:
            p = library.approval_fields(key)
            self.assertEqual(p["scene_approval_status"], "pending")
            self.assertEqual(p["approval_label"], "单品已通过、场景待验收")

    def test_reject_unapproved_scene_promotion(self):
        p = {"id": "bed02", **library.approval_fields("bed02")}
        p["scene_approval_status"] = "approved"
        with self.assertRaises(AssertionError):
            library.validate_approval_fields(p)

    def test_reject_unauthorized_inclusion(self):
        for field, value in (("library_inclusion_authorized", False), ("pure_product_ifc_write_authorized", False),
                             ("formal_ifc_write_authorized", True), ("scene_approval_status", "approved")):
            authority = copy.deepcopy(self.authority)
            authority[field] = value
            with self.subTest(field=field), self.assertRaises(AssertionError):
                library.authorized_scope(self.scope, authority)

    def test_no_partial_or_skipped_products_leak(self):
        authority = copy.deepcopy(self.authority)
        authority["products"].append("marilyn-01")
        with self.assertRaises(AssertionError):
            library.authorized_scope(self.scope, authority)
        with self.assertRaises(AssertionError):
            library.approval_fields("marilyn-01")

    def test_support_is_not_whole_street_sink(self):
        self.assertIn("支撑子构件", library.SINGLE_ACCEPTED["street-h"])
        self.assertNotIn("street", library.authorized_scope(self.scope, self.authority))


if __name__ == "__main__":
    unittest.main()
