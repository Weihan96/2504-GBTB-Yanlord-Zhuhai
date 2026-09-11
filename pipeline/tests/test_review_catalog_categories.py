import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("rules", ROOT / "pipeline/addons/highpoly_review_library/catalog_rules.py")
rules = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rules)


class ReviewCategories(unittest.TestCase):
    def test_full_requires_both_approvals(self):
        self.assertEqual(rules.classify("approved", "approved"), "approved")
        self.assertEqual(rules.classify("approved", "pending"), "partial")

    def test_partial_and_pending_do_not_promote(self):
        self.assertEqual(rules.classify("partially_approved", "pending"), "partial")
        self.assertEqual(rules.classify("pending", "approved"), "pending")
        self.assertEqual(rules.classify("pending", "pending"), "pending")

    def test_skipped_is_explicit(self):
        self.assertEqual(rules.classify("pending", "pending", skipped=True), "skipped")

    def test_unapproved_cannot_claim_approved_payload(self):
        p = dict(single_product_approval_status="pending", scene_approval_status="pending",
                 review_category="pending", category_label="未验收", package_validation="pass",
                 insertion_content="approved_3d_2d")
        with self.assertRaises(AssertionError):
            rules.validate_entry(p, Path, lambda _: "")


if __name__ == "__main__":
    unittest.main()
