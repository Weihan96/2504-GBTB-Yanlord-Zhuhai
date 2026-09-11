import copy
import importlib.util
import json
from pathlib import Path
import hashlib
import struct
import unittest

ROOT = Path(__file__).resolve().parents[2]
RUNTIME = ROOT / 'output/review/approved-product-library/runtime'
spec = importlib.util.spec_from_file_location('preview_rules', ROOT / 'pipeline/addons/highpoly_review_library/catalog_rules.py')
rules = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rules)


class DrawingPreviews(unittest.TestCase):
    def setUp(self):
        self.product = {'ifc_sha256': 'ifc', 'previews': {v: v + '.png' for v in ('plan', 'front', 'side', 'iso')},
                        'preview_evidence': {}}
        for view in ('plan', 'front', 'side'):
            self.product['preview_evidence'][view] = dict(kind='bonsai_body_linework', role=view.title(),
                source_ifc_sha256='ifc', svg=view+'.svg', svg_sha256=view+'.svg', png_sha256=view+'.png')
        self.product['preview_evidence']['iso'] = dict(kind='camera_render_3d', png_sha256='iso.png')

    def validate(self):
        rules.validate_previews(self.product, Path, str)

    def test_directional_linework_and_separate_3d(self):
        self.validate()

    def test_camera_render_cannot_replace_drawing(self):
        self.product['preview_evidence']['plan']['kind'] = 'camera_render_3d'
        with self.assertRaises(AssertionError): self.validate()

    def test_stale_ifc_rejected(self):
        self.product['ifc_sha256'] = 'changed'
        with self.assertRaises(AssertionError): self.validate()

    def test_swapped_front_side_rejected(self):
        self.product['preview_evidence']['front'] = copy.deepcopy(self.product['preview_evidence']['side'])
        with self.assertRaises(AssertionError): self.validate()

    def test_changed_image_rejected(self):
        self.product['previews']['plan'] = 'camera.png'
        with self.assertRaises(AssertionError): self.validate()

    def test_entire_runtime_has_traceable_linework(self):
        catalog = json.loads((RUNTIME / 'catalog.json').read_text())
        self.assertEqual(catalog['preview_contract'], 'plan_front_side_bonsai_linework_iso_camera_render')
        self.assertEqual(len(catalog['products']), 39)
        for product in catalog['products']:
            with self.subTest(product=product['id']):
                rules.validate_previews(product, lambda p: RUNTIME / p,
                    lambda p: hashlib.sha256(p.read_bytes()).hexdigest())

    def test_3d_rasters_are_bounded_thumbnails(self):
        for p in json.loads((RUNTIME / 'catalog.json').read_text())['products']:
            with self.subTest(product=p['id']):
                image = (RUNTIME / p['previews']['iso']).read_bytes()
                self.assertEqual(image[:8], b'\x89PNG\r\n\x1a\n')
                width, height = struct.unpack('>II', image[16:24])
                self.assertLessEqual(max(width, height), 384)
                record = p['preview_evidence']['iso']
                self.assertEqual(record['purpose'], 'asset_thumbnail')
                self.assertEqual((record['width'], record['height']), (width, height))


if __name__ == '__main__':
    unittest.main()
