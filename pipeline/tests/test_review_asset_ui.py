import ast
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import sys
import unittest
from unittest.mock import patch

ADDON = Path(__file__).resolve().parents[1] / 'addons/highpoly_review_library'


class Layout:
    def __init__(self): self.calls = []
    def row(self, **kw): return self
    def column(self, **kw): return self
    def label(self, **kw): self.calls.append(('label', kw))
    def separator(self): pass
    def prop(self, *args, **kw): self.calls.append(('prop', kw))
    def operator(self, op, **kw):
        self.calls.append(('operator', op))
        props = SimpleNamespace()
        self.calls.append(('operator_props', props))
        return props
    def template_asset_view(self, *args, **kw): self.calls.append(('grid', kw))
    def template_icon(self, **kw): self.calls.append(('image', kw))


class AssetPopup(unittest.TestCase):
    def setUp(self):
        bpy = SimpleNamespace(types=SimpleNamespace(Operator=object, Panel=object))
        with patch.dict(sys.modules, {'bpy': bpy}):
            spec = importlib.util.spec_from_file_location('asset_ui', ADDON / 'asset_ui.py')
            self.ui = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(self.ui)
        self.product = dict(id='chair', name='A long product name that must remain readable in a popup',
            source_label='基于原始高模几何生成的简化图纸表达', scene_approval_status='pending',
            approval_label='尚未验收', category_label='未验收', insertion_content='review_body_views')
        self.lib = SimpleNamespace(_catalog={'products': [self.product]}, _preview={}, entry=lambda key: self.product)
        self.ui.library = lambda: self.lib
        class Scene(dict): pass
        scene = Scene()
        scene.review_library = SimpleNamespace(product='chair', show_catalog_settings=False)
        self.context = SimpleNamespace(scene=scene)

    def test_main_panel_no_long_details_or_preview_grid(self):
        panel = self.ui.REVIEWLIB_PT_Library()
        panel.layout = Layout()
        panel.draw(self.context)
        calls = panel.layout.calls
        self.assertTrue(any(c[0] == 'grid' for c in calls))
        self.assertIn(('operator', 'review_library.details'), calls)
        self.assertNotIn(('operator', 'review_library.image'), calls)
        self.assertFalse(any(c[0] == 'image' for c in calls))
        self.assertFalse(any(c[0] == 'label' and c[1]['text'] == self.product['source_label'] for c in calls))
        self.assertEqual(panel.bl_category, 'Assets')

    def test_popup_contains_only_three_drawing_buttons(self):
        layout = Layout()
        self.ui.draw_details(layout, self.product, self.lib)
        self.assertEqual(layout.calls.count(('operator', 'review_library.image')), 3)
        self.assertEqual([c[1].view for c in layout.calls if c[0] == 'operator_props'], ['plan', 'front', 'side'])
        self.assertIn(('label', {'text': '尚未验收', 'icon': 'TIME'}), layout.calls)

    def test_popup_uses_native_popup_not_dialog_confirmation(self):
        result = []
        self.context.window_manager = SimpleNamespace(invoke_popup=lambda op, **kw: result.append(kw))
        self.ui.REVIEWLIB_OT_Details().invoke(self.context, None)
        self.assertEqual(result, [{'width': 360}])

    def test_loading_draws_only_filled_triangles_and_restores_gpu_state(self):
        source = ast.parse((ADDON / 'loading_feedback.py').read_text())
        cls = next(n for n in source.body if isinstance(n, ast.ClassDef))
        method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == 'draw')
        calls = []
        state = {'depth': 'LESS_EQUAL', 'blend': 'NONE', 'mask': True}
        gpu = SimpleNamespace(state=SimpleNamespace(
            depth_test_get=lambda: state['depth'], blend_get=lambda: state['blend'], depth_mask_get=lambda: state['mask'],
            depth_test_set=lambda v: state.update(depth=v), blend_set=lambda v: state.update(blend=v), depth_mask_set=lambda v: state.update(mask=v)))
        scope = {'gpu': gpu, 'bpy': SimpleNamespace(context=SimpleNamespace(area='view')),
            'batch_for_shader': lambda shader, mode, *a, **kw: SimpleNamespace(draw=lambda s: calls.append(mode))}
        exec(compile(ast.Module(body=[method], type_ignores=[]), '<draw>', 'exec'), scope)
        obj = SimpleNamespace(preview=SimpleNamespace(point=(0,0,0), area='view',
            shader=SimpleNamespace(bind=lambda: None, uniform_float=lambda *a: None)), fill_vertices=lambda: [(0,0,0)]*8)
        scope['draw'](obj)
        self.assertEqual(calls, ['TRIS'])
        self.assertEqual(state, {'depth': 'LESS_EQUAL', 'blend': 'NONE', 'mask': True})

    def test_wireframe_preserved_during_loading_and_drag_does_not_open_popup(self):
        native = (ADDON / 'native_assets.py').read_text()
        loading = (ADDON / 'loading_feedback.py').read_text()
        self.assertIn('self.outline_visible = True', native)
        self.assertIn('if not self.outline_visible', native)
        self.assertNotIn('self.preview.outline_visible = False', loading)
        drag = native.split('class REVIEWLIB_OT_DragIFC')[1].split('class REVIEWLIB_OT_SelectAsset')[0]
        self.assertNotIn('open_details', drag)

    def test_loading_enter_keeps_fixed_blue_bounds(self):
        source = ast.parse((ADDON / 'loading_feedback.py').read_text())
        cls = next(n for n in source.body if isinstance(n, ast.ClassDef))
        method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == '__enter__')
        scope = {'_active': None, 'bpy': SimpleNamespace(types=SimpleNamespace(
            SpaceView3D=SimpleNamespace(draw_handler_add=lambda *args: 'fill-handler')))}
        exec(compile(ast.Module(body=[method], type_ignores=[]), '<enter>', 'exec'), scope)
        preview = SimpleNamespace(outline_visible=True, handle='blue-handler')
        obj = SimpleNamespace(preview=preview, draw=lambda: None)
        self.assertIs(scope['__enter__'](obj), obj)
        self.assertTrue(preview.outline_visible)
        self.assertEqual(preview.handle, 'blue-handler')
        self.assertEqual(obj.handle, 'fill-handler')

    def test_asset_press_does_not_steal_drag_but_release_opens_details(self):
        source = ast.parse((ADDON / 'native_assets.py').read_text())
        cls = next(n for n in source.body if isinstance(n, ast.ClassDef)
                   and n.name == 'REVIEWLIB_OT_SelectAsset')
        calls = []
        scope = {'bpy': SimpleNamespace(types=SimpleNamespace(Operator=object)),
                 '__package__': 'drag_test',
                 'asset_product': lambda ctx: {'id': 'chair'}}
        exec(compile(ast.Module(body=[cls], type_ignores=[]), '<select>', 'exec'), scope)
        operator = scope['REVIEWLIB_OT_SelectAsset']()
        popup = SimpleNamespace(open_details=lambda ctx: calls.append('popup'))
        with patch.dict(sys.modules, {'drag_test.asset_ui': popup}):
            for value in ('PRESS', 'CLICK_DRAG', 'DOUBLE_CLICK'):
                result = operator.invoke(self.context, SimpleNamespace(type='LEFTMOUSE', value=value))
                self.assertEqual(result, {'FINISHED'})
                self.assertEqual(calls, [])
            operator.invoke(self.context, SimpleNamespace(type='LEFTMOUSE', value='RELEASE'))
            self.assertEqual(calls, ['popup'])
            operator.execute(self.context)
            self.assertEqual(calls, ['popup'])
        self.assertEqual(self.context.scene.review_library.product, 'chair')


if __name__ == '__main__': unittest.main()
