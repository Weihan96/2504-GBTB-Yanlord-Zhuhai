"""Transient viewport feedback, driven by completed import stages, never time."""
import itertools
import bpy
import gpu
from mathutils import Vector
from gpu_extras.batch import batch_for_shader

_active = None


def notify(product_id, fraction, label):
    if _active is not None and _active.product_id == product_id:
        _active.advance(fraction, label)


class LoadingFeedback:
    def __init__(self, preview, product_id, context):
        self.preview = preview
        self.product_id = product_id
        self.context = context
        self.fraction = 0.0
        self.history = []
        self.handle = None

    def __enter__(self):
        global _active
        if _active is not None:
            raise RuntimeError("已有单品正在加载")
        # Keep the fixed blue placement bounds during loading. Only the
        # rising green fill is borderless; DragIFC.finish closes both overlays.
        self.handle = bpy.types.SpaceView3D.draw_handler_add(self.draw, (), "WINDOW", "POST_VIEW")
        _active = self
        return self

    def advance(self, fraction, label):
        if not self.fraction <= fraction <= 1.0:
            raise ValueError("加载阶段必须单调递增且不超过完成状态")
        self.fraction = fraction
        self.history.append((fraction, label))
        self.context.workspace.status_text_set(label + " · 加载阶段进度 · 完成后 Ctrl+S 保存")
        self.redraw()

    def redraw(self):
        if bpy.app.background:
            return
        self.preview.area.tag_redraw()
        # Bonsai mutation remains one synchronous, undoable transaction. Force
        # a draw at stage boundaries; never run scene edits on a worker thread.
        with self.context.temp_override(area=self.preview.area):
            bpy.ops.wm.redraw_timer(type="DRAW_WIN_SWAP", iterations=2)

    def fill_vertices(self):
        low, high = self.preview.bounds
        top = low[2] + (high[2] - low[2]) * self.fraction
        return [Vector((high[0] if x else low[0], high[1] if y else low[1], top if z else low[2]))
                + self.preview.point for x, y, z in itertools.product((0, 1), repeat=3)]

    def draw(self):
        if self.preview.point is None or bpy.context.area != self.preview.area:
            return
        vertices = self.fill_vertices()
        triangles = [(0, 2, 3), (0, 3, 1), (4, 5, 7), (4, 7, 6),
                     (0, 1, 5), (0, 5, 4), (2, 6, 7), (2, 7, 3),
                     (0, 4, 6), (0, 6, 2), (1, 3, 7), (1, 7, 5)]
        shader = self.preview.shader
        depth, blend = gpu.state.depth_test_get(), gpu.state.blend_get()
        mask = gpu.state.depth_mask_get()
        try:
            gpu.state.depth_test_set("NONE")
            gpu.state.depth_mask_set(False)
            gpu.state.blend_set("ALPHA")
            shader.bind()
            shader.uniform_float("color", (0.08, 0.85, 0.3, 0.16))
            batch_for_shader(shader, "TRIS", {"pos": vertices}, indices=triangles).draw(shader)
        finally:
            gpu.state.depth_mask_set(mask)
            gpu.state.depth_test_set(depth)
            gpu.state.blend_set(blend)

    def __exit__(self, *exc):
        global _active
        if self.handle is not None:
            bpy.types.SpaceView3D.draw_handler_remove(self.handle, "WINDOW")
            self.handle = None
        _active = None
        self.context.workspace.status_text_set(None)
        self.redraw()
