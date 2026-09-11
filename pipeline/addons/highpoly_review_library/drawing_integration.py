"""Select one Body per library product inside native Bonsai Create Drawing.

No IFC mutation, temporary annotation, extra save operator or provider patch.
Other products keep the installed Bonsai context selection unchanged.
"""
import numpy as np
import ifcopenshell.geom
import ifcopenshell.util.placement
from .body_views import views, validate, choose_role

_original = None
_original_bisect = None
selection_log = []
bisect_log = []


def serialize(self, ifc, tree, contexts, context_type, drawing_elements, target_view):
    if context_type == "body":
        # print_all can reuse this operator for another camera. Cut filtering
        # belongs only to the Body pass of the current drawing.
        self._review_selected_bodies = []
    directional = {e: views(e) for e in drawing_elements if getattr(e, "Representation", None)}
    directional = {e: value for e, value in directional.items() if value}
    _original(self, ifc, tree, contexts, context_type, set(drawing_elements) - set(directional), target_view)
    # An independent annotation pass must never draw a library Body again.
    if context_type != "body":
        return
    camera_matrix = ifcopenshell.util.placement.get_local_placement(self.camera_element.ObjectPlacement)
    for product, roles in directional.items():
        validate(product)
        matrix = ifcopenshell.util.placement.get_local_placement(product.ObjectPlacement)
        local_normal = np.linalg.solve(matrix[:3, :3], camera_matrix[:3, 2])
        role = choose_role(roles, local_normal, target_view)
        rep = roles[role][0] if role else next(r for r in product.Representation.Representations if r.ContextOfItems.TargetView == "MODEL_VIEW")
        settings = ifcopenshell.geom.settings()
        settings.set("dimensionality", ifcopenshell.ifcopenshell_wrapper.CURVES_SURFACES_AND_SOLIDS)
        settings.set("iterator-output", ifcopenshell.ifcopenshell_wrapper.NATIVE)
        if target_view in ("PLAN_VIEW", "REFLECTED_PLAN_VIEW"):
            settings.set("model-offset", (0., 0., .002 if target_view == "PLAN_VIEW" else -.002))
        shape = ifcopenshell.geom.create_shape(settings, product, rep)
        if shape is None:
            raise RuntimeError(f"Body 出图失败：{product.GlobalId} / {role}")
        self.serialiser.write(shape)
        tree.add_element(shape)
        record = {"global_id": product.GlobalId, "target_view": target_view,
                  "role": role or "Model", "representation_id": rep.id(),
                  "old_body_or_annotation_emitted": False}
        selection_log.append(record)
        if not hasattr(self, "_review_selected_bodies"):
            self._review_selected_bodies = []
        self._review_selected_bodies.append(record)


def bisect(self, context, root):
    # Native BISECT independently walks visible Blender meshes. Exclude only
    # the products whose directional Body was selected for this exact drawing;
    # otherwise their old high-poly mesh would leak back in through cut lines.
    selected = {r["global_id"] for r in getattr(self, "_review_selected_bodies", ()) if r["role"] != "Model"}
    import bonsai.tool as tool
    excluded = [o.name for o in context.visible_objects if tool.Ifc.get_entity(o) and tool.Ifc.get_entity(o).GlobalId in selected]
    bisect_log.append({"drawing_guid": self.camera_element.GlobalId, "excluded_mesh_objects": excluded})
    class FilteredContext:
        def __getattr__(self, name):
            if name == "visible_objects":
                return [o for o in context.visible_objects if not (tool.Ifc.get_entity(o) and tool.Ifc.get_entity(o).GlobalId in selected)]
            return getattr(context, name)
    return _original_bisect(self, FilteredContext(), root)


def register():
    global _original, _original_bisect
    from bonsai.bim.module.drawing.operator import CreateDrawing
    if _original is not None:
        return
    import inspect
    expected = ("self", "ifc", "tree", "contexts", "context_type", "drawing_elements", "target_view")
    if tuple(inspect.signature(CreateDrawing.serialize_contexts_elements).parameters) != expected:
        raise RuntimeError("Bonsai 出图接口已变化，停止安装单品 Body 适配")
    _original = CreateDrawing.serialize_contexts_elements
    _original_bisect = CreateDrawing.generate_bisect_linework
    CreateDrawing.serialize_contexts_elements = serialize
    CreateDrawing.generate_bisect_linework = bisect


def unregister():
    global _original, _original_bisect
    if _original is not None:
        from bonsai.bim.module.drawing.operator import CreateDrawing
        if CreateDrawing.serialize_contexts_elements is serialize:
            CreateDrawing.serialize_contexts_elements = _original
        if CreateDrawing.generate_bisect_linework is bisect:
            CreateDrawing.generate_bisect_linework = _original_bisect
        _original = None
        _original_bisect = None
