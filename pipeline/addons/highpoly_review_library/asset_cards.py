"""Transient image/identity cards for Blender's native drag grid; no model cache."""
import bpy
from bpy.app.handlers import persistent

OWNER = "highpoly_review_library.card.v1"
KEY = "review_library_card"
_owned = set()


def product_for_id(data, library):
    if data is None or not isinstance(data, bpy.types.Collection):
        return None
    if data.get(KEY) != OWNER or data.as_pointer() not in _owned:
        return None
    if data.get("catalog") != str(library._catalog_path):
        return None
    return library.entry(data.get("product_id", ""))


def clear():
    for collection in list(bpy.data.collections):
        if collection.get(KEY) == OWNER:
            # A user may have linked/edited a card through another Blender UI.
            # Never delete collections or objects that have gained real content.
            if collection.objects or collection.children or collection.users > int(collection.use_fake_user):
                collection.asset_clear()
                collection.pop(KEY, None)
            else:
                bpy.data.collections.remove(collection)
    _owned.clear()


def rebuild(library):
    clear()
    for product in library._catalog["products"]:
        card = bpy.data.collections.new("IFC · " + product["name"])
        card[KEY] = OWNER
        card["product_id"] = product["id"]
        card["catalog"] = str(library._catalog_path)
        card.asset_mark()
        card.asset_data.description = product["name"] + "\n" + product["approval_label"] + " · 从 IFC 拖入；Ctrl+S 保存"
        card.asset_data.tags.new(product["category_label"])
        with bpy.context.temp_override(id=card):
            bpy.ops.ed.lib_id_load_custom_preview(filepath=str(library.resolve_path(product["previews"]["iso"])))
        _owned.add(card.as_pointer())


@persistent
def after_load(_):
    # IFC fresh-session loading may reset all Blender IDs. Recreate metadata only.
    def refresh():
        import sys
        library = sys.modules.get(__package__)
        if library and hasattr(bpy.types.Scene, "review_library"):
            rebuild(library)
        return None
    bpy.app.timers.register(refresh, first_interval=.1)


def register(library):
    rebuild(library)
    bpy.app.handlers.load_post.append(after_load)
    bpy.app.handlers.undo_post.append(after_load)
    bpy.app.handlers.redo_post.append(after_load)


def unregister():
    for handlers in (bpy.app.handlers.load_post, bpy.app.handlers.undo_post, bpy.app.handlers.redo_post):
        if after_load in handlers:
            handlers.remove(after_load)
    clear()
