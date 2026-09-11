"""In-memory IFC insertion. Source identity is provenance, never instance identity."""
import hashlib
from pathlib import Path
import numpy as np
import ifcopenshell
import ifcopenshell.api.project
import ifcopenshell.api.geometry
import ifcopenshell.api.spatial
import ifcopenshell.api.pset
import ifcopenshell.util.unit

# IfcOpenShell 0.8.4 append_asset can retain borrowed source data. Like Bonsai's
# IfcStore.library_file, keep libraries alive for the destination file lifetime.
# Store on the destination itself, also surviving an add-on module reload.


def prepare_source(path, expected_hash, source_guid, target, require_drawings=True):
    with open(path, "rb") as stream:
        if hashlib.sha256(stream.read()).hexdigest() != expected_hash:
            raise ValueError("单品 IFC 哈希不匹配，停止载入")
    source = ifcopenshell.open(str(path))
    if source.schema != target.schema:
        raise ValueError("单品与项目 IFC schema 不一致，请先转换项目")
    product = source.by_guid(source_guid)
    from .body_views import views, validate
    if views(product):
        validate(product)
    if len([e for e in source.by_type("IfcElement") if not e.is_a("IfcOpeningElement")]) != 1 or not product.is_a("IfcElement"):
        raise ValueError("只允许已验证的单品 IFC")
    assert {r.RelatedOpeningElement for r in getattr(product, "HasOpenings", ())} == set(source.by_type("IfcOpeningElement")), "开孔依赖不属于目标单品"
    if require_drawings and len(product.Representation.Representations) < 4:
        raise ValueError("单品缺少 Body 或已批准的三视图")
    if not require_drawings and any(r.RepresentationIdentifier != "Body" for r in product.Representation.Representations):
        raise ValueError("待审候选只允许接入 Body，不能夹带未批准的二维表达")
    source_scale = ifcopenshell.util.unit.calculate_unit_scale(source)
    target_scale = ifcopenshell.util.unit.calculate_unit_scale(target)
    if source_scale != target_scale:
        unit = ifcopenshell.util.unit.get_project_unit(target, "LENGTHUNIT")
        name = (getattr(unit, "Prefix", None) or "") + unit.Name
        source = ifcopenshell.util.unit.convert_file_length_units(source, name)
        product = source.by_guid(source_guid)
    # Library-relative style dependencies must not be resolved against the
    # destination project directory. Only this private source copy is changed.
    for style in source.by_type("IfcExternallyDefinedSurfaceStyle"):
        location = style.Location
        if location and "://" not in location and not Path(location).is_absolute():
            resolved = (Path(path).parent / location).resolve()
            if not resolved.is_file():
                raise ValueError(f"缺少单品材质依赖：{resolved}")
            style.Location = str(resolved)
    # The disk package keeps all original IDs. Only this private in-memory copy
    # gets fresh roots, preventing GUID/type/property collisions on repeat drops.
    for root in source.by_type("IfcRoot"):
        root.GlobalId = ifcopenshell.guid.new()
    return source, product


def insert_prepared(target, source, product, placement_project_m, container, provenance):
    """Caller owns the Bonsai transaction. This function never writes a file."""
    if container is None or not container.is_a("IfcSpatialElement"):
        raise ValueError("请先在 Bonsai 选择默认空间容器（楼层）")
    matrix = np.asarray(placement_project_m, dtype=float)
    if matrix.shape != (4, 4) or not np.isfinite(matrix).all():
        raise ValueError("无效的项目坐标变换")
    target.__dict__.setdefault("_review_library_sources", []).append(source)
    element = ifcopenshell.api.project.append_asset(
        target, library=source, element=product, assume_asset_uniqueness_by_name=False)
    ifcopenshell.api.spatial.assign_container(target, products=[element], relating_structure=container)
    ifcopenshell.api.geometry.edit_object_placement(target, product=element, matrix=matrix, is_si=True)
    pset = ifcopenshell.api.pset.add_pset(target, product=element, name="ReviewLibrarySource")
    ifcopenshell.api.pset.edit_pset(target, pset=pset, properties=provenance)
    return element
