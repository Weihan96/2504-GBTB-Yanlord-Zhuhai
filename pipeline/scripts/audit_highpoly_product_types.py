#!/usr/bin/env python3
"""Inventory high-complexity IFC product types without rendering the whole model."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import ifcopenshell

from falper_sorgente_linework import EXPECTED, ROOT, relative, sha256, write_json


DEFAULT_INPUT = ROOT / "2504 GBTB Yanlord Zhuhai.ifc"
DEFAULT_OUTPUT = ROOT / "pipeline/decisions/highpoly-product-inventory.json"
PRODUCT_CLASSES = {
    "IfcBuildingElementProxy",
    "IfcElectricAppliance",
    "IfcFurniture",
    "IfcSanitaryTerminal",
    "IfcWasteTerminal",
}
SYSTEM_CLASSES = {"IfcFlowSegment", "IfcPipeSegment"}
MANUAL_FOLDERS = {
    "BS01": "falper-sorgente",
    "Geberit 146.140": "geberit-146-140",
    "Geberit Duofix Sigma": "geberit-duofix-sigma-224-212",
    "Libelle": "libelle",
}


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return slug or "unnamed-product"


def product_type(product) -> tuple[str, str | None, str]:
    relations = tuple(getattr(product, "IsTypedBy", ()) or ())
    if relations:
        item = relations[0].RelatingType
        return item.Name or item.GlobalId, item.GlobalId, item.is_a()
    return product.Name or product.GlobalId, None, product.is_a()


def item_complexity(item, seen: set[int]) -> tuple[int, int, set[str]]:
    if item is None or item.id() in seen:
        return 0, 0, set()
    seen.add(item.id())
    kinds = {item.is_a()}
    if item.is_a("IfcMappedItem"):
        vertices = faces = 0
        for child in item.MappingSource.MappedRepresentation.Items:
            child_vertices, child_faces, child_kinds = item_complexity(child, seen)
            vertices += child_vertices
            faces += child_faces
            kinds.update(child_kinds)
        return vertices, faces, kinds
    if item.is_a("IfcTriangulatedFaceSet"):
        return len(item.Coordinates.CoordList), len(item.CoordIndex), kinds
    if item.is_a("IfcPolygonalFaceSet"):
        return len(item.Coordinates.CoordList), len(item.Faces), kinds
    if item.is_a("IfcFacetedBrep"):
        points: set[int] = set()
        faces = list(item.Outer.CfsFaces)
        for face in faces:
            for bound in face.Bounds:
                points.update(point.id() for point in bound.Bound.Polygon)
        return len(points), len(faces), kinds
    if item.is_a("IfcShellBasedSurfaceModel"):
        points: set[int] = set()
        faces = []
        for shell in item.SbsmBoundary:
            faces.extend(shell.CfsFaces)
            for face in shell.CfsFaces:
                for bound in face.Bounds:
                    points.update(point.id() for point in bound.Bound.Polygon)
        return len(points), len(faces), kinds
    return 0, 0, kinds


def status_for(name: str) -> tuple[str, str]:
    if name == "BS01":
        return "completed", "official_native_dwg"
    if name == "Geberit 146.140":
        return "review_ready_pending_approval", "official_native_dwg_archived"
    if name == "Geberit 154.446.KS.1":
        return "review_ready_pending_approval", "official_native_dwg_archived"
    if name == "Geberit Duofix Sigma":
        return "review_ready_pending_approval", "exact_project_selected_224_212_00_2_official_native_dwg_archived"
    if name == "Gessi316 54294":
        return "review_ready_pending_approval", "official_identity_archived_exact_cad_area_pro_not_acquired_geometry_derived_proxy"
    if name == "Gessi316 54146":
        return "review_ready_pending_approval", "exact_official_54146_catalogue_identity_archived_area_pro_native_cad_not_acquired_geometry_derived_proxy"
    if name == "Gessi316 54145":
        return "review_ready_pending_approval", "exact_official_54145_wall_mounted_catalogue_identity_archived_area_pro_native_cad_not_acquired_geometry_derived_proxy"
    if name == "HIMA01":
        return "review_ready_pending_approval", "geometry_derived_simplified_proxy_official_dwg_not_acquired"
    if name == "BED02":
        return "review_ready_pending_approval", "official_baxter_identity_archived_exact_cad_login_not_acquired_geometry_derived_proxy"
    if name == "BED01":
        return "review_ready_pending_approval", "official_baxter_identity_archived_manufacturer_dimension_conflict_exact_cad_login_not_acquired_geometry_derived_proxy"
    if name == "sxb010":
        return "review_ready_pending_approval", "owner_confirmed_hunter_douglas_25mm_family_direction_official_pdf_only_geometry_derived_semantic_views"
    if name == "CHA02":
        return "review_ready_pending_approval", "exact_roda_orson_002_identity_reserved_area_cad_not_acquired_geometry_derived_proxy"
    if name == "CHA01":
        return "review_ready_pending_approval", "exact_baxter_colette_armchair_57x60x73_vector_references_archived_native_cad_login_not_acquired_geometry_derived_proxy"
    if name == "SIS04":
        return "review_ready_pending_approval", "exact_molteni_sistema_7_wall_unit_identity_near_name_dwg_excluded_geometry_derived_proxy"
    if name == "Geberit 154.154.00":
        return "review_ready_pending_approval", "exact_154_154_00_1_identity_native_dwg_404_official_eps_identity_only_geometry_derived_proxy"
    if name == "Geberit 154.154.00.1.F":
        return "review_ready_pending_approval", "exact_parent_154_154_00_1_official_page_and_EPS_archived_project_dot_F_flange_is_not_separate_manufacturer_article_native_component_CAD_not_published_geometry_derived_proxy"
    if name == "505 UP V1.LP.S":
        return "review_ready_pending_approval", "official_molteni_505_up_native_family_dwg_archived_exact_project_configuration_not_matched_geometry_derived_proxy"
    if name == "BST03":
        return "review_ready_pending_approval", "exact_baxter_stone_left_drawer_45_identity_vector_pdf_archived_native_cad_login_not_acquired_geometry_derived_proxy"
    if name == "MiamiSoft E09":
        return "review_ready_pending_approval", "exact_baxter_miami_soft_e09_right_terminal_vector_references_archived_native_cad_login_not_acquired_geometry_derived_proxy"
    if name == "MiamiSoft E07":
        return "review_ready_pending_approval", "exact_baxter_miami_soft_e07_left_dormeuse_vector_references_archived_native_cad_login_not_acquired_geometry_derived_proxy"
    if name == "MiamiSoft F03":
        return "review_ready_pending_approval", "exact_baxter_miami_soft_f03_pouf_vector_references_archived_native_cad_login_not_acquired_geometry_derived_proxy"
    if name == "MiamiSoft H01":
        return "review_ready_pending_approval", "exact_baxter_miami_soft_h01_flexible_cushion_vector_references_archived_native_cad_login_not_acquired_geometry_derived_proxy"
    if name == "MiamiSoft I03":
        return "review_ready_pending_approval", "exact_baxter_miami_soft_i03_roll_vector_references_archived_native_cad_login_not_acquired_geometry_derived_proxy_ifc_description_conflict_recorded"
    if name == "Marilyn 01":
        return "review_ready_pending_approval", "exact_baxter_marilyn_bergere_native_dwg_and_3ds_archived_vector_sources_cross_checked_ifc_description_height_conflict_recorded"
    if name == "Marilyn 02":
        return "review_ready_pending_approval", "exact_baxter_marilyn_pouf_80x62x45_native_dwg_and_3ds_archived_vector_sources_cross_checked"
    if name == "TRAP01":
        return "review_ready_pending_approval", "exact_geberit_151_116_11_1_native_dwg_archived_adjustable_default_not_force_fitted_configured_body_proxy"
    if name == "FAU02":
        return "review_ready_pending_approval", "official_Falper_Cilindro_GH2_native_2d_technical_and_3d_DWG_archived_as_nearest_family_candidate_exact_project_match_rejected_geometry_derived_proxy"
    if name == "Libelle":
        return "excluded_other_task_do_not_touch", "out_of_scope_by_user_instruction"
    return "pending", "official_cad_research_required"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--minimum-faces", type=int, default=1000)
    parser.add_argument("--minimum-vertices", type=int, default=2000)
    args = parser.parse_args()

    input_path = args.input.resolve()
    if sha256(input_path) != EXPECTED["formal_ifc"]:
        raise RuntimeError("formal IFC hash mismatch")
    model = ifcopenshell.open(input_path)
    grouped: dict[tuple[str, str | None, str], dict] = defaultdict(
        lambda: {
            "ifc_classes": set(),
            "instance_global_ids": [],
            "maximum_vertex_count": 0,
            "maximum_face_count": 0,
            "representative_global_id": None,
            "geometry_item_types": set(),
        }
    )
    for product in model.by_type("IfcProduct"):
        if not getattr(product, "Representation", None):
            continue
        key = product_type(product)
        record = grouped[key]
        record["ifc_classes"].add(product.is_a())
        record["instance_global_ids"].append(product.GlobalId)
        vertices = faces = 0
        kinds: set[str] = set()
        for representation in product.Representation.Representations:
            for item in representation.Items:
                item_vertices, item_faces, item_kinds = item_complexity(item, set())
                vertices += item_vertices
                faces += item_faces
                kinds.update(item_kinds)
        record["geometry_item_types"].update(kinds)
        if (
            faces > record["maximum_face_count"]
            or vertices > record["maximum_vertex_count"]
        ):
            record["maximum_vertex_count"] = vertices
            record["maximum_face_count"] = faces
            record["representative_global_id"] = product.GlobalId

    products = []
    system_geometry = []
    for (name, type_global_id, type_class), record in grouped.items():
        if (
            record["maximum_face_count"] < args.minimum_faces
            and record["maximum_vertex_count"] < args.minimum_vertices
        ):
            continue
        classes = sorted(record["ifc_classes"])
        primary_class = classes[0]
        payload = {
            "type_name": name,
            "type_global_id": type_global_id,
            "type_class": type_class,
            "ifc_classes": classes,
            "instance_count": len(record["instance_global_ids"]),
            "instance_global_ids": sorted(record["instance_global_ids"]),
            "representative_global_id": record["representative_global_id"],
            "maximum_vertex_count": record["maximum_vertex_count"],
            "maximum_face_count": record["maximum_face_count"],
            "geometry_item_types": sorted(record["geometry_item_types"]),
        }
        if primary_class in PRODUCT_CLASSES:
            status, source_strategy = status_for(name)
            payload.update(
                {
                    "folder": f"output/review/highpoly-types/{MANUAL_FOLDERS.get(name, slugify(name))}",
                    "status": status,
                    "source_strategy": source_strategy,
                }
            )
            products.append(payload)
        elif primary_class in SYSTEM_CLASSES:
            payload.update(
                {
                    "status": "listed_not_product_cad_target",
                    "source_strategy": "system_geometry_not_manufacturer_product_folder",
                }
            )
            system_geometry.append(payload)

    key = lambda item: (item["maximum_face_count"], item["maximum_vertex_count"])
    products.sort(key=key, reverse=True)
    system_geometry.sort(key=key, reverse=True)
    output = args.output.resolve()
    payload = {
        "schema_version": 1,
        "formal_ifc": relative(input_path),
        "formal_ifc_sha256": EXPECTED["formal_ifc"],
        "inventory_method": "IFC representation topology only; no whole-model rendering or create_shape calls",
        "thresholds": {
            "minimum_polygonal_faces": args.minimum_faces,
            "minimum_coordinates": args.minimum_vertices,
        },
        "product_scope_classes": sorted(PRODUCT_CLASSES),
        "product_count": len(products),
        "completed_product_count": sum(item["status"] == "completed" for item in products),
        "in_progress_product_count": sum(item["status"] == "in_progress" for item in products),
        "review_ready_pending_approval_product_count": sum(
            item["status"] == "review_ready_pending_approval" for item in products
        ),
        "excluded_other_task_count": sum(
            item["status"] == "excluded_other_task_do_not_touch" for item in products
        ),
        "products": products,
        "system_geometry_count": len(system_geometry),
        "system_geometry": system_geometry,
    }
    write_json(output, payload)
    print(output)


if __name__ == "__main__":
    main()
