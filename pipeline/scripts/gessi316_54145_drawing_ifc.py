#!/usr/bin/env python3
"""Approval-gated derived IFC writer for the approved Gessi316 54145 review line."""

from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "official Gessi exact 54145 G000 family reference; not a project shop drawing"
PRODUCT_PAGE = "https://areapro.gessi.com/en/product/54145"
PRODUCT_API = "https://g-ecatalogue-be-prod-we.azurewebsites.net/public/product/GetProductDetails?country=it&language=en&productCode=54145"
DWG_ZIP = "https://gessistorage.blob.core.windows.net/zwa/GPF5414500000G000_arc.zip"
TECHNICAL_PDF = "https://gessistorage.blob.core.windows.net/zc4/GPF5414500000G000_1.pdf"
COLLECTION_CATALOGUE = "https://gwebassets.gessi.com/strapi-uploads/assets/MAGAZINE_GESSI_316_2026_d8f9138389.pdf"
SOURCE_DWG_SHA256 = "9978b68468a61875acb0736aab45a62a98efadfd6be61d347e7ecaa94e08fd09"
TECHNICAL_DRAWING_NUMBER = "GPF5414500000G000"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54145"
SOURCE_DIR = PRODUCT_DIR / "official-source"
COLLECTION_CATALOGUE_ARCHIVE = SOURCE_DIR / "MAGAZINE_GESSI_316_2026.pdf"
EXPECTED_PATH_COUNTS = {"plan": 22, "front": 112, "side": 109}
ORIGINAL_PATH_COUNTS = {"plan": 22, "front": 615, "side": 653}

shared.REPRESENTATIVE_GLOBAL_ID = "3jT4sCgpHC98VSIUdGUNYH"
shared.IFC_TYPE_NAME = "Gessi316 54145"
shared.PROFILE_KEY = "gessi316-54145"
shared.SCOPE = SCOPE
shared.SOURCE_KIND = "native_dwg_review_simplification"
shared.SOURCE_LABEL_ZH = "基于官方54145 G000原生DWG轮廓的简化蓝线审核表达"
shared.OFFICIAL_CAD_USED = True
shared.OFFICIAL_CAD_GEOMETRY_INCLUDED = True
shared.CLOSE_REPRESENTATION_PATHS = False
shared.REPRESENTATIONS = {
    "plan": ("Gessi54145Plan", "PLAN_VIEW"),
    "front": ("Gessi54145Front", "ELEVATION_VIEW"),
    "side": ("Gessi54145Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/gessi316-54145-drawing-approval.json"
shared.PSET_NAME = "Pset_Gessi31654145DrawingSource"
shared.DOCUMENT_ID_PREFIX = "GESSI316-54145-"

shared_require_approval = shared.require_approval


def require_approval(approval: dict, manifest_path: Path) -> None:
    shared_require_approval(approval, manifest_path)
    if approval.get("formal_authoritative_ifc_write_allowed") is not False:
        raise RuntimeError("approval gate rejected IFC write: formal_authoritative_ifc_write_allowed must remain false")


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "54145"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or candidate.get("official_cad_used") is not True
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
        or candidate.get("unaltered_official_cad_used_as_review_representation") is not False
        or candidate.get("original_official_cad_evidence_preserved") is not True
    ):
        raise RuntimeError("pending Gessi316 54145 candidate source gate failed")
    product_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    cross_check = access.get("dimension_cross_check", {})
    if (
        access.get("resolved_article_number") != "54145"
        or access.get("resolved_configuration") != "G000"
        or access.get("scope") != SCOPE
        or product_cad.get("authentication_required") is not False
        or product_cad.get("acquired") is not True
        or product_cad.get("exact_project_configuration_match") is not True
        or product_cad.get("native_dwg", {}).get("sha256") != SOURCE_DWG_SHA256
        or cross_check.get("maximum_ifc_projection_delta_mm") != 0.181335
        or cross_check.get("pass") is not True
        or drawing_source.get("source_kind") != "native_dwg"
        or drawing_source.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or drawing_source.get("official_cad_used") is not True
        or drawing_source.get("third_party_cad_used") is not False
    ):
        raise RuntimeError("Gessi316 54145 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        audit = item.get("handle_line_texture_simplification", {})
        if (
            item.get("source_kind") != shared.SOURCE_KIND
            or item.get("source_dwg_sha256") != SOURCE_DWG_SHA256
            or item.get("unaltered_official_dwg") is not False
            or item.get("original_official_native_dwg_path_count") != ORIGINAL_PATH_COUNTS[view]
            or audit.get("envelope_delta_mm") != [0.0, 0.0]
            or audit.get("centre_delta_mm") != [0.0, 0.0]
            or audit.get("installation_axis_preserved") is not True
            or audit.get("wall_anchor_preserved") is not True
            or audit.get("arm_reach_600mm_preserved") is not True
            or audit.get("pass") is not True
        ):
            raise RuntimeError(f"Gessi316 54145 {view} approved review-simplification gate failed")
        paths[view] = item.get("review_simplified_official_outline_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi316 54145 {view} review path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        ("OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Gessi official exact 54145 product page", "Exact manufacturer product identity"),
        ("OFFICIAL-PRODUCT-API", PRODUCT_API, "Gessi public exact 54145 product API", "Exact identity and 300 x 600 x 115 mm nominal dimensions"),
        ("OFFICIAL-NATIVE-DWG-ZIP", DWG_ZIP, "Gessi exact 54145 G000 native 2D DWG ZIP", "Contains GPF5414500000G000_3.dwg, the exact source of approved Plan and Elevation linework"),
        ("OFFICIAL-TECHNICAL-PDF", TECHNICAL_PDF, "Gessi exact 54145 G000 technical drawing PDF", "Vector identity and dimension cross-check for the native DWG"),
        ("OFFICIAL-COLLECTION-CATALOGUE", COLLECTION_CATALOGUE, "Gessi official Gessi316 catalogue", "Item 94 identifies wall-mounted 54145 and distinguishes ceiling-mounted 54146"),
        ("SOURCE-ACCESS-RECORD", relative(access_path), "Gessi316 54145 G000 source access record", f"SHA-256 {sha256(access_path)}; exact native DWG SHA-256 {SOURCE_DWG_SHA256}; G001, 54146 and third-party CAD excluded"),
    )
    identifiers = []
    for suffix, location, name, description in documents:
        identification = f"GESSI316-54145-{suffix}"
        reference = model.create_entity("IfcDocumentReference", Location=location, Identification=identification, Name=name, Description=description, ReferencedDocument=None)
        model.create_entity(
            "IfcRelAssociatesDocument",
            GlobalId=ifcopenshell.guid.new(),
            OwnerHistory=product.OwnerHistory,
            Name=f"{name} association",
            Description=SCOPE,
            RelatedObjects=[product, product_type],
            RelatingDocument=reference,
        )
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "Gessi",
        "Family": "Gessi316 Meccanica",
        "ArticleNumber": "54145",
        "Configuration": "G000",
        "IFCTypeName": shared.IFC_TYPE_NAME,
        "SourceProductPage": PRODUCT_PAGE,
        "SourceProductApi": PRODUCT_API,
        "SourceNativeDwgZip": DWG_ZIP,
        "SourceTechnicalDrawing": TECHNICAL_PDF,
        "SourceCollectionCatalogue": COLLECTION_CATALOGUE,
        "TechnicalDrawingNumber": TECHNICAL_DRAWING_NUMBER,
        "SourceDwgSha256": SOURCE_DWG_SHA256,
        "OfficialProductCadStatus": "public_official_api_exact_54145_g000_native_dwg_acquired",
        "OfficialCadUsed": "true",
        "OfficialVectorEvidenceUsedAsCadGeometry": "true",
        "UnalteredOfficialDwgUsedAsRepresentation": "false",
        "OriginalOfficialDwgEvidencePreserved": "true",
        "ThirdPartyCadUsed": "false",
        "ExcludedVariant": "G001",
        "ExcludedAdjacentProduct": "54146 ceiling-mounted adjustable headshower",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Gessi54145Plan;Gessi54145Front;Gessi54145Side",
        "RepresentationGeometrySource": "review_simplified_official_outline_paths_mm based on GPF5414500000G000_3.dwg",
        "OfficialCadGeometryIncluded": "true",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "CollectionCatalogueArchiveSha256": sha256(COLLECTION_CATALOGUE_ARCHIVE),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialWidthDepthHeightMm": "300;600;115",
        "ProjectBodyLocalXyzMm": "599.818665;299.818832;119.07444",
        "OfficialNativeDwgPlanEnvelopeMm": "600;299.993081",
        "OfficialNativeDwgFrontEnvelopeMm": "600;118.95",
        "OfficialNativeDwgSideEnvelopeMm": "300;118.9206",
        "MaximumIfcProjectionDeltaMm": "0.181335",
        "DimensionCrossCheckPass": "true",
    }
    properties = [
        model.create_entity("IfcPropertySingleValue", Name=name, Description=None, NominalValue=model.create_entity("IfcText", str(value)), Unit=None)
        for name, value in values.items()
    ]
    pset = model.create_entity(
        "IfcPropertySet",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name=shared.PSET_NAME,
        Description="Mechanically verifiable Gessi316 54145 G000 official-outline review simplification and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Gessi316 54145 drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset
shared.require_approval = require_approval


if __name__ == "__main__":
    shared.main()
