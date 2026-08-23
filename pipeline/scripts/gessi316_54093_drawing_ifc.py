#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Gessi316 54093 linework."""

from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "official Gessi exact 54093 G000 family reference; not a project shop drawing"
PRODUCT_PAGE = "https://areapro.gessi.com/en/product/54093"
PRODUCT_API = "https://g-ecatalogue-be-prod-we.azurewebsites.net/public/product/GetProductDetails?country=it&language=en&productCode=54093"
DWG_ZIP = "https://gessistorage.blob.core.windows.net/zwa/GPF5409300000G000_arc.zip"
TECHNICAL_PDF = "https://gessistorage.blob.core.windows.net/zc4/GPF5409300000G000_1.pdf"
SOURCE_DWG_SHA256 = "a9308fc8498d34c8cf2f68fd28aaf90b11e62bfc59424e0a5a3ac7f0d28d5047"
COLLECTION_CATALOGUE = "https://gwebassets.gessi.com/strapi-uploads/assets/MAGAZINE_GESSI_316_2026_d8f9138389.pdf"
BATHROOM_CATALOGUE = "https://gwebassets.gessi.com/strapi-uploads/assets/Gessi_Cataloghi_Bathroom_bf17b260b2.pdf"
TECHNICAL_DRAWING_NUMBER = "GPF5409300000G000"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54093"
SOURCE_DIR = PRODUCT_DIR / "official-source"
PRODUCT_PAGE_ARCHIVE = SOURCE_DIR / "gessi316-official-product-page.html"
COLLECTION_CATALOGUE_ARCHIVE = SOURCE_DIR / "MAGAZINE_GESSI_316_2026.pdf"
BATHROOM_CATALOGUE_ARCHIVE = SOURCE_DIR / "Gessi_Cataloghi_Bathroom.pdf"
EXPECTED_PATH_COUNTS = {"plan": 7, "front": 914, "side": 21}

shared.REPRESENTATIVE_GLOBAL_ID = "1jM_suNMPAw8_nvZejQSnp"
shared.IFC_TYPE_NAME = "Gessi316 54093"
shared.PROFILE_KEY = "gessi316-54093"
shared.SCOPE = SCOPE
shared.SOURCE_KIND = "native_dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54093 G000 原生 DWG 图纸表达"
shared.OFFICIAL_CAD_USED = True
shared.OFFICIAL_CAD_GEOMETRY_INCLUDED = True
shared.CLOSE_REPRESENTATION_PATHS = False
shared.REPRESENTATIONS = {
    "plan": ("Gessi54093Plan", "PLAN_VIEW"),
    "front": ("Gessi54093Front", "ELEVATION_VIEW"),
    "side": ("Gessi54093Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/gessi316-54093-drawing-approval.json"
shared.PSET_NAME = "Pset_Gessi31654093DrawingSource"
shared.DOCUMENT_ID_PREFIX = "GESSI316-54093-"

shared_require_approval = shared.require_approval


def require_approval(approval: dict, manifest_path: Path) -> None:
    shared_require_approval(approval, manifest_path)
    if approval.get("formal_authoritative_ifc_write_allowed") is not False:
        raise RuntimeError(
            "approval gate rejected IFC write: formal_authoritative_ifc_write_allowed must remain false"
        )


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "54093"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or candidate.get("official_cad_used") is not True
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending Gessi316 54093 candidate source gate failed")
    product_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    cross_check = access.get("dimension_cross_check", {})
    if (
        access.get("resolved_article_number") != "54093"
        or access.get("resolved_configuration") != "G000"
        or access.get("scope") != SCOPE
        or product_cad.get("authentication_required") is not False
        or product_cad.get("acquired") is not True
        or product_cad.get("exact_project_configuration_match") is not True
        or product_cad.get("native_dwg", {}).get("sha256") != SOURCE_DWG_SHA256
        or cross_check.get("maximum_official_nominal_delta_mm") != 0.335804
        or cross_check.get("pass") is not True
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or drawing_source.get("official_cad_used") is not True
        or drawing_source.get("third_party_cad_used") is not False
    ):
        raise RuntimeError("Gessi316 54093 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("source_dwg_sha256") != SOURCE_DWG_SHA256:
            raise RuntimeError(f"Gessi316 54093 {view} must use exact G000 native DWG paths")
        paths[view] = item.get("official_native_dwg_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi316 54093 {view} native-DWG path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        ("OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Gessi official exact 54093 product page", "Exact manufacturer product identity"),
        ("OFFICIAL-PRODUCT-API", PRODUCT_API, "Gessi public exact 54093 product API", "Exact identity and 50 x 190 x 273 mm nominal dimensions"),
        ("OFFICIAL-NATIVE-DWG-ZIP", DWG_ZIP, "Gessi exact 54093 G000 native 2D DWG ZIP", "Contains GPF5409300000G000_3.dwg, the exact source of approved Plan and Elevation linework"),
        ("OFFICIAL-TECHNICAL-PDF", TECHNICAL_PDF, "Gessi exact 54093 G000 technical drawing PDF", "Vector identity and dimension cross-check for the native DWG"),
        ("OFFICIAL-BATHROOM-CATALOGUE", BATHROOM_CATALOGUE, "Gessi official Bathroom catalogue", "Identifies 54093 H 230 mm and distinguishes adjacent 54091"),
        ("SOURCE-ACCESS-RECORD", relative(access_path), "Gessi316 54093 G000 source access record", f"SHA-256 {sha256(access_path)}; exact native DWG SHA-256 {SOURCE_DWG_SHA256}; G001/A004 and third-party CAD excluded"),
    )
    identifiers = []
    for suffix, location, name, description in documents:
        identification = f"GESSI316-54093-{suffix}"
        reference = model.create_entity(
            "IfcDocumentReference",
            Location=location,
            Identification=identification,
            Name=name,
            Description=description,
            ReferencedDocument=None,
        )
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
        "ArticleNumber": "54093",
        "IFCTypeName": shared.IFC_TYPE_NAME,
        "SourceProductPage": PRODUCT_PAGE,
        "SourceProductApi": PRODUCT_API,
        "SourceNativeDwgZip": DWG_ZIP,
        "SourceBathroomCatalogue": BATHROOM_CATALOGUE,
        "SourceTechnicalDrawing": TECHNICAL_PDF,
        "TechnicalDrawingNumber": TECHNICAL_DRAWING_NUMBER,
        "SourceDwgSha256": SOURCE_DWG_SHA256,
        "OfficialTechnicalDrawingLocalArchive": "true",
        "OfficialProductCadStatus": "public_official_api_exact_54093_g000_native_dwg_acquired",
        "OfficialCadUsed": "true",
        "OfficialVectorEvidenceUsedAsCadGeometry": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Gessi54093Plan;Gessi54093Front;Gessi54093Side",
        "RepresentationGeometrySource": "official_native_dwg_paths_mm extracted from GPF5409300000G000_3.dwg",
        "OfficialCadGeometryIncluded": "true",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "BathroomCatalogueArchiveSha256": sha256(BATHROOM_CATALOGUE_ARCHIVE),
        "PlanPathCount": "7",
        "FrontPathCount": "914",
        "SidePathCount": "21",
        "OfficialWidthDepthHeightMm": "50;190;273",
        "ProjectBodyLocalXyzMm": "48.974943;190.281433;273.203123",
        "OfficialNativeDwgPlanEnvelopeMm": "49.973479;190.312335",
        "OfficialNativeDwgFrontEnvelopeMm": "50;273.210406",
        "OfficialNativeDwgSideEnvelopeMm": "190.335804;273.202791",
        "MaximumOfficialNominalDeltaMm": "0.335804",
        "DimensionCrossCheckPass": "true",
        "NativeProjectElevationFound": "false",
    }
    properties = [
        model.create_entity(
            "IfcPropertySingleValue",
            Name=name,
            Description=None,
            NominalValue=model.create_entity("IfcText", str(value)),
            Unit=None,
        )
        for name, value in values.items()
    ]
    pset = model.create_entity(
        "IfcPropertySet",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name=shared.PSET_NAME,
        Description="Mechanically verifiable exact Gessi316 54093 G000 native-DWG source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Gessi316 54093 drawing source properties",
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
