#!/usr/bin/env python3
"""Approval-gated derived IFC writer for exact Gessi316 54146 G000 DWG linework."""

from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "official Gessi exact 54146 G000 family reference; not a project shop drawing"
PRODUCT_ROUTE = "https://areapro.gessi.com/en/product/54146"
PRODUCT_API = "https://g-ecatalogue-be-prod-we.azurewebsites.net/public/product/GetProductDetails?country=it&language=en&productCode=54146"
DWG_ZIP = "https://gessistorage.blob.core.windows.net/zwa/GPF5414600000G000_arc.zip"
TECHNICAL_PDF = "https://gessistorage.blob.core.windows.net/zc4/GPF5414600000G000_1.pdf"
CATALOGUE = "https://gwebassets.gessi.com/strapi-uploads/assets/MAGAZINE_GESSI_316_2026_d8f9138389.pdf"
SOURCE_DWG_SHA256 = "c8ddb90f61565d5273a32574f777812bbb1d9ef33e2b98479f1e7c381e69950d"
TECHNICAL_DRAWING_NUMBER = "GPF5414600000G000"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54146"
EXPECTED_PATH_COUNTS = {"plan": 16, "front": 607, "side": 660}

shared.REPRESENTATIVE_GLOBAL_ID = "04DLh1Jk9Dcu9ibcaE0id8"
shared.IFC_TYPE_NAME = "Gessi316 54146"
shared.PROFILE_KEY = "gessi316-54146"
shared.SCOPE = SCOPE
shared.SOURCE_KIND = "native_dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54146 G000 原生 DWG 图纸表达"
shared.OFFICIAL_CAD_USED = True
shared.OFFICIAL_CAD_GEOMETRY_INCLUDED = True
shared.CLOSE_REPRESENTATION_PATHS = False
shared.REPRESENTATIONS = {
    "plan": ("Gessi54146Plan", "PLAN_VIEW"),
    "front": ("Gessi54146Front", "ELEVATION_VIEW"),
    "side": ("Gessi54146Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/gessi316-54146-drawing-approval.json"
shared.PSET_NAME = "Pset_Gessi31654146DrawingSource"
shared.DOCUMENT_ID_PREFIX = "GESSI316-54146-"

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
        or candidate.get("article_number") != "54146"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or candidate.get("official_cad_used") is not True
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending Gessi316 54146 candidate source gate failed")
    product_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    cross_check = access.get("dimension_cross_check", {})
    if (
        access.get("resolved_article_number") != "54146"
        or access.get("resolved_configuration") != "G000"
        or product_cad.get("authentication_required") is not False
        or product_cad.get("acquired") is not True
        or product_cad.get("exact_project_configuration_match") is not True
        or product_cad.get("native_dwg", {}).get("sha256") != SOURCE_DWG_SHA256
        or cross_check.get("maximum_ifc_projection_delta_mm") != 0.878759
        or cross_check.get("pass") is not True
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or drawing_source.get("official_cad_used") is not True
        or drawing_source.get("third_party_cad_used") is not False
        or access.get("scope") != SCOPE
    ):
        raise RuntimeError("Gessi316 54146 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("source_dwg_sha256") != SOURCE_DWG_SHA256:
            raise RuntimeError(f"Gessi316 54146 {view} must use exact G000 native DWG paths")
        paths[view] = item.get("official_native_dwg_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi316 54146 {view} native-DWG path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        (
            "GESSI316-54146-OFFICIAL-PRODUCT-ROUTE",
            PRODUCT_ROUTE,
            "Gessi official exact 54146 product page",
            "Exact manufacturer product identity",
        ),
        (
            "GESSI316-54146-OFFICIAL-PRODUCT-API",
            PRODUCT_API,
            "Gessi public exact 54146 product API",
            "Exact identity and 300 x 300 x 276 mm nominal dimensions",
        ),
        (
            "GESSI316-54146-OFFICIAL-NATIVE-DWG-ZIP",
            DWG_ZIP,
            "Gessi exact 54146 G000 native 2D DWG ZIP",
            "Contains GPF5414600000G000_3.dwg, the exact source of approved Plan and Elevation linework",
        ),
        (
            "GESSI316-54146-OFFICIAL-TECHNICAL-PDF",
            TECHNICAL_PDF,
            "Gessi exact 54146 G000 technical drawing PDF",
            "Vector identity and dimension cross-check for the native DWG",
        ),
        (
            "GESSI316-54146-OFFICIAL-CATALOGUE",
            CATALOGUE,
            "Gessi official Gessi316 2026 catalogue",
            "Item 95 identifies ceiling-mounted 54146 and distinguishes wall-mounted 54145",
        ),
        (
            "GESSI316-54146-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "Gessi 54146 CAD access record",
            f"SHA-256 {sha256(access_path)}; exact native DWG SHA-256 {SOURCE_DWG_SHA256}; G001, 54145 and third-party CAD excluded",
        ),
    )
    identifiers = []
    for identification, location, name, description in documents:
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
        "ArticleNumber": "54146",
        "Configuration": "G000",
        "IFCTypeName": shared.IFC_TYPE_NAME,
        "SourceProductRoute": PRODUCT_ROUTE,
        "SourceProductApi": PRODUCT_API,
        "SourceNativeDwgZip": DWG_ZIP,
        "SourceTechnicalDrawing": TECHNICAL_PDF,
        "SourceCatalogue": CATALOGUE,
        "TechnicalDrawingNumber": TECHNICAL_DRAWING_NUMBER,
        "SourceDwgSha256": SOURCE_DWG_SHA256,
        "OfficialProductCadStatus": "public_official_api_exact_54146_g000_native_dwg_acquired",
        "OfficialCadUsed": "true",
        "OfficialVectorEvidenceUsedAsCadGeometry": "false",
        "ThirdPartyCadUsed": "false",
        "ExcludedVariant": "G001",
        "ExcludedAdjacentProduct": "54145 wall-mounted adjustable headshower",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Gessi54146Plan;Gessi54146Front;Gessi54146Side",
        "RepresentationGeometrySource": "official_native_dwg_paths_mm extracted from GPF5414600000G000_3.dwg",
        "OfficialCadGeometryIncluded": "true",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialWidthDepthHeightMm": "300;300;276",
        "ProjectBodyLocalXyzMm": "300;299.944916;279.371241",
        "OfficialNativeDwgPlanEnvelopeMm": "300;300",
        "OfficialNativeDwgFrontEnvelopeMm": "300;280.25",
        "OfficialNativeDwgSideEnvelopeMm": "300;280.25",
        "MaximumIfcProjectionDeltaMm": "0.878759",
        "DimensionCrossCheckPass": "true",
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
        Description="Mechanically verifiable exact Gessi316 54146 G000 native-DWG source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Gessi316 54146 drawing source properties",
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
