#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Gessi316 54038 linework."""

from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "official Gessi exact 54038 family reference and 54139_54038 article combination; not a project shop drawing"
PRODUCT_PAGE = "https://areapro.gessi.com/en/product/54038"
PRODUCT_API = "https://g-ecatalogue-be-prod-we.azurewebsites.net/public/product/GetProductDetails?country=it&language=en&productCode=54038"
CATALOGUE = "https://gwebassets.gessi.com/strapi-uploads/assets/Gessi_Cataloghi_Bathroom_bf17b260b2.pdf"
DWG_ZIP = "https://gessistorage.blob.core.windows.net/zwa/GPF5403800000G000_arc.zip"
TECHNICAL_PDF = "https://gessistorage.blob.core.windows.net/zc4/GPF5403800000G000_1.pdf"
SOURCE_DWG_SHA256 = "2c56b532bbbb78e5668050d1a7b52a7d545153a85a05efc377c7a651890897f1"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54038"
EXPECTED_PATH_COUNTS = {"plan": 119, "front": 328, "side": 545}

shared.REPRESENTATIVE_GLOBAL_ID = "245NU$zZL0d9tYTVwwBdk$"
shared.IFC_TYPE_NAME = "Gessi316 54038"
shared.PROFILE_KEY = "gessi316-54038"
shared.SCOPE = SCOPE
shared.SOURCE_KIND = "native_dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54038 原生 DWG 图纸表达"
shared.OFFICIAL_CAD_USED = True
shared.OFFICIAL_CAD_GEOMETRY_INCLUDED = True
shared.CLOSE_REPRESENTATION_PATHS = False
shared.REPRESENTATIONS = {
    "plan": ("Gessi54038Plan", "PLAN_VIEW"),
    "front": ("Gessi54038Front", "ELEVATION_VIEW"),
    "side": ("Gessi54038Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/gessi316-54038-drawing-approval.json"
shared.PSET_NAME = "Pset_Gessi31654038DrawingSource"
shared.DOCUMENT_ID_PREFIX = "GESSI316-54038-"

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
        or candidate.get("article_number") != "54139_54038"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or candidate.get("official_cad_used") is not True
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending Gessi candidate source gate failed")
    drawing_source = access.get("drawing_geometry_source", {})
    if (
        access.get("official_product_cad", {}).get("acquired") is not True
        or access.get("official_product_cad", {}).get("exact_project_configuration_match") is not True
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or drawing_source.get("official_cad_used") is not True
        or drawing_source.get("third_party_cad_used") is not False
        or access.get("scope") != SCOPE
    ):
        raise RuntimeError("Gessi official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("source_dwg_sha256") != SOURCE_DWG_SHA256:
            raise RuntimeError(f"Gessi {view} linework must remain exact native DWG geometry")
        paths[view] = item.get("official_native_dwg_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi {view} native-DWG path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        (
            "GESSI316-54038-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Gessi official exact 54038 product page",
            "Exact manufacturer product identity",
        ),
        (
            "GESSI316-54038-OFFICIAL-PRODUCT-API",
            PRODUCT_API,
            "Gessi public exact 54038 product API",
            "Exact identity and 265 x 101 x 297 mm nominal dimensions",
        ),
        (
            "GESSI316-54038-OFFICIAL-NATIVE-DWG-ZIP",
            DWG_ZIP,
            "Gessi exact 54038 native 2D DWG ZIP",
            "Contains GPF5403800000G000_3.dwg, the exact source of approved Plan and Elevation linework",
        ),
        (
            "GESSI316-54038-OFFICIAL-TECHNICAL-PDF",
            TECHNICAL_PDF,
            "Gessi exact 54038 technical drawing PDF",
            "Vector drawing identity and nominal-dimension cross-check for the native DWG",
        ),
        (
            "GESSI316-54038-OFFICIAL-CATALOGUE",
            CATALOGUE,
            "Gessi official bathroom catalogue",
            "Catalogue identifies exact article combination 54139_54038",
        ),
        (
            "GESSI316-54038-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "Gessi exact-product CAD access record",
            f"SHA-256 {sha256(access_path)}; exact native DWG SHA-256 {SOURCE_DWG_SHA256}; no adjacent or third-party CAD used",
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
        "Family": "Gessi316",
        "ArticleNumber": "54038",
        "IFCTypeName": shared.IFC_TYPE_NAME,
        "SourceProductPage": PRODUCT_PAGE,
        "SourceProductApi": PRODUCT_API,
        "SourceCatalogue": CATALOGUE,
        "SourceNativeDwgZip": DWG_ZIP,
        "SourceTechnicalDrawingPdf": TECHNICAL_PDF,
        "SourceDwgSha256": SOURCE_DWG_SHA256,
        "OfficialProductCadStatus": "public_official_api_exact_54038_native_dwg_acquired",
        "OfficialCadUsed": "true",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Gessi54038Plan;Gessi54038Front;Gessi54038Side",
        "RepresentationGeometrySource": "official_native_dwg_paths_mm extracted from GPF5403800000G000_3.dwg",
        "OfficialCadGeometryIncluded": "true",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
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
        Description="Mechanically verifiable exact native-DWG drawing source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Gessi316 54038 drawing source properties",
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
