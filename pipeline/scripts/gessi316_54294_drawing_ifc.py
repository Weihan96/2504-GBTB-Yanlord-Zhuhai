#!/usr/bin/env python3
"""Approval-gated derived IFC writer for exact Gessi 54294 native-DWG linework."""

from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "official Gessi exact 54294 family reference and 45089_54294 article combination; not a project shop drawing"
PRODUCT_PAGE = "https://areapro.gessi.com/en/product/54294"
DWG_ZIP_URL = "https://gessistorage.blob.core.windows.net/zwa/GPF5429400000G000_arc.zip"
TECHNICAL_PDF_URL = "https://gessistorage.blob.core.windows.net/zc4/GPF5429400000G000_1.pdf"
CATALOGUE = "https://gwebassets.gessi.com/strapi-uploads/assets/MAGAZINE_GESSI_316_2026_d8f9138389.pdf"
REGIONAL_CATALOGUE = "https://gwebassets.gessi.com/strapi-uploads/assets/Bathroom_Digest_Cina_Hong_Kong_c0ef1d4cd3.pdf"
SOURCE_DWG_SHA256 = "dfe8bcf7e100c14187c387b75a35e3fe7f718c61595fadf66adfe92ca7648ff4"
SOURCE_ZIP_SHA256 = "fad0d98e94a83bb488b5f0703472862d628759f21700c1a87c3af343f4c1bdb9"
TECHNICAL_PDF_SHA256 = "82058bb27752f27e6fc92029783d67c15fd443befcd393f7d572fd0147a581e3"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/gessi316-54294"
OFFICIAL_SOURCE_DIR = PRODUCT_DIR / "official-source"
SOURCE_DWG = OFFICIAL_SOURCE_DIR / "GPF5429400000G000_3.dwg"
SOURCE_ZIP = OFFICIAL_SOURCE_DIR / "GPF5429400000G000_arc.zip"
TECHNICAL_PDF = OFFICIAL_SOURCE_DIR / "GPF5429400000G000_1.pdf"
SOURCE_REVALIDATION = OFFICIAL_SOURCE_DIR / "official-source-revalidation.json"
PDF_VERIFICATION = OFFICIAL_SOURCE_DIR / "official-pdf-verification.json"
EXPECTED_PATH_COUNTS = {"plan": 1481, "front": 1099, "side": 745}

shared.REPRESENTATIVE_GLOBAL_ID = "2iKOL78$H0N9Yd9$ky3pW4"
shared.IFC_TYPE_NAME = "Gessi316 54294"
shared.PROFILE_KEY = "gessi316-54294"
shared.SOURCE_KIND = "native_dwg"
shared.SOURCE_LABEL_ZH = "Gessi 官方精确型号 54294 原生 DWG 图纸表达"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Gessi54294Plan", "PLAN_VIEW"),
    "front": ("Gessi54294Front", "ELEVATION_VIEW"),
    "side": ("Gessi54294Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = OFFICIAL_SOURCE_DIR / "source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/gessi316-54294-drawing-approval.json"
shared.PSET_NAME = "Pset_Gessi31654294DrawingSource"
shared.DOCUMENT_ID_PREFIX = "GESSI316-45089-54294-"
shared.CLOSE_REPRESENTATION_PATHS = False
shared.OFFICIAL_CAD_USED = True
shared.OFFICIAL_CAD_GEOMETRY_INCLUDED = True

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
        or candidate.get("article_number") != "45089_54294"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or candidate.get("official_cad_used") is not True
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending Gessi candidate source gate failed")
    official = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    identity_policy = access.get("identity_and_geometry_policy", {})
    if (
        official.get("authentication_required") is not False
        or official.get("acquired") is not True
        or official.get("exact_project_configuration_match") is not True
        or official.get("native_dwg", {}).get("sha256") != SOURCE_DWG_SHA256
        or official.get("native_dwg_zip", {}).get("sha256") != SOURCE_ZIP_SHA256
        or official.get("technical_vector_pdf", {}).get("sha256") != TECHNICAL_PDF_SHA256
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("source_dwg_sha256") != SOURCE_DWG_SHA256
        or drawing_source.get("official_cad_used") is not True
        or drawing_source.get("third_party_cad_used") is not False
        or drawing_source.get("adjacent_product_cad_used") is not False
        or identity_policy.get("54294_native_dwg_used_for_three_views") is not True
        or identity_policy.get("45089_companion_dwg_used_as_54294_geometry") is not False
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
        or sha256(SOURCE_DWG) != SOURCE_DWG_SHA256
        or sha256(SOURCE_ZIP) != SOURCE_ZIP_SHA256
        or sha256(TECHNICAL_PDF) != TECHNICAL_PDF_SHA256
    ):
        raise RuntimeError("Gessi official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if (
            item.get("source_kind") != shared.SOURCE_KIND
            or item.get("source_dwg_sha256") != SOURCE_DWG_SHA256
            or item.get("official_native_dwg_path_count") != EXPECTED_PATH_COUNTS[view]
        ):
            raise RuntimeError(f"Gessi {view} official native-DWG identity gate failed")
        paths[view] = item.get("official_native_dwg_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Gessi {view} official native-DWG path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        ("GESSI316-45089-54294-OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Gessi Area Pro exact 54294 product page", "Exact external visible product identity and public attachment surface"),
        ("GESSI316-45089-54294-OFFICIAL-NATIVE-DWG-ZIP", DWG_ZIP_URL, "Gessi exact 54294 official native DWG archive", f"Archived at {relative(SOURCE_ZIP)}; SHA-256 {SOURCE_ZIP_SHA256}"),
        ("GESSI316-45089-54294-OFFICIAL-NATIVE-DWG", relative(SOURCE_DWG), "Gessi exact 54294 extracted native DWG", f"ZIP member GPF5429400000G000_3.dwg; SHA-256 {SOURCE_DWG_SHA256}; source of Plan and Elevation representations"),
        ("GESSI316-45089-54294-OFFICIAL-TECHNICAL-PDF", TECHNICAL_PDF_URL, "Gessi exact 54294 official technical drawing PDF", f"Archived at {relative(TECHNICAL_PDF)}; SHA-256 {TECHNICAL_PDF_SHA256}; identity and dimension cross-check"),
        ("GESSI316-45089-54294-OFFICIAL-CATALOGUE", CATALOGUE, "Gessi official Gessi316 2026 catalogue", "Item 64 identifies the 45089_54294 article combination"),
        ("GESSI316-45089-54294-OFFICIAL-REGIONAL-CATALOGUE", REGIONAL_CATALOGUE, "Gessi official regional bathroom catalogue", "Maps the companion 45089 built-in part to 54294 external visible product options"),
        ("GESSI316-45089-54294-SOURCE-REVALIDATION", relative(SOURCE_REVALIDATION), "Gessi 54294 public-source revalidation", f"SHA-256 {sha256(SOURCE_REVALIDATION)}; public API, download and ZIP-member verification"),
        ("GESSI316-45089-54294-PDF-VERIFICATION", relative(PDF_VERIFICATION), "Gessi 54294 official PDF verification", f"SHA-256 {sha256(PDF_VERIFICATION)}; technical and catalogue cross-checks"),
        ("GESSI316-45089-54294-SOURCE-ACCESS-RECORD", relative(access_path), "Gessi exact 54294 official-source record", f"SHA-256 {sha256(access_path)}; exact native DWG used; no third-party, adjacent-product or companion-45089 geometry used"),
    )
    identifiers = []
    for identification, location, name, description in documents:
        reference = model.create_entity("IfcDocumentReference", Location=location, Identification=identification, Name=name, Description=description, ReferencedDocument=None)
        model.create_entity("IfcRelAssociatesDocument", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=f"{name} association", Description=SCOPE, RelatedObjects=[product, product_type], RelatingDocument=reference)
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "SourceLabelEn": "drawing representation from the exact Gessi 54294 official native DWG",
        "Manufacturer": "Gessi",
        "Family": "Gessi316 Meccanica",
        "ArticleNumber": "45089_54294",
        "ExternalProductCode": "54294",
        "CompanionBuiltInProductCode": "45089",
        "IFCTypeName": shared.IFC_TYPE_NAME,
        "SourceProductPage": PRODUCT_PAGE,
        "SourceNativeDwgZip": DWG_ZIP_URL,
        "SourceNativeDwgZipSha256": SOURCE_ZIP_SHA256,
        "SourceNativeDwg": relative(SOURCE_DWG),
        "SourceNativeDwgSha256": SOURCE_DWG_SHA256,
        "SourceTechnicalPdf": TECHNICAL_PDF_URL,
        "SourceTechnicalPdfSha256": TECHNICAL_PDF_SHA256,
        "SourceCatalogue": CATALOGUE,
        "OfficialCadUsed": "true",
        "ThirdPartyCadUsed": "false",
        "AdjacentProductCadUsed": "false",
        "Companion45089DwgUsedAs54294Geometry": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Gessi54294Plan;Gessi54294Front;Gessi54294Side",
        "RepresentationGeometrySource": "official_native_dwg_paths_mm from GPF5429400000G000_3.dwg",
        "OfficialCadGeometryIncluded": "true",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "SourceRevalidation": relative(SOURCE_REVALIDATION),
        "SourceRevalidationSha256": sha256(SOURCE_REVALIDATION),
        "PdfVerification": relative(PDF_VERIFICATION),
        "PdfVerificationSha256": sha256(PDF_VERIFICATION),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
    }
    properties = [model.create_entity("IfcPropertySingleValue", Name=name, Description=None, NominalValue=model.create_entity("IfcText", str(value)), Unit=None) for name, value in values.items()]
    pset = model.create_entity("IfcPropertySet", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=shared.PSET_NAME, Description="Mechanically verifiable exact native-DWG drawing source and human approval", HasProperties=properties)
    model.create_entity("IfcRelDefinesByProperties", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name="Gessi316 45089_54294 drawing source properties", Description=None, RelatedObjects=[product, product_type], RelatingPropertyDefinition=pset)


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset
shared.require_approval = require_approval


if __name__ == "__main__":
    shared.main()
