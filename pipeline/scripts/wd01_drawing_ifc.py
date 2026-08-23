#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Poliform Pivot + Senzafine / WD01."""

import json

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "official manufacturer Pivot and Senzafine system identity and remote technical-document evidence only; the project 1202.879 x 673.524 x 2390.023 mm arrangement is a custom configuration, not exact official CAD geometry or a project shop drawing"
PIVOT_PAGE = "https://www.poliform.it/en/products/pivot/"
SENZAFINE_PAGE = "https://www.poliform.it/en/products/senzafine-wardrobe/"
PIVOT_TECHNICAL_PDF = "https://www.poliform.it/wp-content/uploads/pdf/278373-pivot-poliform-en-us.pdf"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/wd01"
PAGE_EVIDENCE = PRODUCT_DIR / "official-source/official-product-page-evidence.json"
EXPECTED_PATH_COUNTS = {"plan": 7, "front": 3, "side": 7}

shared.REPRESENTATIVE_GLOBAL_ID = "2mPTt7$nvBSAn28I9NWt6T"
shared.IFC_TYPE_NAME = "WD01"
shared.PROFILE_KEY = "wd01"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Wd01Plan", "PLAN_VIEW"),
    "front": ("Wd01Front", "ELEVATION_VIEW"),
    "side": ("Wd01Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/wd01-drawing-approval.json"
shared.PSET_NAME = "Pset_Wd01DrawingSource"
shared.DOCUMENT_ID_PREFIX = "POLIFORM-PIVOT-SENZAFINE-WD01-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "Poliform Pivot + Senzafine custom wardrobe / project WD01"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending WD01 candidate source gate failed")
    cad = access.get("official_product_cad", {})
    remote = access.get("official_remote_documents", {})
    drawing_source = access.get("drawing_geometry_source", {})
    identity = access.get("project_identity_cross_check", {})
    dimensions = access.get("dimension_cross_check", {})
    evidence = access.get("official_identity_sources", [])
    if (
        cad.get("registration_form_and_captcha_required") is not True
        or cad.get("public_exact_native_dwg_url_located") is not False
        or cad.get("acquired") is not False
        or cad.get("local_cad_files") != []
        or cad.get("exact_project_configuration_match") is not False
        or remote.get("identity_and_system_relationship_observed") is not True
        or remote.get("binary_archived") is not False
        or remote.get("used_as_representation_geometry") is not False
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or identity.get("pass") is not True
        or dimensions.get("exact_dimension_match_claimed") is not False
        or dimensions.get("geometry_scaled_to_catalogue_example") is not False
        or access.get("scope") != SCOPE
        or len(evidence) != 1
        or evidence[0].get("local_path") != relative(PAGE_EVIDENCE)
        or evidence[0].get("sha256") != sha256(PAGE_EVIDENCE)
    ):
        raise RuntimeError("WD01 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"WD01 {view} linework must remain geometry-derived with no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"WD01 {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path):
    documents = (
        (
            "POLIFORM-PIVOT-SENZAFINE-WD01-OFFICIAL-PIVOT-PAGE",
            PIVOT_PAGE,
            "Poliform official Pivot product page",
            "Manufacturer and Pivot system identity only; not the source of drawing geometry",
        ),
        (
            "POLIFORM-PIVOT-SENZAFINE-WD01-OFFICIAL-SENZAFINE-PAGE",
            SENZAFINE_PAGE,
            "Poliform official Senzafine Wardrobe product page",
            "Manufacturer and modular wardrobe-system identity only; not exact WD01 geometry",
        ),
        (
            "POLIFORM-PIVOT-SENZAFINE-WD01-OFFICIAL-TECHNICAL-PDF",
            PIVOT_TECHNICAL_PDF,
            "Poliform official Pivot technical PDF",
            "Remote system relationship and composition evidence; binary retrieval was Cloudflare-blocked and it is not representation geometry",
        ),
        (
            "POLIFORM-PIVOT-SENZAFINE-WD01-OFFICIAL-PAGE-EVIDENCE",
            relative(PAGE_EVIDENCE),
            "Poliform Pivot and Senzafine official-page evidence",
            f"SHA-256 {sha256(PAGE_EVIDENCE)}; identity and gated-resource observations only",
        ),
        (
            "POLIFORM-PIVOT-SENZAFINE-WD01-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "WD01 manufacturer-source access record",
            f"SHA-256 {sha256(access_path)}; exact native CAD not acquired; no third-party CAD used",
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
    access = json.loads(access_path.read_text(encoding="utf-8"))
    dimensions = access["dimension_cross_check"]
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "Poliform",
        "Family": "Pivot + Senzafine Wardrobe",
        "ProjectTypeCode": "WD01",
        "IFCTypeDescription": "Poliform Pivot Senzafine",
        "IFCElementType": "WARDROBE",
        "SourcePivotProductPage": PIVOT_PAGE,
        "SourceSenzafineProductPage": SENZAFINE_PAGE,
        "SourcePivotTechnicalPdf": PIVOT_TECHNICAL_PDF,
        "OfficialNativeCadStatus": "gated_resource_download_no_public_exact_WD01_configuration_asset_located_not_acquired",
        "OfficialRemotePdfBinaryArchived": "false",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "ExactProjectConfigurationProvenByManufacturer": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Wd01Plan;Wd01Front;Wd01Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialProductPageEvidence": relative(PAGE_EVIDENCE),
        "OfficialProductPageEvidenceSha256": sha256(PAGE_EVIDENCE),
        "ProjectIfcBodyBoundsMm": json.dumps(dimensions["project_ifc_body_bounds_mm"]),
        "DimensionReviewStatus": dimensions["status"],
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
        Description="Mechanically verifiable geometry-derived WD01 drawing source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="WD01 drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
