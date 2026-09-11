#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Molteni&C 505 UP V1.LP.S."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


PRODUCT_PAGE = "https://www.molteni.it/en/ap/product/505-up-system"
DWG_COLLECTION = "https://collection.cloudinary.com/molteni/a8d135c32d8cbaa0513b2eae24638fb4"
TECHNICAL_DWG_URL = "https://res.cloudinary.com/molteni/raw/upload/v1766397289/2021_2D_505-UP_Living-Systems_Indoor.dwg?_s=public-apps"
INSPIRING_DWG_URL = "https://res.cloudinary.com/molteni/raw/upload/v1764598453/2021_2D_505-UP-Sistem_Living-Systems_Indoor-Inspiring-Solution.dwg?_s=public-apps"
TECHNICAL_DWG_SHA256 = "dfe4a6bcd655a3ed19343813a4aba8e139b847cafe8ea69b47623a6e1e71add3"
INSPIRING_DWG_SHA256 = "1e7988d5c7f00522736144fbecde1814032e5673af8cd13a47979e4347231646"
SCOPE = "exact manufacturer family and native family CAD; project-specific V1.LP.S composition is not mechanically matched to an official catalogue composition and is not a project shop drawing"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/505-up-v1-lp-s"
PRODUCT_PAGE_ARCHIVE = PRODUCT_DIR / "official-source/molteni-505-up-system-product-page.html"
TECHNICAL_DWG_ARCHIVE = PRODUCT_DIR / "official-source/2021_2D_505-UP_Living-Systems_Indoor.dwg"
INSPIRING_DWG_ARCHIVE = PRODUCT_DIR / "official-source/2021_2D_505-UP-System_Living-Systems_Indoor-Inspiring-Solution.dwg"
ASSET_INDEX = PRODUCT_DIR / "official-source/molteni-505-up-system-dwg-collection-assets.json"
EXPECTED_PATH_COUNTS = {"plan": 24, "front": 281, "side": 93}

shared.REPRESENTATIVE_GLOBAL_ID = "19MpdkWqXC7uhUNhLQgrce"
shared.IFC_TYPE_NAME = "505 UP V1.LP.S"
shared.PROFILE_KEY = "505-up-v1-lp-s"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Molteni505UpPlan", "PLAN_VIEW"),
    "front": ("Molteni505UpFront", "ELEVATION_VIEW"),
    "side": ("Molteni505UpSide", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/505-up-v1-lp-s-drawing-approval.json"
shared.PSET_NAME = "Pset_Molteni505UpDrawingSource"
shared.DOCUMENT_ID_PREFIX = "MOLTENI-505-UP-"
shared.CLOSE_REPRESENTATION_PATHS = False

shared_require_approval = shared.require_approval


def require_approval(approval: dict, manifest_path: Path) -> None:
    """Permit the explicitly authorized derived revision while re-review is pending."""
    normalized = dict(approval)
    if normalized.get("status") == "revision_pending_review":
        if normalized.get("pending_reapproval_views") != ["plan", "front"]:
            raise RuntimeError("505 revision must remain pending for Plan and Front")
        normalized["status"] = "approved"
    shared_require_approval(normalized, manifest_path)


shared.require_approval = require_approval


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "505 UP System / project V1.LP.S"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status")
        not in ("visual_review_pending", "revision_pending_review")
    ):
        raise RuntimeError("pending Molteni 505 UP candidate source gate failed")
    official_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    cross_check = access.get("configuration_cross_check", {})
    identity = access.get("project_identity_evidence", {})
    if (
        access.get("scope") != SCOPE
        or access.get("pass") is not True
        or official_cad.get("acquired") is not True
        or official_cad.get("exact_project_configuration_match") is not False
        or official_cad.get("official_cad_used_as_candidate_geometry") is not False
        or len(official_cad.get("local_cad_files", [])) != 2
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or cross_check.get("exact_match_count") != 0
        or cross_check.get("pass") is not True
        or identity.get("family_identity_status") != "confirmed"
        or identity.get("project_suffix_status") != "unverified_project_or_asset_composition_suffix"
        or sha256(TECHNICAL_DWG_ARCHIVE) != TECHNICAL_DWG_SHA256
        or sha256(INSPIRING_DWG_ARCHIVE) != INSPIRING_DWG_SHA256
    ):
        raise RuntimeError("Molteni 505 UP official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"Molteni 505 UP {view} must contain no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"Molteni 505 UP {view} geometry-derived path count drifted")
    return paths


def related_objects(product, product_type):
    return [product] + ([product_type] if product_type is not None else [])


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        ("MOLTENI-505-UP-OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Molteni&C official 505 UP product page", "Manufacturer, family and designer identity evidence"),
        ("MOLTENI-505-UP-OFFICIAL-DWG-COLLECTION", DWG_COLLECTION, "Molteni&C official 505 UP DWG collection", "Official family CAD access surface; exact project composition not matched"),
        ("MOLTENI-505-UP-OFFICIAL-TECHNICAL-DWG", TECHNICAL_DWG_URL, "Molteni&C official 505 UP family technical DWG", f"Native DWG SHA-256 {TECHNICAL_DWG_SHA256}; family catalogue evidence only; not inserted as project linework"),
        ("MOLTENI-505-UP-OFFICIAL-INSPIRING-DWG", INSPIRING_DWG_URL, "Molteni&C official 505 UP inspiring-solutions DWG", f"Native DWG SHA-256 {INSPIRING_DWG_SHA256}; catalogue composition evidence only; not inserted as project linework"),
        ("MOLTENI-505-UP-OFFICIAL-PRODUCT-PAGE-ARCHIVE", relative(PRODUCT_PAGE_ARCHIVE), "Archived Molteni&C 505 UP product page", f"SHA-256 {sha256(PRODUCT_PAGE_ARCHIVE)}"),
        ("MOLTENI-505-UP-OFFICIAL-DWG-ASSET-INDEX", relative(ASSET_INDEX), "Archived Molteni&C 505 UP DWG collection index", f"SHA-256 {sha256(ASSET_INDEX)}; two official native DWG assets"),
        ("MOLTENI-505-UP-SOURCE-ACCESS-RECORD", relative(access_path), "Molteni&C 505 UP source access record", f"SHA-256 {sha256(access_path)}; exact project configuration match false"),
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
            RelatedObjects=related_objects(product, product_type),
            RelatingDocument=reference,
        )
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "Molteni&C",
        "Family": "505 UP System",
        "Designer": "Nicola Gallizia",
        "ProjectAssetName": "505 UP V1.LP.S",
        "ProjectSuffixStatus": "unverified_project_or_asset_composition_suffix",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceDwgCollection": DWG_COLLECTION,
        "OfficialFamilyCadAcquired": "true",
        "OfficialCadExactProjectConfigurationMatch": "false",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "OfficialTechnicalDwgSha256": TECHNICAL_DWG_SHA256,
        "OfficialInspiringDwgSha256": INSPIRING_DWG_SHA256,
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Molteni505UpPlan;Molteni505UpFront;Molteni505UpSide",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC MODEL_VIEW Body",
        "OfficialCadGeometryIncluded": "false",
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
        Description="Mechanically verifiable native family-CAD evidence, geometry-derived representation source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Molteni 505 UP drawing source properties",
        Description=None,
        RelatedObjects=related_objects(product, product_type),
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
