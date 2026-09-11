#!/usr/bin/env python3
"""Approval-gated derived IFC writer for the project sxb010 Venetian blind."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "owner-confirmed Hunter Douglas 25 mm Venetian blind family and colour direction only; exact SKU, control system, finish code and project shop dimensions pending; not a project shop drawing and not official CAD geometry"
PRODUCT_PAGE = "https://www.hunterdouglas.cn/product/venetian-blind/16mm-25mm-venetian-blinds"
TECHNICAL_BROCHURE = "https://www.hunterdouglasdam.eu/m/39189f4b3bd19d3b/original/Brochure_VenetianWoodBlinds_EN.pdf"
DOWNLOADS_PAGE = "https://www.hunterdouglasarchitectural.eu/en-eu/downloads/"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/sxb010"
PROFILE = PRODUCT_DIR / "profile.json"
PROJECT_CONTEXT_MANIFEST = PRODUCT_DIR / "project-context-manifest.json"
BONSAI_REVIEW_MANIFEST = PRODUCT_DIR / "bonsai-review-manifest.json"
EXPECTED_PATH_COUNTS = {"plan": 5, "front": 51, "side": 55}
NEAR_LINE_MERGE_THRESHOLD_MM = 22.0

shared.REPRESENTATIVE_GLOBAL_ID = "1O9JRXCI56VRUbpuLJy86Z"
shared.IFC_TYPE_NAME = "sxb010"
shared.PROFILE_KEY = "sxb010"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Sxb010Plan", "PLAN_VIEW"),
    "front": ("Sxb010Front", "ELEVATION_VIEW"),
    "side": ("Sxb010Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/sxb010-drawing-approval.json"
shared.PSET_NAME = "Pset_Sxb010DrawingSource"
shared.DOCUMENT_ID_PREFIX = "HUNTER-DOUGLAS-SXB010-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "project sxb010 / Hunter Douglas 25 mm family direction"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
        or candidate.get("near_line_merge_threshold_mm") != NEAR_LINE_MERGE_THRESHOLD_MM
    ):
        raise RuntimeError("pending sxb010 candidate source gate failed")
    official_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    dimensions = access.get("dimension_cross_check", {})
    if (
        access.get("identity_status") != "owner_confirmed_family_direction_exact_sku_finish_and_operation_pending"
        or access.get("owner_direction_evidence", {}).get("stash_operation") != "read_only_git_show_only_no_pop_no_apply"
        or official_cad.get("acquired") is not False
        or official_cad.get("local_cad_files") != []
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or dimensions.get("status") != "family_and_size_envelope_consistent_exact_product_and_shop_geometry_unconfirmed"
        or access.get("scope") != SCOPE
    ):
        raise RuntimeError("sxb010 source-access record gate failed")
    paths = {}
    expected_axes = {"plan": [0, 2], "front": [0, 1], "side": [2, 1]}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if (
            item.get("source_kind") != shared.SOURCE_KIND
            or item.get("official_cad_paths_mm") != []
            or item.get("projection_axes") != expected_axes[view]
            or item.get("line_simplification", {}).get("merge_threshold_mm") != NEAR_LINE_MERGE_THRESHOLD_MM
            or item.get("line_simplification", {}).get("slat_count_before") != 45
            or item.get("line_simplification", {}).get("slat_centerline_count_after") != 45
            or item.get("line_simplification", {}).get("outer_envelope_preserved") is not True
            or item.get("line_simplification", {}).get("pass") is not True
        ):
            raise RuntimeError(f"sxb010 {view} must remain semantic-axis geometry-derived linework")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"sxb010 {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    access = json.loads(access_path.read_text(encoding="utf-8"))
    evidence = {item["kind"]: item for item in access["official_identity_sources"]}
    brochure = evidence["manufacturer_technical_brochure"]
    flyer = evidence["manufacturer_product_flyer"]
    downloads = evidence["manufacturer_download_catalogue_and_filtered_response"]
    owner = access["owner_direction_evidence"]
    documents = (
        (
            "HUNTER-DOUGLAS-SXB010-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Hunter Douglas official China 16/25 mm Venetian-blind page",
            "Manufacturer family, material and standard size-envelope evidence only; exact project SKU and shop geometry unconfirmed",
        ),
        (
            "HUNTER-DOUGLAS-SXB010-OFFICIAL-TECHNICAL-BROCHURE",
            TECHNICAL_BROCHURE,
            "Hunter Douglas official Venetian and Wood Blinds technical brochure",
            f"SHA-256 {brochure['sha256']}; printed pages 2-3 / PDF file pages 4-5; 25 mm family dimensions and generic 1000 mm vector example; not project drawing geometry",
        ),
        (
            "HUNTER-DOUGLAS-SXB010-OFFICIAL-PRODUCT-FLYER",
            flyer["url"],
            "Hunter Douglas official Venetian Blinds product flyer",
            f"SHA-256 {flyer['sha256']}; family and daylight-control identity evidence only; not project drawing geometry",
        ),
        (
            "HUNTER-DOUGLAS-SXB010-OFFICIAL-DOWNLOADS",
            DOWNLOADS_PAGE,
            "Hunter Douglas official filtered Venetian-blind downloads",
            "Verified product relation 100963 returned two PDF documents and no CAD file",
        ),
        (
            "HUNTER-DOUGLAS-SXB010-OFFICIAL-DOWNLOAD-FILTER",
            relative(ROOT / downloads["filtered_response_local_path"]),
            "Archived Hunter Douglas Venetian-blind filtered download response",
            f"SHA-256 {downloads['filtered_response_sha256']}; two PDF entries and no DWG or DXF",
        ),
        (
            "HUNTER-DOUGLAS-SXB010-OWNER-DIRECTION-EVIDENCE",
            owner["source"],
            "Owner-confirmed Hunter Douglas 25 mm family and colour direction",
            f"Git blob {owner['git_blob_sha1']}; SHA-256 {owner['sha256']}; inspected read-only with git show; stash not applied or popped",
        ),
        (
            "HUNTER-DOUGLAS-SXB010-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "sxb010 family direction, CAD availability and semantic-axis record",
            f"SHA-256 {sha256(access_path)}; owner source inspected read-only from stash@{{1}}; no stash pop or apply; no official or third-party CAD used",
        ),
    )
    related_objects = [product] + ([product_type] if product_type is not None else [])
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
            RelatedObjects=related_objects,
            RelatingDocument=reference,
        )
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    access = json.loads(access_path.read_text(encoding="utf-8"))
    evidence = {item["kind"]: item for item in access["official_identity_sources"]}
    dimensions = access["dimension_cross_check"]
    owner = access["owner_direction_evidence"]
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "ManufacturerDirection": "Hunter Douglas",
        "FamilyDirection": "25 mm aluminium Venetian blind",
        "ProjectObjectCode": "sxb010",
        "ProjectObjectGlobalId": shared.REPRESENTATIVE_GLOBAL_ID,
        "ProjectObjectTypeStatus": "untyped legacy IfcBuildingElementProxy",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceTechnicalBrochure": TECHNICAL_BROCHURE,
        "SourceTechnicalBrochureSha256": evidence["manufacturer_technical_brochure"]["sha256"],
        "SourceProductFlyer": evidence["manufacturer_product_flyer"]["url"],
        "SourceProductFlyerSha256": evidence["manufacturer_product_flyer"]["sha256"],
        "SourceDownloadsPage": DOWNLOADS_PAGE,
        "OfficialDownloadFilterSha256": evidence["manufacturer_download_catalogue_and_filtered_response"]["filtered_response_sha256"],
        "OfficialCadStatus": "not_published_on_verified_manufacturer_surfaces",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "ExactSkuConfirmed": "false",
        "FinishCodeConfirmed": "false",
        "ControlSystemConfirmed": "false",
        "ProjectShopDimensionsConfirmed": "false",
        "FormalAuthoritativeIfcWriteAllowed": "false",
        "OwnerDirectionEvidence": owner["source"],
        "OwnerDirectionEvidenceGitBlobSha1": owner["git_blob_sha1"],
        "OwnerDirectionEvidenceSha256": owner["sha256"],
        "OwnerDirectionConfirmedDate": owner["confirmed_date"],
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Sxb010Plan;Sxb010Front;Sxb010Side",
        "RepresentationGeometrySource": "semantic-axis proxy_paths_mm derived from the isolated representative IFC Body",
        "SemanticViewAxes": "plan=XZ;front=XY;side=ZY",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "Profile": relative(PROFILE),
        "ProfileSha256": sha256(PROFILE),
        "ProjectContextManifest": relative(PROJECT_CONTEXT_MANIFEST),
        "ProjectContextManifestSha256": sha256(PROJECT_CONTEXT_MANIFEST),
        "BonsaiReviewManifest": relative(BONSAI_REVIEW_MANIFEST),
        "BonsaiReviewManifestSha256": sha256(BONSAI_REVIEW_MANIFEST),
        "ProjectIfcBodyWidthHeightDepthMm": json.dumps(dimensions["interpreted_project_width_height_depth_mm"]),
        "OfficialManual25mmWidthLimitsMm": json.dumps(dimensions["official_standard_manual_25mm_width_limits_mm"]),
        "OfficialManual25mmHeightLimitsMm": json.dumps(dimensions["official_standard_manual_25mm_height_limits_mm"]),
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
        Description="Mechanically verifiable family direction, geometry-derived semantic views and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="sxb010 drawing source and semantic-view properties",
        Description=None,
        RelatedObjects=[product],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
