#!/usr/bin/env python3
"""Approval-gated derived IFC writer for Molteni Sistema 7 / SIS04."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "exact Molteni&C Sistema 7 Wall Unit family identity and standard dimensions; project IFC configuration is labelled 4 Doors; exact official CAD not acquired; not a project shop drawing"
PRODUCT_PAGE = "https://www.molteni.it/ap/product/sistema-7-wall-unit"
OFFICIAL_CATALOGUE = "https://www.molteni.it/en/download/document/3c1c2c00b23e4dbb8ab7c89fa68936d8008a7a56"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/sis04"
PROFILE = PRODUCT_DIR / "profile.json"
PROJECT_CONTEXT_MANIFEST = PRODUCT_DIR / "project-context-manifest.json"
BONSAI_REVIEW_MANIFEST = PRODUCT_DIR / "bonsai-review-manifest.json"
EXPECTED_PATH_COUNTS = {"plan": 4, "front": 4, "side": 1}

shared.REPRESENTATIVE_GLOBAL_ID = "1MzM8Ms2vFo8KEm503j9w2"
shared.IFC_TYPE_NAME = "SIS04"
shared.PROFILE_KEY = "sis04"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Sis04Plan", "PLAN_VIEW"),
    "front": ("Sis04Front", "ELEVATION_VIEW"),
    "side": ("Sis04Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/sis04-drawing-approval.json"
shared.PSET_NAME = "Pset_Sis04DrawingSource"
shared.DOCUMENT_ID_PREFIX = "MOLTENI-SISTEMA7-WALL-UNIT-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "Sistema 7 Wall Unit 4 Doors / SIS04"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending SIS04 candidate source gate failed")
    official_cad = access.get("official_product_cad", {})
    drawing_source = access.get("drawing_geometry_source", {})
    dimensions = access.get("dimension_cross_check", {})
    if (
        access.get("identity_status") != "exact_manufacturer_family_and_wall_unit_product_confirmed_by_ifc_description_and_official_dimensions"
        or official_cad.get("acquired") is not False
        or official_cad.get("local_cad_files") != []
        or official_cad.get("near_name_sistema_7_doors_dwg_excluded") is not True
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or dimensions.get("status") != "exact_family_and_wall_unit_size_class_consistent_project_depth_is_not_fabrication_dimension"
        or dimensions.get("pass") is not True
        or access.get("scope") != SCOPE
    ):
        raise RuntimeError("SIS04 source-access record gate failed")
    paths = {}
    expected_axes = {"plan": [0, 1], "front": [0, 2], "side": [1, 2]}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if (
            item.get("source_kind") != shared.SOURCE_KIND
            or item.get("official_cad_paths_mm") != []
            or item.get("projection_axes") != expected_axes[view]
        ):
            raise RuntimeError(f"SIS04 {view} must remain geometry-derived linework")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"SIS04 {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    access = json.loads(access_path.read_text(encoding="utf-8"))
    evidence = {item["kind"]: item for item in access["official_identity_sources"]}
    catalogue = evidence["manufacturer_kitchen_catalogue"]
    near_name = evidence["manufacturer_near_name_download_surface_excluded"]
    documents = (
        (
            "MOLTENI-SISTEMA7-WALL-UNIT-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Molteni&C official Sistema 7 Wall Unit product page",
            "Exact manufacturer family and wall-unit identity evidence only; not drawing geometry",
        ),
        (
            "MOLTENI-SISTEMA7-WALL-UNIT-OFFICIAL-KITCHEN-CATALOGUE",
            OFFICIAL_CATALOGUE,
            "Molteni&C official Kitchen Collection catalogue",
            f"SHA-256 {catalogue['sha256']}; printed pages 272-275 / PDF file pages 156-157; identity and nominal dimensions only",
        ),
        (
            "MOLTENI-SISTEMA7-WALL-UNIT-EXCLUDED-NEAR-NAME-DOWNLOAD",
            near_name["url"],
            "Excluded near-name Molteni Sistema 7 Doors download surface",
            f"SHA-256 {near_name['sha256']}; separate full-height door product; its DWG is not used for this Wall Unit",
        ),
        (
            "MOLTENI-SISTEMA7-WALL-UNIT-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "SIS04 identity, dimension and CAD exclusion record",
            f"SHA-256 {sha256(access_path)}; no official or third-party CAD geometry used",
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
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "Molteni&C",
        "Family": "Sistema 7",
        "Product": "Sistema 7 Wall Unit",
        "ProjectTypeCode": "SIS04",
        "IFCTypeDescription": "Sistema 7 Wall Unit 4 Doors",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceProductPageSha256": evidence["manufacturer_product_page"]["sha256"],
        "SourceOfficialCatalogue": OFFICIAL_CATALOGUE,
        "SourceOfficialCatalogueSha256": evidence["manufacturer_kitchen_catalogue"]["sha256"],
        "SourceOfficialCataloguePages": "printed 272-275; PDF file 156-157",
        "OfficialCadStatus": "exact_wall_unit_native_cad_not_acquired",
        "NearNameSistema7DoorsDwgExcluded": "true",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Sis04Plan;Sis04Front;Sis04Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
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
        "ProjectIfcBodyLocalXyzMm": json.dumps(dimensions["project_ifc_body_local_xyz_mm"]),
        "OfficialSelectedWidthHeightDepthMm": json.dumps(dimensions["selected_official_width_height_depth_mm"]),
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
        Description="Mechanically verifiable Molteni identity, geometry-derived drawing source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="SIS04 drawing source and approval properties",
        Description=None,
        RelatedObjects=[product] + ([product_type] if product_type is not None else []),
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
