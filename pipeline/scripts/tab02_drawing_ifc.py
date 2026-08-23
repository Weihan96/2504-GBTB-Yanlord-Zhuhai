#!/usr/bin/env python3
"""Approval-gated derived IFC writer for RODA Bernardo 367 / TAB02."""

import json

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "exact RODA BERNARDO 367 manufacturer identity and public-catalogue dimensions; official 2D/3D downloads require reserved-area authentication and were not acquired; not a project shop drawing and not official CAD geometry"
PRODUCT_PAGE = "https://www.rodaonline.com/en/collections/bernardo/"
CATALOGUE = "https://www.rodaonline.com/wp-content/uploads/RODA_CATALOGUE_2024.pdf"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/tab02"
CATALOGUE_EXTRACT = PRODUCT_DIR / "official-source/RODA-official-Bernardo-367-catalogue-extract.pdf"
PRODUCT_PAGE_SNAPSHOT = PRODUCT_DIR / "official-source/roda-bernardo-official-product-page.html"
EXPECTED_PATH_COUNTS = {"plan": 1, "front": 1, "side": 1}

shared.REPRESENTATIVE_GLOBAL_ID = "2eX84IyLr8_e34nuOoRPLQ"
shared.IFC_TYPE_NAME = "TAB02"
shared.PROFILE_KEY = "tab02"
shared.SCOPE = SCOPE
shared.PRODUCT_PAGE = PRODUCT_PAGE
shared.TECHNICAL_SHEET = CATALOGUE
shared.REPRESENTATIONS = {
    "plan": ("Tab02Plan", "PLAN_VIEW"),
    "front": ("Tab02Front", "ELEVATION_VIEW"),
    "side": ("Tab02Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/tab02-drawing-approval.json"
shared.PSET_NAME = "Pset_Tab02DrawingSource"
shared.DOCUMENT_ID_PREFIX = "RODA-BERNARDO-367-TAB02-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "RODA BERNARDO 367 / TAB02"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending TAB02 candidate source gate failed")
    drawing_source = access.get("drawing_geometry_source", {})
    official_cad = access.get("official_product_cad", {})
    dimension_check = access.get("dimension_cross_check", {})
    if (
        official_cad.get("authentication_required") is not True
        or official_cad.get("acquired") is not False
        or official_cad.get("local_cad_files") != []
        or official_cad.get("third_party_cad_used") is not False
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or dimension_check.get("status") != "exact_model_geometry_consistent_no_scaling_or_fitting"
        or dimension_check.get("maximum_absolute_delta_mm") != 0.0
        or access.get("scope") != SCOPE
        or sha256(CATALOGUE_EXTRACT) != "e4e58ccefd17c9fb50942b992c75d5018c0d30ef2e3f309dab4590c16aff32d1"
        or sha256(PRODUCT_PAGE_SNAPSHOT) != "8f04c1ddaa757b07a0d4cc62b419990fb10e22c073c62b889247a303d96c2318"
    ):
        raise RuntimeError("TAB02 official-source access record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"TAB02 {view} linework must remain geometry-derived with no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"TAB02 {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path):
    documents = (
        (
            "RODA-BERNARDO-367-TAB02-OFFICIAL-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "RODA official Bernardo collection page",
            "Manufacturer, family, designer and reserved-area CAD access evidence; current public page does not expose model 367 drawing geometry",
        ),
        (
            "RODA-BERNARDO-367-TAB02-OFFICIAL-CATALOGUE",
            CATALOGUE,
            "RODA official Catalogue 2024",
            "Exact BERNARDO 367 identity and 500 x 500 x 670 mm dimensions; not the source of geometry-derived linework",
        ),
        (
            "RODA-BERNARDO-367-TAB02-OFFICIAL-CATALOGUE-EXTRACT",
            relative(CATALOGUE_EXTRACT),
            "RODA official Bernardo 367 catalogue evidence extract",
            f"SHA-256 {sha256(CATALOGUE_EXTRACT)}; physical pages 28 and 192 visually checked; identity and dimensional evidence only",
        ),
        (
            "RODA-BERNARDO-367-TAB02-OFFICIAL-PAGE-SNAPSHOT",
            relative(PRODUCT_PAGE_SNAPSHOT),
            "RODA Bernardo official-page snapshot",
            f"SHA-256 {sha256(PRODUCT_PAGE_SNAPSHOT)}; records family identity and authenticated 2D/3D access boundary",
        ),
        (
            "RODA-BERNARDO-367-TAB02-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "RODA Bernardo 367 CAD access record",
            f"SHA-256 {sha256(access_path)}; reserved-area authentication required; native CAD not acquired; no third-party CAD used",
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
        "Manufacturer": "RODA",
        "Family": "Bernardo",
        "Model": "BERNARDO 367",
        "Designer": "Rodolfo Dordoni",
        "ProjectTypeCode": "TAB02",
        "IFCTypeDescription": "Bernardo 367",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceCatalogue": CATALOGUE,
        "Official2D3DStatus": "reserved_area_authentication_required_not_acquired",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Tab02Plan;Tab02Front;Tab02Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialCatalogueExtract": relative(CATALOGUE_EXTRACT),
        "OfficialCatalogueExtractSha256": sha256(CATALOGUE_EXTRACT),
        "OfficialProductPageSnapshot": relative(PRODUCT_PAGE_SNAPSHOT),
        "OfficialProductPageSnapshotSha256": sha256(PRODUCT_PAGE_SNAPSHOT),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "OfficialModelDimensionsMm": json.dumps(dimensions["official_model_width_depth_height_mm"]),
        "ProjectIfcBodyDimensionsMm": json.dumps(dimensions["project_ifc_body_local_xyz_mm"]),
        "BodyMinusOfficialMm": json.dumps(dimensions["body_minus_official_mm"]),
        "DimensionReviewStatus": dimensions["status"],
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
        Description="Mechanically verifiable geometry-derived drawing source and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="TAB02 drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
