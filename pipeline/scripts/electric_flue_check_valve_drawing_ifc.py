#!/usr/bin/env python3
"""Approval-gated derived IFC writer for the unresolved electric flue proxy."""

from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SCOPE = "project IFC generic identity and observable geometry only; not manufacturer product identification, not official CAD geometry, not a system design and not a project shop drawing"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/electric-flue-check-valve"
EXPECTED_PATH_COUNTS = {"plan": 2, "front": 5, "side": 10}

shared.REPRESENTATIVE_GLOBAL_ID = "1faflkXXH6M9cnYPE9Liir"
shared.IFC_TYPE_NAME = "Electric flue check valve"
shared.PROFILE_KEY = "electric-flue-check-valve"
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("ElectricFlueValvePlan", "PLAN_VIEW"),
    "front": ("ElectricFlueValveFront", "ELEVATION_VIEW"),
    "side": ("ElectricFlueValveSide", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/electric-flue-check-valve-drawing-approval.json"
shared.PSET_NAME = "Pset_ElectricFlueCheckValveDrawingSource"
shared.DOCUMENT_ID_PREFIX = "ELECTRIC-FLUE-VALVE-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "unresolved"
        or candidate.get("source_kind") != shared.SOURCE_KIND
        or candidate.get("source_label_zh") != shared.SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending electric flue proxy source gate failed")
    drawing_source = access.get("drawing_geometry_source", {})
    if (
        access.get("manufacturer") != "unresolved_from_project_IFC"
        or access.get("resolved_article_number") is not None
        or access.get("external_research_boundary", {}).get("exact_manufacturer_or_model_resolved") is not False
        or access.get("official_product_cad", {}).get("acquired") is not False
        or drawing_source.get("source_kind") != shared.SOURCE_KIND
        or drawing_source.get("official_cad_used") is not False
        or drawing_source.get("third_party_cad_used") is not False
        or access.get("scope") != SCOPE
    ):
        raise RuntimeError("electric flue proxy unresolved-source access gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != shared.SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"electric flue {view} must remain geometry-derived with zero official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"electric flue {view} geometry-derived path-count gate failed")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    documents = (
        (
            "ELECTRIC-FLUE-VALVE-M401-PROJECT-REVIEW",
            "pipeline/decisions/m401-existing-review.csv",
            "Existing M401 project review",
            "Project evidence preserving the generic named proxy and prohibiting unsupported product, connection, control and access inference",
        ),
        (
            "ELECTRIC-FLUE-VALVE-ELEC-PROJECT-REVIEW",
            "pipeline/decisions/elec-existing-review.csv",
            "Existing ELEC project review",
            "Project evidence that identity, mounting face, access direction and the purposes of two openings remain unresolved",
        ),
        (
            "ELECTRIC-FLUE-VALVE-RCP-PROJECT-REVIEW",
            "pipeline/decisions/rcp1-existing-review.csv",
            "Existing RCP project review",
            "Project evidence that the proxy has no IFC type, system or ports and must not be automatically converted",
        ),
        (
            "ELECTRIC-FLUE-VALVE-R04-PX-PROJECTION",
            "drawings/elevations/native/EL-03-07-R04-PX.svg",
            "Existing R04 +X model projection",
            "Project-context projection evidence only; not manufacturer CAD",
        ),
        (
            "ELECTRIC-FLUE-VALVE-R04-NY-PROJECTION",
            "drawings/elevations/native/EL-03-08-R04-NY.svg",
            "Existing R04 -Y model projection",
            "Project-context projection evidence only; not manufacturer CAD",
        ),
        (
            "ELECTRIC-FLUE-VALVE-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "Electric flue proxy source and research boundary",
            f"SHA-256 {sha256(access_path)}; manufacturer/model unresolved; exact official CAD not attributable; no third-party CAD used",
        ),
    )
    related = [product] if product_type is None else [product, product_type]
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
            RelatedObjects=related,
            RelatingDocument=reference,
        )
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    values = {
        "SourceKind": shared.SOURCE_KIND,
        "SourceLabelZh": shared.SOURCE_LABEL_ZH,
        "Manufacturer": "unresolved_from_project_IFC",
        "Family": "unresolved_from_project_IFC",
        "ArticleNumber": "not_provided",
        "IFCObjectName": shared.IFC_TYPE_NAME,
        "ProductIdentityStatus": "manufacturer_and_model_unresolved",
        "OfficialProductCadStatus": "no_exact_official_CAD_attributable",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "ElectricFlueValvePlan;ElectricFlueValveFront;ElectricFlueValveSide",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
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
        Description="Mechanically verifiable geometry-derived drawing source, unresolved product identity and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Electric flue proxy drawing source properties",
        Description=None,
        RelatedObjects=[product],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
