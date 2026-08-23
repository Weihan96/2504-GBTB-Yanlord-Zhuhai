#!/usr/bin/env python3
"""Approval-gated derived-IFC writer for project FAU02 linework."""

import json
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
SCOPE = "official Falper Cilindro GH2 nearest family candidate only; project FAU02 exact article is not mechanically established; GH2 CAD is not used as representation and is not a project shop drawing"
PRODUCT_PAGE = "https://falper.it/rubinetteria-cilindro/"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/fau02"
EXPECTED_PATH_COUNTS = {"plan": 31, "front": 1, "side": 3}
SOURCE_FILES = {
    "GH2-2D-DWG": "Falper-Cilindro-GH2-2D.dwg",
    "GH2-3D-DWG": "Falper-Cilindro-GH2-3D.dwg",
    "GH2-TECHNICAL-PDF": "Falper-Cilindro-GH2-technical-sheet.pdf",
    "GH2-NATIVE-DWG-REFERENCE-SVG": "Falper-Cilindro-GH2-native-dwg-reference.svg",
    "OFFICIAL-PRODUCT-PAGE-SNAPSHOT": "falper-cilindro-official-product-page.html",
}

shared.REPRESENTATIVE_GLOBAL_ID = "36ZX3QPyD7SvlXsDKMP8rY"
shared.IFC_TYPE_NAME = "FAU02"
shared.PROFILE_KEY = "fau02"
shared.SOURCE_KIND = SOURCE_KIND
shared.SOURCE_LABEL_ZH = SOURCE_LABEL_ZH
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("FAU02Plan", "PLAN_VIEW"),
    "front": ("FAU02Front", "ELEVATION_VIEW"),
    "side": ("FAU02Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/fau02-drawing-approval.json"
shared.PSET_NAME = "Pset_FAU02DrawingSource"
shared.DOCUMENT_ID_PREFIX = "FALPER-FAU02-"


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("source_label_zh") != SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending FAU02 candidate source gate failed")
    official = access.get("official_product_cad", {})
    drawing = access.get("drawing_geometry_source", {})
    nearest = access.get("nearest_official_family_candidate", {})
    analysis_path = PRODUCT_DIR / "official-source/cad-match-analysis.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    if (
        access.get("scope") != SCOPE
        or access.get("pass") is not True
        or official.get("acquired") is not True
        or official.get("nearest_candidate_article") != "GH2"
        or official.get("exact_project_configuration_match") is not False
        or official.get("used_as_representation") is not False
        or official.get("third_party_cad_used") is not False
        or nearest.get("exact_project_configuration_match") is not False
        or drawing.get("source_kind") != SOURCE_KIND
        or drawing.get("official_cad_used") is not False
        or drawing.get("third_party_cad_used") is not False
        or analysis.get("exact_project_configuration_match") is not False
        or analysis.get("official_CAD_used_as_representation") is not False
    ):
        raise RuntimeError("FAU02 official-source mismatch gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if item.get("source_kind") != SOURCE_KIND or item.get("official_cad_paths_mm") != []:
            raise RuntimeError(f"FAU02 {view} must use Body-derived paths with no official CAD paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"FAU02 {view} geometry-derived path count drifted")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    source = PRODUCT_DIR / "official-source"
    match_analysis = source / "cad-match-analysis.json"
    documents = [
        (
            "FALPER-FAU02-OFFICIAL-CILINDRO-PRODUCT-PAGE",
            PRODUCT_PAGE,
            "Falper official Cilindro product page",
            "Manufacturer family evidence and official download surface; not representation geometry",
        ),
        (
            "FALPER-FAU02-SOURCE-ACCESS-RECORD",
            relative(access_path),
            "FAU02 official-source access record",
            f"SHA-256 {sha256(access_path)}; GH2 nearest candidate and exact-project mismatch",
        ),
        (
            "FALPER-FAU02-CAD-MATCH-ANALYSIS",
            relative(match_analysis),
            "FAU02 to Falper GH2 mechanical comparison",
            f"SHA-256 {sha256(match_analysis)}; exact project article rejected by independent dimensions",
        ),
    ]
    for suffix, filename in SOURCE_FILES.items():
        path = source / filename
        documents.append((
            f"FALPER-FAU02-{suffix}",
            relative(path),
            f"Falper official nearest-candidate evidence: {filename}",
            f"SHA-256 {sha256(path)}; GH2 family comparison only; not representation geometry",
        ))
    identifiers = []
    for identification, location, name, description in documents:
        reference = model.create_entity(
            "IfcDocumentReference",
            Location=location,
            Identification=identification,
            Name=name,
            Description=f"{description}; {SCOPE}",
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
    source = PRODUCT_DIR / "official-source"
    analysis_path = source / "cad-match-analysis.json"
    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    values = {
        "SourceKind": SOURCE_KIND,
        "SourceLabelZh": SOURCE_LABEL_ZH,
        "Manufacturer": "Falper",
        "OfficialFamily": "Cilindro",
        "NearestOfficialCandidate": "GH2",
        "ProjectIFCTypeName": shared.IFC_TYPE_NAME,
        "SourceProductPage": PRODUCT_PAGE,
        "OfficialCadAcquired": "true",
        "OfficialCadExactProjectMatch": "false",
        "OfficialCadUsed": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "FAU02Plan;FAU02Front;FAU02Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "CadMatchAnalysis": relative(analysis_path),
        "CadMatchAnalysisSha256": sha256(analysis_path),
        "OfficialGH2Native2DDwgSha256": sha256(source / SOURCE_FILES["GH2-2D-DWG"]),
        "OfficialGH2Native3DDwgSha256": sha256(source / SOURCE_FILES["GH2-3D-DWG"]),
        "IndependentDimensionDeltas": json.dumps(analysis["mechanical_deltas"], sort_keys=True),
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
        Description="Mechanically verifiable Body-derived drawing source, official GH2 mismatch evidence and human approval",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="FAU02 drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


if __name__ == "__main__":
    shared.main()
