#!/usr/bin/env python3
"""Approval-gated writer for official Geberit 224.212.00.2 drawings."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid
import ifcopenshell.util.element

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json
from geberit_duofix_sigma_224_212_linework import (
    ARTICLE,
    BOX_LABEL_SHA256,
    CATALOGUE_EXTRACT_SHA256,
    EXPECTED,
    EXPECTED_PATH_COUNTS,
    IFC_TYPE_NAME,
    LOCAL_BOX_LABEL,
    LOCAL_CATALOGUE_EXTRACT,
    OFFICIAL_CATALOGUE,
    SCOPE,
)


REPRESENTATIVE_GLOBAL_ID = "3hgNkx97vCTOC2eewCpMNk"
PROFILE_KEY = "geberit-duofix-sigma-224-212"
REQUIRED_VIEWS = {"plan", "front", "side"}
REPRESENTATIONS = {
    "plan": ("Duofix224212Plan", "PLAN_VIEW"),
    "front": ("Duofix224212Front", "ELEVATION_VIEW"),
    "side": ("Duofix224212Side", "ELEVATION_VIEW"),
}
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/geberit-duofix-sigma-224-212"
DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
DEFAULT_LINEWORK = PRODUCT_DIR / "official-native-dwg-linework.json"
DEFAULT_APPROVAL = ROOT / "pipeline/decisions/geberit-duofix-sigma-224-212-drawing-approval.json"


def require_approval(approval: dict, manifest_path: Path) -> None:
    errors = []
    if approval.get("status") != "approved":
        errors.append("status must be approved")
    if approval.get("derived_ifc_write_allowed") is not True:
        errors.append("derived_ifc_write_allowed must be true")
    if approval.get("formal_authoritative_ifc_write_allowed") is not False:
        errors.append("formal_authoritative_ifc_write_allowed must remain false")
    if set(approval.get("approved_views", [])) != REQUIRED_VIEWS:
        errors.append("approved_views must be exactly plan, front and side")
    if not str(approval.get("reviewer") or "").strip():
        errors.append("reviewer is required")
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(approval.get("review_date") or "")):
        errors.append("review_date must be YYYY-MM-DD")
    if approval.get("candidate_manifest_sha256") != sha256(manifest_path):
        errors.append("candidate_manifest_sha256 does not match the reviewed manifest")
    if approval.get("profile_key") != PROFILE_KEY or approval.get("article_number") != ARTICLE:
        errors.append(f"approval identity must be {PROFILE_KEY} / {ARTICLE}")
    if approval.get("scope") != SCOPE:
        errors.append("approval scope must remain the official exact-article family-reference scope")
    if not str(approval.get("approval_evidence") or "").strip():
        errors.append("approval_evidence is required")
    if errors:
        raise RuntimeError("approval gate rejected IFC write: " + "; ".join(errors))


def representation_context(model, identifier: str, target_view: str):
    parent = next(
        context
        for context in model.by_type("IfcGeometricRepresentationContext", include_subtypes=False)
        if context.ContextType == "Model"
    )
    return model.create_entity(
        "IfcGeometricRepresentationSubContext",
        ContextIdentifier=identifier,
        ContextType="Model",
        ParentContext=parent,
        TargetScale=None,
        TargetView=target_view,
        UserDefinedTargetView=None,
    )


def curve_representation(model, context, identifier: str, view: str, paths):
    polylines = []
    for path in paths:
        if len(path) < 2:
            continue
        points = []
        for first, second in path:
            if view == "plan":
                coordinates = (float(first), float(second), 0.0)
            elif view == "front":
                coordinates = (float(first), 0.0, float(second))
            else:
                coordinates = (0.0, float(first), float(second))
            points.append(model.create_entity("IfcCartesianPoint", Coordinates=coordinates))
        polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"approved {view} native-DWG path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    return model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier=identifier,
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )


def official_paths(candidate: dict, linework: dict) -> dict:
    if (
        candidate.get("profile_key") != PROFILE_KEY
        or candidate.get("source_kind") != "native_dwg"
        or candidate.get("official_cad_used") is not True
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("blue_line_present") is not True
        or candidate.get("representation_geometry_source")
        != "official_native_dwg_paths_mm"
        or candidate.get("official_overlay_source_kind") != "native_dwg"
        or candidate.get("official_article_number") != ARTICLE
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending Duofix candidate identity gate failed")
    if (
        linework.get("source_kind") != "native_dwg"
        or linework.get("article_number") != ARTICLE
        or linework.get("scope") != SCOPE
        or not linework.get("pass")
    ):
        raise RuntimeError("official Duofix native-DWG identity gate failed")
    paths = {}
    for view in REQUIRED_VIEWS:
        candidate_view = candidate.get("views", {}).get(view, {})
        linework_view = linework.get("views", {}).get(view, {})
        code = linework_view.get("native_dwg_code")
        if (
            candidate_view.get("native_dwg_code") != code
            or candidate_view.get("native_dwg_sha256") != EXPECTED[code]
            or linework_view.get("source_dwg_sha256") != EXPECTED[code]
        ):
            raise RuntimeError(f"official Duofix {view} native-DWG hash gate failed")
        paths[view] = candidate_view.get("official_native_dwg_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"official Duofix {view} path-count gate failed")
    return paths


def add_document_associations(model, product, product_type):
    documents = []
    for code, view in (("G", "plan"), ("A", "front elevation"), ("L", "side elevation")):
        documents.append((
            f"GEBERIT-224-212-{code}-NATIVE-DWG",
            f"output/review/highpoly-types/geberit-duofix-sigma-224-212/official-source/{ARTICLE}_{code}.dwg",
            f"Geberit official {ARTICLE} {code}.dwg ({view})",
            f"SHA-256 {EXPECTED[code]}; {SCOPE}",
        ))
    documents.append((
        "GEBERIT-224-212-P-NATIVE-DWG-IDENTITY",
        f"output/review/highpoly-types/geberit-duofix-sigma-224-212/official-source/{ARTICLE}_P.dwg",
        f"Geberit official {ARTICLE} P.dwg (3D identity reference)",
        f"SHA-256 {EXPECTED['P']}; official 3D identity evidence only; never substituted for Plan or Elevation linework",
    ))
    documents.extend((
        (
            "GEBERIT-224-212-OFFICIAL-CATALOGUE-EXTRACT",
            relative(LOCAL_CATALOGUE_EXTRACT),
            "Geberit official 224.212.00.2 catalogue extract",
            f"SHA-256 {CATALOGUE_EXTRACT_SHA256}; official catalogue {OFFICIAL_CATALOGUE}",
        ),
        (
            "GEBERIT-224-212-PROJECT-RECEIVED-BOX-LABEL",
            relative(LOCAL_BOX_LABEL),
            "Project received-box label for Geberit 224.212.00.2",
            f"SHA-256 {BOX_LABEL_SHA256}; identity evidence for one of two IFC family instances, not drawing geometry",
        ),
    ))
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


def add_source_pset(model, product, product_type, approval, manifest_hash, linework_path):
    values = {
        "SourceKind": "native_dwg",
        "Manufacturer": "Geberit",
        "Family": "Duofix Sigma concealed cistern frame",
        "ArticleNumber": ARTICLE,
        "IFCTypeName": IFC_TYPE_NAME,
        "OfficialCatalogue": OFFICIAL_CATALOGUE,
        "SourceDwgGPath": f"output/review/highpoly-types/geberit-duofix-sigma-224-212/official-source/{ARTICLE}_G.dwg",
        "SourceDwgGSha256": EXPECTED["G"],
        "SourceDwgAPath": f"output/review/highpoly-types/geberit-duofix-sigma-224-212/official-source/{ARTICLE}_A.dwg",
        "SourceDwgASha256": EXPECTED["A"],
        "SourceDwgLPath": f"output/review/highpoly-types/geberit-duofix-sigma-224-212/official-source/{ARTICLE}_L.dwg",
        "SourceDwgLSha256": EXPECTED["L"],
        "SourceDwgPPath": f"output/review/highpoly-types/geberit-duofix-sigma-224-212/official-source/{ARTICLE}_P.dwg",
        "SourceDwgPSha256": EXPECTED["P"],
        "OfficialCatalogueExtractSha256": CATALOGUE_EXTRACT_SHA256,
        "ProjectReceivedBoxLabelSha256": BOX_LABEL_SHA256,
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Duofix224212Plan;Duofix224212Front;Duofix224212Side",
        "RepresentationSourceMapping": "Duofix224212Plan=G;Duofix224212Front=A;Duofix224212Side=L",
        "RepresentationGeometrySource": "official_native_dwg_paths_mm",
        "ProxyGeometryIncluded": "false",
        "NativeDwgPlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "NativeDwgFrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "NativeDwgSidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "NativeDwgLineworkRegister": relative(linework_path),
        "NativeDwgLineworkRegisterSha256": sha256(linework_path),
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
        Name="Pset_GeberitDuofix224212DrawingSource",
        Description="Mechanically verifiable official CAD source and human approval for derived drawing representations",
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="Geberit Duofix 224.212.00.2 drawing source properties",
        Description=None,
        RelatedObjects=[product, product_type],
        RelatingPropertyDefinition=pset,
    )


def representation_path_count(representation) -> int:
    return sum(len(curve_set.Elements) for curve_set in representation.Items)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--official-linework", type=Path, default=DEFAULT_LINEWORK)
    parser.add_argument("--approval", type=Path, default=DEFAULT_APPROVAL)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    formal = args.input.resolve()
    output = args.output.resolve()
    manifest_path = args.manifest.resolve()
    linework_path = args.official_linework.resolve()
    if not args.apply:
        raise RuntimeError("IFC write requires the explicit --apply flag")
    if output == formal:
        raise RuntimeError("formal IFC cannot be the output; write a separate derived IFC")
    if output.exists():
        raise RuntimeError(f"refusing to overwrite existing output: {output}")
    formal_hash = sha256(formal)
    if formal_hash != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    manifest = load_json(manifest_path)
    candidate = load_json(args.candidate.resolve())
    linework = load_json(linework_path)
    approval = load_json(args.approval.resolve())
    require_approval(approval, manifest_path)
    if (
        manifest.get("formal_ifc_sha256") != formal_hash
        or manifest.get("formal_ifc_bytes_unchanged") is not True
        or manifest.get("source_kind") != "native_dwg"
        or manifest.get("official_cad_used") is not True
        or manifest.get("representation_geometry_source")
        != "official_native_dwg_paths_mm"
        or manifest.get("review_status") != "visual_review_pending"
        or manifest.get("approved_for_drawing_ifc") is not False
        or manifest.get("representative_global_id") != REPRESENTATIVE_GLOBAL_ID
        or manifest.get("bonsai_review", {}).get("mode") != "actual_bonsai_ifc_body_camera_render"
        or manifest.get("project_context", {}).get("pass") is not True
    ):
        raise RuntimeError("review manifest is not the expected pending Duofix candidate")
    paths = official_paths(candidate, linework)
    model = ifcopenshell.open(formal)
    product = model.by_guid(REPRESENTATIVE_GLOBAL_ID)
    if product is None:
        raise RuntimeError("Duofix representative product is missing")
    product_type = next(relation.RelatingType for relation in product.IsTypedBy)
    if product_type.Name != IFC_TYPE_NAME:
        raise RuntimeError("Duofix representative type identity drifted")
    existing_identifiers = {item.RepresentationIdentifier for item in product.Representation.Representations}
    if any(identifier in existing_identifiers for identifier, _ in REPRESENTATIONS.values()):
        raise RuntimeError("Duofix derived drawing representations already exist")
    representations = list(product.Representation.Representations)
    for view, (identifier, target_view) in REPRESENTATIONS.items():
        context = representation_context(model, identifier, target_view)
        representations.append(curve_representation(model, context, identifier, view, paths[view]))
    product.Representation.Representations = representations
    document_ids = add_document_associations(model, product, product_type)
    add_source_pset(model, product, product_type, approval, sha256(manifest_path), linework_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    model.write(output)
    if sha256(formal) != formal_hash:
        output.unlink(missing_ok=True)
        raise RuntimeError("formal IFC bytes changed during derived write")
    derived = ifcopenshell.open(output)
    derived_product = derived.by_guid(REPRESENTATIVE_GLOBAL_ID)
    derived_type = next(relation.RelatingType for relation in derived_product.IsTypedBy)
    representation_counts = {}
    for view, (identifier, _) in REPRESENTATIONS.items():
        representation = next((item for item in derived_product.Representation.Representations if item.RepresentationIdentifier == identifier), None)
        if representation is None:
            output.unlink(missing_ok=True)
            raise RuntimeError(f"derived IFC lost {identifier}")
        representation_counts[view] = representation_path_count(representation)
    if representation_counts != EXPECTED_PATH_COUNTS:
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC native-DWG path counts drifted")
    document_relations = [
        relation
        for relation in derived.get_inverse(derived_product)
        if relation.is_a("IfcRelAssociatesDocument")
        and relation.RelatingDocument.is_a("IfcDocumentReference")
        and (relation.RelatingDocument.Identification or "").startswith("GEBERIT-224-212-")
    ]
    actual_document_ids = sorted(relation.RelatingDocument.Identification for relation in document_relations)
    if actual_document_ids != document_ids or any(derived_type not in relation.RelatedObjects for relation in document_relations):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC Duofix document associations failed post-write verification")
    source_pset = ifcopenshell.util.element.get_pset(derived_product, "Pset_GeberitDuofix224212DrawingSource") or {}
    expected_properties = {
        "SourceKind": "native_dwg",
        "ArticleNumber": ARTICLE,
        "SourceDwgGSha256": EXPECTED["G"],
        "SourceDwgASha256": EXPECTED["A"],
        "SourceDwgLSha256": EXPECTED["L"],
        "SourceDwgPSha256": EXPECTED["P"],
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovedCandidateManifestSha256": sha256(manifest_path),
        "RepresentationGeometrySource": "official_native_dwg_paths_mm",
        "ProxyGeometryIncluded": "false",
    }
    if any(source_pset.get(name) != value for name, value in expected_properties.items()):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC Duofix source property post-write gate failed")
    result = {
        "schema_version": 1,
        "profile_key": PROFILE_KEY,
        "article_number": ARTICLE,
        "representative_global_id": REPRESENTATIVE_GLOBAL_ID,
        "formal_ifc": relative(formal),
        "formal_ifc_sha256": formal_hash,
        "formal_ifc_bytes_unchanged": True,
        "derived_ifc": relative(output),
        "derived_ifc_sha256": sha256(output),
        "candidate_manifest": relative(manifest_path),
        "candidate_manifest_sha256": sha256(manifest_path),
        "approval_record": relative(args.approval.resolve()),
        "approval_record_sha256": sha256(args.approval.resolve()),
        "approval": {
            "status": approval["status"],
            "reviewer": approval["reviewer"],
            "review_date": approval["review_date"],
            "approved_views": approval["approved_views"],
            "scope": approval["scope"],
        },
        "representations": {view: identifier for view, (identifier, _) in REPRESENTATIONS.items()},
        "representation_path_counts": representation_counts,
        "representation_geometry_source": source_pset["RepresentationGeometrySource"],
        "proxy_geometry_included": False,
        "official_linework_register": relative(linework_path),
        "official_linework_register_sha256": sha256(linework_path),
        "source_document_associations": document_ids,
        "source_property_set": "Pset_GeberitDuofix224212DrawingSource",
        "source_kind": source_pset["SourceKind"],
        "pass": True,
    }
    if args.report:
        write_json(args.report.resolve(), result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
