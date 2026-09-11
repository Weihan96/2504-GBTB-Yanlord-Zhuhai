#!/usr/bin/env python3
"""Approval-gated derived-IFC writer for Baxter Marilyn 01 native DWG views."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid
import ifcopenshell.util.element

from falper_sorgente_linework import ROOT, load_json, relative, sha256, write_json


REPRESENTATIVE_GLOBAL_ID = "3l2Ji4k2H9oOTDGWYUq7uV"
IFC_TYPE_NAME = "Marilyn 01"
IFC_TYPE_DESCRIPTION = "Bergère armchair with swivel base W86D100H78"
PROFILE_KEY = "marilyn-01"
SOURCE_KIND = "native_dwg"
SOURCE_LABEL_ZH = "基于 Baxter 精确型号原生 DWG 的官方图纸表达"
SCOPE = "exact Baxter Marilyn bergere 86 x 100 x 94 cm family CAD reference; not a project shop drawing"
FORMAL_SHA256 = "7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c"
DWG_SHA256 = "cb9825eecab7c78ee28e7668e933c2fe413395a63bf75f3afab8cc4959864724"
EXPECTED_PATH_COUNTS = {"plan": 44, "front": 96, "side": 78}
REQUIRED_VIEWS = set(EXPECTED_PATH_COUNTS)
REPRESENTATIONS = {
    "plan": ("Marilyn01Plan", "PLAN_VIEW"),
    "front": ("Marilyn01Front", "ELEVATION_VIEW"),
    "side": ("Marilyn01Side", "ELEVATION_VIEW"),
}
PRODUCT_PAGE = "https://www.baxter.it/en/products/marilyn-sofas-and-armchairs"
NATIVE_ZIP_URL = "https://dam.baxter.it/asset/9b02ac3a-7166-40df-9acc-1e8a2ad46705/Baxter_Marilyn_Armchair_2D_3D.zip"
TECHNICAL_SHEET = "https://productsbook.baxter.it/product-pdf/Marilyn_divani-e-poltrone_TechnicalSheet.pdf?code=MARI&lang=eng&sector=divani-e-poltrone&kind=indoor"
MATTE_SVG = "https://productsbook.baxter.it/models/measurements/MARIPBMN86.svg"
GLOSSY_SVG = "https://productsbook.baxter.it/models/measurements/MARIPBML86.svg"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/marilyn-01"
DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
DEFAULT_LINEWORK = PRODUCT_DIR / "official-native-dwg-linework.json"
DEFAULT_ACCESS = PRODUCT_DIR / "official-source/source-access-record.json"
DEFAULT_APPROVAL = ROOT / "pipeline/decisions/marilyn-01-drawing-approval.json"
PSET_NAME = "Pset_Marilyn01DrawingSource"
DOCUMENT_ID_PREFIX = "BAXTER-MARILYN-01-"


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
    if approval.get("profile_key") != PROFILE_KEY or approval.get("ifc_type_name") != IFC_TYPE_NAME:
        errors.append(f"approval identity must be {PROFILE_KEY} / {IFC_TYPE_NAME}")
    if approval.get("scope") != SCOPE:
        errors.append("approval scope must preserve the exact family-reference limitation")
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
            coordinates = (
                (float(first), float(second), 0.0)
                if view == "plan"
                else (float(first), 0.0, float(second))
                if view == "front"
                else (0.0, float(first), float(second))
            )
            points.append(model.create_entity("IfcCartesianPoint", Coordinates=coordinates))
        polylines.append(model.create_entity("IfcPolyline", Points=points))
    if len(polylines) != EXPECTED_PATH_COUNTS[view]:
        raise RuntimeError(f"approved Marilyn {view} native-DWG path count drifted")
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=polylines)
    return model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier=identifier,
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )


def official_paths(candidate: dict, linework: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != PROFILE_KEY
        or candidate.get("representative_global_id") != REPRESENTATIVE_GLOBAL_ID
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("source_label_zh") != SOURCE_LABEL_ZH
        or candidate.get("official_cad_used") is not True
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
        or candidate.get("source_dwg_sha256") != DWG_SHA256
    ):
        raise RuntimeError("pending Marilyn candidate source gate failed")
    if (
        linework.get("source_kind") != SOURCE_KIND
        or linework.get("source_dwg_sha256") != DWG_SHA256
        or linework.get("pass") is not True
        or access.get("native_cad_selection", {}).get("source_dwg_sha256") != DWG_SHA256
        or access.get("native_cad_selection", {}).get("three_view_path_counts") != EXPECTED_PATH_COUNTS
        or access.get("native_cad_selection", {}).get("official_cad_used") is not True
        or access.get("native_cad_selection", {}).get("third_party_cad_used") is not False
        or access.get("dimension_cross_check", {}).get("pass") is not True
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
    ):
        raise RuntimeError("Marilyn native-DWG/source record gate failed")
    paths = {}
    for view in REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        paths[view] = item.get("official_native_dwg_paths_mm", [])
        if (
            item.get("source_kind") != SOURCE_KIND
            or item.get("source_dwg_sha256") != DWG_SHA256
            or item.get("third_party_cad_used") is not False
            or item.get("geometry_scaled_to_match_ifc") is not False
            or len(paths[view]) != EXPECTED_PATH_COUNTS[view]
        ):
            raise RuntimeError(f"Marilyn {view} native-DWG candidate gate failed")
    return paths


def document_records(access_path: Path, linework_path: Path):
    source = PRODUCT_DIR / "official-source"
    return (
        ("BAXTER-MARILYN-01-OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Baxter official Marilyn product page", "Current bergere identity, designer, dimensions and native package link"),
        ("BAXTER-MARILYN-01-OFFICIAL-PRODUCT-PAGE-ARCHIVE", relative(source / "baxter-marilyn-product-page.html"), "Archived Baxter Marilyn product page", f"SHA-256 {sha256(source / 'baxter-marilyn-product-page.html')}"),
        ("BAXTER-MARILYN-01-OFFICIAL-NATIVE-ZIP", NATIVE_ZIP_URL, "Baxter official Marilyn native 2D/3D package", "Manufacturer package URL"),
        ("BAXTER-MARILYN-01-OFFICIAL-NATIVE-ZIP-ARCHIVE", relative(source / "Baxter_Marilyn_Armchair_2D_3D.zip"), "Archived Baxter Marilyn native 2D/3D package", f"SHA-256 {sha256(source / 'Baxter_Marilyn_Armchair_2D_3D.zip')}"),
        ("BAXTER-MARILYN-01-NATIVE-DWG", relative(source / "Marilyn_Abaco.dwg"), "Baxter official Marilyn_Abaco.dwg", f"SHA-256 {DWG_SHA256}; authoritative Plan/Front/Side linework"),
        ("BAXTER-MARILYN-01-EXACT-3DS", relative(source / "Marilyn_bergere_86x100xh94.3ds"), "Baxter exact Marilyn bergere 3DS", f"SHA-256 {sha256(source / 'Marilyn_bergere_86x100xh94.3ds')}; identity evidence only"),
        ("BAXTER-MARILYN-01-OFFICIAL-TECHNICAL-SHEET", TECHNICAL_SHEET, "Baxter official Marilyn technical sheet", "Page 13 confirms both 86 x 100 x 94 cm bergere finishes"),
        ("BAXTER-MARILYN-01-OFFICIAL-TECHNICAL-SHEET-ARCHIVE", relative(source / "Baxter_Marilyn_current-technical-sheet.pdf"), "Archived Baxter Marilyn technical sheet", f"SHA-256 {sha256(source / 'Baxter_Marilyn_current-technical-sheet.pdf')}"),
        ("BAXTER-MARILYN-01-MATTE-MEASUREMENT-SVG", MATTE_SVG, "Baxter official matte bergere measurement SVG", f"Archived SHA-256 {sha256(source / 'MARIPBMN86.svg')}"),
        ("BAXTER-MARILYN-01-GLOSSY-MEASUREMENT-SVG", GLOSSY_SVG, "Baxter official glossy bergere measurement SVG", f"Archived SHA-256 {sha256(source / 'MARIPBML86.svg')}"),
        ("BAXTER-MARILYN-01-SOURCE-ACCESS-RECORD", relative(access_path), "Baxter Marilyn source access record", f"SHA-256 {sha256(access_path)}"),
        ("BAXTER-MARILYN-01-NATIVE-LINEWORK-REGISTER", relative(linework_path), "Baxter Marilyn native-DWG linework register", f"SHA-256 {sha256(linework_path)}; OCS-corrected path counts 44/96/78"),
    )


def add_document_associations(model, product, product_type, access_path: Path, linework_path: Path):
    identifiers = []
    for identification, location, name, description in document_records(access_path, linework_path):
        reference = model.create_entity("IfcDocumentReference", Location=location, Identification=identification, Name=name, Description=f"{description}; {SCOPE}", ReferencedDocument=None)
        model.create_entity("IfcRelAssociatesDocument", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=f"{name} association", Description=SCOPE, RelatedObjects=[product, product_type], RelatingDocument=reference)
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval: dict, manifest_path: Path, candidate_path: Path, access_path: Path, linework_path: Path):
    access = load_json(access_path)
    values = {
        "SourceKind": SOURCE_KIND,
        "SourceLabelZh": SOURCE_LABEL_ZH,
        "Manufacturer": "Baxter",
        "Family": "Marilyn",
        "Designer": "Draga & Aurel",
        "ModelCode": "Marilyn 01",
        "ResolvedOfficialVariant": access["resolved_variant"],
        "IFCTypeDescription": IFC_TYPE_DESCRIPTION,
        "IfcDescriptionConflict": "The 78 cm description is stale; official native files, current technical sheet and Body identify the 94 cm bergere",
        "SourceProductPage": PRODUCT_PAGE,
        "SourceNativePackage": NATIVE_ZIP_URL,
        "SourceTechnicalSheet": TECHNICAL_SHEET,
        "SourceDwgPath": relative(PRODUCT_DIR / "official-source/Marilyn_Abaco.dwg"),
        "SourceDwgSha256": DWG_SHA256,
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": sha256(manifest_path),
        "ApprovedRepresentationIdentifiers": "Marilyn01Plan;Marilyn01Front;Marilyn01Side",
        "RepresentationGeometrySource": "official_native_dwg_paths_mm",
        "OfficialCadGeometryIncluded": "true",
        "ProxyGeometryIncluded": "false",
        "ThirdPartyCadUsed": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "NativeDwgLineworkRegister": relative(linework_path),
        "NativeDwgLineworkRegisterSha256": sha256(linework_path),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "SideOrientationTransform": "mirror_x_for_ifc_yz_direction",
        "GeometryScaledToMatchIfc": "false",
        "OfficialNominalWidthDepthHeightMm": json.dumps(access["dimension_cross_check"]["official_nominal_width_depth_height_mm"]),
        "ProjectIfcBodyLocalXYZMm": json.dumps(access["dimension_cross_check"]["project_ifc_body_local_xyz_mm"]),
    }
    properties = [model.create_entity("IfcPropertySingleValue", Name=name, Description=None, NominalValue=model.create_entity("IfcText", str(value)), Unit=None) for name, value in values.items()]
    pset = model.create_entity("IfcPropertySet", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=PSET_NAME, Description="Mechanically verifiable Baxter native-DWG source and human approval", HasProperties=properties)
    model.create_entity("IfcRelDefinesByProperties", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name="Marilyn 01 drawing source properties", Description=None, RelatedObjects=[product, product_type], RelatingPropertyDefinition=pset)


def representation_path_count(representation) -> int:
    return sum(len(curve_set.Elements) for curve_set in representation.Items)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--official-linework", type=Path, default=DEFAULT_LINEWORK)
    parser.add_argument("--source-access-record", type=Path, default=DEFAULT_ACCESS)
    parser.add_argument("--approval", type=Path, default=DEFAULT_APPROVAL)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    formal, output = args.input.resolve(), args.output.resolve()
    manifest_path, candidate_path = args.manifest.resolve(), args.candidate.resolve()
    linework_path, access_path = args.official_linework.resolve(), args.source_access_record.resolve()
    if not args.apply:
        raise RuntimeError("IFC write requires the explicit --apply flag")
    if output == formal:
        raise RuntimeError("formal IFC cannot be the output; write a separate derived IFC")
    if output.exists():
        raise RuntimeError(f"refusing to overwrite existing output: {output}")
    formal_hash = sha256(formal)
    if formal_hash != FORMAL_SHA256:
        raise RuntimeError("formal IFC hash mismatch")
    manifest, candidate = load_json(manifest_path), load_json(candidate_path)
    linework, access = load_json(linework_path), load_json(access_path)
    approval = load_json(args.approval.resolve())
    require_approval(approval, manifest_path)
    if (
        manifest.get("formal_ifc_sha256") != formal_hash
        or manifest.get("formal_ifc_bytes_unchanged") is not True
        or manifest.get("review_status") != "visual_review_pending"
        or manifest.get("approved_for_drawing_ifc") is not False
        or manifest.get("representative_global_id") != REPRESENTATIVE_GLOBAL_ID
        or manifest.get("drawing_source", {}).get("source_kind") != SOURCE_KIND
        or manifest.get("project_context", {}).get("pass") is not True
        or manifest.get("project_context", {}).get("manifest_sha256") != sha256(PRODUCT_DIR / "project-context-manifest.json")
        or manifest.get("bonsai_review", {}).get("mode") != "actual_bonsai_ifc_body_camera_render"
        or manifest.get("bonsai_review", {}).get("manifest_sha256") != sha256(PRODUCT_DIR / "bonsai-review-manifest.json")
    ):
        raise RuntimeError("review manifest is not the expected pending Marilyn candidate")
    paths = official_paths(candidate, linework, access)
    model = ifcopenshell.open(formal)
    product = model.by_guid(REPRESENTATIVE_GLOBAL_ID)
    if product is None:
        raise RuntimeError("Marilyn representative product is missing")
    product_type = next(relation.RelatingType for relation in product.IsTypedBy)
    if product_type.Name != IFC_TYPE_NAME or product_type.Description != IFC_TYPE_DESCRIPTION:
        raise RuntimeError("Marilyn representative type identity drifted")
    existing = {representation.RepresentationIdentifier for representation in product.Representation.Representations}
    if any(identifier in existing for identifier, _ in REPRESENTATIONS.values()):
        raise RuntimeError("Marilyn derived drawing representations already exist")
    representations = list(product.Representation.Representations)
    for view, (identifier, target_view) in REPRESENTATIONS.items():
        representations.append(curve_representation(model, representation_context(model, identifier, target_view), identifier, view, paths[view]))
    product.Representation.Representations = representations
    document_ids = add_document_associations(model, product, product_type, access_path, linework_path)
    add_source_pset(model, product, product_type, approval, manifest_path, candidate_path, access_path, linework_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    model.write(output)
    if sha256(formal) != formal_hash:
        output.unlink(missing_ok=True)
        raise RuntimeError("formal IFC bytes changed during derived write")
    derived = ifcopenshell.open(output)
    derived_product = derived.by_guid(REPRESENTATIVE_GLOBAL_ID)
    derived_type = next(relation.RelatingType for relation in derived_product.IsTypedBy)
    counts = {}
    for view, (identifier, _) in REPRESENTATIONS.items():
        representation = next((item for item in derived_product.Representation.Representations if item.RepresentationIdentifier == identifier), None)
        if representation is None:
            output.unlink(missing_ok=True)
            raise RuntimeError(f"derived IFC lost {identifier}")
        counts[view] = representation_path_count(representation)
    if counts != EXPECTED_PATH_COUNTS:
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC native-DWG path counts drifted")
    relations = [relation for relation in derived.get_inverse(derived_product) if relation.is_a("IfcRelAssociatesDocument") and relation.RelatingDocument.is_a("IfcDocumentReference") and (relation.RelatingDocument.Identification or "").startswith(DOCUMENT_ID_PREFIX)]
    actual_ids = sorted(relation.RelatingDocument.Identification for relation in relations)
    if actual_ids != document_ids or any(derived_type not in relation.RelatedObjects for relation in relations):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC Marilyn document associations failed post-write verification")
    pset = ifcopenshell.util.element.get_pset(derived_product, PSET_NAME) or {}
    expected = {"SourceKind": SOURCE_KIND, "SourceDwgSha256": DWG_SHA256, "EvidenceScope": SCOPE, "ReviewStatus": "APPROVED", "Reviewer": approval["reviewer"], "ReviewDate": approval["review_date"], "ApprovedCandidateManifestSha256": sha256(manifest_path), "RepresentationGeometrySource": "official_native_dwg_paths_mm", "OfficialCadGeometryIncluded": "true", "ProxyGeometryIncluded": "false", "ThirdPartyCadUsed": "false"}
    if any(pset.get(name) != value for name, value in expected.items()):
        output.unlink(missing_ok=True)
        raise RuntimeError("derived IFC Marilyn source-property gate failed")
    result = {
        "schema_version": 1,
        "profile_key": PROFILE_KEY,
        "ifc_type_name": IFC_TYPE_NAME,
        "representative_global_id": REPRESENTATIVE_GLOBAL_ID,
        "formal_ifc": relative(formal),
        "formal_ifc_sha256": formal_hash,
        "formal_ifc_bytes_unchanged": True,
        "derived_ifc": relative(output),
        "derived_ifc_sha256": sha256(output),
        "candidate_manifest": relative(manifest_path),
        "candidate_manifest_sha256": sha256(manifest_path),
        "representations": {view: identifier for view, (identifier, _) in REPRESENTATIONS.items()},
        "representation_path_counts": counts,
        "representation_geometry_source": pset["RepresentationGeometrySource"],
        "official_cad_geometry_included": True,
        "proxy_geometry_included": False,
        "official_linework_register": relative(linework_path),
        "official_linework_register_sha256": sha256(linework_path),
        "source_document_associations": document_ids,
        "source_property_set": PSET_NAME,
        "source_kind": pset["SourceKind"],
        "source_label_zh": pset["SourceLabelZh"],
        "pass": True,
    }
    if args.report:
        write_json(args.report.resolve(), result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
