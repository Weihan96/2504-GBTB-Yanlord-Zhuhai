#!/usr/bin/env python3
"""Approval-gated derived-IFC writer for configured Geberit TRAP01 views."""

import json
import sys
from pathlib import Path

import ifcopenshell
import ifcopenshell.guid
import ifcopenshell.util.element

from falper_sorgente_linework import ROOT, relative, sha256
import hima01_drawing_ifc as shared


SOURCE_KIND = "geometry_derived_simplified_proxy"
SOURCE_LABEL_ZH = "基于原始高模几何生成的简化图纸表达"
SCOPE = "exact Geberit 151.116.11.1 adjustable family reference; project instance is a shortened installation configuration; not a project shop drawing"
PRODUCT_PAGE = "https://catalog.geberit.us/en-US/product/PRO_185224"
PRODUCT_DIR = ROOT / "output/review/highpoly-types/trap01"
EXPECTED_PATH_COUNTS = {"plan": 16, "front": 1, "side": 8}
DWG_FILES = {
    "A": "151.116.11.1_A.dwg",
    "G": "151.116.11.1_G.dwg",
    "L": "151.116.11.1_L.dwg",
    "P": "151.116.11.1_P.dwg",
}
EVIDENCE_FILES = {
    "PRODUCT-PAGE-ARCHIVE": "geberit-PRO_185224-product-page.html",
    "PRODUCT-DATA-SHEET": "Geberit-PRO_185224-product-data-sheet.pdf",
    "INSTALLATION-INSTRUCTIONS": "966.798.00.0-installation-instructions.pdf",
    "MAINTENANCE-MANUAL": "969.459.00.0-maintenance-manual.pdf",
}

shared.REPRESENTATIVE_GLOBAL_ID = "2Ak2ma0lvBEA49UpplzUqi"
shared.IFC_TYPE_NAME = "TRAP01"
shared.PROFILE_KEY = "trap01"
shared.SOURCE_KIND = SOURCE_KIND
shared.SOURCE_LABEL_ZH = SOURCE_LABEL_ZH
shared.SCOPE = SCOPE
shared.REPRESENTATIONS = {
    "plan": ("Trap01Plan", "PLAN_VIEW"),
    "front": ("Trap01Front", "ELEVATION_VIEW"),
    "side": ("Trap01Side", "ELEVATION_VIEW"),
}
shared.EXPECTED_PATH_COUNTS = EXPECTED_PATH_COUNTS
shared.PRODUCT_DIR = PRODUCT_DIR
shared.DEFAULT_MANIFEST = PRODUCT_DIR / "manifest.json"
shared.DEFAULT_CANDIDATE = PRODUCT_DIR / "candidate-representations.json"
shared.DEFAULT_ACCESS_RECORD = PRODUCT_DIR / "official-source/source-access-record.json"
shared.DEFAULT_APPROVAL = ROOT / "pipeline/decisions/trap01-drawing-approval.json"
shared.PSET_NAME = "Pset_Trap01DrawingSource"
shared.DOCUMENT_ID_PREFIX = "GEBERIT-151-116-11-1-TRAP01-"

COMPONENT_GROUP_NAME = "TRAP01 adjustable drawing components"
COMPONENT_PSET = "Pset_Trap01DrawingComponent"
CONFIGURATION_PSET = "Pset_Trap01AdjustableConfiguration"
BLUE_SOLID_SEMANTICS = "official native DWG linework retained within the current configured range"
BLUE_DASHED_SEMANTICS = "official adjustable-family excess beyond the current configured endpoint; reference only, not installed project geometry"
COMPONENTS = {
    "fixed_body": {
        "name": "TRAP01 fixed trap body",
        "adjustable": False,
        "axis": "FIXED_BODY_SPINE",
        "start_mm": [53.471272, 0.0, -80.499954],
        "end_mm": [53.471272, 0.0, 41.0],
        "installed_length_mm": 121.499954,
        "minimum_length_mm": 121.499954,
        "maximum_length_mm": 121.499954,
        "official_default_overflow_mm": 0.0,
    },
    "horizontal_adjustable": {
        "name": "TRAP01 horizontal adjustable outlet",
        "adjustable": True,
        "axis": "LOCAL_X",
        "start_mm": [109.288309, 0.0, 94.0],
        "end_mm": [222.690975, 0.0, 94.0],
        "installed_length_mm": 113.402666,
        "minimum_length_mm": 0.0,
        "maximum_length_mm": 252.0,
        "official_default_overflow_mm": 159.52933,
    },
    "vertical_adjustable": {
        "name": "TRAP01 vertical adjustable dip tube",
        "adjustable": True,
        "axis": "LOCAL_Z",
        "start_mm": [53.471272, 0.0, -80.499954],
        "end_mm": [53.471272, 0.0, 111.440025],
        "installed_length_mm": 191.939979,
        "minimum_length_mm": 85.0,
        "maximum_length_mm": 334.0,
        "official_default_overflow_mm": 208.00002,
    },
}


def candidate_paths(candidate: dict, access: dict) -> dict:
    if (
        candidate.get("profile_key") != shared.PROFILE_KEY
        or candidate.get("representative_global_id") != shared.REPRESENTATIVE_GLOBAL_ID
        or candidate.get("article_number") != "151.116.11.1"
        or candidate.get("source_kind") != SOURCE_KIND
        or candidate.get("source_label_zh") != SOURCE_LABEL_ZH
        or candidate.get("official_cad_acquired") is not True
        or candidate.get("official_cad_used") is not False
        or candidate.get("official_cad_used_as_representation") is not False
        or candidate.get("third_party_cad_used") is not False
        or candidate.get("formal_ifc_write_allowed") is not False
        or candidate.get("review_status") != "visual_review_pending"
    ):
        raise RuntimeError("pending TRAP01 candidate source gate failed")
    native = access.get("native_cad_selection", {})
    configuration = access.get("configuration_cross_check", {})
    article = access.get("article_resolution", {})
    if (
        access.get("resolved_article") != "151.116.11.1"
        or access.get("scope") != SCOPE
        or access.get("pass") is not True
        or native.get("source_kind") != "native_dwg"
        or native.get("official_cad_acquired") is not True
        or native.get("third_party_cad_used") is not False
        or native.get("pass") is not True
        or configuration.get("project_instance_is_shortened_configuration") is not True
        or configuration.get("official_default_family_paths_used_as_project_representation") is not False
        or configuration.get("geometry_scaled_or_stretched_to_match") is not False
        or article.get("selected_article_pass") is not True
    ):
        raise RuntimeError("TRAP01 official-source/configuration record gate failed")
    paths = {}
    for view in shared.REQUIRED_VIEWS:
        item = candidate.get("views", {}).get(view, {})
        if (
            item.get("source_kind") != SOURCE_KIND
            or item.get("official_cad_paths_mm") != []
            or item.get("official_family_paths_used_as_project_representation") is not False
            or item.get("geometry_scaled_or_stretched") is not False
        ):
            raise RuntimeError(f"TRAP01 {view} must use only configured Body-derived proxy paths")
        paths[view] = item.get("proxy_paths_mm", [])
        if len(paths[view]) != EXPECTED_PATH_COUNTS[view]:
            raise RuntimeError(f"TRAP01 {view} configured proxy path count drifted")
    return paths


def add_document_associations(model, product, product_type, access_path: Path):
    source = PRODUCT_DIR / "official-source"
    documents = [
        ("OFFICIAL-PRODUCT-PAGE", PRODUCT_PAGE, "Geberit official product page", "Exact family and article identity, attributes and native CAD links"),
        ("SOURCE-ACCESS-RECORD", relative(access_path), "TRAP01 official-source access record", f"SHA-256 {sha256(access_path)}; article selection and configured-instance limitation"),
        ("NATIVE-DWG-LINEWORK", relative(PRODUCT_DIR / "official-native-dwg-linework.json"), "TRAP01 native-DWG linework register", f"SHA-256 {sha256(PRODUCT_DIR / 'official-native-dwg-linework.json')}; family reference only"),
    ]
    for suffix, filename in DWG_FILES.items():
        path = source / filename
        documents.append((f"OFFICIAL-{suffix}-DWG", relative(path), f"Geberit official {suffix} native DWG", f"SHA-256 {sha256(path)}; adjustable-family reference only; not representation geometry"))
    for suffix, filename in EVIDENCE_FILES.items():
        path = source / filename
        documents.append((suffix, relative(path), f"Geberit official {filename}", f"SHA-256 {sha256(path)}; identity and configuration evidence only"))
    identifiers = []
    for suffix, location, name, description in documents:
        identification = f"GEBERIT-151-116-11-1-TRAP01-{suffix}"
        reference = model.create_entity("IfcDocumentReference", Location=location, Identification=identification, Name=name, Description=f"{description}; {SCOPE}", ReferencedDocument=None)
        model.create_entity("IfcRelAssociatesDocument", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=f"{name} association", Description=SCOPE, RelatedObjects=[product, product_type], RelatingDocument=reference)
        identifiers.append(identification)
    return sorted(identifiers)


def add_source_pset(model, product, product_type, approval, manifest_hash, candidate_path, access_path):
    access = json.loads(access_path.read_text(encoding="utf-8"))
    values = {
        "SourceKind": SOURCE_KIND,
        "SourceLabelZh": SOURCE_LABEL_ZH,
        "Manufacturer": "Geberit",
        "Family": access["family"],
        "ArticleNumber": "151.116.11.1",
        "ProjectIFCTypeName": shared.IFC_TYPE_NAME,
        "SourceProductPage": PRODUCT_PAGE,
        "OfficialCadAcquired": "true",
        "OfficialCadUsed": "false",
        "OfficialCadRole": "exact article identity, fixed-width cross-check and 1:1 solid/dashed adjustable review overlay; excluded from installed representation geometry",
        "ProjectConfiguration": "shortened installed configuration represented from isolated IFC Body",
        "GeometryScaledOrStretched": "false",
        "ThirdPartyCadUsed": "false",
        "EvidenceScope": SCOPE,
        "ReviewStatus": "APPROVED",
        "Reviewer": approval["reviewer"],
        "ReviewDate": approval["review_date"],
        "ApprovalEvidence": approval["approval_evidence"],
        "ApprovedCandidateManifestSha256": manifest_hash,
        "ApprovedRepresentationIdentifiers": "Trap01Plan;Trap01Front;Trap01Side",
        "RepresentationGeometrySource": "proxy_paths_mm derived from the configured isolated representative IFC Body",
        "OfficialCadGeometryIncluded": "false",
        "CandidateRepresentations": relative(candidate_path),
        "CandidateRepresentationsSha256": sha256(candidate_path),
        "SourceAccessRecord": relative(access_path),
        "SourceAccessRecordSha256": sha256(access_path),
        "OfficialNativeDwgSha256": json.dumps(access["native_cad_selection"]["sha256"], sort_keys=True),
        "PlanPathCount": str(EXPECTED_PATH_COUNTS["plan"]),
        "FrontPathCount": str(EXPECTED_PATH_COUNTS["front"]),
        "SidePathCount": str(EXPECTED_PATH_COUNTS["side"]),
        "ProjectIFCBodyLocalXYZMm": json.dumps(access["article_resolution"]["project_ifc_body_local_xyz_mm"]),
        "SelectedArticleFixedWidthDeltaMm": str(access["article_resolution"]["selected_fixed_width_absolute_delta_mm"]),
        "ExcludedAlternativeArticle": access["article_resolution"]["excluded_alternative"],
        "DrawingComponentModel": "fixed_body + horizontal_adjustable + vertical_adjustable",
        "HorizontalInstalledLengthMm": str(COMPONENTS["horizontal_adjustable"]["installed_length_mm"]),
        "HorizontalAdjustmentRangeMm": "0-252",
        "VerticalInstalledHMm": str(COMPONENTS["vertical_adjustable"]["installed_length_mm"]),
        "VerticalAdjustmentRangeMm": "85-334",
        "BlueSolidSemantics": BLUE_SOLID_SEMANTICS,
        "BlueDashedSemantics": BLUE_DASHED_SEMANTICS,
        "ReferenceExtensionDisplayedInSceneByDefault": "false",
    }
    properties = [model.create_entity("IfcPropertySingleValue", Name=name, Description=None, NominalValue=model.create_entity("IfcText", str(value)), Unit=None) for name, value in values.items()]
    pset = model.create_entity("IfcPropertySet", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name=shared.PSET_NAME, Description="Mechanically verifiable configured Body drawing source, official CAD evidence and human approval", HasProperties=properties)
    model.create_entity("IfcRelDefinesByProperties", GlobalId=ifcopenshell.guid.new(), OwnerHistory=product.OwnerHistory, Name="TRAP01 drawing source properties", Description=None, RelatedObjects=[product, product_type], RelatingPropertyDefinition=pset)


shared.candidate_paths = candidate_paths
shared.add_document_associations = add_document_associations
shared.add_source_pset = add_source_pset


def add_text_pset(model, owner_history, related_objects, name, description, values):
    properties = [
        model.create_entity(
            "IfcPropertySingleValue",
            Name=key,
            Description=None,
            NominalValue=model.create_entity("IfcText", str(value)),
            Unit=None,
        )
        for key, value in values.items()
    ]
    pset = model.create_entity(
        "IfcPropertySet",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=owner_history,
        Name=name,
        Description=description,
        HasProperties=properties,
    )
    model.create_entity(
        "IfcRelDefinesByProperties",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=owner_history,
        Name=f"{name} assignment",
        Description=None,
        RelatedObjects=related_objects,
        RelatingPropertyDefinition=pset,
    )
    return pset


def add_axis_representation(model, context, component):
    points = [
        model.create_entity("IfcCartesianPoint", Coordinates=tuple(float(value) for value in coordinates))
        for coordinates in (component["start_mm"], component["end_mm"])
    ]
    axis = model.create_entity("IfcPolyline", Points=points)
    curve_set = model.create_entity("IfcGeometricCurveSet", Elements=[axis])
    colour = model.create_entity(
        "IfcColourRgb", Name="TRAP01 current configuration black", Red=17.0 / 255.0, Green=24.0 / 255.0, Blue=32.0 / 255.0
    )
    style = model.create_entity(
        "IfcCurveStyle",
        Name="TRAP01 adjustable component semantic axis",
        CurveFont=None,
        CurveWidth=model.create_entity("IfcPositiveLengthMeasure", 0.35),
        CurveColour=colour,
        ModelOrDraughting=True,
    )
    model.create_entity("IfcStyledItem", Item=curve_set, Styles=[style], Name="Current installed component axis")
    representation = model.create_entity(
        "IfcShapeRepresentation",
        ContextOfItems=context,
        RepresentationIdentifier="Axis",
        RepresentationType="GeometricCurveSet",
        Items=[curve_set],
    )
    return model.create_entity("IfcProductDefinitionShape", Representations=[representation])


def persist_adjustable_components(output: Path, report_path=None):
    model = ifcopenshell.open(output)
    product = model.by_guid(shared.REPRESENTATIVE_GLOBAL_ID)
    if product is None:
        raise RuntimeError("TRAP01 derived product missing before component persistence")
    existing = [item for item in model.by_type("IfcAnnotation") if item.ObjectType == "DRAWING_COMPONENT"]
    if existing:
        raise RuntimeError("TRAP01 drawing components already exist")
    parent = next(
        context for context in model.by_type("IfcGeometricRepresentationContext", include_subtypes=False)
        if context.ContextType == "Model"
    )
    context = model.create_entity(
        "IfcGeometricRepresentationSubContext",
        ContextIdentifier="Trap01AdjustmentAxis",
        ContextType="Model",
        ParentContext=parent,
        TargetScale=None,
        TargetView="MODEL_VIEW",
        UserDefinedTargetView=None,
    )
    annotations = []
    component_records = []
    for role, component in COMPONENTS.items():
        annotation = model.create_entity(
            "IfcAnnotation",
            GlobalId=ifcopenshell.guid.new(),
            OwnerHistory=product.OwnerHistory,
            Name=component["name"],
            Description=(
                f"TRAP01 identifiable drawing component; role={role}; black/current installed configuration. "
                f"{BLUE_DASHED_SEMANTICS}."
            ),
            ObjectType="DRAWING_COMPONENT",
            ObjectPlacement=product.ObjectPlacement,
            Representation=add_axis_representation(model, context, component),
        )
        values = {
            "ComponentRole": role,
            "Adjustable": str(component["adjustable"]).lower(),
            "AdjustmentAxis": component["axis"],
            "InterfacePointMm": json.dumps(component["start_mm"]),
            "InstalledEndpointMm": json.dumps(component["end_mm"]),
            "InstalledLengthMm": component["installed_length_mm"],
            "MinimumLengthMm": component["minimum_length_mm"],
            "MaximumLengthMm": component["maximum_length_mm"],
            "OfficialDefaultOverflowMm": component["official_default_overflow_mm"],
            "CurrentInstalledLineSemantics": "black current project configuration",
            "BlueSolidSemantics": BLUE_SOLID_SEMANTICS,
            "BlueDashedSemantics": BLUE_DASHED_SEMANTICS,
            "DisplayedInSceneByDefault": "true",
            "ReferenceExtensionDisplayedInSceneByDefault": "false",
            "ParameterUpdateStrategy": "update InstalledLengthMm and InstalledEndpointMm, then regenerate dependent Plan/Front/Side drawing representations",
        }
        add_text_pset(
            model,
            product.OwnerHistory,
            [annotation],
            COMPONENT_PSET,
            "TRAP01 product-level adjustable drawing component parameters",
            values,
        )
        annotations.append(annotation)
        component_records.append({"role": role, "global_id": annotation.GlobalId, **component})
    group = model.create_entity(
        "IfcGroup",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name=COMPONENT_GROUP_NAME,
        Description="Product-level adjustable diagram decomposition; current installed geometry remains on the TRAP01 product",
        ObjectType="ADJUSTABLE_DRAWING_COMPONENT_GROUP",
    )
    model.create_entity(
        "IfcRelAssignsToGroup",
        GlobalId=ifcopenshell.guid.new(),
        OwnerHistory=product.OwnerHistory,
        Name="TRAP01 adjustable drawing component membership",
        Description=BLUE_DASHED_SEMANTICS,
        RelatedObjects=[product, *annotations],
        RelatingGroup=group,
    )
    add_text_pset(
        model,
        product.OwnerHistory,
        [product, group],
        CONFIGURATION_PSET,
        "TRAP01 approved current configuration and independent horizontal/vertical adjustment parameters",
        {
            "ComponentGroupGlobalId": group.GlobalId,
            "FixedBodyGlobalId": component_records[0]["global_id"],
            "HorizontalAdjustableGlobalId": component_records[1]["global_id"],
            "VerticalAdjustableGlobalId": component_records[2]["global_id"],
            "HorizontalInstalledLengthMm": COMPONENTS["horizontal_adjustable"]["installed_length_mm"],
            "HorizontalMinimumMm": 0.0,
            "HorizontalMaximumMm": 252.0,
            "VerticalInstalledHMm": COMPONENTS["vertical_adjustable"]["installed_length_mm"],
            "VerticalMinimumHMm": 85.0,
            "VerticalMaximumHMm": 334.0,
            "BlueDashedSemantics": BLUE_DASHED_SEMANTICS,
            "ReferenceExtensionDisplayedInSceneByDefault": "false",
        },
    )
    temporary = output.with_suffix(".ifc.components-next")
    model.write(temporary)
    temporary.replace(output)

    reopened = ifcopenshell.open(output)
    persisted = [item for item in reopened.by_type("IfcAnnotation") if item.ObjectType == "DRAWING_COMPONENT"]
    roles = {}
    for annotation in persisted:
        pset = ifcopenshell.util.element.get_pset(annotation, COMPONENT_PSET) or {}
        roles[pset.get("ComponentRole")] = {
            "global_id": annotation.GlobalId,
            "installed_length_mm": float(pset["InstalledLengthMm"]),
            "axis": pset["AdjustmentAxis"],
            "axis_representation_count": len(annotation.Representation.Representations),
        }
    if set(roles) != set(COMPONENTS) or any(item["axis_representation_count"] != 1 for item in roles.values()):
        raise RuntimeError(f"TRAP01 component persistence verification failed: {roles}")
    if report_path and report_path.is_file():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        report.update({
            "derived_ifc_sha256": sha256(output),
            "adjustable_component_group": {"global_id": group.GlobalId, "name": COMPONENT_GROUP_NAME},
            "adjustable_components": roles,
            "adjustment_parameters_mm": {
                "horizontal": {"installed": 113.402666, "minimum": 0.0, "maximum": 252.0, "official_default_overflow": 159.52933},
                "vertical_h": {"installed": 191.939979, "minimum": 85.0, "maximum": 334.0, "official_default_overflow": 208.00002},
            },
            "blue_solid_semantics": BLUE_SOLID_SEMANTICS,
            "blue_dashed_semantics": BLUE_DASHED_SEMANTICS,
            "reference_extension_displayed_in_scene_by_default": False,
            "component_model_pass": True,
        })
        report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def approved_main():
    shared.main()
    parser = shared.argparse.ArgumentParser(add_help=False)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args, _ = parser.parse_known_args(sys.argv[1:])
    persist_adjustable_components(args.output.resolve(), args.report.resolve() if args.report else None)


if __name__ == "__main__":
    approved_main()
