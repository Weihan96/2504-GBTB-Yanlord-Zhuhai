#!/usr/bin/env python3
"""Run a small native Bonsai elevation batch in an isolated Blender process.

Usage (Blender arguments before ``--`` are omitted here)::

    blender --background --python this_file.py -- work.ifc 05 06

Blender 4.5 can corrupt its Freestyle render depsgraph after repeated Drawing
renders and a project reload.  This runner deliberately loads one checkpoint
IFC, creates a small batch through Bonsai, writes an atomic successor file,
verifies it, and lets the process exit without reloading the project.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import bpy
import ifcopenshell
from bonsai import tool


PROJECT_ROOT = Path(__file__).resolve().parents[2]
GENERATOR = PROJECT_ROOT / "pipeline/scripts/int1_bonsai_elevation.py"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    arguments = sys.argv[sys.argv.index("--") + 1 :]
    preview_only = "--preview-only" in arguments
    if preview_only:
        arguments.remove("--preview-only")
    formal_ifc = None
    if "--formal-ifc" in arguments:
        index = arguments.index("--formal-ifc")
        try:
            formal_ifc = Path(arguments[index + 1]).resolve()
        except IndexError as exc:
            raise SystemExit("--formal-ifc requires a path") from exc
        del arguments[index : index + 2]
    if len(arguments) < 2:
        raise SystemExit(
            "expected: [--preview-only] checkpoint.ifc VIEW_ID [VIEW_ID ...]"
        )
    checkpoint = Path(arguments[0]).resolve()
    view_ids = [view_id.zfill(2) for view_id in arguments[1:]]
    if not checkpoint.is_file():
        raise FileNotFoundError(checkpoint)
    if formal_ifc is not None:
        if not formal_ifc.is_file():
            raise FileNotFoundError(formal_ifc)
        os.environ["INT1_BONSAI_FORMAL_IFC"] = str(formal_ifc)

    os.chdir(PROJECT_ROOT)
    before_sha256 = sha256(checkpoint)
    result = bpy.ops.bim.load_project(
        filepath=str(checkpoint),
        should_start_fresh_session=True,
        use_detailed_tooltip=True,
    )
    if result != {"FINISHED"} or not tool.Ifc.get():
        raise RuntimeError(f"failed to load checkpoint: {result}")

    spec = importlib.util.spec_from_file_location(
        f"int1_bonsai_elevation_headless_{'_'.join(view_ids)}", GENERATOR
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    reports = module.run(view_ids=view_ids)

    if preview_only:
        print(
            "NATIVE_BONSAI_PREVIEW="
            + json.dumps(
                {
                    "checkpoint": str(checkpoint),
                    "checkpoint_sha256": before_sha256,
                    "formal_ifc": str(formal_ifc or checkpoint),
                    "formal_ifc_sha256": sha256(formal_ifc or checkpoint),
                    "drawing_source_is_derived": bool(
                        formal_ifc and formal_ifc != checkpoint
                    ),
                    "view_ids": view_ids,
                    "drawing_names": sorted(
                        report["drawing"]["name"] for report in reports
                    ),
                    "output_svgs": [report["svg"] for report in reports],
                    "checkpoint_write_allowed": False,
                    "pass": True,
                },
                ensure_ascii=False,
            )
        )
        return

    # Checkpoints live under build/int1 but are promoted to the project root.
    # Keep SVG document references relative to the final formal IFC location.
    for report in reports:
        drawing = tool.Ifc.get().by_guid(report["drawing"]["global_id"])
        relations = [
            inverse
            for inverse in tool.Ifc.get().get_inverse(drawing)
            if inverse.is_a("IfcRelAssociatesDocument")
        ]
        if len(relations) != 1:
            raise RuntimeError(f"{drawing.Name}: expected one document relation")
        relations[0].RelatingDocument.Location = report["svg"]

    # Freestyle performs some dependency-graph cleanup after the operator
    # returns.  Give it a quiet boundary before serialising the IFC.
    time.sleep(2)
    temporary = checkpoint.with_name(checkpoint.name + ".next")
    tool.Ifc.get().write(str(temporary))
    verified = ifcopenshell.open(temporary)
    drawings = [
        drawing
        for drawing in verified.by_type("IfcAnnotation")
        if getattr(drawing, "ObjectType", None) == "DRAWING"
    ]
    generated_names = {report["drawing"]["name"] for report in reports}
    missing = generated_names - {drawing.Name for drawing in drawings}
    if missing:
        raise RuntimeError(f"checkpoint verification missed drawings: {sorted(missing)}")
    for name in generated_names:
        drawing = next(item for item in drawings if item.Name == name)
        relation = next(
            inverse
            for inverse in verified.get_inverse(drawing)
            if inverse.is_a("IfcRelAssociatesDocument")
        )
        if not relation.RelatingDocument.Location.startswith("drawings/"):
            raise RuntimeError(f"{name}: document location is not promotion-safe")
    os.replace(temporary, checkpoint)

    print(
        "NATIVE_BONSAI_BATCH="
        + json.dumps(
            {
                "checkpoint": str(checkpoint),
                "before_sha256": before_sha256,
                "after_sha256": sha256(checkpoint),
                "view_ids": view_ids,
                "drawing_names": sorted(generated_names),
                "native_drawing_count": len(drawings),
                "native_elevation_count": sum(
                    (drawing.Name or "").startswith("EL-") for drawing in drawings
                ),
                "pass": True,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
