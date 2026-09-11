"""Approval is independent of technical package validation and library inclusion."""
CATEGORIES = {"approved": "已验收", "partial": "部分验收", "pending": "未验收", "skipped": "已跳过"}


def classify(single_status, scene_status, skipped=False):
    if skipped:
        return "skipped"
    if single_status == "approved" and scene_status == "approved":
        return "approved"
    if single_status in ("approved", "partially_approved"):
        return "partial"
    return "pending"


def validate_entry(p, resolve, digest):
    assert p["review_category"] == classify(p["single_product_approval_status"], p["scene_approval_status"], p.get("skipped", False))
    assert p["category_label"] == CATEGORIES[p["review_category"]]
    assert p["package_validation"] == "pass"
    assert p["insertion_content"] in ("approved_3d_2d", "body_only", "review_body_views")
    if p["insertion_content"] == "review_body_views":
        assert p["representation_contract"] == "one_3d_body_three_directional_body_views"
    if p["insertion_content"] == "approved_3d_2d":
        assert p["single_product_approval_status"] == "approved"
    assert resolve(p["ifc_path"]).is_file()
    assert digest(resolve(p["ifc_path"])) == p["ifc_sha256"]
    assert p["formal_ifc_write_authorized"] is False
    record = p["approval_record"]
    assert digest(resolve(record["path"])) == record["sha256"]
    for view in ("plan", "front", "side", "iso"):
        assert resolve(p["previews"][view]).is_file()
    if "preview_evidence" in p:
        validate_previews(p, resolve, digest)


def validate_previews(p, resolve, digest):
    """A drawing button must show linework from this exact IFC, not a camera render."""
    evidence = p["preview_evidence"]
    for view in ("plan", "front", "side"):
        record = evidence[view]
        assert record["kind"] == "bonsai_body_linework"
        assert record["source_ifc_sha256"] == p["ifc_sha256"]
        assert record["role"] == view.title()
        assert resolve(record["svg"]).suffix.lower() == ".svg"
        assert digest(resolve(record["svg"])) == record["svg_sha256"]
        assert digest(resolve(p["previews"][view])) == record["png_sha256"]
    assert evidence["iso"]["kind"] == "camera_render_3d"
    assert digest(resolve(p["previews"]["iso"])) == evidence["iso"]["png_sha256"]
