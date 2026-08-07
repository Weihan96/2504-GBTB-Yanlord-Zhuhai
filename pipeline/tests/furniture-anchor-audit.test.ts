import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/furniture_anchor_audit.py");

test("explicit loose furniture semantics become no-write controlled exceptions", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("furniture_anchor_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(json.dumps({key: module.classify_furniture(key) for key in ["BED", "CHAIR", "SOFA", "TABLE", "SHELF", "USERDEFINED", None]}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const records = JSON.parse(result.stdout.toString());
  for (const key of ["BED", "CHAIR", "SOFA", "TABLE"]) {
    expect(records[key].installation_role).toBe("loose_furniture");
    expect(records[key].review_required).toBe(false);
    expect(records[key].automatic_write_allowed).toBe(false);
  }
  for (const key of ["SHELF", "USERDEFINED", "null"]) {
    expect(records[key].installation_role).toBe(
      "fixed_or_loose_requires_review",
    );
    expect(records[key].review_required).toBe(true);
    expect(records[key].automatic_write_allowed).toBe(false);
  }
});

test("project role decisions classify fixed joinery and loose bedside tables", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("furniture_anchor_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
fixed={"selector_kind":"type_name","selector_value":"SB03","installation_role":"fixed_furniture","basis":"cabinet side panel","confidence":1.0,"human_review_required":"no","status":"implemented"}
loose={"selector_kind":"type_name","selector_value":"BST01","installation_role":"loose_furniture","basis":"bedside table","confidence":0.99,"human_review_required":"no","status":"implemented"}
print(json.dumps({"fixed":module.classify_furniture("USERDEFINED",fixed),"loose":module.classify_furniture("USERDEFINED",loose)}))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const records = JSON.parse(result.stdout.toString());
  expect(records.fixed.installation_role).toBe("fixed_furniture");
  expect(records.fixed.normalization_disposition).toBe(
    "fixed_installation_anchor_required",
  );
  expect(records.fixed.review_required).toBe(false);
  expect(records.loose.installation_role).toBe("loose_furniture");
  expect(records.loose.normalization_disposition).toBe(
    "controlled_exception_keep_complex_origin",
  );
});

test("product identity stays separate from installation role", () => {
  const source = `
import importlib.util, json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(${JSON.stringify(modulePath)}).parent))
spec = importlib.util.spec_from_file_location("furniture_anchor_audit", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
products=[{"selector_kind":"type_name","selector_value":"Libelle","manufacturer":"Baxter","product_name":"Libelle","product_variant":"","product_category":"modular_bookcase","intended_use":"wine and cocktail tools","source_url":"https://example.test/libelle","identity_basis":"user confirmed","confidence":1.0,"human_review_required":"no","status":"confirmed"}]
matched=module.find_selector_record("gid", "Libelle", products)
print(json.dumps(module.product_identity(matched)))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  const identity = JSON.parse(result.stdout.toString());
  expect(identity.manufacturer).toBe("Baxter");
  expect(identity.product_name).toBe("Libelle");
  expect(identity.intended_use).toBe("wine and cocktail tools");
  expect(identity).not.toHaveProperty("installation_role");
});
