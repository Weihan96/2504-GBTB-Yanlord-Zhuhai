import { expect, test } from "bun:test";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "../..");
const modulePath = resolve(root, "pipeline/scripts/render_svg_pdf.py");

test("SVG PDF renderer parses a single custom-size page", () => {
  const source = `
import importlib.util, json, pathlib, sys
spec = importlib.util.spec_from_file_location("render_svg_pdf", ${JSON.stringify(modulePath)})
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
print(json.dumps(module.parse_pdfinfo("Pages:           1\\nPage size:       1416.96 x 1134 pts (custom)\\n")))
`;
  const result = Bun.spawnSync(["python3", "-c", source], { cwd: root });
  expect(result.exitCode).toBe(0);
  expect(JSON.parse(result.stdout.toString())).toEqual({
    pages: 1,
    width_points: 1416.96,
    height_points: 1134,
  });
});
