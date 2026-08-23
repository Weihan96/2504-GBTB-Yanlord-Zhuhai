# Geberit CleanLine installation set 154.154.00.1

项目 IFC 类型 `Geberit 154.154.00` 的单件代表实例为 `14EazrLgP8whYZWY_yCuKy`。本文件夹只生成该实例，不渲染整个模型。

- `plan.svg`、`front.svg`、`side.svg`：灰线为实际 IFC Body 网格投影，黑线为“基于原始高模几何生成的简化图纸表达”。
- `project-context-*.svg`：把同一候选以 1:50 统一比例放回项目平面和 R12 立面，保留墙体、家具及周边设备。
- `bonsai-camera-*.png`：Bonsai 实际打开隔离 IFC、保持 `Body` representation、创建四台相机并执行 Blender Render 的结果。
- `official-source/`：归档精确官方产品页与三份官方 EPS。EPS 仅用作型号、形态和标准尺寸证据，不是 CAD representation 线稿。
- `official-source/official-source-revalidation.json`：重新请求官方产品页、文章 API、四个标准 DWG 端点和三份 EPS，并记录 HTTP 状态、SHA-256 与本地 EPS 字节比对。

精确官方页面的文章条目仍显示 `cadDrawings` 未提供，标准 `A/G/L/P` 原生 DWG URL 均返回 404；匿名文章 API 当前返回 401，这只记录为访问状态，不用来推断 CAD 是否存在。因此本候选没有蓝色产品 CAD 线。项目平面中已有的蓝色给排水管线仍属于原项目图层，不能解释为产品官方 CAD。

当前状态为待人眼审核。审批记录保持 `pending` 时，写入器即使带 `--apply` 也必须拒绝；正式权威 IFC 永远不是允许的输出目标。
