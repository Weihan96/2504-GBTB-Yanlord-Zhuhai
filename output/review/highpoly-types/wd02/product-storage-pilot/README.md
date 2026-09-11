# WD02 单品 IFC 保存方式试行

本目录只持久保存一份单品 IFC：`WD02-product.ifc`。它包含 WD02 原始 Body、已批准的三视图表达、图纸相机/注释设置、必要空间层级、属性、类型及样式依赖。不是完整项目副本。

保留代表 GlobalId `3yuXF4PtnBHgIXDlmXFJ$7`，单位为毫米。原点仍在项目 `(5897, 896, 48) mm`，原始旋转不变；没有移到世界原点或重新摆正。

来源为“基于原始高模几何生成的简化图纸表达”，没有官方 CAD 蓝线。原批准记录、官方身份资料及其哈希不变。

## 文件职责

- `WD02-product.ifc`：唯一的新持久 IFC；独立打开、保存及重载验证。
- `WD02-SCENE-PLAN.svg`、`WD02-SCENE-FRONT.svg`、`WD02-SCENE-SIDE.svg`：真实 Bonsai Create Drawing 场景出图。PNG 只是 SVG 的可视预览。
- `package-manifest.json`：产品身份、来源、批准记录哈希、坐标、接入规则和输出入口。
- `validation.json`：保存/重载、坐标、实体数量、正式 IFC 前后哈希、独立校验及已批准线稿对照。
- `cleanup-proposal.json`：旧完整审核副本、旧完整场景 Blend 和临时目录的待清理清单。只是清单，未删除。

## 场景工作方式

1. 先核验正式 IFC 基线 SHA-256。
2. 创建唯一临时目录，将正式 IFC 复制到其中；不会把正式 IFC 作为可写工作文件载入。
3. 从单品包读取已批准表达，匹配项目既有 GlobalId；Body、单位和放置矩阵必须相符。只接入二维表达及其依赖，不追加第二个家具实例。
4. 再接入一次，检查没有新实体，证明操作可重复。
5. 真实 Bonsai Create Drawing 出三张 SVG；仅保存并重载临时项目副本，验证模型和线的坐标。
6. 检查 SVG 和预览，逐图对照此前批准的场景黑线。正式 IFC 前后哈希必须相同。
7. 验证完重新打开单品 IFC；完整项目临时副本只留作待确认清理对象，不作为长期审核资产。

本试行不变更正式 IFC，也不扩大任何产品批准范围。旧文件继续保留，不能仅凭本说明删除它们。

## 已知边界

目前是 WD02 试行，未批量迁移。带嵌套子构件、开洞或连接端口的产品需要单独定义依赖范围，转换器会停止而不是静默丢失依赖。单品相机中的项目 Include IDs 是场景配方，不要求其引用的墙体也装入单品 IFC；只有临时接回项目后才生成场景图。

旧 IFC 的 Plan 上下文缺少必填 WCS，新单品包补为显式二维恒等坐标系；旧应用记录缺开发者，新包使用“源文件未注明”的明确占位，不伪造作者。原项目 Body 和已批准线稿均未改变。正式项目中既有的 `OD_Textures/Materials.blend` 缺失单独记为源文件问题，不影响本次二维出图。

## 复验

代码使用 Python 的 IfcOpenShell 完成 IFC 运算，用 Bun 调用现有公共 Bonsai bridge；没有修改或另建 Provider。

```sh
python3 -m unittest discover -s pipeline/tests -p test_review_product_package.py
python3 pipeline/scripts/verify_wd02_product_storage.py
```

首次提取使用 `pipeline/scripts/pilot_wd02_product_storage.py`；若单品 IFC 已存在会拒绝覆盖。场景流程必须在本任务的 Bonsai 会话中运行 `pilot_wd02_product_storage.scene()`，会再次生成临时目录与本试行场景 SVG。不要运行旧的完整副本生成器来覆盖本包，也不要将临时 IFC 路径当作正式文件路径。
