# SXB010 单品 IFC 包

这是亨特 25 mm 百叶窗帘已批准成果的存储迁移，不是重新设计，也不是官方 CAD。

- `SXB010-product.ifc`：唯一产品 `1O9JRXCI56VRUbpuLJy86Z`，保留全部 2 个 Body 表示、批准的 Plan/Front/Side，以及必要属性、样式和项目坐标依赖。
- IFC 内没有 Annotation、Drawing 相机、场景 Include/Exclude、空间几何或其他产品。
- `scene-recipe.json`：最终三个已批准 session 的相机、二维线的原位置、Drawing 设置和文档依赖；不存 3D Body 或二维线几何。出图线只能来自上述单品 IFC。
- 场景验证在临时正式项目副本中，按 GlobalId 替换/接入同一产品的图纸表达，不能新增第二个 Body。真实 Bonsai Create Drawing 生成 `SXB010-SCENE-*.svg`。
- 场景沿用已批准的蓝色**审核高亮**，不表示官方 CAD 蓝线。单品 SVG 用黑线，来源始终为“基于原始高模几何生成的简化图纸表达”。
- `SXB010-SINGLE-*.svg/.png` 直接从纯 IFC 中的 Approved 三视图表达导出，用于图库，不是场景 SVG。`SXB010-BODY-ISO.png` 复用已批准的实际 Body 保存相机渲染，原始文件和哈希记录在 manifest。
- 保留 45 片百叶；三视图分别 5、51、55 条路径。此迁移不重算简化、不拉伸、不移动项目坐标。
- `validation.json` / `manifest.json` / `handoff.json` 记录持久保存、重载、三视图几何对比及目视结果；只有状态为 complete/pass 才表示验证完成。
- 正式 IFC 必须保持 SHA-256 `7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c`。旧完整副本未删除，`cleanup-proposal.json` 仅是等待用户确认的清单。

运行顺序：外部 Python 执行 `migrate_sxb010.prepare()`；Bonsai launcher 打开纯 IFC，并通过任务专属 9888 桥执行 `migrate_sxb010.scene()`；最后外部 Python 执行 `verify_sxb010.py` 并目视 PNG。运行时场景步骤不读取旧完整审核 IFC。
