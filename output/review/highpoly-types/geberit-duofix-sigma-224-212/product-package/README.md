# Duofix 224.212.00.2 单品包

`DUOFIX-product.ifc` 是主卫已批准实例 `3pvAlH5C14v8uVEJ1LmK8M` 的纯单品文件：保留全部 3 个 Body 表示、ApprovedPlan／ApprovedFront／ApprovedSide 以及原 GlobalId、单位、坐标链和必要属性样式。它不包含场景 Drawing、Annotation、空间几何或其他 WC 部件。

三张蓝线来自最终通过的主卫会话，不使用早期未对齐的代表实例。原 WC_01 聚合关系记录在 validation 中；单品不复制父装配，而接回项目时按 GlobalId 匹配原水箱，保留项目已有聚合成员关系。

`scene-recipe.json` 只记录相机、显示与文档依赖，不存第二份三维或二维线条。运行 `migrate_duofix.scene()` 时，只读取本单品 IFC、该配方和正式项目的临时副本；真实 Bonsai Create Drawing 生成 `DUOFIX-SCENE-*.svg`。`DUOFIX-SINGLE-*.svg` 为相同已生成二维线条的单品重框预览，不是另一套几何算法。

迁移与核验结果见 `validation.json`、`manifest.json` 和 `handoff.json`。旧完整项目副本、官方 DWG 与批准记录均保留；`cleanup-proposal.json` 仅提出待清理项，未经用户确认不得删除。正式 IFC 禁止永久写入，操作前后均检查 SHA-256。

来源：Geberit 精确型号 224.212.00.2 官方 G/A/L 原生 DWG。仅官方型号参考，不宣称是项目施工深化图。官方下载链接及 SHA-256 保留在上一级 `official-source/source-access-record.json`，并索引于 manifest。
