# TAB02 / RODA Bernardo 367

单品已通过、场景待验收。来源：基于原始高模几何生成的简化图纸表达。

`TAB02-product.ifc` 是唯一持久化产品几何来源：原始 Body、原项目放置矩阵、三个已批准的语义边界 Representation。Plan 1 条路径；Front 和 Side 各 4 条路径。两条立柱边线保持开放，接触边界由桌面底边和底座顶边表达。没有取得官方 DWG；官网和目录仅用于产品身份及尺寸核对。

`scene-recipe.json` 只保存相机、样式、分组和放置等场景元数据，不保存产品 Body 或二维曲线。运行时按原 GlobalId 将纯产品的二维表示接到正式项目临时副本，保留项目的 714 个构件及其 Body 和位置。IfcOpenShell 导入的短线合并问题通过从纯 IFC 逐点重建临时 LINEWORK 网格处理，随后仍运行真实 Bonsai Create Drawing。

Front / Side 使用世界 Z 向上的紧凑画幅；背景是浅灰门墙和地面，批准边界保持黑色。SVG 的样式步骤不移动、删除或重画任何线段。三张 SINGLE 图仅从已验证的真实场景 SVG 中裁取产品 LINEWORK，供图库展示。

`validation.json` 记录纯 IFC schema 检查、正式基线与临时场景继承错误对比、逐点批准坐标和 SVG 逐段投影核对、两次接入幂等性、714 个原构件 Body/放置一致、公共 Provider 保存重载以及最终 Fresh IFC 会话。`cleanup-proposal.json` 只是旧副本和临时目录清单，没有执行删除。

本包使用 Bonsai Course Operator 的 embedded-course-index 证据和 bonsai-launcher；公共 bonsai-mcp Provider 未修改。场景 SVG 的技术通过不等于用户场景验收。
