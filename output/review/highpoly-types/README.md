# 高模图纸代理视觉验收

这个目录只用于检查高模简化后的图纸表达，不是施工放样文件。

每次生成只调用一个 IFC 类型的一件代表实例，不生成整套 IFC 场景几何。

- 灰线：该单件高模的原始三角网格投影。
- 黑线：几何投影得到的简化外轮廓与产品类型语义线候选。
- 蓝线：优先来自已校验的官方原生 CAD/DWG，带白色遮罩并始终位于黑色代理线之上；无官方 CAD 时只能标记为“基于原始高模几何生成的简化图纸表达”。

## Geberit Sigma01 / 115.770.11.5

- [可点击审核目录](geberit-115-770/index.html)
- [汇总审核图](geberit-115-770/review-contact-sheet.png)
- [平面 PLAN](geberit-115-770/plan.svg)
- [正立面 FRONT](geberit-115-770/front.svg)
- [侧立面 SIDE](geberit-115-770/side.svg)
- [生成证据 manifest](geberit-115-770/manifest.json)
- [官方 CAD 与哈希](geberit-115-770/official-native-dwg-linework.json)
- [实际 Bonsai 相机渲染证据](geberit-115-770/bonsai-review-manifest.json)
- [含墙体和卫浴设备的完整项目平面](geberit-115-770/project-context-sanitary-plan.svg)
- [R17 正立面](geberit-115-770/project-context-front-elevation.svg)
- [R17 侧立面](geberit-115-770/project-context-side-elevation.svg)

项目设备登记和 Geberit 官方目录把该类型精确锁定为白色亮光双冲 `115.770.11.5`。蓝线来自官网直接列出的 `G/A/L.dwg`；`P.dwg` 只作为 3D 型号身份文件，未冒充二维视图。官方 DWG 为约 `245.12 × 164.23 × 12.01 mm`，旧 IFC Body 为 `254 × 170.18 × 12.45 mm`，两者均保持原比例，最大 `8.88 mm` 差异显式保留。正式 IFC 中错误的 `WCSEAT` 预定义类型只登记为语义问题，未修改；正确产品语义为双冲水面板。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Geberit 146.140 / AquaClean Sela

- [平面 PLAN](geberit-146-140/plan.svg)
- [正立面 FRONT](geberit-146-140/front.svg)
- [侧立面 SIDE](geberit-146-140/side.svg)
- [生成证据 manifest](geberit-146-140/manifest.json)
- [官方 CAD 与哈希](geberit-146-140/official-native-dwg-linework.json)
- [实际 Bonsai 相机渲染证据](geberit-146-140/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](geberit-146-140/project-context-sanitary-plan.svg)
- [含项目上下文的正立面](geberit-146-140/project-context-front-elevation.svg)
- [含项目上下文的侧立面](geberit-146-140/project-context-side-elevation.svg)

两只马桶共用一套类型代理，不按实例重复维护。蓝线来自 Geberit 官方旧型号 `146.140.11.1_G/A/L.dwg`，不是替代型号；三方向与实际 IFC Body 的最大尺寸差为 0.015 mm。当前仍待用户视觉审核，未写派生 IFC。

## Geberit CleanLine50 / 154.446.KS.1

- [平面 PLAN](geberit-154-446-ks-1/plan.svg)
- [正立面 FRONT](geberit-154-446-ks-1/front.svg)
- [侧立面 SIDE](geberit-154-446-ks-1/side.svg)
- [生成证据 manifest](geberit-154-446-ks-1/manifest.json)
- [官方 CAD 与哈希](geberit-154-446-ks-1/official-native-dwg-linework.json)
- [实际 Bonsai 相机渲染证据](geberit-154-446-ks-1/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](geberit-154-446-ks-1/project-context-sanitary-plan.svg)
- [含项目上下文的正立面](geberit-154-446-ks-1/project-context-front-elevation.svg)
- [含项目上下文的侧立面](geberit-154-446-ks-1/project-context-side-elevation.svg)

蓝线来自 Geberit 官方归档旧型号 `154.446.KS.1_G/A/L.dwg`，没有使用替代型号 `154.446.KS.2`。平面完整 DWG 范围与实际 IFC Body 的差值不超过 0.151 mm；侧立面中官方可见槽体宽 53.4 mm，而 IFC Body 还包含隐藏安装构件、总宽 74.99 mm，此差异保留为显式审核项。当前仍待用户视觉审核，未写派生 IFC。

## Geberit Duofix Sigma / 224.212.00.2

- [平面 PLAN](geberit-duofix-sigma-224-212/plan.svg)
- [正立面 FRONT](geberit-duofix-sigma-224-212/front.svg)
- [侧立面 SIDE](geberit-duofix-sigma-224-212/side.svg)
- [生成证据 manifest](geberit-duofix-sigma-224-212/manifest.json)
- [官方 CAD 与哈希](geberit-duofix-sigma-224-212/official-native-dwg-linework.json)
- [实际 Bonsai 相机渲染证据](geberit-duofix-sigma-224-212/bonsai-review-manifest.json)
- [含墙体和卫浴设备的完整项目平面](geberit-duofix-sigma-224-212/project-context-sanitary-plan.svg)
- [含项目上下文的正立面](geberit-duofix-sigma-224-212/project-context-front-elevation.svg)
- [含项目上下文的侧立面](geberit-duofix-sigma-224-212/project-context-side-elevation.svg)

精确型号由项目收货标签、设备登记和 Geberit 官方目录共同确认为 `224.212.00.2`，不使用欧洲替代型号。蓝线来自该精确 SKU 的官方 `G/A/L.dwg`，`P.dwg` 只归档作为 3D 身份证据。官方正立面 `500 × 1215 mm` 与 IFC Body 差值小于 0.5 mm；官方完整安装深度比当前 IFC Body 包络多约 14.8 mm，在验收图中保留原比例差异。当前仍待用户视觉审核，未写派生 IFC。

## Falper Sorgente / BS01

- [平面 PLAN](falper-sorgente/plan.svg)
- [正立面 FRONT](falper-sorgente/front.svg)
- [侧立面 SIDE](falper-sorgente/side.svg)
- [生成证据 manifest](falper-sorgente/manifest.json)

蓝线为 Falper Sorgente WFB 官方原生 DWG；项目保存的 WFA/WFB 官方矢量 PDF 仅用于逐视图机械交叉验证。WFB 属于产品族参考，不是项目加工图。

## Poliform Hima / HIMA01

- [平面 PLAN](hima01/plan.svg)
- [正立面 FRONT](hima01/front.svg)
- [侧立面 SIDE](hima01/side.svg)
- [生成证据 manifest](hima01/manifest.json)
- [官方来源访问记录](hima01/official-source/source-access-record.json)
- [实际 Bonsai 相机渲染证据](hima01/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](hima01/project-context-furniture-plan.svg)
- [含项目上下文的正立面](hima01/project-context-front-elevation.svg)
- [含项目上下文的侧立面](hima01/project-context-side-elevation.svg)

官方产品页已证实 Hima 有 2D DWG，但下载需要厂商注册流程，本项目尚未取得文件。因此三视图不显示蓝线，黑线必须标记为“基于原始高模几何生成的简化图纸表达”；官方页和技术表只是品牌、产品族与名义尺寸的身份证据，不是线稿几何来源。当前仍待用户视觉审核，未写派生 IFC。

## Gessi316 / 54038 two-way built-in shower mixer external parts

- [可点击审核目录](gessi316-54038/index.html)
- [汇总审核图](gessi316-54038/review-contact-sheet.png)
- [平面 PLAN](gessi316-54038/plan.svg)
- [正立面 FRONT](gessi316-54038/front.svg)
- [侧立面 SIDE](gessi316-54038/side.svg)
- [生成证据 manifest](gessi316-54038/manifest.json)
- [官方来源访问记录](gessi316-54038/official-source/source-access-record.json)
- [Gessi Bathroom 官方目录第 17 页预览](gessi316-54038/official-source/Gessi_Cataloghi_Bathroom-page-17.png)
- [实际 Bonsai 相机渲染证据](gessi316-54038/bonsai-review-manifest.json)
- [实际 Bonsai 平面渲染](gessi316-54038/bonsai-camera-plan.png)
- [实际 Bonsai 正立面渲染](gessi316-54038/bonsai-camera-front-elevation.png)
- [实际 Bonsai 侧立面渲染](gessi316-54038/bonsai-camera-side-elevation.png)
- [含墙体和卫浴设备的完整 FFL 平面](gessi316-54038/project-context-ffl-plan.svg)
- [R17 正立面](gessi316-54038/project-context-front-elevation.svg)
- [R17 侧立面](gessi316-54038/project-context-side-elevation.svg)

Gessi 官方 Bathroom 目录在印刷页 30-31（PDF 文件页 17）把 `54139_54038` 明确标为两路暗装淋浴混水器，图 26 显示手持花洒、软管、出水口和控制件；官方零件图 `GPF540380G001G000` 与安装说明 `GIS001820` 又精确覆盖 `54038`。这些 PDF 只用于身份和组件关系核验，不作为 CAD 线稿。精确 Area Pro 路由存在，但本次没有取得原生 DWG/DXF，因此蓝线严格为零，Plan/Front/Side 的 `4/3/6` 条黑线均标记为“基于原始高模几何生成的简化图纸表达”。单件实际 IFC Body 包络为 `265.347 × 98.695 × 661.908 mm`；完整 FFL 平面按正式 IFC 世界坐标叠加，R17 的 NY/NX 原生 Bonsai 立面分别对应局部 Side/Front，三张图固定 1:50 比例、无拉伸并保留墙体和周边设备。实际 IFC Body 另由 Bonsai 加四台正交相机执行 Blender Render。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Gessi316 Meccanica / 45089_54294

- [平面 PLAN](gessi316-54294/plan.svg)
- [正立面 FRONT](gessi316-54294/front.svg)
- [侧立面 SIDE](gessi316-54294/side.svg)
- [生成证据 manifest](gessi316-54294/manifest.json)
- [官方来源访问记录](gessi316-54294/official-source/source-access-record.json)
- [实际 Bonsai 相机渲染证据](gessi316-54294/bonsai-review-manifest.json)
- [含墙体和卫浴设备的完整项目平面](gessi316-54294/project-context-sanitary-plan.svg)
- [含项目上下文的正立面](gessi316-54294/project-context-front-elevation.svg)
- [含项目上下文的侧立面](gessi316-54294/project-context-side-elevation.svg)

项目 IFC 的 `Gessi316 54294` 描述与 Gessi 官方 2026 目录一致，精确组合型号为内置件 `45089` 加 Meccanica 外露件 `54294`，即长嘴三孔暗装面盆龙头、无下水。产品级 CAD 只在需登录的 Gessi Area Pro 提供，本项目尚未取得精确原生文件；因此蓝线严格为零，黑线统一标记为“基于原始高模几何生成的简化图纸表达”，官方网页和目录只作为产品身份来源。当前仍待用户视觉审核，未写派生 IFC。

## Gessi316 / 54146 ceiling-mounted headshower

- [可点击审核目录](gessi316-54146/index.html)
- [汇总审核图](gessi316-54146/review-contact-sheet.png)
- [平面 PLAN](gessi316-54146/plan.svg)
- [正立面 FRONT](gessi316-54146/front.svg)
- [侧立面 SIDE](gessi316-54146/side.svg)
- [生成证据 manifest](gessi316-54146/manifest.json)
- [官方来源访问记录](gessi316-54146/official-source/source-access-record.json)
- [官方原生 DWG 线稿提取记录](gessi316-54146/official-native-dwg-linework.json)
- [官方 G000 技术 PDF](gessi316-54146/official-source/GPF5414600000G000_1.pdf)
- [DWG/PDF 机械交叉验证](gessi316-54146/official-source/official-pdf-verification.json)
- [官方目录第 36-37 页视觉核验](gessi316-54146/official-source/catalog-page-19-render.png)
- [实际 Bonsai 相机渲染证据](gessi316-54146/bonsai-review-manifest.json)
- [实际 Bonsai 平面渲染](gessi316-54146/bonsai-camera-plan.png)
- [实际 Bonsai 正立面渲染](gessi316-54146/bonsai-camera-front-elevation.png)
- [实际 Bonsai 侧立面渲染](gessi316-54146/bonsai-camera-side-elevation.png)
- [含墙体和家具的完整 FFL 平面](gessi316-54146/project-context-ffl-plan.svg)
- [R12 正立面](gessi316-54146/project-context-front-elevation.svg)
- [R12 侧立面](gessi316-54146/project-context-side-elevation.svg)

Gessi 公共附件 API 已提供精确 `54146 G000` 原生 2D DWG 与技术 PDF；G001、壁装 `54145` 和其他相邻产品均被门禁排除。原生 DWG 的 Plan/Front/Side 分别含 `16/607/660` 条蓝线，三视图包络为 `300.000 × 300.000`、`300.000 × 280.250`、`300.000 × 280.250 mm`；与项目单件 IFC Body `300.000 × 299.945 × 279.371 mm` 的最大投影差为 `0.878759 mm`，未拉伸几何。平面和 R12 两张立面保留原项目墙体、家具与相邻产品，官方蓝线置顶并带白色遮罩；实际 IFC Body 仍由 Bonsai 保存的四台正交相机执行 Blender Render。当前待用户视觉审核，审批文件保持 pending，未写派生 IFC；候选审核包已独立提交，可从新工作树继续比对。

## Gessi316 / 54145 wall-mounted headshower

- [可点击审核目录](gessi316-54145/index.html)
- [汇总审核图](gessi316-54145/review-contact-sheet.png)
- [平面 PLAN](gessi316-54145/plan.svg)
- [正立面 FRONT](gessi316-54145/front.svg)
- [侧立面 SIDE](gessi316-54145/side.svg)
- [生成证据 manifest](gessi316-54145/manifest.json)
- [官方来源访问记录](gessi316-54145/official-source/source-access-record.json)
- [官方原生 DWG 线稿提取记录](gessi316-54145/official-native-dwg-linework.json)
- [官方 G000 技术 PDF](gessi316-54145/official-source/GPF5414500000G000_1.pdf)
- [DWG/PDF 机械交叉验证](gessi316-54145/official-source/official-pdf-verification.json)
- [官方目录第 36-37 页视觉核验](gessi316-54145/official-source/catalog-page-19-render.png)
- [实际 Bonsai 相机渲染证据](gessi316-54145/bonsai-review-manifest.json)
- [实际 Bonsai 平面渲染](gessi316-54145/bonsai-camera-plan.png)
- [实际 Bonsai 正立面渲染](gessi316-54145/bonsai-camera-front-elevation.png)
- [实际 Bonsai 侧立面渲染](gessi316-54145/bonsai-camera-side-elevation.png)
- [含墙体和家具的完整 FFL 平面](gessi316-54145/project-context-ffl-plan.svg)
- [R17 正立面](gessi316-54145/project-context-front-elevation.svg)
- [R17 侧立面](gessi316-54145/project-context-side-elevation.svg)

Gessi 公共附件 API 已提供精确 `54145 G000` 原生 2D DWG 与技术 PDF；G001 和相邻吊顶型号 `54146` 均被门禁排除。原生 DWG 的 Plan/Front/Side 分别含 `22/615/653` 条蓝线，三视图包络为 `600.000 × 299.993`、`600.000 × 118.950`、`300.000 × 118.921 mm`；与项目单件 IFC Body `599.819 × 299.819 × 119.074 mm` 的最大投影差为 `0.181335 mm`，未拉伸几何。平面和 R17 两张立面保留原项目墙体、家具与相邻产品，官方蓝线置顶并带白色遮罩；实际 IFC Body 仍由 Bonsai 保存的四台正交相机执行 Blender Render。当前待用户视觉审核，审批文件保持 pending，未写派生 IFC；候选审核包已独立提交，可从新工作树继续比对。

## Gessi316 Meccanica / 54093 high counter basin spout

- [可点击审核目录](gessi316-54093/index.html)
- [汇总审核图](gessi316-54093/review-contact-sheet.png)
- [平面 PLAN](gessi316-54093/plan.svg)
- [正立面 FRONT](gessi316-54093/front.svg)
- [侧立面 SIDE](gessi316-54093/side.svg)
- [生成证据 manifest](gessi316-54093/manifest.json)
- [官方来源访问记录](gessi316-54093/official-source/source-access-record.json)
- [Gessi316 2026 官方目录第 18 页预览](gessi316-54093/official-source/MAGAZINE_GESSI_316_2026-page-18.png)
- [Gessi Bathroom 官方目录第 17 页预览](gessi316-54093/official-source/Gessi_Cataloghi_Bathroom-page-17.png)
- [实际 Bonsai 相机渲染证据](gessi316-54093/bonsai-review-manifest.json)
- [实际 Bonsai 平面渲染](gessi316-54093/bonsai-camera-plan.png)
- [实际 Bonsai 正立面渲染](gessi316-54093/bonsai-camera-front-elevation.png)
- [实际 Bonsai 侧立面渲染](gessi316-54093/bonsai-camera-side-elevation.png)
- [含墙体和卫浴设备的完整项目平面](gessi316-54093/project-context-sanitary-plan.svg)

Gessi 两份官方目录把 `54093` 精确锁定为高位台面面盆出水嘴，索引技术图 `GPF5409300000G001` 给出的主体直径、水平伸出和总高为 `49.2125×190.5×273.05 mm`；实际 IFC Body 为 `48.974943×190.281433×273.203123 mm`，三个轴的差值均小于 `0.24 mm`。技术图当前无法下载，原生 CAD 仍需 Area Pro 登录且未取得，所以 Plan/Front/Side 的 `1/2/1` 条黑线全部标记为“基于原始高模几何生成的简化图纸表达”，蓝线严格为零。项目 Sanitary Plan 中存在同一 GlobalId 的真实投影，审核副本保留墙体和卫浴设备并用白色遮罩把黑色代理置于最上层；项目没有任何原生立面 SVG 包含该 GlobalId，因此没有伪造 Front/Side 项目上下文。实际 IFC Body 另由 Bonsai 加四台正交相机执行 Blender Render。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Baxter Viktor / BED02

- [平面 PLAN](bed02/plan.svg)
- [正立面 FRONT](bed02/front.svg)
- [侧立面 SIDE](bed02/side.svg)
- [生成证据 manifest](bed02/manifest.json)
- [官方来源访问记录](bed02/official-source/source-access-record.json)
- [实际 Bonsai 相机渲染证据](bed02/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](bed02/project-context-furniture-plan.svg)

项目 IFC 的 `BED02` 描述为 `VIKTOR 162x234xh106`，Baxter 官方产品页和技术表确认 Viktor 的 160×200 cm 排骨架版本外形为 `1720×2340×1060 mm`，实际 IFC Body 包络为 `1675.711×2404.018×1107.146 mm`。三组尺寸差异原样保留，候选没有通过非等比缩放去拟合官方尺寸。官方 2D、3D、BIM 下载需要 Baxter 登录，本项目未取得精确原生 CAD，也未使用第三方替代，因此蓝线严格为零，黑线标记为“基于原始高模几何生成的简化图纸表达”。项目中只有家具平面图包含该床的有效上下文，未伪造不存在的项目立面。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Casablanca / BED01

- [平面 PLAN](bed01/plan.svg)
- [正立面 FRONT](bed01/front.svg)
- [侧立面 SIDE](bed01/side.svg)
- [生成证据 manifest](bed01/manifest.json)
- [官方来源访问记录](bed01/official-source/source-access-record.json)
- [实际 Bonsai 相机渲染证据](bed01/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](bed01/project-context-furniture-plan.svg)

项目 IFC 的 `BED01` 描述为 `Casablanca 180`。Baxter 当前官方产品页把 180×200 cm 排骨架版本标为 `2100×2520×900 mm`，官方技术表第 4 页则标为 `2200×2520×900 mm`，实际 IFC Body 包络为 `2264.167×2538.802×947.268 mm`。两个厂商来源之间的 100 mm 宽度冲突以及高模差异均原样保留；候选没有非等比拟合。官方 2D、3D、BIM 下载需要登录，本项目未取得精确原生 CAD，也未使用第三方替代，因此蓝线严格为零。项目家具平面以 1:50 统一比例叠加实际 Body 代理，并记录旧项目轮廓约 87×42 mm 的包络差。当前仍待用户视觉审核，未写派生 IFC。

## Hunter Douglas 25 mm 百叶帘方向 / sxb010

- [平面 PLAN](sxb010/plan.svg)
- [正立面 FRONT](sxb010/front.svg)
- [侧立面 SIDE](sxb010/side.svg)
- [生成证据 manifest](sxb010/manifest.json)
- [官方来源与 owner direction 访问记录](sxb010/official-source/source-access-record.json)
- [实际 Bonsai 相机渲染证据](sxb010/bonsai-review-manifest.json)
- [含墙体、家具和立面索引的完整 FFL 平面](sxb010/project-context-ffl-plan.svg)

项目旧代理 `sxb010` 是未绑定 IfcType 的 `IfcBuildingElementProxy`。只读检查 `stash@{1}` 的 owner 窗饰确认文件后，可确认 Hunter Douglas 25 mm 铝合金百叶帘产品族和哑光石墨灰／灰黑颜色方向；精确 SKU、色号、操作系统、操作侧与现场复尺尚未确认，stash 未 pop 或 apply。官方中国产品页和官方技术手册证明该高模的 `1812×1596 mm` 宽高落在标准手动 25 mm 产品范围内，但官方下载筛选只返回两份 PDF，没有 DWG；手册中的 1000 mm 矢量示例也不是项目加工图，因此蓝线严格为零。旧高模局部轴为 X=宽、Y=高、Z=深，候选明确使用 Plan=XZ、Front=XY、Side=ZY；完整项目图只存在可机械对应的平面投影，EL-P02 中的薄截面没有被冒充完整立面。当前仍待用户视觉审核，未写派生 IFC。

## Molteni&C Sistema 7 Wall Unit / SIS04

- [平面 PLAN](sis04/plan.svg)
- [正立面 FRONT](sis04/front.svg)
- [侧立面 SIDE](sis04/side.svg)
- [生成证据 manifest](sis04/manifest.json)
- [官方来源访问记录](sis04/official-source/source-access-record.json)
- [实际 Bonsai 相机渲染证据](sis04/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](sis04/project-context-furniture-plan.svg)
- [R04 厨房正立面](sis04/project-context-r04-front-elevation.svg)
- [R04 厨房侧立面](sis04/project-context-r04-side-elevation.svg)

项目 IFC 类型描述为 `Sistema 7 Wall Unit 4 Doors`。Molteni&C 官方 Kitchen Collection 的印刷页 272–275（PDF 文件页 156–157）确认该产品和 `980/1960 × 722/1083 × 331 mm` 标准规格；实际 IFC Body 局部包络为 `1961.975 × 369 × 721.991 mm`，宽高匹配 1960×722 mm 规格，增加的 38 mm 深度只登记为项目高模差异。公开厂家表面未取得这件 Wall Unit 的原生 CAD；近似名称 `Sistema 7 Doors` 的登录区 DWG 属于另一件全高门产品，已明确列为反例并禁止替代，因此蓝线严格为零。三视图黑线标记为“基于原始高模几何生成的简化图纸表达”，项目平面和 R04 立面保留周边墙体与家具；当前仍待用户视觉审核，未写派生 IFC。

## Geberit CleanLine installation set / 154.154.00.1

- [平面 PLAN](geberit-154-154-00/plan.svg)
- [正立面 FRONT](geberit-154-154-00/front.svg)
- [侧立面 SIDE](geberit-154-154-00/side.svg)
- [生成证据 manifest](geberit-154-154-00/manifest.json)
- [官方来源访问记录](geberit-154-154-00/official-source/source-access-record.json)
- [实际 Bonsai 相机渲染证据](geberit-154-154-00/bonsai-review-manifest.json)
- [含墙体和卫浴设备的完整项目平面](geberit-154-154-00/project-context-sanitary-plan.svg)
- [R12 正立面](geberit-154-154-00/project-context-r12-front-elevation.svg)
- [R12 侧立面](geberit-154-154-00/project-context-r12-side-elevation.svg)

Geberit 官方精确页面把项目类型解析为 `154.154.00.1`，标准尺寸 `L1=358、B=130、H=90 mm` 与实际 IFC Body 的 `361.700×130×90.000 mm` 在明确容差内一致。官网内嵌 `cadDrawings` 为空，标准 `A/G/L/P` 原生 DWG URL 均返回 404；官方 EPS 只归档为身份、形态和尺寸证据，未转作 CAD representation。因此三视图蓝色产品 CAD 线严格为零，黑线标记为“基于原始高模几何生成的简化图纸表达”。项目平面中原有的蓝色给排水线属于项目图层，不是产品 CAD。当前仍待用户视觉审核，未写派生 IFC。

## Geberit CleanLine / project 154.154.00.1.F flange component

- [可点击审核目录](geberit-154-154-00-1-f/index.html)
- [汇总审核图](geberit-154-154-00-1-f/review-contact-sheet.png)
- [平面 PLAN](geberit-154-154-00-1-f/plan.svg)
- [正立面 FRONT](geberit-154-154-00-1-f/front.svg)
- [侧立面 SIDE](geberit-154-154-00-1-f/side.svg)
- [生成证据 manifest](geberit-154-154-00-1-f/manifest.json)
- [官方来源访问记录](geberit-154-154-00-1-f/official-source/source-access-record.json)
- [实际 Bonsai 相机渲染证据](geberit-154-154-00-1-f/bonsai-review-manifest.json)
- [含墙体和卫浴设备的完整项目平面](geberit-154-154-00-1-f/project-context-sanitary-plan.svg)
- [R12 正立面](geberit-154-154-00-1-f/project-context-r12-front-elevation.svg)
- [R12 侧立面](geberit-154-154-00-1-f/project-context-r12-side-elevation.svg)

项目 IFC 把官方 `154.154.00.1` 安装套件拆成 parent、`.C` cover 和 `.F` flange 三个协同构件；`.F` 是项目建模分解标签，不是 Geberit 另行发布的商品号。官方产品页与三份 EPS 已归档，但它们画的是完整父套件，不能冒充独立 flange 的蓝色 CAD。该单件实际 Body 包络为 `88×210×104.999998 mm`，Plan/Front/Side 的 `2/7/6` 条黑线均标记为“基于原始高模几何生成的简化图纸表达”。项目平面投影与 Body 一致；R12 两张原生立面分别裁掉约 47 mm 和 87 mm 的隐藏构件高度，审核副本以固定 1:50 比例恢复完整 Body，并保留墙体、地坪、排水组件和周边图纸上下文。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Baxter Stone / BST03

- [平面 PLAN](bst03/plan.svg)
- [正立面 FRONT](bst03/front.svg)
- [侧立面 SIDE](bst03/side.svg)
- [生成证据 manifest](bst03/manifest.json)
- [官方来源访问记录](bst03/official-source/source-access-record.json)
- [Baxter 官方技术表第 10 页预览](bst03/official-source/Baxter_Stone_technical-sheet-page-10-preview.png)
- [实际 Bonsai 相机渲染证据](bst03/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](bst03/project-context-furniture-plan.svg)

项目 `BST03` 已由产品登记、IFC 描述和 Baxter 官方页面共同锁定为 Stone 系列 45 cm 左抽屉独立床头柜。官方技术表第 10 页包含 45 cm 左/右版本的矢量平、立、侧示意，名义尺寸 `450×450×460 mm`；实际 IFC Body 为约 `462.461×462.480×466.219 mm`，皮革、圆角和把手包络差异原样保留。原生 2D/3D/BIM 下载需要 Baxter 登录且尚未取得，因此蓝色产品 CAD 线严格为零，技术表只作为身份、抽屉方向和名义尺寸证据。项目中只有家具平面包含这件活动家具的可机械对应投影，未伪造项目上下文立面。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Beside / BST02

- [可点击审核目录](bst02/index.html)
- [汇总审核图](bst02/review-contact-sheet.png)
- [平面 PLAN](bst02/plan.svg)
- [正立面 FRONT](bst02/front.svg)
- [侧立面 SIDE](bst02/side.svg)
- [生成证据 manifest](bst02/manifest.json)
- [官方来源访问记录](bst02/official-source/source-access-record.json)
- [Baxter 官方公开尺寸 SVG](bst02/official-source/BESICOBS55.svg)
- [Baxter 官方技术表第 7 页](bst02/official-source/Baxter_Beside_TechnicalSheet-page-7.png)
- [实际 Bonsai 相机渲染证据](bst02/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](bst02/project-context-furniture-plan.svg)
- [R09 原生 Side 立面上下文](bst02/project-context-side-elevation.svg)

项目 `BST02` 的 IFC 类型描述和项目产品登记均指向 Baxter Beside `55×58×35 cm`。Baxter 官方产品页、技术表和公开尺寸 SVG 的 Side/Front/Plan 均已归档，但登录区里的原生 2D/3D/BIM 未取得；公开 PDF/SVG 只作为身份、名义尺寸和视图方向证据，不被冒充 native CAD 蓝线。实际 IFC MODEL_VIEW Body 为 `560×577.889×250 mm`，平面宽深与官方尺寸在 15 mm 容差内，实际高度则少 `100 mm`。候选不拉伸几何，并将高度差明确保留为人眼审核项。项目家具平面和 R09 NY 原生立面分别机械对应 Plan 与 Side；项目没有该实例的 Front 立面，因此未伪造。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Molteni&C 505 UP System / project 505 UP V1.LP.S

- [平面 PLAN](505-up-v1-lp-s/plan.svg)
- [正立面 FRONT](505-up-v1-lp-s/front.svg)
- [侧立面 SIDE](505-up-v1-lp-s/side.svg)
- [生成证据 manifest](505-up-v1-lp-s/manifest.json)
- [官方来源访问记录](505-up-v1-lp-s/official-source/source-access-record.json)
- [官方 2021 技术库 CAD SVG](505-up-v1-lp-s/official-source/molteni-505-up-technical-library-full-preview.svg)
- [官方 Inspiring Solutions CAD SVG](505-up-v1-lp-s/official-source/molteni-505-up-inspiring-solution-full-preview.svg)
- [实际 Bonsai 相机渲染证据](505-up-v1-lp-s/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](505-up-v1-lp-s/project-context-furniture-plan.svg)

Molteni&C 官方 505 UP 页面及其公开资源集合已归档两份原生 DWG，文件哈希和转换后的完整官方 CAD SVG 均有记录。但官方目录中所有宽 2912 mm 的现成组合，都不能同时匹配项目高模 `2912×436.5×2380 mm` 的左侧竖向格栅、大型中央开口、右侧展示龛和连续底柜。因此官方 DWG 只作为产品族与目录组合证据，不被重排或冒充项目组合蓝线；三视图黑线明确标记为“基于原始高模几何生成的简化图纸表达”。项目只有家具平面包含该 GlobalId，完整项目平面按 1:50 等比叠加，未伪造项目上下文立面。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Miami Soft E09 / right terminal module

- [平面 PLAN](miamisoft-e09/plan.svg)
- [正立面 FRONT](miamisoft-e09/front.svg)
- [侧立面 SIDE](miamisoft-e09/side.svg)
- [生成证据 manifest](miamisoft-e09/manifest.json)
- [官方来源访问记录](miamisoft-e09/official-source/source-access-record.json)
- [Baxter 官方技术表第 6 页预览](miamisoft-e09/official-source/Baxter_MiamiSoft_technical-sheet-page-6-preview.png)
- [Baxter 官方 E09 测量 SVG](miamisoft-e09/official-source/MIAMSOESE09D.svg)
- [实际 Bonsai 相机渲染证据](miamisoft-e09/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](miamisoft-e09/project-context-furniture-plan.svg)
- [R20 +Y 完整项目侧立面](miamisoft-e09/project-context-r20-side-elevation.svg)

项目类型、描述与 Baxter 官方页面、矢量技术表和 E09 测量 SVG 共同锁定为 `E09 - dx/r` 右端模块，名义尺寸 `1300×1080×700/800 mm`。实际 IFC Body 包络为约 `1343.173×1144.754×797.480 mm`，软包外轮廓差异在明确的 70 mm 容差内且未进行非等比拉伸。原生 2D/3D/BIM 下载需登录且尚未取得，因此蓝色 CAD 线严格为零，官方 PDF/SVG 只作为身份、左右手和尺寸证据；黑线统一标记为“基于原始高模几何生成的简化图纸表达”。项目家具平面与 R20 +Y 侧立面均保留墙体和周边家具，并用白色遮罩把代理置于最上层；R20 的 ±X 投影存在遮挡或裁切，未伪造成完整正立面。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Miami Soft E07 / left dormeuse

- [平面 PLAN](miamisoft-e07/plan.svg)
- [正立面 FRONT](miamisoft-e07/front.svg)
- [侧立面 SIDE](miamisoft-e07/side.svg)
- [生成证据 manifest](miamisoft-e07/manifest.json)
- [官方来源访问记录](miamisoft-e07/official-source/source-access-record.json)
- [Baxter 官方技术表第 7 页预览](miamisoft-e07/official-source/Baxter_MiamiSoft_technical-sheet-page-7-preview.png)
- [Baxter 官方 E07 测量 SVG](miamisoft-e07/official-source/MIAMSOESE07S.svg)
- [实际 Bonsai 相机渲染证据](miamisoft-e07/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](miamisoft-e07/project-context-furniture-plan.svg)
- [R20 +Y 完整项目侧立面](miamisoft-e07/project-context-r20-side-elevation.svg)

项目类型描述和 Baxter 官方第 7 页技术表共同锁定为 `E07 - sx/l` 左 dormeuse，名义尺寸 `1500×1700×700/800 mm`；不得与 `E06 - dx/r` 右手版本混用。实际 IFC Body 包络约 `1543.831×1710.486×797.485 mm`，最大差值约 43.831 mm，在明确的 50 mm 软包容差内且未进行非等比拟合。原生 2D/3D/BIM 下载需登录且尚未取得，因此蓝色 CAD 线严格为零，官方 PDF/SVG 只作为身份、左右手和尺寸证据。项目家具平面和 R20 +Y 侧立面保留墙体与周边家具，白色遮罩后的黑色代理位于最上层；R20 的 ±X 产品投影不完整，未伪造成正立面上下文。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Miami Soft F03 / pouf

- [平面 PLAN](miamisoft-f03/plan.svg)
- [正立面 FRONT](miamisoft-f03/front.svg)
- [侧立面 SIDE](miamisoft-f03/side.svg)
- [生成证据 manifest](miamisoft-f03/manifest.json)
- [官方来源访问记录](miamisoft-f03/official-source/source-access-record.json)
- [Baxter 官方技术表第 7 页预览](miamisoft-f03/official-source/Baxter_MiamiSoft_technical-sheet-page-7-preview.png)
- [Baxter 官方 F03 测量 SVG](miamisoft-f03/official-source/MIAMSOPSF03C.svg)
- [实际 Bonsai 相机渲染证据](miamisoft-f03/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](miamisoft-f03/project-context-furniture-plan.svg)
- [R20 +Y 完整项目侧立面](miamisoft-f03/project-context-r20-side-elevation.svg)

项目 IFC 描述、Baxter 官方第 7 页技术表和独立测量 SVG 共同锁定为 `F03` pouf，名义尺寸 `1700×1080×400 mm`。实际 IFC Body 包络约 `1747.539×1103.160×464.939 mm`；其中高度比名义值大约 64.939 mm，作为软垫鼓起后的真实包络明确保留，并未压缩到 400 mm。原生 2D/3D/BIM 下载需登录且尚未取得，因此蓝色 CAD 线严格为零，官方 PDF/SVG 只作为身份和尺寸证据。项目家具平面和 R20 +Y 侧立面保留墙体与周边家具；+X 投影高度被遮挡，未伪造成完整正立面上下文。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Miami Soft H01 / cushion

- [平面 PLAN](miamisoft-h01/plan.svg)
- [正立面 FRONT](miamisoft-h01/front.svg)
- [侧立面 SIDE](miamisoft-h01/side.svg)
- [生成证据 manifest](miamisoft-h01/manifest.json)
- [官方来源访问记录](miamisoft-h01/official-source/source-access-record.json)
- [Baxter 官方技术表第 8 页预览](miamisoft-h01/official-source/Baxter_MiamiSoft_technical-sheet-page-8-preview.png)
- [Baxter 官方 H01 测量 SVG](miamisoft-h01/official-source/MIAMSOCUSH01.svg)
- [实际 Bonsai 相机渲染证据](miamisoft-h01/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](miamisoft-h01/project-context-furniture-plan.svg)
- [R20 +Y 完整项目侧立面](miamisoft-h01/project-context-r20-side-elevation.svg)

项目 IFC 类型、描述、Baxter 官方第 8 页技术表和独立测量 SVG 共同锁定为 `H01` 80×80 cm 柔性靠垫。项目有两个共享该类型和映射 Body 的实例；候选和 Bonsai 仅隔离一个代表实例。代表 Body 局部包络约 `809.757×564.976×403.964 mm`，靠垫在沙发上立放并受压，官方 800×800 mm 是柔性面尺寸而非刚性三维盒；仅 809.757 mm 宽度用于机械尺寸核验，候选没有把高模拉平或缩放成目录方形。原生 2D/3D/BIM 下载需登录且尚未取得，因此蓝色 CAD 线严格为零，官方 PDF/SVG 只作为身份和尺寸证据。上下文平面发现旧 `PLAN_VIEW` 比实际 Body 投影高约 150 mm，审核副本只隐藏代表实例旧投影并以统一比例 Body 代理替换；第二个 H01、墙体和其他家具均保留。R20 +Y 侧立面也保留完整家具环境。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Miami Soft I03 / roll cushion

- [平面 PLAN](miamisoft-i03/plan.svg)
- [正立面 FRONT](miamisoft-i03/front.svg)
- [侧立面 SIDE](miamisoft-i03/side.svg)
- [生成证据 manifest](miamisoft-i03/manifest.json)
- [官方来源访问记录](miamisoft-i03/official-source/source-access-record.json)
- [Baxter 官方技术表第 8 页预览](miamisoft-i03/official-source/Baxter_MiamiSoft_technical-sheet-page-8-preview.png)
- [Baxter 官方 I03 测量 SVG](miamisoft-i03/official-source/MIAMSORUSI03.svg)
- [实际 Bonsai 相机渲染证据](miamisoft-i03/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](miamisoft-i03/project-context-furniture-plan.svg)
- [R20 +Y 完整项目侧立面](miamisoft-i03/project-context-r20-side-elevation.svg)

项目 IFC 类型名是 `MiamiSoft I03`，实际 Body 包络约 `1711.750×274.262×272.809 mm`，与 Baxter 官方 I03 `1700×Ø270 mm` 在 15 mm 容差内一致。项目 IFC 的 Description 却写成 `Roll 108 x Ø27 cm`；108 cm 是官方相邻 I01 的长度，已作为明确数据冲突记录，未修改正式 IFC，也没有据此把 I03 缩短。原生 2D/3D/BIM 下载需登录且尚未取得，因此蓝色 CAD 线严格为零，官方第 8 页矢量技术表和 I03 测量 SVG 只用于锁定型号、长度和直径。候选黑线统一标记为“基于原始高模几何生成的简化图纸表达”。项目家具平面和 R20 +Y 侧立面均保留墙体和周边家具；+X 投影被裁切，未伪造成完整正立面上下文。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Marilyn / project Marilyn 01 bergere

- [平面 PLAN](marilyn-01/plan.svg)
- [正立面 FRONT](marilyn-01/front.svg)
- [侧立面 SIDE](marilyn-01/side.svg)
- [生成证据 manifest](marilyn-01/manifest.json)
- [官方来源访问记录](marilyn-01/official-source/source-access-record.json)
- [官方原生 DWG 线稿登记](marilyn-01/official-native-dwg-linework.json)
- [Baxter 官方技术表第 13 页预览](marilyn-01/official-source/Baxter_Marilyn_technical-sheet-page-13-preview.png)
- [实际 Bonsai 相机渲染证据](marilyn-01/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](marilyn-01/project-context-furniture-plan.svg)
- [R22 完整项目正立面](marilyn-01/project-context-r22-front-elevation.svg)
- [R22 完整项目侧立面](marilyn-01/project-context-r22-side-elevation.svg)

项目 `Marilyn 01` 已由 Baxter 官方产品页、当前技术表第 13 页、精确测量 SVG、原生 2D/3D ZIP，以及包内 `Marilyn_bergere_86x100xh94.3ds` 共同锁定为 `860×1000×940 mm` bergere。原生 `Marilyn_Abaco.dwg` 直接解析得到 Plan/Front/Side 的 `24/68/54` 条路径，蓝线全部来自该 DWG；失败的通用 DWG→SVG 转换仅保留为反例。实际 IFC Body 约 `877.175×1013.912×970.477 mm`，在 35 mm 软包包络容差内一致。IFC Description 中的 `H78` 是过时的标准扶手椅高度，未用于型号选择、缩放或修改正式 IFC。项目家具平面和两张 R22 立面保留墙体与周边家具，并以白色遮罩把官方蓝线置于最上层；审核裁切只隐藏遮挡产品的定位圆圈或残差诊断框。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Marilyn / project Marilyn 02 pouf

- [平面 PLAN](marilyn-02/plan.svg)
- [正立面 FRONT](marilyn-02/front.svg)
- [侧立面 SIDE](marilyn-02/side.svg)
- [生成证据 manifest](marilyn-02/manifest.json)
- [官方来源访问记录](marilyn-02/official-source/source-access-record.json)
- [官方原生 DWG 线稿登记](marilyn-02/official-native-dwg-linework.json)
- [实际 Bonsai 相机渲染证据](marilyn-02/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](marilyn-02/project-context-furniture-plan.svg)
- [R22 完整项目正立面](marilyn-02/project-context-r22-front-elevation.svg)
- [R22 完整项目侧立面](marilyn-02/project-context-r22-side-elevation.svg)

项目 `Marilyn 02` 由 Baxter 官方产品页、精确测量 SVG、原生 2D/3D ZIP 和包内 `Marilyn_pouf_80x62xh45.3ds` 共同锁定为 `800×620×450 mm` 旋转底座 pouf。原生 `Marilyn_Abaco.dwg` 中精确 pouf 区域直接解析出 Plan/Front/Side 的 `10/52/56` 条路径；蓝线没有从 Marilyn 01 bergere 借用，也没有按项目 Body 非等比拟合。实际 IFC Body 约 `810.615×591.949×456.582 mm`，与官方 DWG 各视图最大包络差约 `28.677 mm`，在明确的 35 mm 软包容差内。项目家具平面与两张 R22 立面保留墙体和周边家具，官方开放路径不会被错误闭合成长对角线；蓝线和白色遮罩位于最上层。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Geberit TRAP01 / 151.116.11.1

- [平面 PLAN](trap01/plan.svg)
- [正立面 FRONT](trap01/front.svg)
- [侧立面 SIDE](trap01/side.svg)
- [生成证据 manifest](trap01/manifest.json)
- [官方来源访问记录](trap01/official-source/source-access-record.json)
- [官方原生 DWG 线稿登记](trap01/official-native-dwg-linework.json)
- [实际 Bonsai 相机渲染证据](trap01/bonsai-review-manifest.json)
- [含墙体、台盆及周边设备的完整项目给排水平面](trap01/project-context-sanitary-plan.svg)

项目 `TRAP01` 由 IFC Body 的固定宽度 `76.540974 mm` 机械锁定为 Geberit `151.116.11.1 / d32`：与官方原生 CAD 的 `76.541707 mm` 仅差 `0.000733 mm`，并排除了 `151.117.11.1 / d40`。官方 `G/A/L/P` DWG 已完整归档，Plan/Front/Side 分别直接解析出 `92/77/117` 条蓝色路径；但官方 CAD 是可调式产品的默认伸展状态，而项目实例是缩短后的安装配置，因此蓝线只在独立右侧面板作为型号和可调包络参考，绝不缩放、拉伸或写入项目 representation。候选黑线统一标记为“基于项目内原始高模几何生成的配置化简化图纸表达”，Plan/Front/Side 为 `16/1/8` 条路径。项目给排水平面保留墙体、台盆和周边设备；原投影受台盆遮挡，完整候选轮廓通过 853 条精确边向量投票定位后置于最上层。项目没有包含该 GlobalId 的原生立面图，因此没有伪造立面上下文。当前仍待用户视觉审核，未写派生 IFC。

## Baxter Colette / project CHA01 armchair

- [平面 PLAN](cha01/plan.svg)
- [正立面 FRONT](cha01/front.svg)
- [侧立面 SIDE](cha01/side.svg)
- [生成证据 manifest](cha01/manifest.json)
- [官方来源访问记录](cha01/official-source/source-access-record.json)
- [Baxter 官方当前技术表第 37 页预览](cha01/official-source/Baxter_Colette_current-technical-sheet-page-37-preview.png)
- [Baxter 官方稳定技术表第 6 页预览](cha01/official-source/Baxter_Colette_technical-sheet-page-6-preview.png)
- [Baxter 官方 COLEPOCO57 测量 SVG](cha01/official-source/COLEPOCO57.svg)
- [实际 Bonsai 相机渲染证据](cha01/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](cha01/project-context-furniture-plan.svg)
- [R07 完整项目正立面](cha01/project-context-r07-front-elevation.svg)
- [R07 完整项目侧立面](cha01/project-context-r07-side-elevation.svg)

项目类型描述、Colette 材质名、Baxter 官方产品页、两份矢量技术表和精确测量 SVG 共同锁定为 `570×600×730 mm` Colette 扶手椅。实际 IFC Body 约 `568.466×593.483×734.465 mm`，最大差值约 `6.517 mm`，无需缩放或拟合。Baxter 原生 2D/3D/BIM 下载需要登录且尚未取得，因此蓝色 CAD 线严格为零，官方 PDF/SVG 只作为身份和尺寸证据；黑线统一标记为“基于原始高模几何生成的简化图纸表达”。项目有两个 CHA01 实例，家具平面和两张 R07 立面均保留墙体与周边家具，以 1:50 等比黑线和白色遮罩替换旧投影；旧平面投影宽度不足约 `131.169 mm`，正立面高度裁切约 `27.465 mm`，侧立面与实际 Body 一致。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Falper floor-mounted basin spout / project FAU02

- [平面 PLAN](fau02/plan.svg)
- [正立面 FRONT](fau02/front.svg)
- [侧立面 SIDE](fau02/side.svg)
- [生成证据 manifest](fau02/manifest.json)
- [官方来源访问记录](fau02/official-source/source-access-record.json)
- [与 GH2 的机械匹配分析](fau02/official-source/cad-match-analysis.json)
- [Falper GH2 官方原生 DWG 解析参考 SVG（不用于 representation）](fau02/official-source/Falper-Cilindro-GH2-native-dwg-reference.svg)
- [实际 Bonsai 相机渲染证据](fau02/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](fau02/project-context-furniture-plan.svg)
- [含墙体和家具的完整项目正立面](fau02/project-context-front-elevation.svg)
- [含墙体和家具的完整项目侧立面](fau02/project-context-side-elevation.svg)

项目类型描述为 `Falper SORGENTE Faucet`，实际 Body 的短水平段加斜出水嘴最接近 Falper Cilindro `GH2`。已从 Falper 官方 ZIP 归档 GH2 的原生 2D/技术 DWG、3D DWG 和技术 PDF，并确认 2D 与技术 DWG 字节完全一致；但项目高模的高度、杆径、底座直径分别约为 `1093.286 / 19.330 / 98.492 mm`，官方 GH2 为 `1128 / 22 / 115 mm`，三个独立尺寸均未通过精确型号门。因此官方 GH2 CAD 只作为最近产品族候选和反证，不作为蓝线或 IFC representation；Plan/Front/Side 的 `31/1/3` 条黑线统一标记为“基于原始高模几何生成的简化图纸表达”。完整项目平面与两张立面保留实际 IFC 投影及台盆、墙体、家具的真实遮挡，未强行把完整轮廓置顶。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Poliform Senzafine / project WD03 wardrobe

- [平面 PLAN](wd03/plan.svg)
- [正立面 FRONT](wd03/front.svg)
- [侧立面 SIDE](wd03/side.svg)
- [生成证据 manifest](wd03/manifest.json)
- [官方来源访问记录](wd03/official-source/source-access-record.json)
- [官方产品页观察与访问边界](wd03/official-source/official-product-page-evidence.json)
- [实际 Bonsai 相机渲染证据](wd03/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](wd03/project-context-furniture-plan.svg)
- [R14 侧向项目立面](wd03/project-context-side-elevation.svg)

项目 IFC 类型名为 `WD03`，Description 为 `Poliform SENZAFINE`，ElementType 为 `WARDROBE`；Poliform 官方产品页证明 Senzafine 是可按建筑空间定制的模块化衣柜系统。项目只有一个 WD03 实例，实际 IFC Body 包络约 `1312.964×585.639×2389.500 mm`。厂家公开页面的资源下载区要求个人资料登记和 CAPTCHA，且公开文本没有给出精确 WD03 或该项目配置的原生 CAD URL；未提交表单、未绕过验证、未使用第三方 CAD。因此 Plan/Front/Side 的 `3/2/6` 条黑线统一标记为“基于原始高模几何生成的简化图纸表达”，蓝线严格为零，也没有把目录组合缩放成项目配置。项目家具平面包含实际 WD03 投影；R14 原生立面只保留同一 GlobalId 的诊断包围框，其尺寸对应侧向投影，审核副本已移除残差诊断颜色并按原生 1:50 比例叠加黑线。没有包含 WD03 投影的项目正立面，因此未伪造正立面上下文。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Poliform Pivot + Senzafine / project WD01 custom wardrobe

- [可点击审核目录](wd01/index.html)
- [汇总审核图](wd01/review-contact-sheet.png)
- [平面 PLAN](wd01/plan.svg)
- [正立面 FRONT](wd01/front.svg)
- [侧立面 SIDE](wd01/side.svg)
- [生成证据 manifest](wd01/manifest.json)
- [官方来源访问记录](wd01/official-source/source-access-record.json)
- [官方页面与技术文档观察记录](wd01/official-source/official-product-page-evidence.json)
- [实际 Bonsai 四相机渲染证据](wd01/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](wd01/project-context-furniture-plan.svg)
- [R09 +X 完整项目正立面](wd01/project-context-front-elevation.svg)
- [R09 +Y 完整项目侧立面](wd01/project-context-side-elevation.svg)

项目 IFC 类型 `WD01` 的 Description 为 `Poliform Pivot Senzafine`，ElementType 为 `WARDROBE`；项目产品登记明确把它定义为 Senzafine 衣柜与墙装 Pivot 旋转门组合。Poliform 官方 Pivot、Senzafine 产品页和官方 Pivot 技术文档共同证明两个系统可集成，并发布“between cabinet and wall”等组合，但项目 Body 包络约 `1202.879×673.524×2390.023 mm`，属于项目定制配置。官网原生资源需要提交个人资料和 CAPTCHA，公开文本没有精确 WD01 原生 CAD URL；直接获取官方 PDF 二进制也被 Cloudflare 403 阻止，未绕过验证、未使用第三方 CAD。因此 Plan/Front/Side 的 `7/3/7` 条黑线统一标记为“基于原始高模几何生成的简化图纸表达”，蓝线严格为零。Furniture Plan 与 R09 +X/+Y 原生立面均包含同一 GlobalId 的实际投影组；三张审核副本按原生 1:50 等比叠加并保留墙体、柜体与周边家具。+Y 原投影侧宽较完整 Body 少约 14 mm，审核叠加保持完整 Body 并居中，没有拉伸。实际 Bonsai 文件只导入这一件 Body，保存 Plan/Front/Side/Iso 四个正交相机和 Blender Render 结果。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Poliform Senzafine / project WD02 glass wardrobe

- [平面 PLAN](wd02/plan.svg)
- [正立面 FRONT](wd02/front.svg)
- [侧立面 SIDE](wd02/side.svg)
- [生成证据 manifest](wd02/manifest.json)
- [官方来源访问记录](wd02/official-source/source-access-record.json)
- [官方产品页观察与访问边界](wd02/official-source/official-product-page-evidence.json)
- [实际 Bonsai 相机渲染证据](wd02/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](wd02/project-context-furniture-plan.svg)
- [R22 +X 完整项目正立面](wd02/project-context-front-elevation.svg)
- [R22 -Y 完整项目侧立面](wd02/project-context-side-elevation.svg)

项目 IFC 类型名为 `WD02`，Description 为 `Poliform Glass Wardrobe`，ElementType 为 `WARDROBE`；项目家具产品登记与来源登记将其明确映射到 Poliform Senzafine 模块化衣柜系统。项目只有一个 WD02 实例，实际 IFC Body 局部包络约 `600.000×628.406×2392.753 mm`。厂家公开资源仍需要个人资料登记和 CAPTCHA，且没有公开精确 WD02 玻璃衣柜配置的原生 CAD URL；未提交表单、未绕过验证、未使用第三方 CAD。因此 Plan/Front/Side 的 `4/3/3` 条黑线统一标记为“基于原始高模几何生成的简化图纸表达”，蓝线严格为零。Furniture Plan、R22 +X 与 R22 -Y 原生 SVG 均包含同一 GlobalId 的实际投影组；三张审核副本按原生 1:50 刚体等比叠加，保留墙体和周边家具，并仅移除审核用残差诊断颜色。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## antoniolupi Street / project STREET-H sink-holder subcomponent

- [平面 PLAN](street-h/plan.svg)
- [正立面 FRONT](street-h/front.svg)
- [侧立面 SIDE](street-h/side.svg)
- [生成证据 manifest](street-h/manifest.json)
- [官方来源访问记录](street-h/official-source/source-access-record.json)
- [官方原生 DXF 机械解析记录](street-h/official-native-dxf-linework.json)
- [antoniolupi 官方 AL_Street.dxf](street-h/official-source/AL_Street.dxf)
- [官方 Street 技术图 PDF](street-h/official-source/ANTONIOLUPI-official-Street-technical.pdf)
- [实际 Bonsai 相机渲染证据](street-h/bonsai-review-manifest.json)
- [含墙体和卫浴设备的完整项目平面](street-h/project-context-sanitary-plan.svg)

项目登记把两个 `STREET-H` / `Sink holder` 子部件关联到 BATHM 的 antoniolupi Street 国内定制台盆，审核只隔离代表实例 `2ajpw0I9n1dBypfISg3ejX`。官方 `AL_Street.dxf` 和技术图中的精确父产品族簇为 `STREET240 prof. 40 + STREET4054 prof. 40`：DXF 平面、立面分别是 `108×40`、`108×25` 图形单位，官方技术图第 1 页把同一完整台盆顶明确为最小 `1080×400×250 mm`。项目 STREET-H Body 只有 `300×150×100 mm`，且是重复的内部子部件，不是该完整目录台盆顶。因此官方 DXF 已归档并通过 PDF 尺寸交叉验证，但被机械门拒绝作为 STREET-H representation；Plan/Front/Side 的 `1/1/5` 条黑线统一标记为“基于原始高模几何生成的简化图纸表达”，蓝线严格为零。项目原图只有 Sanitary Plan 包含代表实例的实际投影组，没有包含该 GlobalId 的原生项目立面，未伪造立面上下文。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## antoniolupi Street / project STREET domestic-custom washbasin top

- [可点击审核目录](street/index.html)
- [汇总审核图](street/review-contact-sheet.png)
- [平面 PLAN](street/plan.svg)
- [正立面 FRONT](street/front.svg)
- [侧立面 SIDE](street/side.svg)
- [生成证据 manifest](street/manifest.json)
- [官方来源访问记录](street/official-source/source-access-record.json)
- [官方原生 DXF 配置审计](street/official-native-dxf-configuration-audit.json)
- [antoniolupi 官方 AL_Street.dxf](street/official-source/AL_Street.dxf)
- [实际 Bonsai 四相机渲染证据](street/bonsai-review-manifest.json)
- [含墙体和卫浴设备的完整项目平面](street/project-context-sanitary-plan.svg)

项目登记把 `STREET` 标记为 antoniolupi Street 国内定制一体台盆，隔离的实际 IFC Body 包络为 `1000×470×250 mm`。官方 `AL_Street.dxf` 已逐配置机械扫描：IFC Description 对应的 `street240 prof. 40 + street4054 prof. 40` 标准簇为 `1080×400×250 mm`，同深度最近的 `street147 prof. 47 + street4754 prof. 40` 标准簇为 `1080×470×250 mm`；所有 470 mm 深簇均没有 1000 mm 宽版本。因此官方家族 DXF 仅作为反证和来源关联归档，不缩放、不转写为蓝线，也不冒充国内定制项目图。Plan/Front/Side 的 `3/1/2` 条黑线统一标记为“基于原始高模几何生成的简化图纸表达”。Sanitary Plan 审核图保留墙体、洁具和周边项目元素，按原生 1:50 等比叠加完整 Body 平面；项目原生立面 SVG 没有该 GlobalId，所以未伪造项目立面上下文。实际 Bonsai 文件只导入这一件 IFC Body，并保存 Plan/Front/Side/Iso 四个正交相机和 Blender Render 结果。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## RODA Bernardo 367 / project TAB02 side table

- [平面 PLAN](tab02/plan.svg)
- [正立面 FRONT](tab02/front.svg)
- [侧立面 SIDE](tab02/side.svg)
- [生成证据 manifest](tab02/manifest.json)
- [官方来源访问记录](tab02/official-source/source-access-record.json)
- [RODA 官方 2024 目录 Bernardo 367 证据页](tab02/official-source/RODA-official-Bernardo-367-catalogue-extract.pdf)
- [实际 Bonsai 相机渲染证据](tab02/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](tab02/project-context-furniture-plan.svg)

项目 IFC 类型描述为 `Bernardo 367`，材质样式包含 `Bernardo_serizzo` 与 `Bernardo Metal milk`；RODA 官方 2024 目录明确列出 BERNARDO 367 方形边桌、Rodolfo Dordoni 设计、尺寸 `500×500×670 mm` 和 `Ø300×40 mm` Sempione 白石底座。隔离的 IFC Body 包络同样为 `500×500×670 mm`，三个轴零差值，无需缩放或拟合。RODA 当前 Bernardo 产品页仍把 2D/3D 下载导向登录区，且公开页面已只列 353/354；未创建账号、未取得原生 CAD、未使用第三方文件。因此 Plan/Front/Side 均为 1 条黑色“基于原始高模几何生成的简化图纸表达”，蓝线严格为零。项目 Furniture Plan 中存在同一 GlobalId 的真实投影，按原生 1:50 比例精确叠加并保留墙体和其他家具；没有任何项目原生立面 SVG 包含 TAB02，因此未伪造立面上下文。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Baxter Ninfea / BST01（开口侧未解析）

- [可点击审核目录](bst01/index.html)
- [汇总审核图](bst01/review-contact-sheet.png)
- [平面 PLAN](bst01/plan.svg)
- [正立面 FRONT](bst01/front.svg)
- [侧立面 SIDE](bst01/side.svg)
- [生成证据 manifest](bst01/manifest.json)
- [官方来源访问记录](bst01/official-source/source-access-record.json)
- [Baxter 官方右开版本公开矢量证据](bst01/official-source/NINFCLE40D-right-opening.svg)
- [Baxter 官方左开版本公开矢量证据](bst01/official-source/NINFCLE40S-left-opening.svg)
- [实际 Bonsai 相机渲染证据](bst01/bonsai-review-manifest.json)
- [含墙体和家具的完整项目平面](bst01/project-context-furniture-plan.svg)
- [R09 原生 Front 立面上下文](bst01/project-context-front-elevation.svg)
- [R09 原生 Side 立面上下文](bst01/project-context-side-elevation.svg)

项目 IFC 的 `BST01` 描述为 `Nifea Comodino diam42xh45`；Baxter 官方名称是 Ninfea，公开产品页与技术表确认床头柜名义尺寸为 `420×420×450 mm`，并分别发布右开和左开版本。项目类型、Body 与实际 Bonsai 渲染均不能可靠确定开口侧，因此两份官方 R/L 矢量只作为产品族、尺寸和变体存在性的证据，未选择、拼接或转写为蓝色 CAD representation。原生 2D/3D/BIM 下载需要 Baxter 登录且尚未取得，所以三视图蓝线严格为零，黑线标记为“基于原始高模几何生成的简化图纸表达”。实际 IFC Body 包络约 `399.031×399.032×450 mm`，与官方名义平面包络各差约 21 mm，几何保持原尺寸。完整家具平面和 R09 两张原生立面均保留墙体与周边家具，白色遮罩后的黑色代理位于最上层。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## Electric flue check valve / unresolved project proxy

- [可点击审核目录](electric-flue-check-valve/index.html)
- [汇总审核图](electric-flue-check-valve/review-contact-sheet.png)
- [平面 PLAN](electric-flue-check-valve/plan.svg)
- [正立面 FRONT](electric-flue-check-valve/front.svg)
- [侧立面 SIDE](electric-flue-check-valve/side.svg)
- [生成证据 manifest](electric-flue-check-valve/manifest.json)
- [来源与检索边界记录](electric-flue-check-valve/official-source/source-access-record.json)
- [实际 Bonsai 四相机渲染证据](electric-flue-check-valve/bonsai-review-manifest.json)
- [R04 +X 完整项目立面](electric-flue-check-valve/project-context-r04-px-elevation.svg)
- [R04 −Y 完整项目立面](electric-flue-check-valve/project-context-r04-ny-elevation.svg)

正式 IFC 只有一个未类型化 `IfcBuildingElementProxy`，名称为 `Electric flue check valve`，没有厂商、型号、系统、端口、材料或产品文档关联。隔离 Body 的局部包络约 `210×240×247.5 mm`，形态包含球形阀体、两端套筒、方形安装板与执行器盒；该外形只支持泛化分组，不能证明具体产品。Novy 官方 `906271` 与专利 `CN209977425U` 已登记为被拒绝的近似产品/技术类型参考，均未冒充项目身份或 CAD。Plan/Front/Side 的 `2/5/10` 条黑线统一标记为“基于原始高模几何生成的简化图纸表达”，蓝线严格为零。R04 +X 与 −Y 原生立面含同一 GlobalId 的真实投影，按原生 1:50 比例叠加，包围框差分均为零，并保留墙体、柜体和周边构件。`E303` 的 PX009 只作为定位包络，未冒充平面投影。实际 Bonsai 文件只渲染这一件 Body 并保存 Plan/Front/Side/Iso 四个正交相机。当前仍待用户视觉审核，审批文件保持 pending，未写派生 IFC。

## 放行规则

在用户确认三个方向的 silhouette 和必要语义线之前，`manifest.json` 中的 `approved_for_drawing_ifc` 必须保持 `false`，不得替换已发布立面图的图纸代理。
