# IFC → 施工图批量流水线

本目录定义滨海湾项目的可重复编译流程。目标是把已确认的 IFC 几何与语义，经过标准化信息验证、人审决策和机械 QA，生成可追溯的施工图候选及发布包。

## 产品形态

- **用户工作**：编译、审查、确认、发布施工图，而不是逐张手工补图。
- **操作入口**：`bun run pipeline:*` 命令和 `build/reports` 中的统一报告。
- **人眼验收入口**：Blender + Bonsai；仅用于需要空间定位或视觉判断的对象。
- **候选区**：`build/candidates`，可随时重建，绝不覆盖正式图纸。
- **发布区**：`releases/<version>`；只有全部发布门通过后才可生成。

## 标准数据边界

| 数据 | 责任 | 是否为源数据 |
| --- | --- | --- |
| `2504 GBTB Yanlord Zhuhai.ifc` | 已确认几何和语义 | 是 |
| `pipeline/ids/*.ids` | buildingSMART IDS 机器信息要求 | 是 |
| `pipeline/standards/*.md` | IDS 无法表达的几何、模度、锚点和验收规则 | 是 |
| `pipeline/decisions/*.csv` | IDS 无法表达的选择、推断依据和人审状态 | 是 |
| `pipeline/project.json` | 项目路径、基线和出图登记配置 | 是 |
| `build/**` | 快照、候选、QA 报告 | 否，可重建 |
| `releases/**` | 通过发布门的不可变交付包 | 是 |

不使用自造的 `contract.yaml`。自动推断必须进入 CSV 审查队列，至少记录对象 GUID、依据、置信度、是否需要人审和状态。

## 安全不变量

1. 只有一个 IFC 写入器；快照、IDS 与图纸 QA 均只读。
2. 默认命令不写 IFC，不调用 `save_ifc_file`，不改变 Blender 场景。
3. IFC 修改必须在独立批次完成，并在同一轮保存、reload/rebuild、磁盘验证和视觉验收。
4. `build/` 中的候选不得作为正式图纸或反向覆盖源 IFC。
5. 未通过人审的空间、门窗或施工语义不得自动写回 IFC。

## 命令

```bash
bun run pipeline:snapshot
bun run pipeline:check
bun run pipeline:gate
bun run pipeline:coordinate-audit
bun run pipeline:geometry-alignment-audit
bun run pipeline:geometry-diff-head
bun run pipeline:c003-wall-candidates
bun run pipeline:c003-group-candidate -- --help
bun run pipeline:coordinate-noise-candidate
bun run pipeline:direction-noise-candidate
bun run pipeline:normalization-exit-audit
bun run pipeline:a105-flooring-delete-candidate
bun test
```

- `pipeline:snapshot`：读取 IFC 和现有 SVG，生成源哈希与结构快照。
- `pipeline:check`：生成 QA JSON/Markdown；允许已知阻塞存在，便于持续盘点。
- `pipeline:gate`：执行同一套检查；存在阻塞时返回非零退出码，用于发布门。
- `pipeline:coordinate-audit`：只读统计对象 origin、IFC 长度数值和非整数分布；当前 origin 复核阈值为 `0.1 mm`。
- `pipeline:geometry-alignment-audit`：只读提取墙体世界坐标几何，检查近共面边和墙体接缝；默认判定门槛 `0.1 mm`，候选搜索窗口 `1.0 mm`，且墙体必须至少有 `100 mm` 垂直重叠才会比较平面缝隙，避免把上下楼层投影误判为连接。报告会把超差墙对按连通关系归并为 `C003-G*` 人审组，审核人不需要阅读散乱记录。可直接调用脚本并传入 `--tolerance-mm <毫米>`、`--search-window-mm <毫米>`、`--minimum-vertical-overlap-mm <毫米>`。
- `pipeline:geometry-diff-head`：在不创建 worktree 或临时 IFC 的情况下直接读取 Git `HEAD` 中的 IFC，与当前工作区逐墙比较 ObjectPlacement 和世界坐标几何。它会区分“没有变化”“只改 Placement 但世界几何保持”“世界几何发生变化”。
- `pipeline:c003-wall-candidates`：只在 `build/c003-wall-candidates/` 生成客卫飘窗墙的三种归整副本，比较原点重设、刚性平移和正交归整。每种副本都会恢复宿主 Opening 的原世界 placement，并记录墙体 Hausdorff 差、共面与接缝结果；不会覆盖源 IFC，也不会原位修改共享 Profile。
- `pipeline:c003-group-candidate`：按传入的 IFC、组 ID、对齐报告和 `pipeline/decisions/c003-wall-group-transforms.csv` 生成只读候选；检查已确认墙对是否全部达到 `0.1 mm`、拓扑和 Root GlobalId 是否保持、是否出现非预期对象变化，并让 20 m 飘窗上下裁切体保持世界坐标不变。已写入批次的候选与 QA 报告保留在 `build/c003-wall-groups/`，不得对当前 IFC 重复应用历史变换。
- `pipeline:space-semantics-candidate`：生成只读 Space 语义候选；把 IFC4 `IfcSpace` 从普通构件 containment 改为 Storey aggregation，并按已确认 Grid 语义边界写入 `Qto_SpaceBaseQuantities.GrossFloorArea`。候选门会验证 22 个 Space 世界几何不变、既有 Root GlobalId 全部保留、`NetFloorArea` 仍为空。
- `pipeline:parametric-json-noise`：清理 Bonsai `BBIM_Door/BBIM_Window.Data` JSON 中不超过 `0.01 mm` 的整数尾差；比例等非近整数数值保持不变，并验证全部门窗世界几何、schema、实体数和 Root GlobalId 不变。
- `pipeline:direction-noise-candidate`：只在 `build/candidates/` 生成方向尾差副本。它只归整由单个 `IfcAxis2Placement3D.Axis/RefDirection` 独占、归一化后距正交轴不超过 `0.000001` 的三维 `IfcDirection`；方向公差是无量纲比例，不与 `0.1 mm` 几何公差混用。候选必须证明只有目标 Direction 实体改变、全部带表示的 IfcProduct 世界顶点移动不超过传入的 `0.1 mm`、网格数量、schema、实体数和 Root GlobalId 不变。Direction `#43042` 控制直接切割窗 `1KBBoRNWb1Svpz9__rmYKy` 的 Opening；`#2081575/#2081576` 控制直接切割排烟止回阀 `1faflkXXH6M9cnYPE9Liir` 的 Opening；`#2101029/#2101030` 控制该止回阀自身坐标系。归零会跨过布尔三角化阈值，因此这 5 个 Direction 作为明确列出的拓扑敏感例外保留。
- `pipeline:normalization-exit-audit`：读取当前 IFC、15 步标准、决策 CSV 及当前哈希的派生审计，输出 15 步 `pass/in_progress/pending/blocked` 状态；过期报告会硬失败。9 类受控例外或非整数设计值必须在 `p0-review.csv` 中保持唯一的 `implemented` 记录；5 个 Direction ID 直接从该 CSV 读取，不再依赖命令行硬编码。剩余原点只有在类别专用报告与对应决策的精确对象集合一致时才能从独立队列扣除；当前支持 19 件活动复杂家具、4 个插座接触保护例外、2 个 PVC110 组合对象，以及满足“单一宿主、无 Filling、直接相对宿主 placement”的 114 个从属 Opening。另有 81 个没有唯一安装锚点的对象，必须由 8 条 `c003-origin-responsibility` 记录无重复、无漏项地精确移交 A104、A105、WFIN、PLUM、ELEC、RCP1、INT1、DET1；移交不是受控例外，也不授权当前 IFC 写入。报告只写入 `build/`，不成为新的 PM 真相源。`pipeline:normalization-exit-gate` 使用同一证据，但在仍有未分配原点、未决 C003 决定、例外/移交登记失配或未提交 IFC 批次时返回非零退出码。
- `pipeline:a105-flooring-delete-candidate`：只在 `build/candidates/` 生成 A-105 候选；把 18 个本层 `TerrazzoMosaicTile` Covering 写为 `FLOORING`，并按明确的人审动作删除 4 个 `-1020..-1000 mm` 参考副本。每个参考副本各宿主 1 个无填充线性地漏 Opening；删除候选会逐一验证宿主关系、相对 placement 和无填充状态，并把这 4 个附属 Opening 纳入预期删除边界。`pipeline:a105-flooring-keep-candidate` 保留全部参考对象。两个候选均通过 8 个本层地砖宿主的 `IfcOpeningElement` Boolean 切割体把两樘线性地漏可见留口从约 `14.58 mm` 对称调整为 `14.6 mm`，不移动地砖 ObjectPlacement、地漏或留口中心线；同时验证本层 18 个对象的变化不超过 `0.1 mm`、共享 RepresentationMap 保留、非目标产品无超差变化及 Root GlobalId 变化符合动作。
- `pipeline:construction-surface-audit`：只读盘点 Covering、Slab、Furniture 与 SanitaryTerminal 的世界包围盒、朝向、容器、材料、原点整数残差和最近顶点距离，并生成 Blender 人审队列。水平对象通过世界 XY 投影与 IfcSpace 轮廓求交，记录完整覆盖率及主要房间，不按包围盒中心猜房间。宽泛队列进一步按可观察证据拆成子组，例如全高竖向表面、高位竖向/水平条带、抬高水平表面、FFL/RFL Slab、服务区/房间家具和卫生洁具；子组名称只描述几何位置与既有容器，不自动赋予地面、天花、墙面、屋面或固定安装语义。已有 `PredefinedType`、整数原点或原点恰好落在顶点上，也不能单独授权锚点写入，必须由对象类别专用候选证明施工锚点。
- `pipeline:remaining-origin-audit`：以可传入的 `--tolerance-mm` 只读分类仍超过门槛的产品原点，分别列出无宿主门、已填充/未填充 Opening、施工表面、家具、机电设备及次要建筑构件。原点残差不等于世界几何错误；所有队列均保持 `automatic_write_allowed=false`，必须由对象类别专用候选证明施工锚点后才能写入。
- `pipeline:furniture-anchor-audit`：只读区分家具安装角色。显式 `IfcFurnitureType.PredefinedType=BED/CHAIR/SOFA/TABLE` 直接登记为可移动复杂家具；其余对象读取 `pipeline/decisions/furniture-installation-role.csv`，按 exact type name 或 GlobalId 应用已记录的安装角色、依据、置信度和人审要求。品牌、产品名、型号变体、用途和官方来源独立记录在 `pipeline/decisions/furniture-product-register.csv`，不与“固定/活动”安装角色混成同一备注。当前 89 件家具已分成 19 件活动家具受控例外与 70 件固定安装构件；70/70 固定构件已通过类别专用候选完成整数安装锚点。
- `pipeline:furniture-group-origin-candidate`：读取 `pipeline/decisions/c003-fixed-furniture-anchor-targets.csv` 和明确的接触余量表，把固定家具按安装接触组生成候选。锚点只能使用现有物理顶点、物理边中点或已机械证明位于原几何表面的整数点；复杂网格不逐点取整。门同时检查目标 placement、锚点到面、全部家具内部接触、墙/板/饰面接触、新相交、非目标产品世界几何、schema、实体数和 Root GlobalId。`--tolerance-mm` 与 `--maximum-axis-shift-mm` 分别控制验收精度和每轴最大刚体位移。
- `pipeline:light-fixture-anchor-candidate`：为统一 `IfcLightFixtureType=RA.LP / DIRECTIONSOURCE` 的灯具生成安装中心候选。只有外形为 `55×55×89 mm`、ObjectPlacement 精确位于 XY 中心、Z 已为整数且每轴吸附不超过 `0.5 mm` 的实例才进入目标集。候选保持灯具内部中心关系，检查全部灯具、天花/楼板/固定 Proxy/梁/墙的接触与相交集合、非目标产品世界几何、schema、实体数和 Root GlobalId；不能加入实体碰撞树的 Curve/CAD Proxy 会显式登记，且源/候选排除集合必须一致。
- `pipeline:socket-anchor-candidate`：盘点 11 个 `SOC01/SOC04/SOCF04` IfcElectricAppliance，只将有统一外形与既有安装面中心证据的 7 个安全对象吸附到整数毫米。4 个会破坏 CAB-A、WW20 或 WAL200 安装接触的精确 GlobalId 由 `COORD-SOCKET-C003-CONTROLLED` 保护且不写入；3 个 SOCF04 只允许关闭与岛台间精确记录的安装余量。候选门检查原点、内部/上下文接触、新增/丢失相交、非目标世界几何以及 schema/实体/Root GlobalId。
- `pipeline:electric-appliance-anchor-options`：对 `AC700 / WD01 / OV01 / HD01` 比较“刚体吸附旧原点”与“纯原点重设”。前者会破坏 HD01 与 WW20 的接触；后者证明全部世界几何可保持，但只有原点仍落在自身网格上的对象才能写入，不允许只为数字好看而保留无意义锚点。
- `pipeline:electric-appliance-hd-origin-reset-candidate` 与 `pipeline:stacked-appliance-origin-reset-candidate`：分别把 HD01 油烟机、WD01 饮水机和 OV01 烤箱的 ObjectPlacement 重设到自身现有底面整数点。局部 MappingTarget 施加逆变换，设备不移动；候选必须证明锚点到网格不超过 `0.1 mm`、全模型世界几何不变、实体与 Root 集合保持。
- `pipeline:ac700-anchor-candidate`：为 AC700 / RPIZ-22FSLN5QD/P 生成上安装面和下检修面两个候选。两者均先把原点无损重设到现有面中心，再使设备每轴不超过 `0.5 mm` 吸附到整数锚点。候选检查目标点到移动后几何、预期世界位移、接触/相交、非目标变化及 schema/实体/Root。项目已明确按安装便利设计，因此固定设备优先采用直接对应吊顶或安装支承关系的上安装面；下检修面保留为检修净空语义。
- `pipeline:flush-plate-anchor-candidate`：针对两块类型名为 `Geberit 115.770`、几何为约 `12.445×254×170.180 mm` 的墙装冲水面板，选择朝向管井墙一侧的安装面中心并吸附到整数毫米。候选检查产品身份、面中心到网格、每轴最大 `0.5 mm`、全部产品世界几何、接触/相交、schema、实体数和 Root GlobalId；本命令不修改其当前错误的 `WCSEAT` 语义。
- `pipeline:opening-anchor-audit`：只读区分未填充 Opening 的独立定位责任。只有恰有一个宿主、没有 Filling、且 `ObjectPlacement.PlacementRelTo` 直接等于宿主 `ObjectPlacement` 的 Opening，才能登记为“随宿主验收”；它不得独立吸附整数，也不重复占用独立红点队列。共享 placement、非直接父子或其他关系形态继续作为独立 Boolean/锚点审核对象。
- `pipeline:opening-group-origin-reset-candidate`：针对多个 `IfcOpeningElement` 共用同一 ObjectPlacement 与 Tessellation 的裁切组生成纯原点候选。目标必须是既有裁切盒表面内的整数点；共享表示只做一次逆变换，组内全部 Opening 同时移动原点，世界切割几何保持。当前两个组覆盖 5 个 Opening 和 5 个 IfcCovering 宿主；命令本身只生成候选，这 5 个 Opening 已随批准原子批次正式写入。
- `pipeline:c003-shared-surface-host-origin-reset-candidate`：针对上述 5 个 IfcCovering 宿主生成专用原子候选。脚本先保存 7 个直属 Opening 的 4 个唯一 LocalPlacement 组，再重设宿主，最后每组只恢复一次；门检查 5 个宿主锚点、Opening 世界 placement、共享 placement/representation 身份、全部产品世界几何、墙体共面/接缝、schema、实体数和 Root GlobalId。该批已正式写入。
- `pipeline:covering-origin-reset-candidate`：对整数几何锚点审计已证明的 10 个 IfcCovering 生成纯原点候选，不写地面、墙面或天花语义。Axis/Body 施加逆变换，19 个直属 Opening 保持世界 placement；候选同时检查全模型产品几何、0.1 mm 墙体关系、实体增量和 Root GlobalId。
- `pipeline:c003-safe-origin-batch-candidate`：把两个共享 Opening 组与 10 个 Covering 串成一次候选原子批次，避免正式写入时重复保存/reload。当前直接对正式 IFC 的机械比较覆盖 39 个相关对象：15 个只改变原点且世界几何保持，24 个完全不变。
- `pipeline:c003-surface-edge-origin-reset-candidate`：在上述原子候选后追加 2 个竖向 Covering 与 1 个 Slab；锚点分别位于既有 20 mm 模度斜边和 100 mm 模度顶边。该阶段机械比较覆盖 46 个相关对象：18 个只改变原点且世界几何保持，28 个完全不变。
- `pipeline:c003-construction-surface-origin-reset-candidate`：在前述原子链后追加 45 个 IfcCovering 与 9 个 IfcSlab；目标均为现有轴向三角面内的整数点。五个与共享 Opening 组耦合的 Covering 宿主被显式排除。最终 72 对象批次比较 138 个相关对象：72 个只改变 ObjectPlacement 且世界几何保持，66 个完全不变；该批写入时正式 IFC 与批准候选逐字节一致，随后只清理写入产生的 5 个正交 Direction 尾差；全模型最大世界顶点变化 `0.0 mm`，墙体 183/183 共面与 165/165 接头满足 `0.1 mm`。
- `pipeline:integer-geometry-anchor-audit`：只读扫描原点超过 `0.1 mm` 的产品，查找已经存在于世界网格上的整数顶点、物理轴向/非轴向边整数点；仅对 IfcSlab/IfcCovering 允许轴向三角面内整数点回退，共面三角化对角线不会被误作构件边。报告按 `100/50/30/20/1 mm` 模度、IFC 类别和审核子组汇总，但始终保持 `automatic_write_allowed=0`：几何点存在不等于该点符合对象类别的安装锚点，仍需类别专用规则和世界几何不动候选。
- `pipeline:a103-wall-semantics-candidate`：只在 `build/candidates/` 生成用户确认后的 A-103 墙体语义候选。84 面最终保留墙写入 `Status=EXISTING`，其中 64 面橙色 Aircrete 为 `LoadBearing=false`、20 面蓝灰 Concrete 为 `LoadBearing=true`；4 面新建墙写入 `Status=NEW`。厨房移门门垛 `2Ca2tGerPBj8kvUTjlUtDl` 的完成面厚度从 `153.886884 mm` 归整为 `154 mm`，固定 `Y=-3678 mm` 一侧，只把另一侧移动 `0.113116 mm`。候选必须证明已调整门洞墙 `0hKdvAZkn1TejLgJhK_vDp`、其 Opening `1YxMx6s0r3ZPPohkRKXWbl`、厨房门垛 Opening 和全部非目标产品世界几何不变，并保持 183/183 共面和 165/165 接头满足 `0.1 mm`。
- `pipeline:a103-candidate`：从当前磁盘 IFC 和 Bonsai 刷新的 `Wall Plan.svg` 只读生成 A-103 SVG 候选、88 面墙后台登记表和机械 QA。它动态读取 20 条 GridAxis（强制包含 Grid 08），检查两方向尺寸链闭合，标注 4 面新建墙和 5 个有宿主门洞，并对生成文字做碰撞检查。图例直接读取已经写入 IFC 的墙体状态：绿色为 `NEW`，橙色为 `EXISTING + LoadBearing=false`，蓝灰为 `EXISTING + LoadBearing=true`；语义与材料分组不一致时硬失败。3 樘无宿主门继续移交 A-104。随后保留矢量导出单页 `500×400 mm` PDF，并用 Poppler 检查页数/纸张、回渲 PNG 到已忽略的 `tmp/` 供目检；PDF 稳定输出到 `output/pdf/`。
- `pipeline:a102-demolition-candidate`：只读校验交付基础 DWG、拆除情况 PDF 和正式 IFC 的 SHA-256。脚本从 PDF 图面提取 11 个黄色轮廓“已经拆除”矩形和 2 个红色轮廓“计划拆除”矩形，排除图例；以 14 个 Grid 锚点拟合 PDF/项目坐标，最大残差不得超过 `0.5 mm`，再生成 13 个位于 `50 mm` 控制网格上的待审墙段。每条记录保存原始矩形、原始世界坐标、候选边界、依据、置信度和 `review_required=yes`，但 `formal_ifc_write_allowed=no`。命令在已忽略的 `build/candidates/2504-GBTB-a102-demolition.ifc` 中新增 13 面真实 `IfcWall`，写入确定的 GlobalId、`Status=DEMOLISH`、D01–D13 Tag、来源/置信度/审核状态，以及 `Qto_WallBaseQuantities.Length/Width/Height`；13 面墙的世界包围盒必须在 `0.1 mm` 内等于登记候选，原 IFC 的 88 墙、Opening、门、窗必须全部保持世界几何。Blender 直接加载候选 IFC：D01–D11 墙为红色、D12–D13 墙为紫红色，白色顶面标签显示 `编号 长×厚`；不存在旧实体边框，标签不承担几何。Solid、X-Ray 关闭、`show_in_front=0`。D13 显式保护已调整门洞墙 `0hKdvAZkn1TejLgJhK_vDp` 及 Opening `1YxMx6s0r3ZPPohkRKXWbl`。
- `pipeline:a102-demolition-postwrite`：在 D01–D13 获得“平面位置、方向和长度大差不差”的近似人审结论后，机械核对正式 IFC 与批准候选、Git 写前基线。正式对象使用 `Status=DEMOLISH`，审核属性使用 `ReviewStatus=CONFIRMED_APPROXIMATE`、`FormalIfcWriteAllowed=true`，并明确 `POSITION_DIRECTION_LENGTH_APPROXIMATE_NOT_SURVEY_GRADE`；不得把图纸配准结果解释为现场测量。门限为可传入的 `0.1 mm`：13 面拆除墙包围盒及 QTO 必须符合登记值，Git 基线中的 284 个墙/Opening/门/窗不得发生世界几何变化，已调整门洞墙 `0hKdvAZkn1TejLgJhK_vDp` 与 Opening `1YxMx6s0r3ZPPohkRKXWbl` 必须保持 `0.0 mm`。
- `pipeline:beam-anchor-candidate`：读取 `pipeline/decisions/c003-beam-anchor-targets.csv` 中 6 根梁的明确整数毫米锚点，只在 `build/candidates/` 生成刚体平移副本。候选允许梁及其直属 Opening 随宿主移动，但要求目标 placement 在 `0.1 mm` 内、每根梁移动不超过 `1 mm`、非目标产品几何不变、拓扑/schema/实体数/Root GlobalId 保持。
- `pipeline:beam-relationship-audit`：只读比较待人审梁在当前 IFC 与候选 IFC 中对 `IfcBeam/IfcWall/IfcSlab` 的结构面关系。脚本分别检查相对面接触与同侧面齐平，记录候选是改善、保持还是把原本不超过 `0.1 mm` 的关系变成超差；搜索窗和最小交叉投影均可传入。出现回退时只说明梁不能单独刚体移动，必须联动相邻构件或改用纯原点重设，不授权写入。
- `pipeline:beam-origin-reset-candidate`：为不能安全刚体移动的梁生成“纯原点重设”副本。候选把 ObjectPlacement 放到现有 Body 顶面的整数 Grid 交点，并在梁局部 Axis/Body 中施加逆变换；直属 Opening 恢复原世界 placement。写入门要求目标点到实际三角面不超过 `0.1 mm`，全模型产品世界几何、schema、实体数和 Root GlobalId 均保持。
- `pipeline:slab-origin-reset-candidate`：对整数几何锚点审计确认的 7 块 IfcSlab 生成纯原点候选。目标均为现有物理轴向边上的 `100 mm` Grid 点；ObjectPlacement 重设时对 Axis/Body 施加逆变换并恢复直属 Opening 的世界 placement。候选必须验证锚点确实落在三角面上且全模型产品世界几何不变，不依赖尚未确认的 FFL/RFL 施工语义。
- `pipeline:structural-origin-reset-candidate`：把已经分别通过的 3 根梁和 7 块板合并到同一当前哈希候选，作为下一次正式 IFC 写入的唯一结构原点批次。它避免依次应用两个旧基线候选；同一份报告验证 10 个目标 placement、13 个直属 Opening、锚点到面距离、全模型产品几何、预计新增表示实体、schema 和 Root GlobalId。
- `pipeline:ceiling-proxy-origin-reset-candidate`：为 15 个名称和几何均明确的固定天花/灯槽 IfcBuildingElementProxy 生成纯原点候选。目标来自现有物理轴向边整数点；Tessellation 的独占点表及可选 Box 表示做局部逆变换，不移动世界几何。`AC Diffuser Embeded` 不因位于天花区而纳入此批，它继续按机电安装锚点审核。
- `pipeline:fixed-origin-reset-candidate`：把已分别通过的 3 根梁、7 块板和 15 个固定天花/灯槽对象合并到同一当前哈希原子候选，作为正式 IFC 写入的唯一依据；统一复核 25 个 placement、13 个直属 Opening、锚点到面距离、全模型产品世界几何、预期实体增量、schema 和 Root GlobalId。
- `pipeline:baseboard-origin-reset-candidate`：为剩余 2 条已明确为 `SKIRTINGBOARD` 的踢脚线生成纯原点候选；目标点位于现有底部物理轴向边的 `100 mm` Grid 点，Axis/Body 局部表示做逆变换，世界几何不得移动。
- `pipeline:fixed-surface-service-origin-reset-candidate`：为 39 个已有整数几何锚点的固定对象生成纯原点候选，包括 5 个固定 Proxy、6 个施工表面、4 个空调设备、8 个灯具、4 段排水管和 12 个洁具。Tessellation、SweptSolid、Clipping 和独占 `IfcMappedItem.MappingTarget` 均做局部逆变换；Curve2D Axis 目标保持原世界标高。候选不得把备选几何点解释为连接尺寸或安装定位语义。
- `pipeline:covering-semantics-candidate`：只读生成踢脚线语义候选；只有名称明确以 `Baseboard` 开头、几何高度 `90-110 mm`、厚度不超过 `25 mm` 且长度至少 `100 mm` 的 IfcCovering 才候选写为 `PredefinedType=SKIRTINGBOARD`。候选必须保持全部 Covering 世界几何、Root GlobalId 和实体数量不变。
- `pipeline:space-placement-recovery`：只读修复候选，用已验证 Git IFC 基线按 GlobalId 恢复异常归零的 IfcSpace ObjectPlacement。基线目标必须是整数毫米平移和单位旋转；候选必须与基线 Space 世界几何在传入的 `0.1 mm` 内一致，并证明当前 IFC 中只有 22 个目标 Space 的世界几何发生变化，Space 语义、IFC4 和 Root GlobalId 保持。
- `pipeline:coordinate-noise-candidate`：只在 `build/candidates/` 生成副本，将距离整数不超过 `0.01 mm` 的 `IfcLengthMeasure` 尾数写为精确整数；不会覆盖源 IFC。
- 坐标工具不提供“全部原点取整”或“施工几何取整”写入命令；超过复核阈值的 origin 和墙体尺寸必须先确定几何锚点与联动对象，再经过全模型机械差分、关系检查和受控 IFC 写入。Blender 只审核整体位置、拓扑和语义，不负责判断亚毫米对齐。
- IDS/IfcTester 继续负责 IFC 信息要求，例如实体类型、属性、分类、材料、关系和属性值；墙体共面、接缝间距、ObjectPlacement 与世界几何差分属于几何 QA，由本项目的只读 IfcOpenShell 脚本负责，不写入 IDS。
- `bun test`：运行解析器和机械检查的单元测试。

生成文件：

```text
build/
  snapshots/ifc-snapshot.json
  reports/qa-report.json
  reports/qa-report.md
```

## 当前 P0 顺序

1. Space 合并和楼梯删除已进入 `f216251` 基线；继续关闭剩余 `PROVISIONAL_GRID_CELL`、楼层关系和基础数量。
2. A-103 房间标注、墙体完成面、墙厚和门洞定位候选。
3. A-104 门窗编号、洞口、开启方向、房间邻接及门窗表。
4. A-105 地坪类型、材料、边界、完成面标高及现场复核项。

## 发布门

1. **Source**：工作目录、IFC 路径、schema 和 SHA-256 明确。
2. **IFC integrity**：STEP 编号及 GUID 唯一，关系和数量满足要求。
3. **IDS**：P0 信息要求通过 buildingSMART IDS 验证器。
4. **Decision**：所有需要人审的推断均为 `confirmed` 或 `rejected`。
5. **Drawing**：图号、比例、版本、编号、尺寸闭合和文件引用完整。
6. **Visual**：空白、裁切、文字/标注碰撞和对象遮挡检查通过。
7. **Release**：IFC、SVG、PDF、报告和哈希写入同一 release manifest。
