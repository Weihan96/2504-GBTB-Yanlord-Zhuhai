# IFC Assets v0.8.3 · IFC 直接导入

v0.8.3：加载时保留固定蓝色位置框，绿色填充自身不描边；完成或取消后由原拖放生命周期清理。用户已通过的拖放及点击查看不变。本轮 43 项测试通过；`loading-bounds-blue-green.png` 是在现有验收窗口绘制的 65% 显示样式测试，不是真实导入进度。`loading-bounds-validation.json` 核验更新前后 IFC 内存/磁盘、对象、选择及相机状态相同；没有重启或保存用户场景。

v0.8.2：修复详情弹窗在鼠标按下时抢走拖动事件的问题。点击松手才打开详情；按住拖动直接走原生 IFC 放置，取消或释放回侧栏不会插入。原生鼠标事件队列验证覆盖点击、Esc 取消、侧栏取消和实际 IFC 插入，保留四个 Body 表达且文件不自动保存；见 `drag-gesture-validation.json`。重启验收窗口时连接端口若被其他任务占用，会选择空闲端口并记录 PID/port，不操作其他任务。

v0.8.1：弹窗只保留一排 Plan / Front / Side。3D 渲染图原本已被原生资产卡片复用，因此移除弹窗中的重复入口，并将 runtime 内 39 张 iso PNG 压至最长边 384 像素、128 色。体积从 30,368,323 降至 1,752,393 字节；原始大图保留在历史审核资料，来源路径与哈希见 `thumbnail-compression-validation.json`。117 张线稿 PNG/SVG、IFC、批准状态不变。`compress_asset_thumbnails.mjs` 总是从保留的原图生成缩略图，避免重复有损缩放。

2026-09-10 UI 更新：N 侧栏改名 **Assets**。点击原生资产卡片并松手后打开详情弹窗，显示产品名、验收状态、来源、Plan / Front / Side；顶部 ⓘ 为当前项入口。主面板下方不再保留重复详情和预览。拖动仍直接放置；蓝色定位包围盒仅在拖动时显示，释放后的加载反馈改为无描边的绿色半透明体积。未改 IFC 保存边界、模型几何或批准记录。

2026-09-10：已解除单品资产 `.blend` 依赖。原生缩略图网格现在使用仅在内存中的空卡片，模型只从 IFC 导入。39 个磁盘资产缓存经核验后移到可恢复归档；摆放数据、全部 IFC、预览和材质库保留。新验证见 `ifc-card-release-validation.json`，最终运行文件哈希见 `portable-runtime-manifest.json`。原 v0.6 缓存加载报告和迁移报告保留为历史证据。

保留包围盒预览、绿色加载阶段反馈、原生保存／重载和撤销。卡片没有几何、不进入 IFC；非本库本地资产即使在原生网格出现，也不能触发本库 IFC 插入。本轮核验了真实网格截图、拖放处理器的合成事件、39+1 实例保存重载和两个样本的六张 SVG；不把处理器测试冒充物理鼠标验收。当前测试窗口是本任务 PID 33367 的空白 IFC，未替换 Libelle 的窗口。

2026-09-09：当前运行入口已迁移到 `runtime/catalog.json`。将整个 `runtime/` 目录一起复制即可使用；运行其中的 `open_library.py`，无需旧审核目录。材质库保留 18 个必要材质并内嵌 63 张纹理。没有增加保存按钮：插入留在内存，Ctrl+S 由 Bonsai 保存 IFC。

[打开自包含运行包](</Users/jiaxinchen/.codex/worktrees/3a80/2504 GBTB Yanlord Zhuhai/output/review/approved-product-library/runtime/README.md>)

全部 39 项在操作系统禁止访问原仓库及本机 Poliigon 目录的隔离环境完成加载、插入、保存重载和真实 Bonsai Create Drawing，共 117 张 SVG；已检查全部栅格预览。这是运行验证，不改变 8 项已验收、11 项部分验收、20 项未验收的设计状态。

[隔离出图验证](</Users/jiaxinchen/.codex/worktrees/3a80/2504 GBTB Yanlord Zhuhai/output/review/approved-product-library/portable-isolation-validation.json>)

337 个旧 IFC、Blend、备份和旧流程脚本共 12,748,971,342 字节已移到仓库外的可恢复归档，未永久销毁。421 个来源／批准文件和正式 IFC 通过哈希保护核验。恢复时依照归档中的 `recovery-manifest.json` 逐文件还原，不覆盖新文件。

[清理与恢复清单](</Users/jiaxinchen/.codex/worktrees/3a80/2504 GBTB Yanlord Zhuhai/output/review/approved-product-library/portable-cleanup-result.json>)

本批次开始前的 1,732 项暂存仍是保护基线；新运行包、代码修改和清理删除均未暂存。不要直接提交旧索引：它仍含旧大文件。本轮没有提交、改写历史或操作 stash；仓库外归档和既有 Git 对象仍占磁盘空间。

`all-review-catalog.json`、旧 manifest、旧验证 JSON 和原审核 SVG/PNG 保留为历史证据，其中路径可能指向已归档文件，不再用于运行。`portable-build.json` 记录迁移输入，材质后续变化见 `portable-materials.json`；最终文件哈希以 `portable-runtime-manifest.json` 为准。`build_portable_review_library.py` 是一次性迁移脚本，原输入已归档，不应作为日常重建入口。

## v0.7.1 三视图预览修复

全部 39 个单品的 Plan / Front / Side 改为当前 IFC 经真实 Bonsai Create Drawing 生成的线稿预览，3D 仍用原渲染图。117 张 SVG 与出图记录和单品 IFC 哈希逐一核对，PNG 仅作 Blender Image Editor 的显示缓存；对应 SVG 一并放入自包含 runtime。使用默认黑色打印样式，不伪装成官方蓝线对照资料，不改变设计批准状态。

`pipeline/scripts/repair_review_drawing_previews.mjs` 从现有 `portable-validation/<id>/result.json` 验证过的 SVG 更新预览；IFC 哈希变化时拒绝使用旧图，必须先重新出图。原先的相机/审核图片保留在历史源目录。`drawing-preview-validation.json` 记录逐图来源。插件加载新预览时核验方向及 IFC/SVG/PNG 哈希，防止再次将照片当线稿。只刷新预览数据即可更新已有窗口，不重载 IFC 或替换 RNA。

## 以下为 v0.6 历史验收说明

2026-09-09：用户已验收本轮 Body 表达接入和绿色加载反馈，并授权将当前工作树全部暂存。验收记录为 `library-v06-approval.json`。这不改变各单品原有的设计验收类别，也不授权写入正式 IFC 或删除旧文件。下文的自动化验证记录保留执行当时状态，后续人眼验收以该批准记录为准。

## 表达与验收范围

全部 39 项采用四个 Body：三维 MODEL_VIEW、Plan PLAN_VIEW、Front / Side ELEVATION_VIEW。新单品包移除了旧平立面 Body 和产品 Annotation 表达，保留三维 Body。验收状态仍是 8 项已验收、11 项部分验收、20 项未验收；Libelle 不属于本任务。

15 个原三视图包的已有线稿几何保持不变；其余 24 项接入当前候选线稿，用于继续审核，不代表获批。逐视图来源记录在 body-normalization-validation.json。当前运行包使用 runtime/ifc/<id>.ifc；旧 all-review-catalog.json 只作迁移来源记录。

每包保留原 GlobalId、单位、项目坐标、3D 几何和必要依赖；止回阀还保留两个必要开孔。每包一个非开孔产品，不含 Drawing 相机或 IfcAnnotation 实体。官方 DWG、下载链接、哈希、批准记录及场景配方保留；旧 IFC 的可恢复位置见上面的清理清单。

## Body 与方向

曲线采用产品局部三维坐标，使用维数 3 的 Model / Body 上下文，RepresentationType 保留曲线类型。课程 109000 展示 Plan / Body / PLAN_VIEW / Annotation2D，支持“产品二维表达属于 Body”。Annotation2D 这个 RepresentationType 不等于 Annotation 上下文或 IfcAnnotation 实体；本次没有把三维坐标错误标成二维坐标。

Bonsai 0.8.4 原生只按 TargetView 选上下文，不能区分同为 ELEVATION_VIEW 的 Front 与 Side。本插件在真实 Bonsai Create Drawing 的表达选择处增加适配：

- 每套视图的产品局部方向以 IfcShapeAspect 持久记录；插入时合并上下文、保存及重载后仍可辨认。
- Plan 选 Plan；正交立面将相机方向换算到产品局部，选 Front 或 Side。旋转实例不依赖全局 X/Y 或 STEP 编号。
- 一张标准视图只选一套 Body，不进入独立 Annotation 通道；BISECT 网格剖切路径也排除该产品，防止又叠回高模细节。普通构件保留原生逻辑。
- 反方向使用同一平面的反投影，不冒充独立获批的背面设计。偏离标准方向超过 5°、剖面或其他视图回落到三维 Body。
- 未修改安装目录的 Bonsai 或桥接插件。其他会话出图也须加载本库插件；原生 Bonsai 单独运行仍不能区分两套立面。

测试 SVG 使用 Bonsai 默认黑色打印样式；源 IFC 已有线稿样式不改写。蓝色官方来源的对照 SVG 仍是独立审核资料。缩略图不是本轮逐项重新生成的场景 SVG，不据此宣称场景获批。

## 在 Blender 中验收

先用 Bonsai 打开目标 IFC、选择默认楼层，再加载仓库脚本；不自动插入、不保存偏好或 IFC。

[加载单品库](</Users/jiaxinchen/.codex/worktrees/3a80/2504 GBTB Yanlord Zhuhai/pipeline/scripts/activate_approved_library.py>)

1. 3D 视图按 N，进入“单品库”。网格显示全部 39 项及分类，可搜索和查看三视图 / 3D 资料。
2. 从本面板网格拖进三维视图，蓝色实际尺寸包围盒跟随落点。释放后进入真实 IFC 插入事务；右键、Esc 或视口外释放取消。
3. 松手后包围盒固定，绿色半透明体积按真实载入阶段自下向上填充，事务完成后到达全高。它是阶段进度，不是耗时预测；较长几何计算期间保持上一阶段。成功或异常均清理绘制回调。
4. 每次拖入新建 GlobalId，只对应一个 Blender 模型对象；来源 GlobalId / 哈希存于 ReviewLibrarySource。允许重复放置，不覆盖现有构件。
5. 拖入只改变内存；Ctrl+S 由 Bonsai 原生 bim.save_project 保存。没有额外写入按钮，没有自动保存。

普通 Asset Browser 默认 Collection 拖放只处理 Blender 数据，请使用“单品库”网格。runtime/native-assets/ 只保留摆放数据，权威模型为 runtime/ifc/。旧只读阵列已归档；内部阵列工具仍可直接读取 IFC，不再持久保存 Blend 副本。

环境为 Blender 4.5.3 LTS、Bonsai / IfcOpenShell 0.8.4。升级旧版本须先保存自己的工作，再重新打开，不在资产面板显示期间热替换 RNA。安装包 highpoly_review_library.zip，也可直接运行仓库脚本。

## 验证与边界

- body-normalization-validation.json：39 项四个 Body、零 Annotation、3D 与坐标保持、保存重载。
- body-view-pilot-validation.json：水龙头 provider 保存重载后，真实三视图及栅格检查。
- body-view-*-validation.json：代表样本与 BISECT 的真实 SVG 测试。技术验证不等于用户批准。
- native-integration-validation.json：39 项加一次重复插入、唯一构件/对象、单位、放置、四个 Body、原生保存及重载。
- library-v06-validation.json：本轮汇总及保护状态。旧 v0.4/v0.5 记录不冒充新版鼠标验收。

正式 IFC 未写入，SHA-256 为 7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c。测试只保存可清理隔离副本；正式场景接入应按原 GlobalId 匹配替换审核表达，不叠加第二个模型。

实施期间以 1,639 项暂存为保护基线，新代码和结果保持未暂存，部分与暂存文件重叠；当时索引 SHA-256 为 2d4f1c5ae18f3415abe7ce9c38b118bee3b7efd00c76890b821d5b96a2687f4e。用户完成验收后授权全部暂存，保护基线与本轮修改一并纳入索引。没有提交、stash 操作或旧文件删除。

cleanup-audit.json 是历史候选盘点；此次用户授权的清理实际结果以 portable-cleanup-result.json 为准。
