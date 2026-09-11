# IFC Assets v0.8.3 · IFC 直接导入

侧栏标签为 **Assets**。点击资产缩略图并松手后弹出名称、验收状态、来源与一排 Plan / Front / Side 三个线稿预览按钮；顶部 ⓘ 也可打开当前资产详情。按下时只选择，避免弹窗截断拖动。3D 图只作为网格 Thumbnail，不在弹窗重复显示。详情不再堆放在网格下方。拖动保持原生放置行为，不弹出详情；Esc 取消或释放回侧栏均不插入。拖动定位时显示蓝色包围盒；释放进入加载后蓝框固定保留，框内无描边的绿色半透明填充自下而上增长，完成或取消后清理。加载进度仍由实际导入阶段驱动。

Plan / Front / Side 按钮显示当前单品 IFC 的 Bonsai 单品 SVG 线稿预览。线稿 PNG 保持 1600 × 1200，对应 SVG 同目录保存；3D Thumbnail 则限制最长边 384 像素，以 128 色 PNG 存储，39 张共约 1.75 MB。`catalog.json` 的 `preview_evidence` 校验方向、IFC、SVG 和 PNG 哈希，并将 iso 的用途明确标为 asset_thumbnail。默认打印样式为黑线，不等同于官方蓝线对照图；未验收模型的线稿仍是候选，不会因预览修复而获批。点击图片不会切换场景对象的 representation，也不会重新出图或保存 IFC。

本目录整体复制到其他位置即可使用；需要 Blender 4.5.3 LTS、Bonsai / IfcOpenShell 0.8.4。不要只移动其中的 IFC 或插件子目录。

先用 Bonsai 打开目标 IFC，选择默认楼层，再在 Blender 中运行同目录的 `open_library.py`。按 N → Assets，从资产网格拖入。

39 个单品均含三维、Plan、Front、Side 四个 Body。出图时保持图库插件加载，由它按相机方向选择对应 Body。绿色加载反馈按真实阶段显示。插入仅改变内存，Ctrl+S 使用 Bonsai 原生保存，没有额外写入按钮。

内容：`catalog.json`、`ifc/`、`previews/`、`native-assets/placements.json`、`approvals/`、`addon/`。不再保存任何单品资产 `.blend`。材质用的 `ifc/Materials.blend` 内嵌必要纹理，仍须保留，不依赖原项目 Textures 或本机 Poliigon 目录。批准记录中的历史路径是来源证据，不是运行时依赖。

缩略图使用仅在内存中的空 Collection 资产卡片，包含图片和产品 ID，没有模型对象、网格或 IFC 实体，也不链接到场景。原生资产网格保留拖放手感；放置后直接从对应 IFC 生成一个真实 Bonsai 构件。卡片在刷新、文件加载及撤销／重做后重建，不保存 `.blend` 或用户偏好。若当前文件另有 Collection 资产，原生本地网格可能同时显示它们；非本库卡片不能触发 IFC 导入。

首次升级请在新的 Blender 会话加载，不对旧版 RNA 热替换。旧 v0.6 缓存检查和一次性迁移脚本是历史证据，不适用于 v0.7。`build_native_review_assets.py` 现在只重建摆放数据，不再写资产 `.blend`。

验收类别仍为 8 项已验收、11 项部分验收、20 项未验收；技术验证不改变设计验收状态。源 IFC 保留原 GlobalId 与项目坐标，拖入的新实例生成新 GlobalId，并记录来源。

预览图保留原审核资料；Bonsai 默认黑色打印样式的真实出图验证在上级 `portable-validation/`，并非每项新的设计批准。
