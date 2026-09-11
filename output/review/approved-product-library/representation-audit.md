# 单品表示核查与清理候选

状态：39 项核查完成；IFC 修正尚未执行。所有 IFC、批准状态及当前 Blender 会话保持不变。本报告不是修正版 IFC 的验收证明。

## 核查结论

- 39 份活动单品 IFC 均只有一个实际产品（电动止回阀另有两个必要开孔），均没有 IfcAnnotation 实体。截图中的 Annotation 是表示上下文名称，不是额外实体数量。
- 15 份包含已批准三视图：其中 13 份放在 Annotation 上下文；WD03 与 Falper 使用原生 Drawing 不选取的自定义上下文。14 份同时保留旧 Body 平面／立面。
- 24 份是 Body-only 候选，不能当作已完成三视图接入。不能从“representation 数量足够”推断接入正确。
- 没有发现与 MODEL_VIEW 的 Items 前向图指纹完全相同的旧 Body 表示；不能把所有旧平面／立面称作重复高模。旧版本与批准版本共存的问题仍需解决。
- 安装版 Bonsai 的上下文选择函数已直接执行；15 份的 Front/Side 都使用 ELEVATION_VIEW，原生选择器不按相机方向区分。Gessi 54294 的内存实验确认：只将 Annotation 改名 Body，会同时交给几何迭代器旧立面、Front、Side 三组几何。没有生成新 SVG，因此不声称已观察到最终图纸重影。

## 方法与需要确认的边界

Bonsai Course Operator 课程 109000 的 00:27 支持产品平面使用 Plan / Body / PLAN_VIEW；Annotation2D 是几何类型，不是 Annotation 上下文。证据为 embedded-course-index，未声称直接打开课程私有截图。安装版源码与 IFC4 资料用于交叉核对。

建议保留原始 3D Body 和已批准三视图；用批准的 Plan 替换旧 Plan Body，而不是追加一套。Front/Side 必须独立保留视向语义，不能改成两个无法区分的 ELEVATION_VIEW Body。

需要人工确认的架构选择：继续保留独立官方二维图，并让图库在出图时按相机方向选用正确线稿；或改成由统一简化三维几何生成立面（不能保证保留已批准官方二维画法）。推荐前者，但不能声称它是现成的原生 Bonsai 自动行为。Bonsai Course Operator 对分歧架构要求人工审阅，因此本轮不批量改写。

实现前者还必须处理 append_asset 的上下文合并：当前版本按 ContextType/ContextIdentifier/TargetView 合并，单独增添 Front/Side 名称不足以保留方向。不得修改公共 Provider 或增设保存按钮；仍由 Bonsai Ctrl+S 保存。

修正验收：先 Gessi 54294 试行；原点、单位、3D Body、已批准线稿和样式不变；保存重载后真实 Create Drawing 的 Plan/Front/Side SVG 必须分别验证，没有旧线或另一视向混入；再批量覆盖 39 项对应规则。正式 IFC 前后 SHA-256 必须不变。

## 全部单品

| 产品 | 表示数 | 旧 Body 非模型视图 | 已批准线稿表示 | 类别 |
|---|---:|---:|---:|---|
| hima01 | 6 | 2 | 3 | partial |
| geberit-duofix-sigma-224-212 | 6 | 2 | 3 | approved |
| geberit-154-446-ks-1 | 6 | 2 | 3 | partial |
| gessi316-54294 | 6 | 2 | 3 | approved |
| bed02 | 6 | 2 | 3 | partial |
| bed01 | 6 | 2 | 3 | approved |
| sxb010 | 5 | 1 | 3 | approved |
| cha02 | 2 | 1 | 0 | pending |
| sis04 | 3 | 2 | 0 | pending |
| 505-up-v1-lp-s | 2 | 1 | 0 | partial |
| miamisoft-e09 | 5 | 1 | 3 | approved |
| miamisoft-e07 | 2 | 1 | 0 | pending |
| marilyn-01 | 2 | 1 | 0 | partial |
| miamisoft-f03 | 2 | 1 | 0 | pending |
| trap01 | 5 | 1 | 3 | approved |
| marilyn-02 | 2 | 1 | 0 | pending |
| cha01 | 2 | 1 | 0 | pending |
| miamisoft-h01 | 2 | 1 | 0 | pending |
| falper-sorgente | 5 | 1 | 3 | partial |
| geberit-154-154-00 | 2 | 1 | 0 | pending |
| bst03 | 2 | 1 | 0 | pending |
| fau02 | 2 | 1 | 0 | pending |
| gessi316-54146 | 2 | 1 | 0 | partial |
| wd03 | 5 | 1 | 3 | approved |
| miamisoft-i03 | 2 | 1 | 0 | pending |
| geberit-146-140 | 2 | 1 | 0 | partial |
| gessi316-54145 | 5 | 1 | 3 | partial |
| street-h | 5 | 1 | 3 | partial |
| wd02 | 5 | 1 | 3 | approved |
| tab02 | 4 | 0 | 3 | partial |
| bst02 | 2 | 1 | 0 | pending |
| street | 2 | 1 | 0 | pending |
| gessi316-54093 | 1 | 0 | 0 | pending |
| geberit-154-154-00-1-f | 2 | 1 | 0 | pending |
| geberit-115-770 | 2 | 1 | 0 | pending |
| wd01 | 2 | 1 | 0 | pending |
| gessi316-54038 | 2 | 1 | 0 | pending |
| electric-flue-check-valve | 1 | 0 | 0 | pending |
| bst01 | 2 | 1 | 0 | pending |

## 清理清单（未授权删除）

这里只盘点 output/review 下 IFC 与 Blender 文件；不是对仓库其他内容的删除授权。字面路径引用检查不能证明运行时可删除。所有 SVG/PNG 验收证据、官方 CAD 原件及下载链接/哈希、批准记录、scene-recipe 和必要材质依赖保留。

| 分类 | 文件数 | 大小（十进制 GB） |
|---|---:|---:|
| blend_cache_candidate_requires_dependency_and_session_check | 63 | 10.310 |
| full_scene_legacy_candidate_blocked | 18 | 1.582 |
| isolated_legacy_candidate_requires_dependency_and_render_check | 37 | 0.048 |
| keep_active_library_ifc | 39 | 0.154 |
| keep_active_material_or_library_asset | 49 | 0.142 |
| keep_full_scene_pending_product_or_unresolved_scope | 4 | 0.352 |
| keep_small_multi_element_fixture_dependencies_unresolved | 4 | 0.007 |

### 优先复查的旧完整场景 IFC

以下每份含 714 个 IfcElement，与正式项目的元素数量相同；已有对应单品包，但表示修正与替代出图验证未完成，暂不可立即删除。

| 文件 | MB |
|---|---:|
| output/review/highpoly-types/bed01/Baxter-Casablanca-BED01-derived-drawing.ifc | 87.7 |
| output/review/highpoly-types/geberit-duofix-sigma-224-212/Geberit-Duofix-Sigma-224-212-derived-drawing.ifc | 88.6 |
| output/review/highpoly-types/geberit-duofix-sigma-224-212/bonsai-drawings/main-bathroom/duofix-main-bath-front-session.ifc | 89.4 |
| output/review/highpoly-types/geberit-duofix-sigma-224-212/bonsai-drawings/main-bathroom/duofix-main-bath-left-side-session.ifc | 88.9 |
| output/review/highpoly-types/geberit-duofix-sigma-224-212/bonsai-drawings/main-bathroom/duofix-main-bath-plan-session.ifc | 89.0 |
| output/review/highpoly-types/gessi316-54294/bonsai-drawings/GESSI316-54294-MAIN-BATH-FRONT-session.ifc | 88.2 |
| output/review/highpoly-types/gessi316-54294/bonsai-drawings/GESSI316-54294-MAIN-BATH-PLAN-session.ifc | 87.4 |
| output/review/highpoly-types/gessi316-54294/bonsai-drawings/GESSI316-54294-MAIN-BATH-SIDE-session.ifc | 87.4 |
| output/review/highpoly-types/miamisoft-e09/Baxter-Miami-Soft-E09-derived-drawing.ifc | 87.7 |
| output/review/highpoly-types/sxb010/bonsai-drawings/dining-bay/sxb010-dining-bay-front-session.ifc | 87.4 |
| output/review/highpoly-types/sxb010/bonsai-drawings/dining-bay/sxb010-dining-bay-plan-session.ifc | 87.4 |
| output/review/highpoly-types/sxb010/bonsai-drawings/dining-bay/sxb010-dining-bay-side-session.ifc | 87.4 |
| output/review/highpoly-types/sxb010/sxb010-derived-drawing.ifc | 87.4 |
| output/review/highpoly-types/trap01/Geberit-151.116.11.1-TRAP01-derived-drawing.ifc | 87.5 |
| output/review/highpoly-types/trap01/Geberit-151.116.11.1-TRAP01-official-detail-v3.ifc | 87.9 |
| output/review/highpoly-types/trap01/Geberit-151.116.11.1-TRAP01-official-detail-v4.ifc | 88.0 |
| output/review/highpoly-types/wd02/Poliform-Senzafine-WD02-derived-drawing.ifc | 87.4 |
| output/review/highpoly-types/wd03/Poliform-Senzafine-WD03-derived-drawing.ifc | 87.4 |

其余 .blend/.blend1 包括审核相机、待审产品和旧阵列，需要逐个确认替代产物及是否正被使用；不应整目录清除。详见 cleanup-audit.json 的逐文件哈希、分类和引用清单。

正式 IFC SHA-256：`7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c`，操作前后不变。此次未操作 stash、未修改正式 IFC、未删除旧审核文件。
