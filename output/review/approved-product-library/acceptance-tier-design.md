# 图库验收状态分层

## Product Shape

- User Job / Workflow Mode：浏览已批准的单品，加载独立 IFC 展示副本，继续复审尚未验收的场景。
- Product-Best Shape / Experience Architecture：一个单品库入口，每项缩略图、来源及验收状态同时可见；场景技术验证不能替代用户批准。
- Surface Ownership：Blender「单品库」N 面板；计数在目录顶部，所选产品的状态在产品名称与预览之间。
- Information Architecture：单品批准、场景批准、迁移技术验证分别存储；收录授权单独保留证据。
- What Must Be Separate：正式 IFC 写入权限、单品包写入权限、展示副本操作；前者未授权。
- Old Surfaces To Remove/Avoid：不再用笼统的“已验收”或同一个勾号表示所有产品均已完成场景验收。

## Acceptance Contract

- Capability：两类产品均可搜索、看四张预览和载入 XY 阵列。
- Placement：同一 N 面板显示总数、场景已验收/待验收计数，所选产品显示完整状态。
- Negative Placement：不得将新增七项移入 whole_product_accepted，不能把迁移 pass 变成 scene approved。
- Workflow：查看预览时仍保留所选产品与状态，不为待审产品另造重复入口。
- State/Result：原八项场景 approved，新增七项场景 pending；尚未技术交付者进入迁移 pending，不提前显示为可用资产。
- Evidence：单位测试拒绝缺授权或批准状态提升；真实 Blender 窄面板截图覆盖两种状态；搜索、预览、载入、保存重载与 IFC 哈希验证通过。

## Implementation Path

- Product-best version：完整 15 项统一图库、明确状态、独立 IFC 和真实缩略图。
- Safer incremental version：产品代理逐项完成技术验证后收录，保留原八项 IFC 与批准记录。
- Prototype to test：选一个场景待验收产品检查窄面板状态与四张预览，再验证全量阵列。
- Open Questions：无新增授权问题；七项场景最终仍由用户验收。
