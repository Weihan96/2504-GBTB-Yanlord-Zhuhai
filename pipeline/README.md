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
bun run pipeline:coordinate-noise-candidate
bun test
```

- `pipeline:snapshot`：读取 IFC 和现有 SVG，生成源哈希与结构快照。
- `pipeline:check`：生成 QA JSON/Markdown；允许已知阻塞存在，便于持续盘点。
- `pipeline:gate`：执行同一套检查；存在阻塞时返回非零退出码，用于发布门。
- `pipeline:coordinate-audit`：只读统计对象 origin、IFC 长度数值和非整数分布；当前 origin 人审阈值为 `0.1 mm`。
- `pipeline:coordinate-noise-candidate`：只在 `build/candidates/` 生成副本，将距离整数不超过 `0.01 mm` 的 `IfcLengthMeasure` 尾数写为精确整数；不会覆盖源 IFC。
- 坐标工具不提供“全部原点取整”或“施工几何取整”写入命令；超过人审阈值的 origin 和墙体尺寸必须经过 Blender 审核与受控 IFC 写入批次。
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
