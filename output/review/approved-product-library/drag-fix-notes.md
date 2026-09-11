# Assets 0.8.2 拖放修复

原生资产列表在鼠标按下时激活行；0.8.1 的选择回调直接打开模态详情弹窗，截断随后拖动。已改为按下仅选择，完成点击的 RELEASE/CLICK 才打开详情；拖动继续交给原生 drag_operator，不增加保存入口。`execute` 本身不弹窗。

验证使用 Blender 4.5.3 的真实窗口鼠标事件队列，未直接调用插入函数替代拖放。`drag-gesture-validation.json` 记录了完整事件：点击弹窗、Esc 取消、拖回侧栏取消、成功插入 CHA01。两种取消后 IFC 内存哈希均与操作前相同；成功插入增加一个产品，保留四个 Body 表达。另一次 BST02 实际拖放记录了 0、0.05、0.25、0.5、0.65、0.85、0.95、1.0 的真实阶段，`drag-actual-loading.png` 是真实导入时的无描边绿色反馈，不是合成进度。42 项单元测试通过。

测试期间一次卸载/重注册资产选择 Operator 导致 Blender 原生 UIList 保留失效 RNA 指针而崩溃（堆栈落在 panel draw 的 template_asset_view）。没有保存数据；随后改用完整冷启动，交付窗口没有热替换 RNA。旧截图中的“双击查看”是中间测试文案，最终保留单击松手查看。

启动时发现 9881 已被另一个橱柜任务占用；bootstrap 现在探测可用端口并输出 PID/port。没有关闭或修改那个任务。

最终冷启动：PID 68207，任务 owner `d63e2a037659c572260e`，bridge 54299；Assets 0.8.2、39 项、IFC 元素 0、无活动拖动、未保存 Blend、无 Bonsai Blend 警告。Scene 为 `🔍 逐个审核高模图纸候选`；显示单位 MILLIMETERS，与验收 IFC 声明的毫米一致。两次测试插入仅在测试内存中，最终窗口从未变更的空验收 IFC 重开。

正式 IFC SHA-256：`7a50b87e8f48a7c2bfbaa9f1dabdfca7155fa684b2325c66aa1ae15c8ab0a25c`。验收 fixture SHA-256：`f8b249ccd1f8bc9c9e3d58d2b878e27818d8a8237631e89806a05210f8daa75c`。两者操作前后相同；421 个保护来源文件核验不变。安装 zip、runtime addon 镜像和 manifest 已同步。暂存区哈希仍为 `fa9b7329175207683069cb18cc2f5015e1e7043ac52e64f882e95720168c53ea`；本轮未暂存、未提交、未操作 stash。
