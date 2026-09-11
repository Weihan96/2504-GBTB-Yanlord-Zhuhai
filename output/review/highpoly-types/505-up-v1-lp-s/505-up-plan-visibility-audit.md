# 505 UP Plan 深度可见性诊断

本轮只做诊断，不改最终 Plan、不写派生 IFC、不刷新 Bonsai Drawing。Front 与 Side 已冻结。

## 结论

上一轮把局部 `z≈0` 的 76 mm 板误认作顶盖；实例放置矩阵反转了局部 Z，所以它们实际位于世界高程 0–76 mm。真正顶面为 A（2380 mm）。当前黑线 P01 把 A/C/D/E 合成外轮廓，P20 只补了 A/C；D/E 及 C/E 的深度区域没有完整语义分割，P02–P05 更是没有产品语义的微小闭环，因此当前 Plan 不可批准。

## 诊断图

![组件分色](./505-up-plan-component-visibility-diagnostic.svg)

![深度图](./505-up-plan-depth-visibility-diagnostic.svg)

红色虚线为实际可见区域边界，黑线为当前候选。红线未被黑线覆盖的位置就是缺失的语义分割。
