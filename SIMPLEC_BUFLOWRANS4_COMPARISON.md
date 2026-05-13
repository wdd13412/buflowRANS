# `buflowRANS4.f90` 与当前 SIMPLEC 的真实对比报告

## 1. 文件定位结论

这次已经确认真正的 `buflowRANS4.f90` 可以读取。它不在当前本地 `work` 工作树中，而在 GitHub 远程分支：

- `github/main:buflowRANS4.f90`
- `github/codex/verify-simplec-solver-correctness-yothce:buflowRANS4.f90`

两个分支中的 `buflowRANS4.f90` 是同一个 blob：`d9520fa040f8faa86affbeb3c1bbbee0885d31c6`，文件大小 `187615 bytes`，共 `5183` 行。

> 更正：之前用本地历史 commit `b17a5cb` 近似代表 `buflowRANS4.f90` 是不严谨的。后续对比必须以 `github/main:buflowRANS4.f90` 为准。

## 2. 实测对比结论

我用当前汽车网格编译并运行了真正的 `buflowRANS4.f90` 到 SIMPLEC 350 步：

- `Umax = 17.000 m/s`
- `dp = 200.000 Pa`
- 入口边界映射为 inlet，`symmetry` 映射为 wall

这说明用户判断正确：真正的 `buflowRANS4.f90` 不会出现当前版本的 `Umax≈55.4 m/s`。

但是，这个 `Umax=17` 和 `dp=200` 不是纯方程自然收敛得到的，而是由 `simplec_apply_low_mach_bounds` 每步强制控制得到的：

1. 全场表压先减均值；
2. 若全场压差超过 `SIMPLE_MAX_PRESSURE_RANGE = 200 Pa`，则整体缩放压力；
3. 若任意单元速度超过 `SIMPLE_MAX_SPEED = 17 m/s`，则按速度模长整体缩放该单元速度；
4. 随后用 clipped 后的速度重建所有内部面通量和边界面通量。

因此 `buflowRANS4.f90` 的速度云图能向流场扩散，不能简单解释为“当前真实 SIMPLEC 方程已经正确而后续某个小改动破坏了它”。更准确地说：`buflowRANS4.f90` 用一个强 bounded 后处理同时完成了数值稳定、极值控制和通量重建。

## 3. 与当前版本的关键差异

### 3.1 最大差异：`simplec_apply_low_mach_bounds`

`buflowRANS4.f90` 在每次压力修正后调用：

```fortran
call simplec_apply_low_mach_bounds(mesh, ss, boundaryConditions, U_init)
```

当前版本不再调用这个函数，而是改为：

```fortran
call simple_enforce_inlet_velocity(mesh, ss, boundaryConditions, U_init)
```

这正好解释了两个现象：

- `buflowRANS4.f90` 没有 `Umax≈55`，因为速度被硬限制到 `17 m/s`；
- 当前版本暴露 `Umax≈55`，因为取消了全场速度 clipping，只保留入口缓冲约束。

### 3.2 压力/速度梯度模型不同

`buflowRANS4.f90` 的 SIMPLEC 仍采用 x 方向压力梯度模型：

- 动量方程中 `gP = 0` 后只赋值 `gP(:,1)`；
- Rhie-Chow 压力插值中也只使用 x 方向压力梯度；
- 压力修正速度更新中也只修正 x 方向。

当前版本改成了完整三维压力梯度和三维速度修正。这从 CFD 理论上更完整，但它也让横向速度和局部近壁加速真正暴露出来，不再被 x-only 模型隐藏。

### 3.3 松弛系数不同

`buflowRANS4.f90`：

- `SIMPLE_ALPHA_U = 0.5`
- `SIMPLE_ALPHA_P = 0.002`

当前版本：

- `SIMPLE_ALPHA_U = 0.20`
- `SIMPLE_ALPHA_P = 0.02`，并叠加动压尺度单步限制

这意味着当前版本的动量预测更保守，但压力修正框架更复杂；它不是 `buflowRANS4.f90` 的简单延续。

### 3.4 当前版本新增了方程层稳定化

当前版本相对于 `buflowRANS4.f90` 新增了：

- 伪瞬态动量对角项；
- `aP-H` / floor 的压力修正 D 系数；
- `SIMPLEC_RC_DAMPING`；
- 非正交显式修正；
- 出口压力弱参考；
- 入口缓冲固定速度；
- 大量诊断输出。

这些修改的目的都是让去掉全场 clipping 后仍能运行，但实测说明目前还没有同时满足：

- 不用 clipping；
- `Umax` 接近 11–17 m/s；
- 车身速度扰动能自然扩散到流场；
- 压力/速度云图有物理意义。

## 4. 当前最可能的问题重点

基于真正的 `buflowRANS4.f90` 对比，问题重点不应再归因于单独的 `simple_apply_wall_functions`。更准确的结论是：

1. **`buflowRANS4.f90` 的好云图很大程度上依赖 `simplec_apply_low_mach_bounds` 的全场限幅和面通量重建。** 这能压住 `Umax`，也能让通量场每步跟 clipped 后的速度一致。
2. **当前版本取消全场限幅后，真实方程层稳定化不足。** 因此局部近壁/尾迹速度峰值升到约 `55 m/s`。
3. **当前版本的“速度只集中第一层”与“Umax 偏高”是同一个问题的两面：** 既要让扰动扩散，又要避免局部峰值，需要更物理的动量扩散/湍流/壁面剪切处理，而不是只调边界名。
4. **如果短期目标是复现 `buflowRANS4.f90` 云图形态，必须恢复或等价替代它的 bounded 通量重建机制。** 如果长期目标是高精度 CFD，则不能简单恢复全场 clipping，而应把这个机制替换为物理壁函数、受限对流格式和一致的压力-速度耦合。

## 5. 建议的下一步验证顺序

为了避免继续“无头苍蝇”式修改，下一步建议做三个受控实验，每次只改一个因素：

1. **当前版本临时恢复 `simplec_apply_low_mach_bounds`**：验证是否立刻恢复 `Umax=17/dp=200` 和类似 `buflowRANS4.f90` 的云图扩散。
2. **在 `buflowRANS4.f90` 中关闭 `simplec_apply_low_mach_bounds`**：验证它是否也会出现 `Umax` 飙升或云图变坏。
3. **只保留通量重建、不做速度/压力 clipping**：区分问题到底来自“bounded 限幅”还是“每步按速度重建面通量”。

这三个实验完成后，才能确定应当恢复哪一部分、替换哪一部分，而不是继续整体换 SIMPLEC 方法。
