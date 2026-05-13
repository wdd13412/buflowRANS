# SIMPLEC 云图回归受控实验开关

## 目的

为了确认 `buflowRANS4.f90` 与当前版本速度云图差异到底来自哪一段代码，本轮在 `buflowRANS.f90` 顶部新增了一个编译期实验开关：

```fortran
integer, parameter :: SIMPLEC_EXPERIMENT_MODE = SIMPLEC_EXPERIMENT_BUFLOW4_BOUNDS
```

修改这个值后重新编译即可切换实验，不需要反复手工改多个子程序。

## 模式说明

| 模式 | 参数值 | 作用 | 主要验证点 |
|---|---:|---|---|
| `SIMPLEC_EXPERIMENT_CURRENT` | `0` | 当前方程层求解器：完整三维压力梯度，无全场 clipping，只保留入口缓冲。 | 检查不靠全场限幅时真实 SIMPLEC 方程是否仍会出现 `Umax` 偏高/云图不扩散。 |
| `SIMPLEC_EXPERIMENT_BUFLOW4_BOUNDS` | `1` | buflowRANS4 对照实验：x 方向压力梯度 + 全场压力/速度 bounded 后处理 + 面通量重建。 | 检查 `buflowRANS4.f90` 的云图行为是否主要来自 bounded 后处理。当前默认是这个模式。 |
| `SIMPLEC_EXPERIMENT_FLUX_REBUILD_ONLY` | `2` | 当前方程层求解器 + 入口缓冲 + 只重建面通量，不裁剪压力/速度。 | 区分“面通量重建”与“全场 clipping”到底哪个对速度云图影响更大。 |

## 当前默认实验

当前默认：

```fortran
integer, parameter :: SIMPLEC_EXPERIMENT_MODE = SIMPLEC_EXPERIMENT_BUFLOW4_BOUNDS
```

也就是说，你直接编译运行时，会做 `buflowRANS4` 风格的后处理对照实验：

1. SIMPLEC 动量/压力修正中的压力梯度切到 x-only；
2. 每步压力修正后执行 `simplec_apply_low_mach_bounds`；
3. 如果速度超过 `17 m/s` 则按速度模长缩放；
4. 如果压力范围超过 `200 Pa` 才缩放压力；
5. 最后用后处理后的速度重建内部/边界面通量。

注意：这不是完全回到 `buflowRANS4.f90`，因为当前主方程里仍保留 PCG、伪瞬态对角、`aP-H`、Rhie-Chow damping、非正交项等改动。因此本模式用于定位，不代表最终高精度 CFD 方案。

## 建议你跑图的顺序

建议每次只改一行 `SIMPLEC_EXPERIMENT_MODE` 并重新编译：

1. 先跑默认 `SIMPLEC_EXPERIMENT_BUFLOW4_BOUNDS`：如果速度云图明显接近 `buflowRANS4.f90`，说明 bounded 后处理/通量重建影响最大。
2. 再改成 `SIMPLEC_EXPERIMENT_FLUX_REBUILD_ONLY`：如果云图仍好，但 `Umax` 不被硬压到 17，说明关键是面通量重建；如果云图又坏，说明 clipping 本身是主要因素。
3. 最后改成 `SIMPLEC_EXPERIMENT_CURRENT`：作为当前无 clipping 方程层版本的基线。

## 运行命令

```bash
gfortran -O2 -o buflow_run buflowRANS.f90 run_parameter.f90
./buflow_run
```

如果只想自动跑到 350 步并看数值范围：

```bash
./check_simplec_350.py --exe ./buflow_run --target-iter 350 --timeout 300 --skip-diagnostics-check --max-umax 17.1 --max-dp 60
```

默认模式下我这里跑到第 350 步的结果是：

- `Umax = 17.0 m/s`
- `dp ≈ 34 Pa`

这说明当前代码的 bounded 后处理能压住速度，但因为当前方程层已经不同于 `buflowRANS4.f90`，压力范围没有自然增长到 200 Pa；这正是下一步需要通过模式 2/0 区分的重点。

## 入口/出口压力带说明

根据你跑图的反馈，默认 `SIMPLEC_EXPERIMENT_BUFLOW4_BOUNDS` 已能恢复较正常的速度扩散，但入口会出现高压带、出口会出现低压带。这个现象说明 `buflowRANS4` 风格的 bounded 后处理虽然能重建较平滑的速度/通量场，但没有给入口/出口远场压力足够强的表压参考。

本轮在 `simplec_apply_low_mach_bounds` 内加入了 `simplec_enforce_farfield_pressure`：只把 inlet/outlet patch 相邻 owner cells 的表压设为 `0 Pa`，即大气表压。它不会全场缩放压力，也不会改车身壁面压力结构；目的只是验证入口/出口高低压带是否来自远场压力参考缺失。

如果这一步后速度云图仍保持较正常，同时入口/出口压力带消失，则说明下一步要把入口/出口远场压力条件写成方程级边界条件，而不是只做后处理。
