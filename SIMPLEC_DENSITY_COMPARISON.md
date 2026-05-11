# SOLVER_DENSITY_BASED 与 SOLVER_SIMPLE/SIMPLEC 简要对照报告

## 1. 本次结论

- 不再改动已经确认需要保留的两点：密度基 JST 压力传感器与 SIMPLEC 压力耦合均只使用流向（x）压力梯度；`symmetry` 在 SIMPLEC 边界名映射中与密度基保持一致，按 wall 处理。
- 当前 OpenFOAM 边界名映射不是“入口/壁面搞反”的主要嫌疑：`airfoil` 映射为 wall，`inlet` 映射为 inlet，`outlet` 映射为 outlet，`symmetry` 映射为 wall。运行时会打印映射，便于直接检查。
- 若 SIMPLEC 云图表现为压力极值主要出现在入口/出口、壁面压力差异很弱，更大的嫌疑是 SIMPLEC 当前的低 Mach bounded 后处理（尤其全场 pressure-range scaling）和压力基边界模型过强地限制了压力动态范围，而不是 OpenFOAM boundary 读取顺序反了。

## 2. 两种求解器的大方向差异

### 2.1 求解变量不同

- `SOLVER_DENSITY_BASED` 直接推进守恒变量：密度、动量、总能量以及湍流守恒量；压力和温度由状态方程/守恒量恢复。
- `SOLVER_SIMPLE`/SIMPLEC 直接求原始变量层面的压力修正与速度修正；目前低 Mach 处理使用参考密度，压力以 gauge pressure 形式参与修正。

### 2.2 压力的角色不同

- 密度基方法中，压力进入欧拉通量、能量方程和人工耗散传感器；压力波/声速会影响谱半径和耗散。
- SIMPLEC 中，压力主要通过动量源项、Rhie-Chow 面通量和压力修正方程来强制连续性；压力本身不由能量方程推进。

### 2.3 边界条件用法不同

- 密度基边界通过 `updateBoundaryConditions_RANS` 分发到 wall/empty/inlet/outlet，边界面上写入 RANS 通量。
- SIMPLEC 边界直接参与动量方程、压力修正方程和面质量通量：入口固定速度/通量，出口固定压力修正参考，wall/symmetry 目前按 wall 处理为零法向通量。

### 2.4 数值耗散与稳定化不同

- 密度基方法依赖中央差分通量 + JST 人工耗散，并用局部时间步推进。
- SIMPLEC 依赖压力修正、动量方程松弛、PCG 压力修正求解、伪瞬态对角项，以及当前低 Mach bounded 后处理。

### 2.5 云图差异的主要来源

- 密度基方法的壁面压力来自守恒方程、能量方程和边界通量的耦合，通常更容易形成壁面停滞/绕流压力分布。
- SIMPLEC 当前为了避免汽车网格发散，对压力范围和速度幅值做了硬限制；这会保证 `dp≈200 Pa` 与 `Umax≤17 m/s`，但也可能压平壁面压力差异，使入口/出口边界附近成为全场极值位置。

## 3. 后续建议

1. 保留当前 x 方向压力梯度和 symmetry-as-wall 设置，不再回退。
2. 若要让 SIMPLEC 云图更有物理意义，下一步应优先检查/替换 `simplec_apply_low_mach_bounds` 中的全场压力缩放，改成更局部、更物理的压力参考处理，而不是全场线性压缩压力范围。
3. 入口/出口边界建议继续对照密度基通量边界：入口速度方向与面法向符号、出口压力参考、wall 零通量三者应分别打印质量通量积分进行验证。
