# Core Simulation Change Review

## Bottom line

目前没有看到会让轨迹积分本身明显错误或失控的致命数值改动。相反，RK45、自由传播边界处理、非有限数保护、单轨迹超时保护等部分，相比上游 `temken/DaMaSCUS-SUN` 更稳健。

真正需要警惕的不是“轨迹解算器坏了”，而是下面三类改动已经改变了物理统计量的定义或底层散射模型。因此，如果目标是“得到一个自洽的当前模型数密度分布”，本分支仍然可能是正确的；但如果目标是“严格复现上游 capture/reflection/evaporation 定义和数值结果”，当前实现不能直接视为等价。

## High-risk changes

### 1. Captured 的统计定义已经改变

- `Trajectory_Result::Particle_Captured()` 在 `src/Simulation_Trajectory.cpp` 中只检查最终态是否低于局域逃逸速度。
- 但真正用于统计 `number_of_captured_particles` 的条件在 `src/Data_Generation.cpp` 中是 `trajectory.bincount.is_captured`。
- `bincount.is_captured` 的含义不是“最终束缚”，而是“轨迹在演化过程中曾经进入过负能量区”。
- 在 capture mode 下，代码还会在一旦检测到负能量后提前终止轨迹。

这件事本身不会让轨迹积分失真，但会直接改变下列物理量的含义：

- captured fraction
- evaporation sample
- captured / not captured 的 bincount
- 任何依赖“最终束缚态样本”的数密度分布

结论：

- 如果你的目标是研究“首次进入束缚区后的后续蒸发时间”和相关时间权重统计，这个定义可能是有意且合理的。
- 如果你的目标是复现上游 `DaMaSCUS-SUN` 的 capture 定义，那么这是一个高风险偏离，不能直接拿当前结果与上游表格或论文数值逐项对比。

### 2. 散射核已经不是上游那套实现

本地实现中：

- 总散射率在 `src/Solar_Model.cpp` 中按 `n * sigma_total * v_rel` 计算。
- 末态速度在 `src/Simulation_Trajectory.cpp` 中通过“采样靶粒子热速度 + 采样散射角 + 两体弹性散射运动学”构造。

上游实现中：

- 使用 `q` 和角度采样的 formalism 来构造散射过程。
- 保留了 medium effects 和 `q_cutoff_parameter` 相关接口。

这意味着当前分支不再是对上游散射核的逐行等价实现，而是一个重新参数化后的散射模型。

风险判断需要分情况：

- 对于接触型、短程、角分布良性的模型，这套实现仍然可能是自洽的，未必会导致错误轨迹。
- 对于依赖 medium effects、`q` cutoff、或者 IR 敏感的长程相互作用，这个改动是高风险的，因为模型空间已经和上游不同。

当前 `Dark_Photon` 实现也说明了这一点：

- 本地 `src/Dark_Photon.cpp` 中 `Sigma_Total_Electron/Nucleus` 与角分布接口可以支持 `Contact` 和 `General` 形式。
- 对 `Long-Range` / `Electric-Dipole` 这类 IR 敏感情况，代码会直接报错退出，而不是像上游那样继续走更通用的 `q` 采样路径。

结论：

- 如果当前项目主要跑的是接触型或等效短程模型，这不是“立刻致命”的错误。
- 如果希望宣称与上游 dark photon / medium-effect 结果数值一致，这个改动不能跳过重新标定与对照测试。

### 3. Solar model 的 target set 和接口定义已偏离上游

当前 `Solar_Model` 与上游相比有两个明显变化：

- `target_isotopes` 现在包含全部 isotope，而不是上游默认采用的较小子集。
- `include/Solar_Model.hpp` 的公开构造函数已经简化为无参数版本，不再暴露 `medium_effects` 和 `q_cutoff_parameter` 这两个上游模型开关。

这会带来两类后果：

- 总散射率的归一化和径向分布都会改变。
- 当前结果在模型定义上已经不等于上游默认设置，即使代码没有数值 bug，也不能把结果解释成“原版 DaMaSCUS-SUN 输出”。

这同样不是轨迹求解器的致命错误，但它是一个高风险的 physics-definition drift。

## Likely non-fatal and probably beneficial changes

下面这些改动目前看更像是数值稳健性增强，而不是物理错误来源：

- RK45 的 Fehlberg 系数修正为 `2197 / 4104`，避免使用明显可疑的 `2197 / 4101`。
- RK45 各 stage 直接使用 `solar_model.Mass(r_i)`，而不是对整步只使用单个常数质量近似。
- 针对 `NaN` / `Inf` 的步长收缩与最小步长保护，避免内层自适应循环卡死。
- 外层自由传播阶段增加最大步长限制，减少越过边界和一步跨太远的风险。
- 单轨迹 wall-clock 超时保护，避免单个病态轨迹把整个 MPI rank 卡住。

这些改动更像是“让轨迹更容易算完，而且不至于在病态参数下挂死”，目前没有看到它们会系统性破坏物理分布的证据。

## Practical judgment

### 可以相对放心的结论

- 当前分支没有暴露出“轨迹积分器明显错误”这一类致命问题。
- 如果目标是得到当前模型定义下自洽的束缚/蒸发样本及其数密度分布，那么最大的不确定性不在 RK45，而在物理定义是否与你想研究的对象一致。

### 不能直接宣称的结论

- 不能直接宣称当前结果与上游 `DaMaSCUS-SUN` 的 capture fraction、reflection fraction、evaporation sample 或对应数密度分布完全等价。
- 不能在没有额外 benchmark 的情况下，把当前结果当作“原版模型只是做了工程优化”的结果来引用。

## What to validate next

要判断“虽然改动很多，但最终数密度分布仍然正确”是否成立，最关键的不是继续盯源码，而是做侧对侧 benchmark：

1. 在同一组 `m_DM`、`sigma`、初始条件和 solar model 设定下，对比本地与上游的总散射率 `Gamma(r, v)`。
2. 对比反射分支和束缚分支的分类比例，而不是只看最终总数。
3. 对比 `last scattering radius`、`deepest radius`、`bound time`、`evaporation time` 的分布。
4. 对比最终关心的径向数密度 `n(r)` 是否只是在归一化上不同，还是形状也已经改变。

只有当第 1 到 4 项都说明差别可控时，才可以说“虽然实现和上游不同，但轨迹结果仍然是正确的”。

## Current validation status

已完成的本地验证：

- `build/tests/test_Simulation_Trajectory` 已在 sandbox 中运行并通过，结果为 13 / 13。

需要明确的是：

- 这些单元测试说明当前分支内部行为自洽。
- 它们并不证明当前实现与上游 physics model 等价。

## Scatter kernel: likely reason for the rewrite

这里能从代码结构推断出实现动机，但不能从现有仓库直接证明作者当时的主观意图。

从实现形态看，当前分支把散射过程改写成了下面这套统一接口：

- 先由 `Solar_Model::Total_DM_Scattering_Rate()` 给出总散射率。
- 再由 `Trajectory_Simulator::Sample_Target()` 选择电子或某个 nucleus。
- 然后用 `Sample_Target_Velocity()` 采样热靶速度。
- 最后调用 `DM->Sample_Scattering_Angle_Nucleus/Electron()` 给出角分布，并在两体弹性散射运动学中构造末态速度。

这套写法的直接好处是：

- 轨迹模拟器只依赖于 obscura 的通用粒子接口，不再强耦合上游 `Scattering_Rates.hpp` 那套专门的 `q`-sampling 实现。
- SI、SD、Dark photon 都可以走同一个轨迹引擎，只需要各自提供总截面和角分布。
- 更容易把蒸发时间、capture 标记和在线 bincount 与同一条轨迹主循环集成到一起。

因此，更合理的表述不是“作者随意改了散射核”，而是：

- 当前分支把散射过程重写成了一个以总散射率和角分布为中心的统一 Monte Carlo 方案。

这件事本身并不自动等于错误，但它确实让当前模型不再与上游的 `q`-sampling 版本逐项等价。

## Impact for the current config file

当前 [bin/config_Lingyu.cfg](bin/config_Lingyu.cfg) 对应的关键设置是：

- `run_mode = Parameter point`
- `capture_mode = false`
- `DM_interaction = SD`
- `DM_light = true`
- `DM_cross_section_nucleon = 1e-36 cm^2`
- `DM_cross_section_electron = 1e-80 cm^2`
- `interpolation_points = 1000`
- 主程序中轨迹传播边界设置为 `2 rSun`

对这份配置，散射核改写的实际影响比我之前讨论 Dark photon 时要小很多，原因如下：

### 1. 当前根本没有走 Dark photon 那条分支

配置里是 `SD`，所以实际构造的是 obscura 的 `DM_Particle_SD`，不是 `DM_Particle_Dark_Photon`。

因此以下问题对当前配置不是一阶效应：

- medium effects 缺失
- `q_cutoff_parameter` 缺失
- Dark photon 的 IR 敏感长程散射 formalism 被替换

这些差异对将来跑 Dark photon 很重要，但对你现在这份 SD 配置不是主导项。

### 2. 对当前 SD 配置，角分布本身非常简单

在 obscura 的 `DM_Particle_SD` 实现里：

- 当 `DM_light = true` 时，核散射角是各向同性的，`cos(theta)` 在 `[-1, 1]` 上均匀采样。
- `Sigma_Total_Nucleus()` 对速度没有额外复杂依赖。

这意味着对于当前配置，散射核改写后最核心的动力学输入其实就是：

- 哪个 radius 上发生散射
- 散射前 DM 速度是多少
- 选中了哪个有自旋的核同位素
- 热靶速度抽样给末态速度带来了多少展宽

换句话说，当前配置下真正需要担心的不是“有没有丢掉复杂的 `q`-physics”，而是下面两件更实际的事：

- 总散射率 `Gamma(r, v)` 的空间分布是否和你希望的 SD 模型一致
- 散射后轨迹被分到 captured / not captured 的规则是否符合你的物理定义

### 3. 当前 SD 配置下，电子散射几乎可以忽略

配置里 `DM_cross_section_electron = 1e-80 cm^2`，几乎等于关掉电子散射。

所以对当前 bincount 来说，主导项几乎全部来自核散射，而不是电子散射。

### 4. 当前 Solar model 的“全 isotope 集合”改动，对 SD 的影响比对 SI 小

这是因为 obscura 的 SD 截面实现里：

- 自旋为 0 的 isotope，`Sigma_Total_Nucleus()` 直接就是 0。

因此，虽然当前 `Solar_Model` 把所有 isotope 都放进了 `target_isotopes`，但其中大量 spinless isotope 对当前 SD 配置其实不会贡献散射率。

所以对当前配置而言：

- “target set 扩大”仍然可能影响结果，
- 但影响主要来自那些额外被纳入的 spinful isotope，
- 不会像 SI 那样对所有 isotope 都全面增重。

## What controls the bincount in the current code

如果完全只关注 bincount，那么最重要的是下面这几类因素。

### 1. 轨迹分类规则

这是当前实现里最重要的一条。

最终输出文件里的 bincount 被拆成两套：

- captured_dt_hist / captured_v2dt_hist
- not_captured_dt_hist / not_captured_v2dt_hist

而这两套不是按“最终态是否束缚”来分，而是按：

- `trajectory.bincount.is_captured` 是否为真

也就是：

- 只要轨迹在演化过程中某一步出现过 `E < 0`，它就会被算进 captured 那一边。

这件事对 bincount 的影响非常直接：

- 它不会改变“所有轨迹合起来的总停留时间”，
- 但会显著改变 captured 与 not captured 两侧各自的形状和归一化。

对你现在的分析，如果你关心的是束缚态数密度，这一条比“散射角采样公式的细枝末节”更重要。

### 2. 总散射率 Gamma(r, v)

这决定粒子在什么半径、以多大频率发生散射。

而 `Gamma(r, v)` 本身受以下因素控制：

- DM 类型：当前是 SD
- 核散射截面：当前是 `1e-36 cm^2`
- 电子散射截面：当前几乎可忽略
- Solar model 中各核的 number density
- 只有 spinful isotope 才对当前 SD 模型贡献非零散射率
- `interpolation_points = 1000` 时使用的散射率插值表精度

如果你看到 bincount 形状变了，第一怀疑对象通常应该是 `Gamma(r, v)` 的径向分布，而不是 RK45 本身。

### 3. 散射后的运动学更新

每次散射之后，末态速度不是任意给的，而是由下面几件事共同决定：

- 被选中的 target mass
- 热靶速度采样
- 当前 DM 入射速度
- SD 角分布采样

这些量决定散射后粒子是：

- 更容易向内掉进太阳深处
- 还是更容易被弹回外层
- 以及它在不同半径壳层里停留多久

因此它们会直接影响：

- 哪些 bin 的 `dt_hist` 被加权更多
- 哪些 bin 的 `v2dt_hist` 被加权更多

### 4. 外边界和 binning 几何

当前代码中：

- `NUM_BINS = 2000`
- `BIN_MAX_KM = 2 * R_sun`
- 半径 bin 是线性等宽的
- `r < 0` 或 `r >= 2 R_sun` 的步不会计入 bincount

这意味着：

- bincount 只描述 `0` 到 `2 R_sun` 范围内的停留时间
- 轨迹一旦主要在 `2 R_sun` 外演化，对最终 bincount 就不再贡献

而主程序也把轨迹传播边界设成了 `2 rSun`。所以对当前输出而言：

- `2 rSun` 既是 propagation stopping surface，
- 也是 bincount 的硬截断范围。

如果你未来想从 bincount 反推出更大范围的数密度，这个边界本身就是一个模型假设，不是无害细节。

### 5. 自适应步长与时间积分

bincount 是逐 RK45 步累加的：

- 每一步把上一个位置的 `r` 和 `v^2` 用本步的 `dt` 加权
- 加到对应的径向 bin 里

所以理论上它近似的是时间积分；数值上它会受下面这些因素影响：

- RK45 容差
- 太阳外的自由传播最大步长限制
- 太阳内 `time_step_max = 0.1 / total_rate`
- 是否发生 NaN/Inf 中止
- 是否发生超高速 `v > v_max` 中止

这些通常影响的是数值精度和极端病态轨迹，不是当前 bincount 的主导物理来源，但它们会影响尾部和少数异常轨迹。

### 6. capture_mode 对 bincount 是决定性的开关

当前配置中 `capture_mode = false`，所以会正常累加并写出 bincount。

如果把它改成 true：

- 轨迹一旦第一次进入负能量区就会提前终止
- 最终不会写普通的 bincount 输出文件

因此：

- capture_mode 不是一个小选项，
- 它会直接改变“你到底在积什么量”。

对于当前配置，这一点还好，因为它现在是 false。

### 7. Snapshot 不是最终 bincount 的主导因素

当前配置里 snapshot 是打开的，但 `max_trajectory_wall_time_sec = 0`，表示不做单轨迹 wall-time 截断。

在这种设置下：

- snapshot 主要影响中间状态文件的写出，
- 不应改变最终 bincount 的物理定义。

真正会改变最终 bincount 的不是 snapshot 本身，而是：

- 如果某条轨迹数值崩溃或被别的 guard 中止，
- 它只会贡献中止前已经累加进去的那部分时间。

### 8. sample size 影响统计噪声，不决定单条轨迹物理

当前 `sample_size = 10000`。

它对 bincount 的作用主要是：

- 减小或放大 Monte Carlo 噪声
- 改变误差条
- 改变少数稀有轨迹对整体 histogram 的波动程度

它不会改变单条轨迹的动力学规律。

## Practical takeaway for the current SD run

如果只看你当前这份配置，那么我会按下面的优先级判断 bincount 是否可信：

1. 首先检查 captured / not captured 的划分是否符合你要的物理定义。
2. 然后检查 `Gamma(r, v)` 在太阳内的径向分布是否合理。
3. 再检查 spinful isotope 集合是否和你想比较的理论模型一致。
4. 最后才去担心 RK45 和散射角采样是不是主导误差源。

对当前 SD、低质量模式、电子散射几乎关闭的配置来说：

- 散射核改写不是最危险的部分；
- bincount 最敏感的仍然是分类规则、总散射率分布和 `2 R_sun` 截断。

## Total scattering rate: comparison with the old code

这里专门比较“总散射率”和各 target 对总率的贡献分率。

### 1. 对当前 SD 配置，公式本身基本没有变

旧代码的 `Total_Scattering_Rate_Nucleus()` 有两条路径：

- 如果 `DM.Is_Sigma_Total_V_Dependent()` 为真，或者打开 medium effects，或者 `qMin > 0`，则积分 `dGamma / dq / dcos(theta)`。
- 否则直接使用

```text
Gamma_i(r, v) = n_i(r) * sigma_i(v) * <v_rel>(T(r), m_i, v)
```

当前代码的 `Solar_Model::DM_Scattering_Rate_Nucleus()` 也是：

```text
Gamma_i(r, v) = n_i(r) * DM.Sigma_Total_Nucleus(i, v, r) * <v_rel>(T(r), m_i, v)
```

对当前 [bin/config_Lingyu.cfg](bin/config_Lingyu.cfg)：

- `DM_interaction = SD`
- `DM_light = true`
- `DM_form_factor` 对 SD 无关
- medium effects 在当前代码中不存在，等价于旧代码默认的 off
- `qMin = 0`
- obscura 的 `DM_Particle_SD::Is_Sigma_Total_V_Dependent()` 返回 false

所以，对于当前 SD 配置，旧代码和当前代码都会落到同一个 `n * sigma_total * <v_rel>` 总散射率公式。这里没有发现公式层面的致命差别。

### 2. 真正的差别是 target isotope 集合

旧代码默认只保留：

```text
included_target_indices = {0, 1, 2, 7, 54}
```

对应到太阳模型表，大致是：

- H-1
- He-4
- He-3
- O-16
- Fe-56

但对于 SD 相互作用：

- He-4 自旋为 0，贡献为 0。
- O-16 自旋为 0，贡献为 0。
- Fe-56 自旋为 0，贡献为 0。

因此，在当前 SD 配置下，旧代码有效的核散射总率几乎就是：

```text
Gamma_old_SD ≈ Gamma_H1 + Gamma_He3
```

当前代码的 `Solar_Model` 会纳入 63 个 isotope。对 SD 来说，spinless isotope 仍然贡献 0，但所有 spinful 且 `<sp>, <sn>` 不为零的 isotope 都会进入总率，例如：

- H-1
- He-3
- C-13
- N-15
- O-17
- Na-23
- Mg-25
- Al-27
- Si-29
- P-31
- K-39
- Ti-47 / Ti-49
- V-51
- Mn-55
- Co-59

所以当前代码有效的核散射总率是：

```text
Gamma_current_SD = Gamma_H1 + Gamma_He3 + sum(extra spinful isotopes)
```

这意味着：

- 当前总散射率一定大于或等于旧代码的 SD 总散射率。
- 差别不是来自公式，而是来自额外 spinful isotope 的贡献。
- 由于 H-1 丰度极高，H-1 仍然很可能主导总率。
- 额外 heavy isotope 丰度低，但它们有更大的 reduced mass 和不同 spin matrix element，因此会给总率带来非零修正。

### 3. 电子散射分率对当前配置可忽略

当前配置中：

```text
DM_cross_section_electron = 1e-80 cm^2
```

所以电子散射总率形式上被加入总率，但数值上可以视为关闭。

因此当前配置下的总散射率分率主要是核靶贡献分率：

```text
f_i(r, v) = Gamma_i(r, v) / sum_j Gamma_j(r, v)
```

而不是 electron / nucleus 竞争。

### 4. 对 bincount 的影响

总散射率进入轨迹演化的位置是：

```text
minus_log_xi -= dt * Gamma_total(r, v)
```

所以如果当前 `Gamma_current_SD` 比旧代码更大，则会导致：

- 平均自由传播时间变短。
- 散射次数增加。
- 粒子更早在太阳内部被重新定向或减速/加速。
- captured / not captured 分类比例可能改变。
- `dt_hist` 和 `v2dt_hist` 的径向形状可能改变。

但这个改变的来源应表述为：

```text
额外 spinful isotope 改变了 Gamma_total 和 target sampling fraction
```

而不是：

```text
当前代码用了错误的总散射率公式
```

### 5. 需要进一步精确量化的量

要精确回答“分率差多少”，最应该输出以下诊断表：

```text
r/Rsun, v[km/s], Gamma_total_old, Gamma_total_current,
Gamma_current / Gamma_old,
Gamma_H1 / Gamma_current,
Gamma_He3 / Gamma_current,
sum_extra_spinful / Gamma_current
```

其中旧代码的 SD target set 可以按 `{H-1, He-4, He-3, O-16, Fe-56}` 实现；但因为 He-4/O-16/Fe-56 对 SD 为 0，实际就是 `{H-1, He-3}`。

结论：

- 对当前 SD 配置，总散射率公式和旧代码默认路径一致。
- 主要差别是 target set：旧代码 SD 有效贡献约为 H-1 + He-3，当前代码还包含其他 spinful isotope。
- 这会提高总散射率并改变 target contribution fraction，但是否显著需要用上面的诊断表定量。
