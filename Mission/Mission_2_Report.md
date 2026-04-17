# 任务书 2 执行报告

---

## 任务 4：输出数据梳理

### 运行模式

程序支持两种运行模式，通过配置文件 `run_mode` 参数选择：

### 模式一：`"Parameter point"` —— 单参数点模拟

**入口**：`main.cpp` L48-74

**流程**：
1. 解析配置 → 加载太阳模型 → 预计算散射率插值表
2. `Simulation_Data::Generate_Data()` 运行 MPI 并行模拟循环
3. 每个 rank 独立模拟轨迹直到达到目标捕获数或安全阀上限
4. MPI 归约所有 rank 的直方图和蒸发记录
5. 输出文件 + 可选反射谱

**输出目录**：`{output_dir}/results_{log10(m_DM/GeV)}_{log10(σ_p/cm²)}/`

**输出文件**：

| 文件名 | 内容 | 格式 |
|--------|------|------|
| `captured_bincount.txt` | 所有被捕获轨迹的 r-histogram 和 v²-histogram | 头部元信息 + 2000行，每行：`bin_index Σdt[s] Σ(v²·dt)[km²/s]` |
| `not_captured_bincount.txt` | 所有未捕获轨迹的同上 | 同上 |
| `evaporation_summary.txt` | 被捕获轨迹的蒸发时间列表 | 头部元信息 + 每行：`trajectory_id t_evap[s] truncated(0/1)` |
| `snapshot_{time}s.txt` | （可选）bincount 在指定物理时间的快照 | 同 bincount 格式 |

**文件头部元信息**（每个文件均包含）：
```
# DM_mass_GeV = ...
# DM_sigma_cm2 = ...
# total_trajectories = ...
# captured_particles = ...
# capture_rate = ...
# EARLY_STOP: max_trajectories reached  (仅当触发安全阀时)
```

**终端输出摘要**（`Print_Summary()`）：
- 总模拟轨迹数、平均散射次数
- Free/Reflected/Captured 粒子比例
- 蒸发时间中位数（仅非截断轨迹）
- 轨迹生成速率、捕获速率、总运行时间

### 模式二：`"Parameter scan"` —— 参数空间扫描

**入口**：`main.cpp` L76-101

**流程**：
1. 可选：计算 halo 约束曲线
2. 构建 `(m_DM, σ)` 参数网格
3. Full scan 或 STA 边界跟踪算法
4. 每个网格点调用 `Compute_p_Value()` → 模拟 + 反射谱 → p 值

**输出目录**：`{output_dir}/results/{ID}/`

**输出文件**：

| 文件名 | 内容 |
|--------|------|
| `P_Values_Grid.txt` | p 值网格（cross_sections × masses 矩阵） |
| `Halo_Limit_{CL}.txt` | halo 约束排斥曲线 |

### Bincount 直方图详解

- **分bin方式**：以径向距离 $r$ 为自变量，$r \in [0, 2R_\odot]$，线性等宽 2000 bins
- **bin 宽度**：$\Delta r = 2R_\odot / 2000 \approx 695.7$ km
- **r-histogram** (`dt_hist`)：$\sum \Delta t$，每个 bin 中粒子停留的总时间
- **v²-histogram** (`v2dt_hist`)：$\sum v^2 \cdot \Delta t$，每个 bin 中速度平方加权的停留时间
- **时间权重**：前向差分 $\Delta t_i = t_{i+1} - t_i$，第 $i$ 步的值归入 $r_i$ 所在的 bin
- **按轨迹类型分类**：captured（$E \leq 0$ 出现过）和 not_captured 分别累积

---

## 任务 5：误差分析方案

### 5.1 问题描述

当前 bincount 文件输出的是所有轨迹的**总和**直方图，没有每个 bin 的统计误差信息。需要为每个 bin 添加上下误差，以及捕获率的统计误差。

### 5.2 误差来源分析

1. **采样统计误差**（主要）：有限的轨迹数量导致每个 bin 的 $\sum \Delta t$ 有统计涨落
2. **捕获率 Poisson 误差**：捕获/未捕获是二项分布
3. **轨迹间方差**：不同轨迹对同一 bin 的贡献不同

### 5.3 方案设计

#### 方案 A：在线方差累积（推荐）

在模拟过程中，对每个 bin 同时累积 $\sum x_i$ 和 $\sum x_i^2$（Welford 在线算法），其中 $x_i$ 是第 $i$ 条轨迹对该 bin 的贡献。模拟结束后可直接计算方差和标准误差。

**具体修改**：

1. **新增数据结构**：在 `Simulation_Data` 中，为每个 bin 增加一个 `_sq` 平方和数组：

```cpp
// 在 Data_Generation.hpp 中新增
std::array<double, NUM_BINS> captured_dt_sq_hist;      // Σ (Δt_i)²  per bin per trajectory
std::array<double, NUM_BINS> captured_v2dt_sq_hist;    // Σ (v²Δt_i)² per bin per trajectory
std::array<double, NUM_BINS> not_captured_dt_sq_hist;
std::array<double, NUM_BINS> not_captured_v2dt_sq_hist;
```

2. **累积逻辑**：每条轨迹完成后，其单条轨迹的 `dt_hist[b]` 记为 $x_b$：
   - `captured_dt_hist[b] += x_b`（已有）
   - `captured_dt_sq_hist[b] += x_b * x_b`（新增）

3. **误差计算**：设共 $N$ 条被捕获轨迹，某 bin 的总和 $S = \sum x_i$，平方和 $Q = \sum x_i^2$：
   - 均值：$\bar{x} = S / N$
   - 方差：$\text{Var}(x) = Q/N - \bar{x}^2$
   - **总和的标准误差**：$\sigma_S = \sqrt{N \cdot \text{Var}(x)} = \sqrt{Q - S^2/N}$（如果想得到总和的误差）
   - 或者用 $\sigma_{\bar{x}} = \sqrt{\text{Var}(x)/N}$ 表示均值的标准误差

4. **MPI 归约**：平方和数组同样用 `MPI_SUM` 归约

5. **输出格式**扩展为：
```
# bin_index  Sigma_dt  Sigma_v2dt  err_dt  err_v2dt
```

#### 方案 B：Bootstrap 重采样（后处理方式）

保存每条轨迹的单独 bincount，模拟结束后通过 bootstrap 重采样得到置信区间。优点是精确，缺点是内存消耗大（$N \times 2000 \times 2$ 个 double）。

**不推荐**：对于 $N > 10000$ 条轨迹，内存约 300+ MB，且需要额外的 MPI 通信。

#### 捕获率误差

捕获率 $p = N_{\text{cap}} / N_{\text{total}}$ 服从二项分布，标准误差：
$$\sigma_p = \sqrt{\frac{p(1-p)}{N_{\text{total}}}}$$

当 $N_{\text{cap}}$ 较小时，使用 Wilson 区间或 Clopper-Pearson 精确区间更佳。

**实现**：在 `Write_Output_Files()` 和 `Print_Summary()` 中输出：
```
# capture_rate = 0.05230000
# capture_rate_err = 0.00223456  (binomial standard error)
# capture_rate_CI_lower = 0.04812  (95% Clopper-Pearson)
# capture_rate_CI_upper = 0.05689
```

### 5.4 推荐实施步骤

1. 在 `Data_Generation.hpp` 新增 4 个 `std::array<double, NUM_BINS>` 平方和数组
2. 在 `Data_Generation.cpp` 的 `Generate_Data()` 循环中累积平方和
3. 在 `Perform_MPI_Reductions()` 中对平方和做 `MPI_SUM` 归约
4. 在 `Write_Output_Files()` 中计算误差并输出
5. 在 `Print_Summary()` 中输出捕获率误差

---

## 任务 6：代码审视与优化建议

### 6.1 Bug（确认的错误）

#### Bug 1 [严重]：`Compute_p_Value()` 构造函数参数错误
- **位置**：`Parameter_Scan.cpp` L339
- **问题**：`Simulation_Data data_set(sample_size, u_min)` 调用了新版 4 参数构造函数的前两个参数，但第二个参数语义已从 `u_min` 变为 `max_trajectories`。`u_min`（double 速度值）被隐式转换为 `unsigned int max_trajectories`，而实际的 u_min 和 iso_rings 使用了默认值。
- **影响**：Parameter scan 模式下，(1) max_trajectories 被设为错误值（可能为 0 → 无限制），(2) 反射最小速度阈值丢失（u_min=0），(3) iso_rings=1
- **修复**：`Simulation_Data data_set(sample_size, g_max_trajectories, u_min);`

#### Bug 2 [中等]：RK45 Fehlberg 系数错误
- **位置**：`Simulation_Trajectory.cpp` L438-440（以及 `main2.cpp` L136-138）
- **问题**：RK4 解的第四级系数 `2197.0 / 4101.0` 应为 `2197.0 / 4104.0`
- **参考**：标准 Fehlberg RK45 系数 $b_4 = 2197/4104$（Fehlberg 1969；同见 Numerical Recipes, Hairer et al.）
- **影响**：RK4 解的精度略有偏差（相对误差约 $7.3 \times 10^{-4}$），但由于自适应步长控制使用 RK4 与 RK5 的差值估计误差，系统性偏移会被步长调整部分补偿。对最终结果影响可能较小，但应修正。

#### Bug 3 [中等]：最后一步 bincount 未累积
- **位置**：`Simulation_Trajectory.cpp` `Propagate_Freely()` L207-210
- **问题**：注释说"Accumulate the last step"，但实际没有任何代码执行最后一步的累积。在线 bincount 使用前向差分 $\Delta t_i = t_{i+1} - t_i$，循环退出时最后一个点 `(prev_r_km, prev_v2_km2s2)` 没有被写入直方图（因为没有下一步来计算 dt）。
- **影响**：每次 `Propagate_Freely()` 调用损失一个步的 bincount 贡献。对于散射次数多的轨迹（调用多次 Propagate_Freely），累积遗漏可能不可忽略。
- **修复**：在循环退出后，使用上一个 dt 作为最后一步的权重：

```cpp
// After the while loop, accumulate the final point
if(prev_time_sec >= 0.0 && time_steps > 1)
{
    double dt_last = t_now_sec - prev_time_sec;  // reuse last available dt
    Accumulate_Bincount_Step(prev_r_km, prev_v2_km2s2, dt_last);
}
```

但需注意：`t_now_sec` 在循环退出后等于最后一步的时间（已存入 prev_time_sec），而 `dt_last` 应是前一个间隔的 dt。需要额外记录 `prev_dt` 变量。

#### Bug 4 [低]：Snapshot 文件竞态写入
- **位置**：`Simulation_Trajectory.cpp` L160-173
- **问题**：所有 MPI rank 的所有轨迹写入同一个 snapshot 文件路径（无 rank_id/trajectory_id 区分），后者覆盖前者。
- **影响**：Snapshot 功能实际上只保留最后一个跨越阈值的轨迹的快照，且多 rank 并发写入可能导致文件损坏。
- **修复**：文件名中加入 `rank_id` 和 `trajectory_id`，或将 snapshot 改为聚合模式。

#### Bug 5 [低]：散射概率使用错误的 dt
- **位置**：`Simulation_Trajectory.cpp` L185-195
- **问题**：`minus_log_xi -= particle_propagator.time_step * total_rate` 中，`time_step` 此时已被 RK45 自适应步长更新为**下一步**的预期步长，而非当前步实际消耗的时间。
- **影响**：由于自适应步长的连续性，相邻步长差异不大，但在步长剧烈变化时（如进入/离开太阳时）会引入误差。
- **修复**：在 RK45 步之前保存当前 time，之后用 `actual_dt = t_after - t_before` 计算散射概率。

### 6.2 代码问题（非 bug 但需改进）

#### 问题 1：RK45 递归可能栈溢出
- **位置**：`Simulation_Trajectory.cpp` `Runge_Kutta_45_Step()` L461
- **描述**：当误差超出容限时，缩小步长后递归调用自身。如果初始步长远大于所需步长，可能递归很深。
- **建议**：改为 `while` 循环。

#### 问题 2：废弃文件增加维护混淆
- `src/main2.cpp`：草稿代码，包含全局可变状态、硬编码常量、整数除法 bug
- `src/Data_Generation copy.cpp`：旧版本代码
- 这两个文件不参与编译，但在源码目录中会误导开发者。

#### 问题 3：`Propagate_Freely()` 中的 `v_after < 0.0` 检查
- **位置**：L188
- **问题**：速度是通过 `Current_Speed()` 返回的标量（总是 ≥0），在正常情况下不可能为负。只有当 `radius=0` 且 `v_radial` 为负时 `fabs(v_radial)` 仍返回正值。此检查本身无害，但如果真的触发说明有更深层的数值问题。

#### 问题 4：Snapshot 在模拟循环内打开/关闭文件
- 高频 I/O 可能成为性能瓶颈（虽然 snapshot 仅在时间阈值被跨越时触发，每条轨迹最多触发 `time_thresholds.size()` 次，实际影响有限）。

### 6.3 优化建议

| 优先级 | 建议 | 理由 |
|--------|------|------|
| 高 | 修复 `Compute_p_Value()` 构造函数调用 | Parameter scan 模式无法正确运行 |
| 高 | 修复 RK45 系数 `4101→4104` | 数值精度 |
| 高 | 实现最后一步 bincount 累积 | 数据完整性 |
| 中 | 添加每 bin 统计误差（方案 A） | 科学分析必需 |
| 中 | RK45 递归改循环 | 鲁棒性 |
| 中 | 修复散射概率 dt | 物理正确性 |
| 低 | 删除废弃文件 | 代码卫生 |
| 低 | 修复 Snapshot 文件命名 | 功能正确性 |
| 低 | 添加捕获率 Clopper-Pearson 置信区间 | 统计分析 |
