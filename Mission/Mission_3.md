# 任务书 v3

## 前情提要

任务书 v2 完成了核心模拟改动：在线 bincount 累积、捕获率检测、蒸发时间记录、输出文件结构、sample_size 语义变更、max_trajectories 安全阀、snapshot 测试接口。

在任务书 v2 执行报告（`Mission_2_Report.md`）中对代码进行了全面审视，发现以下问题：

**确认的 Bug：**

| 编号 | 严重度 | 位置 | 描述 |
|------|--------|------|------|
| B1 | 严重 | `Parameter_Scan.cpp` L339 | `Compute_p_Value()` 中 `Simulation_Data(sample_size, u_min)` 调用参数错位，`u_min`(double) 被隐式转换为 `max_trajectories`(unsigned int)，实际的最小速度阈值丢失 |
| B2 | 中等 | `Simulation_Trajectory.cpp` L438-440 | RK45 Fehlberg RK4 解系数 `2197.0/4101.0` 应为 `2197.0/4104.0` |
| B3 | 中等 | `Simulation_Trajectory.cpp` `Propagate_Freely()` L207 | 最后一步 bincount 未累积：注释说"accumulate the last step"但无实际代码 |
| B4 | 低 | `Simulation_Trajectory.cpp` L160-173 | Snapshot 文件无 rank_id/trajectory_id，多轨迹覆盖、多 rank 竞态写入 |
| B5 | 低 | `Simulation_Trajectory.cpp` L185-195 | 散射概率使用 RK45 更新后的 time_step（下一步预期值），而非当前步实际消耗的 dt |

**代码质量问题：**

| 编号 | 位置 | 描述 |
|------|------|------|
| Q1 | `Simulation_Trajectory.cpp` RK45 | `Runge_Kutta_45_Step()` 误差过大时递归调用自身，极端情况栈溢出 |
| Q2 | `src/main2.cpp`, `Data_Generation copy.cpp` | 废弃文件不参与编译但在源码目录中增加混淆 |

**功能缺失：**

| 编号 | 描述 |
|------|------|
| F1 | bincount 输出文件没有每 bin 的统计误差 |
| F2 | 捕获率没有置信区间 |

---

## 任务内容

### 阶段一：Bug 修复

以下修复必须在当前代码基础上进行，修复后进行编译验证。

**任务 1. 修复 `Compute_p_Value()` 构造函数调用 [B1]**

在 `Parameter_Scan.cpp` 的 `Compute_p_Value()` 函数中，将：
```cpp
Simulation_Data data_set(sample_size, u_min);
```
修改为：
```cpp
Simulation_Data data_set(sample_size, g_max_trajectories, u_min);
```

这确保：
- `max_trajectories` 使用配置文件中的全局值
- `u_min` 正确传递为反射最小速度阈值
- `iso_rings` 使用默认值 1

**验证**：修改后能正常编译，且 Parameter scan 模式下的 p 值计算结果合理。

**任务 2. 修复 RK45 Fehlberg 系数 [B2]**

在 `Simulation_Trajectory.cpp` 的 `Runge_Kutta_45_Step()` 中，将 RK4 解的三行系数：
```cpp
double radius_4   = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4101.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4101.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
double phi_4      = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4101.0 * k_p[3] - 1.0 / 5.0 * k_p[4];
```
中的 `4101.0` 全部替换为 `4104.0`。

**参考**：标准 Fehlberg RK45 系数 $b_4 = 2197/4104$（E. Fehlberg, NASA TR R-315, 1969）。

**任务 3. 修复最后一步 bincount 遗漏 [B3]**

在 `Simulation_Trajectory.cpp` 的 `Propagate_Freely()` 函数中，需要在循环退出后补充最后一个点的 bincount 累积。

**实现方案**：
1. 在 `Trajectory_Simulator` 类中增加一个 `double prev_dt_sec` 成员变量，记录上一步的 $\Delta t$
2. 在循环内每次计算 `dt_sec` 后同步更新 `prev_dt_sec = dt_sec`
3. 在循环退出后：
```cpp
// After the while loop, accumulate the final point using the last known dt
if(prev_dt_sec > 0.0)
{
    Accumulate_Bincount_Step(prev_r_km, prev_v2_km2s2, prev_dt_sec);
}
```
4. 在 `Simulate()` 函数中初始化 `prev_dt_sec = 0.0`

**任务 4. 修复散射概率使用错误的 dt [B5]**

在 `Propagate_Freely()` 中，散射概率计算应使用当前步的实际时间消耗，而非 RK45 自适应更新后的下一步预期步长。

**实现方案**：在 RK45 步之前记录时间，之后计算实际 dt：
```cpp
double t_before = particle_propagator.Current_Time();
particle_propagator.Runge_Kutta_45_Step(solar_model.Mass(r_before));
double actual_dt = particle_propagator.Current_Time() - t_before;
// ...
// 散射部分使用 actual_dt 而非 particle_propagator.time_step：
minus_log_xi -= actual_dt * total_rate;
// 步长限制仍然使用 time_step（影响下一步）
```

**注意**：步长上限约束 `if(particle_propagator.time_step > time_step_max)` 仍应作用于 `time_step`（影响下一步的大小），但概率计算使用 `actual_dt`。

### 阶段二：RK45 递归改循环 [Q1]

**任务 5. RK45 递归改循环**

将 `Runge_Kutta_45_Step()` 中的递归步长重试改为 `do-while` 循环，避免极端情况下的栈溢出：

```cpp
void Free_Particle_Propagator::Runge_Kutta_45_Step(double mass)
{
    bool accepted = false;
    while(!accepted)
    {
        // ... 计算 k_r, k_v, k_p ...
        // ... 计算 RK4 和 RK5 的结果 ...
        // ... 计算误差和新步长 ...

        if(errors[0] < error_tolerances[0] && errors[1] < error_tolerances[1] && errors[2] < error_tolerances[2])
        {
            time     += time_step;
            radius    = radius_4;
            if(radius < 0.0) radius = 0.0;
            v_radial  = v_radial_4;
            phi       = phi_4;
            time_step = time_step_new;
            accepted  = true;
        }
        else
        {
            time_step = time_step_new;
            // Loop continues with smaller time_step
        }
    }
}
```

### 阶段三：误差分析实现 [F1, F2]

**任务 6. 在线方差累积——数据结构扩展**

在 `Data_Generation.hpp` 的 `Simulation_Data` 类中新增：
```cpp
// Per-bin sum of squares for error estimation
std::array<double, NUM_BINS> captured_dt_sq_hist;      // Σ (per-traj dt)²
std::array<double, NUM_BINS> captured_v2dt_sq_hist;    // Σ (per-traj v²dt)²
std::array<double, NUM_BINS> not_captured_dt_sq_hist;
std::array<double, NUM_BINS> not_captured_v2dt_sq_hist;
```

在构造函数中初始化为 0。

**任务 7. 在线方差累积——累积逻辑**

在 `Data_Generation.cpp` 的 `Generate_Data()` 中，轨迹完成后累积平方和：
```cpp
if(trajectory.bincount.is_captured)
{
    for(int b = 0; b < NUM_BINS; b++)
    {
        captured_dt_hist[b]    += trajectory.bincount.dt_hist[b];
        captured_v2dt_hist[b]  += trajectory.bincount.v2dt_hist[b];
        // 新增：平方和
        captured_dt_sq_hist[b]   += trajectory.bincount.dt_hist[b] * trajectory.bincount.dt_hist[b];
        captured_v2dt_sq_hist[b] += trajectory.bincount.v2dt_hist[b] * trajectory.bincount.v2dt_hist[b];
    }
}
// not_captured 部分同理
```

**任务 8. 在线方差累积——MPI 归约**

在 `Perform_MPI_Reductions()` 中增加平方和数组的归约：
```cpp
MPI_Allreduce(MPI_IN_PLACE, captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, not_captured_dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, not_captured_v2dt_sq_hist.data(), NUM_BINS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
```

**任务 9. 输出文件格式扩展**

修改 `Write_Output_Files()` 的 bincount 输出格式为 5 列：
```
# bin_index  Sigma_dt[s]  Sigma_v2dt[km2/s]  err_dt[s]  err_v2dt[km2/s]
```

其中误差列的计算（以 captured_dt 为例）：
```cpp
double N = number_of_captured_particles;
if(N > 1)
{
    double mean = captured_dt_hist[b] / N;
    double var  = captured_dt_sq_hist[b] / N - mean * mean;
    if(var < 0.0) var = 0.0;  // 防止浮点误差
    double err_sum = sqrt(N * var);  // 总和的标准误差
    // 输出 err_sum
}
```

**任务 10. 捕获率误差输出**

在 `Write_Output_Files()` 和 `Print_Summary()` 中，增加捕获率的统计误差：

1. **二项标准误差**：$\sigma_p = \sqrt{p(1-p)/N}$
2. **Clopper-Pearson 95% 置信区间**：使用 beta 分布分位点
   - 下界：$B(\alpha/2; k, n-k+1)$
   - 上界：$B(1-\alpha/2; k+1, n-k)$
   - 其中 $k = N_{\text{cap}}, n = N_{\text{total}}, \alpha = 0.05$

如果项目中没有 beta 分布的逆 CDF 实现，可使用 Wilson 区间近似：
$$p_{\pm} = \frac{p + z^2/(2n) \pm z\sqrt{p(1-p)/n + z^2/(4n^2)}}{1 + z^2/n}$$
其中 $z = 1.96$（95% CL）。

输出格式（在文件头部）：
```
# capture_rate = 0.05230000
# capture_rate_err = 0.00223456
# capture_rate_CI_95_lower = 0.04812345
# capture_rate_CI_95_upper = 0.05689012
```

### 阶段四：Snapshot 修复与清理

**任务 11. 修复 Snapshot 文件命名 [B4]**

修改 snapshot 文件命名方案，包含 MPI rank 和轨迹编号，避免覆盖和竞态：
```
snapshot_{time}s_rank{rank_id}_traj{trajectory_id}.txt
```

需要：
1. 在 `Trajectory_Simulator` 中传递当前的 `mpi_rank` 和轨迹计数
2. 修改 `Propagate_Freely()` 中的文件名生成逻辑

**任务 12. 清理废弃文件 [Q2]**

删除不参与编译的废弃文件：
- `src/main2.cpp`
- `src/Data_Generation copy.cpp`

或移入一个 `archive/` 目录以保留历史参考。

---

## 工作要求

1. 一切工作过程要保证上下文的限制，在限制内保证任务完成的正确性和效率。
2. 可以使用 Subagent 来提升工作的效率与正确性。
3. 每完成一个任务阶段，进行一次 git commit，commit message 清晰描述改动内容。
4. Bug 修复（阶段一）完成后，必须在 sandbox 上编译验证通过。
5. 误差分析实现（阶段三）完成后，运行一个小规模测试（sample_size=5），确认输出文件中误差列的值合理（非零、非 NaN、误差 < 总和）。
6. 暂时保留反射统计相关代码（`Reflection_Spectrum` 等），留作后续版本处理。
