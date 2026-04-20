# 任务书 v4 — 计算时间 & RK45 步数统计 + 步级切片 Bincount

## 概述

本任务分两阶段实施：
1. **阶段一（功能 B）**：记录每条轨迹的 wall-clock 计算时间 & RK45 步数，输出对数直方图和分类汇总。完成后用户根据步数统计数据决定切片参数。
2. **阶段二（功能 A，待定）**：基于阶段一的数据，确定合理的步数切片间隔后实施步级切片 Bincount。

---

## 功能 B：计算时间 & RK45 步数统计（优先实施）

### B0. Trajectory_Result 新增 RK45 步数字段

**文件**: `include/Simulation_Trajectory.hpp`, `src/Simulation_Trajectory.cpp`

`Trajectory_Result` 新增：
```cpp
unsigned long int total_rk45_steps;  // 该轨迹的总 RK45 步数
```

`Trajectory_Simulator` 新增私有成员：
```cpp
unsigned long int total_rk45_steps_current_traj;  // 累计当前轨迹所有 Propagate_Freely 调用的步数
```

- 在 `Simulate()` 开始时 reset 为 0
- 在每次 `Propagate_Freely()` 结束时累加 `time_steps`
- `Trajectory_Result` 构造函数接收并存储

### B1. 数据采集

**文件**: `src/Data_Generation.cpp` (`Generate_Data` 函数)

每条轨迹前后用 `std::chrono::high_resolution_clock` 计时：
```cpp
auto t0 = std::chrono::high_resolution_clock::now();
Trajectory_Result trajectory = simulator.Simulate(IC, DM, mpi_rank);
auto t1 = std::chrono::high_resolution_clock::now();
double wall_sec = 1e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
```

### B2. 新增数据结构

**文件**: `include/Data_Generation.hpp`

```cpp
// --- 计算时间统计 ---
double total_wall_time_captured = 0.0;
double total_wall_time_not_captured = 0.0;

// Wall-clock 时间对数直方图：10^{-4} ~ 10^{8} 秒，每十倍频程 10 bin，共 120 bin
// 使用 overflow/underflow 计数处理极端值
static constexpr int WALL_TIME_BINS = 120;
static constexpr double WALL_TIME_LOG_MIN = -4.0;   // 0.1 ms
static constexpr double WALL_TIME_LOG_MAX =  8.0;   // ~3 年
std::array<unsigned long int, WALL_TIME_BINS> wall_time_hist_captured;
std::array<unsigned long int, WALL_TIME_BINS> wall_time_hist_not_captured;
unsigned long int wall_time_overflow_captured = 0;    // > 10^8 s
unsigned long int wall_time_overflow_not_captured = 0;

// --- RK45 步数统计 ---
unsigned long int total_rk45_steps_captured = 0;
unsigned long int total_rk45_steps_not_captured = 0;

// RK45 步数对数直方图：10^{0} ~ 10^{14} 步，每十倍频程 10 bin，共 140 bin
static constexpr int STEP_COUNT_BINS = 140;
static constexpr double STEP_COUNT_LOG_MIN = 0.0;   // 1 步
static constexpr double STEP_COUNT_LOG_MAX = 14.0;   // 10^14 步 (匹配 maximum_free_time_steps = 1e12 × 多次散射)
std::array<unsigned long int, STEP_COUNT_BINS> step_count_hist_captured;
std::array<unsigned long int, STEP_COUNT_BINS> step_count_hist_not_captured;
unsigned long int step_count_overflow_captured = 0;
unsigned long int step_count_overflow_not_captured = 0;
```

### B3. 累积逻辑

在主循环中，根据 captured / not_captured 分别累积 wall-clock 时间和步数：
```cpp
double wall_sec = ...;
unsigned long int steps = trajectory.total_rk45_steps;

auto bin_wall = [](double t) -> int {
    if(t <= 0.0) return -1;
    return (int)((log10(t) - WALL_TIME_LOG_MIN) * 10.0);
};
auto bin_step = [](unsigned long int s) -> int {
    if(s == 0) return -1;
    return (int)((log10((double)s) - STEP_COUNT_LOG_MIN) * 10.0);
};

if(trajectory.bincount.is_captured) {
    total_wall_time_captured += wall_sec;
    total_rk45_steps_captured += steps;
    int wb = bin_wall(wall_sec);
    if(wb >= 0 && wb < WALL_TIME_BINS) wall_time_hist_captured[wb]++;
    else if(wb >= WALL_TIME_BINS) wall_time_overflow_captured++;
    int sb = bin_step(steps);
    if(sb >= 0 && sb < STEP_COUNT_BINS) step_count_hist_captured[sb]++;
    else if(sb >= STEP_COUNT_BINS) step_count_overflow_captured++;
} else {
    // 同理 not_captured
}
```

### B4. MPI 归约

**文件**: `src/Data_Generation.cpp` (`Perform_MPI_Reductions`)

```cpp
MPI_Allreduce(MPI_IN_PLACE, &total_wall_time_captured, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, &total_wall_time_not_captured, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, &total_rk45_steps_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, &total_rk45_steps_not_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, wall_time_hist_captured.data(), WALL_TIME_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, wall_time_hist_not_captured.data(), WALL_TIME_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, &wall_time_overflow_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, &wall_time_overflow_not_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, step_count_hist_captured.data(), STEP_COUNT_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, step_count_hist_not_captured.data(), STEP_COUNT_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, &step_count_overflow_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(MPI_IN_PLACE, &step_count_overflow_not_captured, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
```

### B5. 输出文件

**文件**: `src/Data_Generation.cpp` (`Write_Output_Files`)

新增文件 `computation_time_summary.txt`：
```
# DM_mass_GeV = ...
# DM_sigma_cm2 = ...
# total_wall_time_captured_sec = 12345.678
# total_wall_time_not_captured_sec = 67890.123
# mean_wall_time_captured_sec = 1.234
# mean_wall_time_not_captured_sec = 0.045
# total_rk45_steps_captured = 123456789
# total_rk45_steps_not_captured = 987654321
# mean_rk45_steps_captured = 12345.6
# mean_rk45_steps_not_captured = 678.9
#
# [Wall-clock time histogram]
# bin_index  log10_t_lower  log10_t_upper  count_captured  count_not_captured
0   -4.0   -3.9   0   5
1   -3.9   -3.8   2   12
...
119  7.9   8.0   0   0
# overflow_captured = 0
# overflow_not_captured = 0
#
# [RK45 step count histogram]
# bin_index  log10_steps_lower  log10_steps_upper  count_captured  count_not_captured
0   0.0   0.1   0   100
1   0.1   0.2   0   50
...
139  13.9  14.0  0   0
# overflow_captured = 0
# overflow_not_captured = 0
```

### B6. Print_Summary 扩展

**文件**: `src/Data_Generation.cpp` (`Print_Summary`)

新增输出行：
```
Captured wall-clock time:          3h 25m 45s
Not-captured wall-clock time:      18h 51m 30s
Mean captured trajectory time:     1.23 s
Mean not-captured trajectory time: 0.045 s
Total RK45 steps (captured):       123,456,789
Total RK45 steps (not captured):   987,654,321
Mean RK45 steps/traj (captured):   12345
Mean RK45 steps/traj (not captured): 678
```

---

## 功能 A：步级切片 Bincount（待阶段一完成后决定）

功能 A 的方案设计见 Mission_3 讨论。实施参数（`step_snapshot_interval`）需根据功能 B 的输出数据确定。

阶段一完成后，用户将基于以下数据决策：
- 典型捕获轨迹的 RK45 步数分布
- 步数与蒸发时间的对应关系（可从 evaporation_summary.txt 和步数数据交叉分析）
- 步数与 wall-clock 时间的对应关系

---

## 修改文件清单（阶段一）

| 文件 | 修改内容 |
|------|----------|
| `include/Simulation_Trajectory.hpp` | `Trajectory_Result` 加 `total_rk45_steps`；`Trajectory_Simulator` 加步数累计变量 |
| `src/Simulation_Trajectory.cpp` | `Propagate_Freely` 累计步数；`Simulate` 重置和传递步数；构造函数扩展 |
| `include/Data_Generation.hpp` | wall-clock 时间 & RK45 步数统计成员变量 |
| `src/Data_Generation.cpp` | 计时逻辑；步数累积；MPI 归约；输出文件；Print_Summary |

## 实施顺序

1. 修改头文件（数据结构）→ git commit
2. 实现 Simulation_Trajectory 步数统计 → git commit
3. 实现 Data_Generation 计时 + 步数 + MPI + 输出 → git commit
4. sandbox 编译验证 → git commit (fix if needed)
5. **用户决策点**：根据数据决定是否以及如何实施功能 A

---

## 工作要求

1. 一切工作过程要保证上下文的限制，在限制内保证任务完成的正确性和效率。
2. 可以使用 Subagent 来提升工作的效率与正确性。
3. 每完成一个任务阶段，进行一次 git commit，commit message 清晰描述改动内容。
4. 代码修改完成后，必须在 sandbox 上编译验证通过。

---

**请审阅本任务书，确认后开始执行代码修改。**
