# 任务书 v2

## 前情提要

本项目基于 DaMaSCUS-SUN 数值模拟框架，用于模拟暗物质粒子在太阳引力场中的轨迹演化。原始代码的主要流程为：

1. **配置加载**：通过 libconfig++ 读取 `.cfg` 配置文件，获取暗物质质量、截面、样本数等参数。
2. **轨迹模拟**：使用 RK45 自适应步长积分器（`Free_Particle_Propagator`）演化粒子轨迹，记录 `(t, r, vx, vy, vz, E)` 六列数据并写入二进制 `.dat` 文件。
3. **反射谱统计**：原代码以**反射粒子**为核心统计量——`sample_size` 定义为反射数据点的目标数量，通过 KDE 平滑得到反射速度谱，计算暗物质微分通量。
4. **后处理**：Python 脚本 `bincount2.py` 读取轨迹文件，使用梯形法则时间加权统计径向分布直方图（r-bincount）和速度平方分布直方图（v²-bincount）。

**v1 版本审阅意见与决策记录：**

| 议题 | 审阅意见 | 决策 |
|------|---------|------|
| bincount 输出规格 | bin 定义需明确（独立/联合、范围、间距） | 两个独立 1D 直方图，r ∈ [0, 2R☉] 线性 2000 bins |
| 输出文件结构 | "每条轨迹输出"有歧义 | 输出三个汇总文件：捕获粒子 bincount、未捕获粒子 bincount、蒸发时间列表 |
| 流式 bincount 写入 | 梯形法则需前后时间步信息 | 采用相邻两步差值（前向差分）作为时间权重，轨迹结束时将最后一步直接累积 |
| 蒸发时间定义 | 粒子可能未在模拟时间内蒸发完成 | 增加布尔标志标记是否被截断（censored） |
| recording_step_override | 删除后需清理相关自适应校准逻辑 | 同时清理 perihelion detection、Get_Recording_Interval() 查找表等死代码 |
| 删除反射统计 | 涉及代码量巨大 | 暂时保留反射统计相关代码，留作后续任务处理 |
| sample_size 语义 | 全局精确计数需频繁 MPI 同步 | 每个 rank 的本地目标设为 ceil(sample_size / N_ranks)，各 rank 独立终止 |
| 低捕获率场景 | 达到 sample_size 可能需要极长时间 | 增加 `max_trajectories` 上限参数作为安全阀（见任务 8） |

---

## 任务内容

### 阶段一：代码备份与版本管理

**任务 1.** 完整阅览本项目中所有代码，充分理解代码逻辑与架构。

**任务 2.** 使用 git 将代码上传至 GitHub 仓库：
- 远程地址：`git@github.com:Funyday-k/DaMaSCUS-SUN-EVAP.git`
- 用户：Funyday-k
- 编写 `README.md`（项目简介、构建方法、使用说明）和 `.gitignore`（排除 build/、data/、results/、*.dat、*.o 等）

### 阶段二：核心模拟改动

以下所有修改必须在 git 备份完成后进行。

**任务 3. bincount 在线累积**

将原来的"写出完整轨迹文件 + Python 后处理"模式，改为在 C++ 模拟过程中**在线累积 bincount**。

- **输出两个独立的 1D 直方图**（与 `bincount2.py` 算法一致）：
  - **r-histogram**: $\sum \Delta t$ per bin，r ∈ [0, 2R☉]，2000 个线性等宽 bin，bin 宽 ≈ 695.7 km
  - **v²-histogram**: $\sum v^2 \Delta t$ per bin，同样以 r 为自变量分 bin（即统计每个径向 bin 内的 v² 加权时间）
  - 常量：$R_\odot = 6.957 \times 10^5$ km
- **时间加权方法**：使用相邻两步的时间差 $\Delta t_i = t_{i+1} - t_i$ 作为第 $i$ 步的时间权重（前向差分）。轨迹结束时，最后一步使用 $\Delta t_N = t_N - t_{N-1}$（即复用上一步的差分值）。
- **增量写入策略**：每条轨迹在内存中累积 bincount（2000 bins × 2 arrays × 8 bytes ≈ 32 KB），轨迹完成后一次性归类写入到对应的汇总数组中。无需中途写入文件——单条轨迹的 bincount 内存开销可忽略。
- **不再输出任何轨迹 `.dat` 文件**。

**任务 4. 能量检测与捕获率**

- 在每个 RK45 步中计算能量 $E = \frac{1}{2}m_\chi(v^2 - v_{\rm esc}^2(r))$（沿用现有逻辑）。
- 如果轨迹中**任意步**出现 $E \leq 0$，则标记该轨迹为 **captured**。
- 模拟结束后输出：
  - 总模拟粒子数
  - 被捕获粒子数与比例（capture rate）

**任务 5. 蒸发时间记录**

- 对每条被捕获的轨迹，记录：
  - $t_{\rm first}$：首次出现 $E \leq 0$ 的时间
  - $t_{\rm last}$：最后一次出现 $E \leq 0$ 的时间
  - 蒸发时间 $t_{\rm evap} = t_{\rm last} - t_{\rm first}$
  - **截断标志**（bool）：如果轨迹在模拟结束时仍处于 $E \leq 0$ 状态（即最后一步的 $E \leq 0$），标记为 `truncated = true`，表示蒸发时间为下界估计

**任务 6. 输出文件结构**

每个参数点的模拟在指定输出目录下生成以下三个汇总文件（由各 MPI rank 在模拟结束后归约输出）：

| 文件名 | 内容 | 格式 |
|--------|------|------|
| `captured_bincount.txt` | 所有被捕获轨迹的 r-histogram 和 v²-histogram 汇总 | 文本，2000 行，每行：`bin_index  Σdt  Σ(v²·dt)` |
| `not_captured_bincount.txt` | 所有未被捕获轨迹的 r-histogram 和 v²-histogram 汇总 | 同上 |
| `evaporation_summary.txt` | 所有被捕获轨迹的蒸发时间列表 | 文本，每行：`trajectory_id  t_evap[s]  truncated(0/1)` |

文件头部应包含元信息：暗物质质量、截面、总模拟数、捕获数、捕获率。

**任务 7. 删除 recording_step_override 及相关逻辑**

由于 bincount 在线累积需要统计每个 RK45 步，`recording_step_override` 变量及其相关逻辑不再有意义。需要删除：
- 配置文件中的 `recording_step_override` 参数
- `g_recording_step_override` 全局变量及其读取代码
- perihelion detection 自适应校准机制（`step_interval_calibrated`、`perihelion_count` 等）
- `Get_Recording_Interval()` 查找表函数
- 所有基于 `recording_step_interval` 的"每 N 步记录一次"的分支逻辑

**任务 8. sample_size 语义变更与并行终止策略**

- `sample_size` 的新语义：**目标被捕获粒子总数**（而非反射数据点数量）。
- **并行终止策略**：每个 MPI rank 的本地目标设为 $\lceil \texttt{sample\_size} / N_{\rm ranks} \rceil$。各 rank 达到本地目标后独立终止，无需频繁全局同步。最终实际捕获数可能略大于 `sample_size`（至多多出 $N_{\rm ranks} - 1$ 个）。
- **安全阀机制**：在配置文件中新增 `max_trajectories` 参数（默认值：`sample_size × 1000`）。当单个 rank 模拟的轨迹总数达到 $\lceil \texttt{max\_trajectories} / N_{\rm ranks} \rceil$ 时，即使未达到捕获目标也强制终止。终止时在输出文件头部标注 `EARLY_STOP: max_trajectories reached`，并输出当前已获得的统计量。
  - **使用场景**：当捕获率极低（<0.1%）时，达到 `sample_size` 可能需要模拟数百万条轨迹。`max_trajectories` 保证程序不会无限运行。用户可根据预估的捕获率调整此参数：若已知捕获率约为 $p$，建议设置 `max_trajectories ≈ sample_size / p × safety_factor`（safety_factor ≈ 2-5）。

**任务 9. 模拟终端输出**

模拟完成时在标准输出（stdout）打印以下统计摘要：
- 总模拟粒子数
- 被捕获粒子数与比例
- 捕获率（capture rate）
- 蒸发时间中位数（仅限未截断的轨迹）
- 是否触发 `max_trajectories` 安全阀

### 阶段三：测试与验证

**任务 10. 鲁棒性测试接口**

编写一个独立的测试模式（可通过配置文件或命令行参数激活），用于测试长时间模拟粒子的 bincount 分布收敛性：

- 在模拟物理时间达到指定阈值时，输出当前 bincount 的快照（snapshot）。
- 默认时间阈值序列：1000 s、10000 s、100000 s（可在配置中自定义）。
- 每个快照为独立文件，命名格式：`snapshot_{time}s.txt`，内容格式同 `captured_bincount.txt`。
- 快照是当前时刻 bincount 的**拷贝**，不影响主 bincount 的继续累积。
- 此模式用于验证：对于那些模拟时间极长的粒子，其 bincount 分布是否随时间收敛。

---

## 工作要求

1. 一切工作过程要保证上下文的限制，在限制内保证任务完成的正确性和效率。
2. 可以使用 Subagent 来提升工作的效率与正确性。
3. 每完成一个任务阶段，进行一次 git commit，commit message 清晰描述改动内容。
4. 暂时保留反射统计相关代码（`Reflection_Spectrum` 等），留作后续版本处理。
5. 修改代码前应先用现有 `bincount2.py` 的输出作为基准，确保 C++ 在线累积的 bincount 与 Python 后处理结果数值一致（允许浮点误差 < 1e-6 相对误差）。
