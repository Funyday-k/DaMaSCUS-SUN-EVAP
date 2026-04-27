# DaMaSCUS-SUN 项目完整技术文档

## 目录

1. [项目概览与物理背景](#1-项目概览与物理背景)
2. [顶层设计与基本物理思想](#2-顶层设计与基本物理思想)
3. [核心数值模拟方法](#3-核心数值模拟方法)
4. [代码架构与数据管线](#4-代码架构与数据管线)
5. [模拟阶段的详细实现](#5-模拟阶段的详细实现)
6. [后处理数据管线](#6-后处理数据管线)
7. [计算效率优化与物理影响](#7-计算效率优化与物理影响)
8. [与文献的对照分析](#8-与文献的对照分析)
9. [总结](#9-总结)

---

## 1. 项目概览与物理背景

### 1.1 项目名称与定位

**DaMaSCUS-SUN**（**Da**rk **Ma**tter **S**imulation **C**ode for **U**nderground **S**cattering — **S**olar）是一个MPI并行化的暗物质蒙特卡洛模拟框架，专门用于模拟轻暗物质（$m_\chi \lesssim$ 几GeV）在太阳内部的运动、散射、捕获和蒸发过程。

### 1.2 物理动机

在宇宙中，暗物质占总物质-能量密度的约27%。直接探测实验试图通过暗物质与探测器中原子核的弹性散射来发现暗物质。然而，对于质量较轻的暗物质粒子（$m_\chi \lesssim 1 \text{ GeV}$），传统的地下实验灵敏度急剧下降，因为核反冲能量 $E_R \propto m_\chi^2$ 低于探测器阈值。

**太阳作为天然暗物质探测器**提供了一个突破方案：

1. **引力聚焦（Gravitational Focusing）**：太阳的强引力场将银河系晕中的暗物质粒子吸引并加速，增大了散射截面的有效面积。进入太阳的暗物质通量增强因子为：

$$\Gamma_\text{entering} = \pi R_\odot^2 \, n_\chi \langle v_\text{eff} \rangle$$

其中 $n_\chi = \rho_\chi / m_\chi$ 是局域暗物质数密度，$\langle v_\text{eff} \rangle$ 包含晕速度分布和太阳表面逃逸速度 $v_\text{esc} \approx 618 \text{ km/s}$ 的修正。

2. **暗物质捕获（Capture）**：暗物质粒子在太阳内部与核子（主要是氢和氦）或电子发生弹性散射，损失动能。如果散射后速度低于局域逃逸速度 $v < v_\text{esc}(r)$（即总能量变为负），粒子被引力束缚——即**被捕获**。

3. **暗物质蒸发（Evaporation）**：被捕获的轻质暗物质粒子在太阳热浴中经历多次散射，有概率获得足够能量从太阳逃逸。蒸发率 $E_\odot(m_\chi, \sigma)$ 决定了太阳中暗物质的平衡数量。存在一个**蒸发质量** $m_\text{evap}$，低于该质量的暗物质粒子蒸发时标短于太阳年龄（$\sim 4.6 \times 10^9$ 年），无法在太阳中积累。

4. **太阳反射暗物质（Solar Reflected Dark Matter, SRDM）**：部分暗物质粒子在太阳中散射后获得额外动能，以高于初始银晕速度的速度逃逸太阳。这些"被太阳加速"的暗物质粒子到达地球后，其速度谱不同于银晕背景分布，可作为新的探测信号源。参考文献 Emken (2022) 和 Nguyen & Linden (2026) 详细讨论了这一机制。

### 1.3 与文献的关联

本项目的物理基础和数值方法密切关联以下参考文献：

| 文献 | 核心内容 | 在代码中的对应 |
|------|---------|---------------|
| **Garani & Palomares-Ruiz (2025, Erratum)** | 暗物质在太阳中的散射：电子 vs 核子，蒸发率计算公式 | `evaporation_theory.py` 直接实现其Eq. (5.1), (A.23), (B.11), (C.8)-(C.13)，主要实现SD相互作用 |
| **Emken (2022)** | 轻暗物质的太阳反射（重介质子情形），反射谱计算 | `Reflection_Spectrum.cpp`, `Dark_Photon.cpp` 实现暗光子形因子和反射通量 |
| **Busoni et al. (2013)** | 太阳中微子约束的最小暗物质质量——蒸发边界 | 蒸发质量概念，`Parameter_Scan.cpp` 中排斥曲线计算 |
| **Nguyen & Linden (2026)** | 太阳约束自旋相关暗物质-核子散射（蒸发极限以下） | `config_Lingyu.cfg` 中SD相互作用配置，自旋相关截面处理 |

---

## 2. 顶层设计与基本物理思想

### 2.1 基本物理方法：蒙特卡洛轨迹模拟

DaMaSCUS-SUN 的核心方法是**直接蒙特卡洛轨迹模拟**（Direct Monte Carlo Trajectory Simulation, MCTS）：不使用解析近似，而是逐一模拟大量暗物质粒子在太阳引力场和多组分等离子体中的完整运动-散射历史。

$$\text{单粒子轨迹} = \sum_{i=0}^{N_\text{scat}} \left[ \text{自由传播}_i \to \text{弹性散射}_i \right]$$

其中每一步的自由传播遵循Kepler轨道力学（引力 + 太阳质量分布），散射事件为与太阳物质的弹性碰撞。

### 2.2 基本物理假设

1. **弹性散射**：暗物质与靶粒子（核子或电子）的相互作用为弹性散射，通过微分截面 $d\sigma/dq^2$ 参数化。散射角分布由暗物质的形因子（如暗光子模型的 $q$-依赖形因子）和核的Helm形因子共同决定。

2. **太阳标准模型**：太阳内部的密度、温度、化学组成等物理量采用 **AGSS09 标准太阳模型**（Asplund, Grevesse, Sauval & Scott 2009），提供从太阳中心到表面的径向剖面数据。模型包含29种同位素的丰度分布。

3. **暗物质晕模型**：入射暗物质的速度分布取自银河系标准晕模型（Standard Halo Model），通常为截断 Maxwell-Boltzmann 分布，参数为：
   - 局域暗物质密度 $\rho_\chi \approx 0.3 \text{ GeV/cm}^3$
   - 晕速度色散 $v_0 \approx 220 \text{ km/s}$
   - 银河逃逸速度 $v_\text{esc, gal} \approx 544 \text{ km/s}$

4. **引力传播**：太阳外部为纯Kepler（$1/r^2$）引力问题，内部需要考虑太阳质量分布 $M(r)$ 对引力场的修正。

### 2.3 模拟的三种运行模式

```
                    ┌──────────────────────────────┐
                    │     DaMaSCUS-SUN 运行模式     │
                    └──────────┬───────────────────┘
                               │
        ┌───────────────┬───────────────┬───────────────┐
        │               │               │               │
┌───────▼───────┐ ┌─────▼──────┐ ┌──────▼───────┐
│ Parameter     │ │ Capture    │ │ Parameter    │
│ Point         │ │ Mode       │ │ Scan         │
│ 单参数点完整模拟│ │ 快速捕获率  │ │ 参数空间扫描  │
└───────┬───────┘ └─────┬──────┘ └──────┬───────┘
        │               │               │
┌───────▼───────┐ ┌─────▼──────┐ ┌──────▼───────┐
│ 生成bincount   │ │ E<0即停止   │ │ STA/Full扫描 │
│ 蒸发时间统计   │ │ 只打印捕获率 │ │ 计算p值网格  │
│ 可生成反射谱   │ │ 不写输出文件 │ │ 提取排斥曲线 │
└───────────────┘ └────────────┘ └──────────────┘
```

- **"Parameter Point" 模式**：对单一 $(m_\chi, \sigma)$ 参数点进行详细模拟，生成大量轨迹，输出反射速度谱和探测器信号率。
- **Capture Mode 快速捕获率模式**：通过 `capture_mode = true` 或 `run_mode = "Capture"` 启用。轨迹一旦出现总能量 $E < 0$ 就立即终止并记为捕获；不写 `bincount`、蒸发、snapshot 或反射谱输出，只在终端打印捕获率和 95% Wilson 上下误差。
- **"Parameter Scan" 模式**：在二维 $(m_\chi, \sigma)$ 参数空间中扫描，对每个点计算统计检验的p值，最终提取一定置信水平（如90% CL）下的排斥极限曲线。

---

## 3. 核心数值模拟方法

### 3.1 自适应Runge-Kutta-Fehlberg (RK45) 积分器

暗物质粒子在太阳引力场中的自由传播通过求解以下 ODE 系统实现，采用二维极坐标（利用角动量守恒将三维问题化为二维）：

$$\frac{dr}{dt} = v_r, \quad \frac{dv_r}{dt} = \frac{J^2}{r^3} - \frac{G_N M(r)}{r^2}, \quad \frac{d\phi}{dt} = \frac{J}{r^2}$$

其中 $J = |\vec{r} \times \vec{v}|$ 为守恒角动量，$M(r)$ 为太阳在半径 $r$ 以内的包含质量（由AGSS09模型提供并内插）。

**RK45 自适应步长控制（当前代码：`Free_Particle_Propagator`）**：
- 初始步长为 `0.1 s`。
- 采用 Fehlberg 系数的6阶段嵌入式方法，同时给出4阶和5阶近似解。
- 误差估计：$\Delta_i = |y_4^{(i)} - y_5^{(i)}|$，对三个分量分别设置容差：

| 分量 | 容差 |
|------|------|
| 位置 $r$ | 1 km |
| 径向速度 $v_r$ | $10^{-3}$ km/s |
| 角度 $\phi$ | $10^{-7}$ rad |

- 步长调整公式：$\delta t_\text{new} = 0.84 \times (\text{tol}/\text{err})^{1/4} \times \delta t_\text{old}$，并限制每次变化因子在 `[0.1, 4.0]`。
- 若误差为 NaN/Inf，步长强制乘以 `0.1`，避免非有限误差导致步长发散。
- 若内层尝试超过 `2000` 次，或新步长小于 `1.0e-8 s`，代码会强制接受一个最小步长，确保单个 RK45 步不会无限卡死。
- RK45 步长有绝对上限 `1.0e6 s`，防止太阳外弱引力区步长增长到 `inf`。
- 在太阳外传播时，步长还会受单步跨越距离和局域 Kepler 动力学时间限制，避免一步跳过 `2R_\odot` 边界并产生非物理巨大半径。
- **太阳内部**：步长还受散射率约束，$\delta t \leq 0.1 / \Gamma_\text{total}(r, v)$。
- **太阳外部**：无散射事件，RK45 步长只由轨道误差控制；从 `1000 AU` 到 `2R_\odot` 的入射段，以及逃逸后到 `1 AU` 的传播，可用 `Hyperbolic_Kepler_Shift()` 做解析 Kepler 推进。
- **单轨迹 wall-time 保护**：`max_trajectory_wall_time_sec` 默认 `300 s`，每 `256` 个 RK45 步检查一次；超过后中止该轨迹，避免 MPI rank 被病态轨迹长期卡住。设为 `0` 表示不限制。

### 3.2 弹性散射的蒙特卡洛采样

每次散射事件的处理流程：

1. **靶粒子选择**：按各同位素的散射率权重随机选择靶核（或电子）：

$$P(\text{target}_i) = \frac{\Gamma_i(r, v)}{\Gamma_\text{total}(r, v)}$$

2. **热靶速度采样**：靶粒子速度从 Maxwell-Boltzmann 分布中抽样（温度 $T(r)$ 由太阳模型提供），采用 Romano & Walsh (2011) 剔除法（rejection sampling）：

$$f(\vec{v}_\text{target}) \propto v_\text{target}^2 \exp\left(-\frac{m_\text{target} v_\text{target}^2}{2 T(r)}\right)$$

参数 $\kappa = \sqrt{m_\text{target} / 2T}$ 控制分布宽度。

3. **散射角采样**：从微分截面导出的PDF中抽样散射角 $\cos\alpha$。对于暗光子模型：

$$P(\cos\alpha) \propto \frac{m_{A'}^2 (m_{A'}^2 + q_\text{max}^2)}{(2m_{A'}^2 + q_\text{max}^2(1 - \cos\alpha))^2}$$

方位角 $\phi$ 在 $[0, 2\pi)$ 均匀分布。

4. **速度更新**：在质心系中计算散射后的暗物质速度，然后变换回太阳参考系。弹性碰撞运动学保证能量和动量守恒。

### 3.3 散射事件的判定：光学深度方法

暗物质粒子在自由传播过程中是否发生散射，通过对随机光学深度的监控实现：

$$\xi = -\frac{\ln(u_0)}{\Gamma_\text{total}(r, v)}$$

其中 $u_0 \sim \text{Uniform}(0, 1)$ 为随机数。每个积分步中累积光学深度，当 $\xi < 0$ 时触发一次散射事件。

### 3.4 初始条件生成

暗物质粒子的初始条件从银河系晕分布中抽样：

1. **渐近速度 $u$**：从修正后的晕速度PDF中采样：$P(u) \propto f_\text{halo}(v) \times (u^2 + v_\text{esc}^2)/u$
2. **入射角 $\cos\theta$**：给定 $u$，从角度条件PDF中采样，受能量和角动量守恒约束
3. **碰撞参数 $b$**：$b \leq b_\text{max} = R_\odot \sqrt{1 + v_\text{esc}^2/u^2}$（引力聚焦效应）
4. **Kepler解析推进**：将粒子从渐近距离（1000 AU）沿双曲线轨道解析传播至模拟起始半径（$2R_\odot$），利用双曲Kepler轨道元素：
   - 半长轴 $a = -G_N M_\odot / u^2$
   - 半通径 $p = J^2 / (G_N M_\odot)$
   - 离心率 $e = \sqrt{1 + p/a}$

---

## 4. 代码架构与数据管线

### 4.1 系统架构总览

```
┌────────────────────────────────────────────────────────────────────┐
│                    DaMaSCUS-SUN 完整数据管线                       │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  ┌─────────────────────────────────────────┐                      │
│  │          配置层 (Configuration)          │                      │
│  │  config_Lingyu.cfg → Parameter_Scan.cpp │                      │
│  │  - DM质量, 截面, 相互作用类型            │                      │
│  │  - 样本数, 运行模式                      │                      │
│  │  - 暗光子参数 (ε, α_D, m_A')            │                      │
│  └──────────────────┬──────────────────────┘                      │
│                     │                                              │
│  ┌──────────────────▼──────────────────────┐                      │
│  │          太阳模型层 (Solar_Model)        │                      │
│  │  model_agss09.dat → Solar_Model.cpp     │                      │
│  │  - AGSS09 密度/温度/元素丰度剖面         │                      │
│  │  - 29种同位素散射率                      │                      │
│  │  - 2D内插表: Γ(r, v)                    │                      │
│  └──────────────────┬──────────────────────┘                      │
│                     │                                              │
│  ┌──────────────────▼──────────────────────┐                      │
│  │         模拟引擎层 (Trajectory)          │                      │
│  │  Simulation_Trajectory.cpp              │                      │
│  │  + Simulation_Utilities.cpp             │                      │
│  │  - RK45自适应传播                        │                      │
│  │  - 蒙特卡洛散射                          │                      │
│  │  - 在线bincount与捕获判定                 │                      │
│  └──────────────────┬──────────────────────┘                      │
│                ┌────┴────┐                                         │
│           ┌────▼───┐ ┌──▼───┐                                     │
│           │ 捕获    │ │ 反射 │                                     │
│           └────┬───┘ └──┬───┘                                     │
│                │        │                                          │
│  ┌─────────────▼────────▼──────────────────┐                      │
│  │         数据生成层 (Data_Generation)      │                      │
│  │  - MPI AllReduce/AllGather 汇总          │                      │
│  │  - captured/not_captured统计             │                      │
│  │  - Capture Mode快速捕获率                 │                      │
│  └──────────────────┬──────────────────────┘                      │
│                     │                                              │
│  ┌──────────────────▼──────────────────────┐                      │
│  │         反射谱层 (Reflection_Spectrum)    │                      │
│  │  - KDE核密度估计 → 速度分布P(v)          │                      │
│  │  - 进入率/反射率/通量计算                │                      │
│  │  - 微分DM通量 dΦ/dv                     │                      │
│  └──────────────────┬──────────────────────┘                      │
│                     │                                              │
│  ┌──────────────────▼──────────────────────┐                      │
│  │        分析层 (Parameter_Scan / 检测)     │                      │
│  │  - obscura库: 探测器p值计算              │                      │
│  │  - STA边界跟踪 / 全网格扫描             │                      │
│  │  - 排斥极限曲线提取                      │                      │
│  └─────────────────────────────────────────┘                      │
│                                                                    │
│  ┌─────────────────────────────────────────────────────────────┐  │
│  │              后处理Python管线                                │  │
│  │  evaporation.py ─── 二进制轨迹 → 捕获时间分析              │  │
│  │  evaporation_theory.py ─ 理论蒸发率E_⊙计算                 │  │
│  │  bincount2.py ──── 捕获轨迹 → 径向直方图 + <v²> 分布       │  │
│  └─────────────────────────────────────────────────────────────┘  │
└────────────────────────────────────────────────────────────────────┘
```

### 4.2 核心C++模块关系

| 模块 | 文件 | 职责 |
|------|------|------|
| **主入口** | `main.cpp` | 程序入口，MPI初始化，运行模式分发 |
| **配置** | `Parameter_Scan.cpp/hpp` | 配置文件解析（libconfig++），参数空间扫描 |
| **太阳模型** | `Solar_Model.cpp/hpp` | AGSS09数据加载，物理量内插，散射率计算 |
| **轨迹模拟** | `Simulation_Trajectory.cpp/hpp` | RK45积分器，散射处理，轨迹记录 |
| **工具** | `Simulation_Utilities.cpp/hpp` | Event结构体，初始条件生成，Kepler推进 |
| **数据生成** | `Data_Generation.cpp/hpp` | 蒙特卡洛主循环，MPI并行，统计汇总 |
| **反射谱** | `Reflection_Spectrum.cpp/hpp` | KDE速度分布，DM通量计算 |
| **暗光子** | `Dark_Photon.cpp/hpp` | 暗光子相互作用模型，形因子，截面 |

外部依赖：
- **obscura** (v1.0.1)：提供暗物质粒子基类、晕模型、探测器框架、核形因子
- **libphysica**：数学工具（内插、积分、统计、自然单位）
- **MPI**：进程间并行通信
- **Boost**：C++工具支持
- **libconfig++**：配置文件解析

### 4.3 构建系统

```
CMakeLists.txt (根目录)
├── C++11标准
├── MPI Required
├── FetchContent: obscura v1.0.1 (from GitHub)
├── FetchContent: Google Test (测试框架)
├── 静态库: lib_damascus_sun
│   ├── Dark_Photon.cpp
│   ├── Data_Generation.cpp
│   ├── Parameter_Scan.cpp
│   ├── Reflection_Spectrum.cpp
│   ├── Simulation_Trajectory.cpp
│   ├── Simulation_Utilities.cpp
│   └── Solar_Model.cpp
└── 可执行: DaMaSCUS-SUN → 链接 lib_damascus_sun + obscura + MPI
```

构建命令：
```bash
cd build && cmake --build . --config Release && cmake --install .
```

运行命令：
```bash
mpirun -n 32 ./DaMaSCUS-SUN config_Lingyu.cfg
```

---

## 5. 模拟阶段的详细实现

### 5.1 Phase 1: 初始化

**配置文件解析**（`Parameter_Scan.cpp :: Configuration`）：

读取 `.cfg` 文件中的所有参数，包括：
- 暗物质属性：质量 $m_\chi$，自旋，截面 $\sigma_p$（或 $\sigma_e$），相互作用类型（SI/SD/暗光子）
- 运行参数：样本大小，最大轨迹数，输出目录，等反射环数
- 暗光子专属参数：$\epsilon$（动能混合），$\alpha_D$（暗规范耦合），$m_{A'}$（介质子质量），形因子类型
- 可选数值模拟参数：`max_trajectories`（未设置时默认为 `sample_size * 1000`）、`snapshot_enabled`、`snapshot_interval`、`max_trajectory_wall_time_sec`
- 快速捕获率参数：`capture_mode = true`（也可使用 `run_mode = "Capture"`）。该模式只用于估计捕获率，$E < 0$ 后立即停止当前轨迹，不写模拟输出文件。

示例：

```cfg
run_mode = "Parameter point";
capture_mode = true;
sample_size = 1000;
max_trajectories = 100000;
```

**全局配置变量**（`Parameter_Scan.cpp`）：
- `g_top_level_dir`：从配置文件读取的输出目录
- `g_max_trajectories`：从配置文件读取的最大轨迹数安全阀

**太阳模型加载**（`Solar_Model.cpp`）：

从 `data/model_agss09.dat` 加载 AGSS09 标准太阳模型数据（~200个径向点 × 35列），构建以下内插表（1D spline）：

- $M(r)$：包含质量 $[g]$
- $T(r)$：温度 $[K]$
- $\rho(r)$：质量密度 $[g/cm^3]$
- $v_\text{esc}^2(r)$：局域逃逸速度平方
- $n_e(r)$：电子数密度
- $n_i(r)$：29种同位素各自的数密度

**散射率2D内插**（`Solar_Model.cpp :: Interpolate_Total_DM_Scattering_Rate`）：

预计算散射率的二维查找表 $\Gamma_\text{total}(r, v)$，在 $N_r \times N_v$ 网格上（默认 $1000 \times 1000$），速度范围 $[0, 0.75c]$。此预计算通过MPI分发到各进程，大幅减少模拟中的重复散射率运算。

### 5.2 Phase 2: 单粒子轨迹模拟

**核心循环**（`Trajectory_Simulator::Simulate`）：

```
输入: 初始条件 Event(t, r, v)
while (未终止):
    1. 自由传播 Propagate_Freely()
       - RK45 积分 Kepler 方程组
       - 累积光学深度，触发散射或逃逸
       - 每个 RK45 步在线累积径向 bincount
       - 每个 RK45 步计算 E = ½mχ(v² - v_esc²)
    2. if 散射被触发:
       Scatter()
       - 选择靶粒子 (Sample_Target)  
       - 抽样靶热速度 (Sample_Target_Velocity)
       - 弹性碰撞 → 新暗物质速度 (New_DM_Velocity)
    3. 检查终止条件:
       - r > R_max 且 v > v_esc → 反射/自由逃逸
       - 普通模式：E < 0 → 标记为捕获，但继续完整轨迹
       - Capture Mode：E < 0 → 立即终止当前轨迹并计为捕获
       - 散射次数达到 maximum_number_of_scatterings → 强制终止
       - 自由传播步数达到 maximum_free_time_steps → 本段传播终止
       - 速度超过 0.75c、非有限数值或单轨迹 wall-time 超限 → 中止该轨迹
输出: Trajectory_Result {initial, final, N_scat}
```

**在线统计与输出策略**：

当前 C++ 主模拟不再把每条轨迹写成 `.dat` 文件，而是在内存中在线累积径向箱计数。径向范围为 `[0, 2R_\odot)`，共 `2000` 个 bin，bin 宽约 `695.7 km`。每个 RK45 步用前一状态和当前步长累积：

$$\sum \Delta t,\qquad \sum v^2 \Delta t$$

普通模式按 captured / not_captured 分别累积 `bincount` 与每轨迹平方和，用于输出误差估计。Capture Mode 为了快速扫描捕获率，会跳过 `bincount`、蒸发记录、snapshot 和反射谱数据累积，只保留轨迹总数、捕获数和捕获率误差。

### 5.3 Phase 3: 粒子命运分类

模拟结束后，每个粒子被分类为三种结局之一：

| 结局 | 判定条件 | 物理含义 |
|------|---------|---------|
| **自由（Free）** | $N_\text{scat} = 0$，$r > R_\odot$ | 粒子穿越太阳但未散射（光学薄极限） |
| **反射（Reflected）** | $N_\text{scat} \geq 1$，$v > v_\text{esc}$，$r > R_\odot$ | 散射后仍保持正能量并逃逸 |
| **捕获（Captured）** | $E = \frac{1}{2}m_\chi(v^2 - v_\text{esc}^2) < 0$ | 散射后总能量为负，被引力束缚 |

**样本计数逻辑**：`sample_size` 在当前实现中表示目标 captured 数量，而不是固定总入射轨迹数。每个 MPI rank 的目标捕获数为 `ceil(sample_size / N_ranks)`；总轨迹数由实际达到目标捕获数前经历的所有 free / reflected / captured 轨迹共同决定。`max_trajectories` 是安全阀，达到后提前停止并给出 `EARLY STOP` 标记。

### 5.4 Phase 4: MPI数据汇总

**Data_Generation.cpp** 的 MPI 汇总方式：

1. **独立采样**：各 rank 独立生成轨迹，直到达到本 rank 的目标捕获数或最大轨迹数。
2. **全局聚合**：`MPI_Allreduce` 汇总轨迹数、捕获数、bincount、计算时间、RK45步数和 early-stop 标记。
3. **蒸发记录收集**：普通模式用 `MPI_Allgather` / `MPI_Allgatherv` 汇总各 rank 的 `EvaporationRecord`。
4. **Snapshot 合并**：启用 snapshot 时，各 rank 写 checkpoint，任意 rank 可尝试合并已就绪的 snapshot 报告。

### 5.5 Phase 5: 反射谱构建

**Reflection_Spectrum.cpp**：

1. **核密度估计（KDE）**：对所有反射粒子的最终速度样本应用Gaussian核密度估计，得到平滑的速度概率分布 $P(v)$，带边界校正因子 $f_\text{KDE} = 0.75$。

2. **通量计算**：

$$\Phi_\text{reflect}(v) = \frac{N_\text{rings}}{4\pi d^2} \times \Gamma_\text{entering} \times R_\text{reflect} \times P(v)$$

其中 $R_\text{reflect}$ 为反射粒子占总入射粒子的比例，$d$ 为探测距离（默认1 AU）。

3. **等反射环分层**：将太阳反射方向按角度分为 $N$ 个等面积环（$\theta_i = \arccos(\cos\theta_{i-1} - 2/N)$），对每个环独立计算反射谱，以提高角度分辨率。

### 5.6 Phase 6: 统计分析与排斥极限

**Parameter_Scan.cpp**：

- **p值计算**：对每个 $(m_\chi, \sigma)$ 参数点，将模拟反射谱与 obscura 库中的探测器模型（如LUX、XENON等）比对，计算统计检验p值。
- **排斥极限**：在二维参数网格中，标记 $p < p_\text{critical}$（如5%对应95% CL）的区域为"被排斥"。
- **STA算法**（Square Tracing Algorithm）：高效边界跟踪算法，不需评估整个 $N^2$ 网格，仅沿排斥区域边界行走，计算复杂度降至 $O(\text{周长})$。
- **轮廓提取**：在相邻网格点之间通过对数空间线性内插，精确确定排斥边界位置。

---

## 6. 后处理数据管线

### 6.1 蒸发率理论计算 (`evaporation_theory.py`)

基于 Garani & Palomares-Ruiz (2025, Erratum) 的解析公式，直接从太阳模型计算理论蒸发率 $E_\odot(m_\chi, \sigma)$。代码主要实现自旋相关（SD）相互作用的蒸发率计算（`compute_evap_vs_sigma_SD_full()`），包含10种太阳元素靶粒子：H1, He4, C12, N14, O16, Ne, Mg, Si, S, Fe。

还实现了完整的捕获率计算（`capture_rate_SD()`）、湮灭率计算（`annihilation_rate()`）和暗物质粒子数演化（`N_total()`），支持对任意 $(m_\chi, \sigma)$ 参数点计算完整的捕获-蒸发-湮灭平衡状态。

**数值方法**：使用 Gauss-Legendre 求积（速度积分 $n_w=80$ 点，出射速度 $n_v=150$ 点），并对指数因子实施 `np.clip(exponent, -500, 500)` 截断以防止数值溢出。抑制因子 $s(r)$ 中使用超几何函数 $_0F_1$ 与 Bessel 函数回退以保证中等光学深度区域的数值稳定性。

优化：对 $\sigma$-无关部分预计算并缓存（LTE分布、等温分布、单位捕获率），仅对 Knudsen-LTE 混合比例随 $\sigma$ 变化，大幅加速参数空间扫描。

**核心公式 — Eq. (5.1):**

$$E_\odot = \int_0^{R_\odot} 4\pi r^2 dr \sum_i \int_0^{v_e(r)} 4\pi w^2 [n_\chi f_\chi(w, r)] \int_{v_e(r)}^{\infty} R_i^+(w \to v) \, dv \, dw$$

其中：
- $n_\chi f_\chi(w, r)$ 是暗物质在太阳中的相空间分布
- $R_i^+(w \to v)$ 是从速度 $w$ 散射到速度 $v$ 的正向散射率（Eq. A.23）
- $v_e(r)$ 是局域逃逸速度
- 求和遍历所有靶粒子种类 $i$

**暗物质相空间分布**采用 Knudsen 数控制的 LTE-等温过渡（Eq. 4.11）：

$$n_\chi f_\chi = f(K) \cdot n_\text{LTE} f_\text{LTE} + [1 - f(K)] \cdot n_\text{iso} f_\text{iso}$$

其中 Knudsen 数 $K = \ell(0) / r_\chi$（平均自由程与暗物质尺度高度之比），过渡函数 $f(K) = 1/(1 + (K/K_0)^2)$，$K_0 = 0.4$。

**抑制因子 $s(r)$**（Eq. C.8-C.13）修正了光学深度效应和多次散射：

$$s(r) = \eta_\text{ang}(r) \times \eta_\text{mult}(r) \times e^{-\tau(r)}$$

**暗物质温度 $T_\chi$** 通过能量平衡方程（Eq. B.11）数值求根获得：

$$\int \sum_i n_i \frac{m_\chi m_i}{(m_\chi + m_i)^2} \langle v_\text{rel} \rangle (T_\odot - T_\chi) e^{-m_\chi \phi / T_\chi} r^2 dr = 0$$

### 6.2 蒸发时间与径向统计（当前 C++ 在线输出）

普通 `Parameter point` 模式直接在 `Data_Generation.cpp` 中输出三类文件：

- `bincount.txt`：captured 与 not_captured 的径向占据时间 $\sum \Delta t$、速度二阶矩 $\sum v^2\Delta t$，以及逐 bin 误差估计。
- `evaporation_summary.txt`：每条 captured 轨迹的 `t_evap = t_last_negative - t_first_negative` 和 `truncated` 标记。
- `computation_time_summary.txt`：captured / not_captured 轨迹的 wall-clock 时间与 RK45 步数统计。

`capture_rate`、`capture_rate_err` 和 `capture_rate_CI_95_lower/upper` 会写入这些文件头部；Capture Mode 不写这些文件，只在终端打印同类捕获率统计。

### 6.3 Snapshot 输出

若设置 `snapshot_enabled = true`，代码按 wall-clock 间隔输出累计 snapshot。相关参数为：

- `snapshot_interval`：snapshot 间隔，默认 `60 s`。
- `max_trajectory_wall_time_sec`：单条轨迹 wall-clock 上限，默认 `300 s`；用于避免某个 MPI rank 卡在单条轨迹上，导致 snapshot 一直处于 waiting 状态。

Snapshot 会合并各 rank 的当前进度，包括已完成轨迹的 captured / not_captured bincount，以及正在运行轨迹的临时 bincount。Capture Mode 会强制关闭 snapshot，因为该模式只需要终端捕获率。

---

## 7. 计算效率优化与物理影响

### 7.1 散射率2D内插表

**优化内容**：预计算 $\Gamma_\text{total}(r, v)$ 的二维内插表（`Interpolate_Total_DM_Scattering_Rate`），替代每步的逐组分直接计算。

**效率提升**：散射率计算涉及对29种同位素求和，每个需要数密度查询、热平均和截面积分。预计算将 $O(N_\text{iso} \times N_\text{steps})$ 降低至 $O(N_\text{steps})$（仅需一次2D内插查询）。

**物理影响**：内插引入的数值误差通常在百分之几以内。由于散射率变化平滑（位置和速度上均缓慢变化），双线性/样条内插精度足够。然而，在太阳核心区域（$r < 0.1 R_\odot$），温度和密度梯度陡峭，需要较密的径向网格（默认1000点）以保证精度。

### 7.2 在线统计替代逐轨迹文件

**优化内容**：将原始逐轨迹文件记录改为 C++ 内部在线统计。每个 RK45 步只更新当前轨迹的 `TrajectoryBincount`，轨迹结束后再按 captured / not_captured 汇总到全局数组。

**效率提升**：避免了大量小文件 I/O 和后处理扫描，MPI rank 之间只需要在模拟结束或 snapshot 边界合并聚合量。

**物理影响**：径向统计仍按真实 RK45 步长加权，保留 $\sum \Delta t$ 与 $\sum v^2\Delta t$，可直接用于估计被捕获暗物质的径向分布与速度二阶矩。

### 7.3 太阳内外的步长策略分区

**优化内容**：
- **太阳内部**（$r < R_\odot$）：步长受 $0.1/\Gamma_\text{total}$ 约束，确保正确采样散射事件
- **太阳外部**（$r > R_\odot$）：无散射可能，RK45步长可自由增长到大值

**物理影响**：太阳外部为纯 Kepler 问题，数值积分可使用大步长而不损失精度（轨道是光滑解析解）。这显著加速了粒子在太阳外区域的传播（从 $R_\odot$ 到 $2R_\odot$ 的初始/最终传播区域，以及双曲轨道在远距离的演化）。

### 7.4 从固定数组到动态向量的内存管理改进

**优化内容**：将原始的固定大小数组 `int traj_captured[10000]` 改为 `std::vector<int> traj_captured`，使用 `push_back` 动态增长。

**效率与安全影响**：
- 避免缓冲区溢出（当捕获粒子超过10000个时）
- 在捕获粒子少时节省内存
- 不影响物理结果，纯粹的工程改进

### 7.5 在线 bincount 与 Capture Mode 快速统计

**优化内容**：当前主模拟不再依赖逐轨迹 `.dat` 文件保存/删除流程，而是在 C++ 内部在线累积 captured 与 not_captured 的径向 `bincount`、误差平方和、蒸发时间记录和计算时间统计。

当配置 `capture_mode = true` 或 `run_mode = "Capture"` 时，模拟只关心捕获率：任意轨迹第一次满足 $E < 0$ 就立即停止并计为 captured，同时跳过 `bincount`、蒸发记录、snapshot 和反射谱输出。

**效率影响**：普通模式保留完整统计用于后处理；Capture Mode 避免了被捕获粒子后续长时间束缚轨道模拟，适合快速扫描多个参数点的捕获率。

### 7.6 历史后处理脚本

`bincount2.py` 是早期逐轨迹文件工作流的后处理脚本，使用多进程并行读取轨迹文件并生成径向统计。当前主 C++ 模拟已经在线生成 `bincount.txt`，通常不再需要该脚本参与标准流程。

保留该脚本主要用于读取历史数据或交叉检查旧格式输出。

### 7.7 Kepler 解析推进替代数值积分

**优化内容**：在远距离（$>2R_\odot$）使用双曲 Kepler 轨道解析公式 `Hyperbolic_Kepler_Shift()` 代替逐步RK45数值积分。

**物理影响**：完全精确（Kepler问题有精确解析解），同时省去了长距离（从1000 AU到$2R_\odot$，或从$2R_\odot$到1 AU）的数值积分成本。这是合理的，因为太阳外的引力场为纯 $1/r^2$ 形式。

### 7.8 MPI 合并与负载分配

每个 MPI rank 独立生成轨迹，目标捕获数为 `ceil(sample_size / N_ranks)`，最大轨迹数为 `ceil(max_trajectories / N_ranks)`。模拟完成后通过 `MPI_Allreduce` 合并轨迹数、捕获数、bincount、计算时间和 RK45 步数；蒸发记录通过 `MPI_Allgather` / `MPI_Allgatherv` 汇总。

Snapshot 模式额外使用每 rank 的二进制 checkpoint 文件来合并中间进度，避免长时间运行时只能等最终 `MPI_Allreduce`。

---

## 8. 与文献的对照分析

### 8.1 Garani & Palomares-Ruiz (2025 Erratum): 暗物质在太阳中的散射

本文（以及其原始版本 arXiv:1702.02768）是 `evaporation_theory.py` 的直接理论基础。

**对应关系**：
- **Eq. (5.1)**：蒸发率主公式 → `compute_evap_vs_sigma_SD_full()` 的三重积分
- **Eq. (A.23)**：常数截面的正向散射核 $R^+$ → `R_plus_const_vec()` 函数
- **Eq. (B.11)**：暗物质温度能量平衡 → `solve_T_chi()` 使用 Brent 方法求根
- **Eq. (C.8)-(C.13)**：抑制因子（光学深度、多次散射校正）→ `compute_suppression()` 函数
- **Eq. (4.7)**：LTE分布 → `compute_n_LTE()` 函数
- **Eq. (4.10)-(4.12)**：Knudsen数与LTE-等温过渡 → `knudsen_number()`, `f_K()` 函数

散射核 $R^+$ 的表达式中包含 erf 函数对，精确描述了具有 Maxwell-Boltzmann 热分布的靶粒子束与暗物质的散射运动学。代码对指数因子实施了 `np.clip(exponent, -500, 500)` 截断以防止数值溢出——这在极端质量比 $\mu \gg 1$ 或 $\mu \ll 1$ 情况下是必要的。

Gould & Raffelt (1990) 的热扩散系数 $\alpha_0(\mu)$ 查表也被代码完整实现，用于LTE分布的构建。

### 8.2 Emken (2022): 轻暗物质的太阳反射（重介质子情形）

本文的核心概念是**太阳反射暗物质（SRDM）**：暗物质粒子经太阳散射后以加速的速度逃逸，到达地球。

**对应关系**：
- **暗光子模型**：Dark_Photon.cpp 实现的动能混合暗光子 $A'$ 即本文讨论的重介质子模型（$m_{A'} > 0$）
- **形因子**：$F_\text{DM}^2(q) = (q_\text{ref}^2 + m_{A'}^2)^2 / (q^2 + m_{A'}^2)^2$ 对应论文中的一般情形
- **微分DM通量**：Reflection_Spectrum 计算的 $d\Phi/dv$ 对应论文中 Figure 3 的反射速度谱
- **蒙特卡洛方法**：本文指出解析方法对多次散射不精确，DaMaSCUS-SUN 通过直接轨迹模拟避免了这个问题

本代码实现了论文中讨论的四种形因子类型（Contact、Long-Range、General、Electric-Dipole），均对应不同的暗光子质量极限：
- Contact ($m_{A'} \gg q_\text{max}$)：截面与动量转移无关
- Long-Range ($m_{A'} \to 0$)：截面 $\propto 1/q^4$，前向散射主导
- General：过渡区域，物理最丰富

### 8.3 Busoni et al. (2013): 太阳中微子约束的最小暗物质质量

本文引入了**蒸发质量**（evaporation mass）的概念：暗物质粒子质量低于此阈值时，蒸发过程在太阳寿命内将其完全清除，太阳无法有效积累暗物质。

**对应关系**：
- **蒸发边界**：`generate_recording_table.py` 和 `evaporation_theory.py` 计算的 $E_\odot(m_\chi, \sigma)$ 正是用于确定蒸发质量边界
- **捕获-蒸发平衡**：当 $C_\odot \cdot t_\odot \ll E_\odot \cdot t_\odot$ 时，暗物质数量趋于零；平衡条件 $N_\chi \approx C_\odot / E_\odot$
- **记录间隔物理**：查找表中蒸发率极低的区域对应 $m_\chi$ 远高于蒸发质量（粒子被稳定捕获），而蒸发率高的区域对应 $m_\chi$ 接近或低于蒸发质量

### 8.4 Nguyen & Linden (2026): 太阳约束的自旋相关暗物质-核子散射

本文是本项目的直接研究目标之一，讨论太阳如何约束**蒸发极限以下**的自旋相关（SD）暗物质-核子散射。

**对应关系**：
- **SD相互作用**：`config_Lingyu.cfg` 中 `DM_interaction = "SD"` 配置对应本文的自旋相关分析
- **质子耦合**：`DM_relative_couplings = (1.0, 0.0)` 表示纯质子耦合，正是论文的焦点（SD-proton）
- **轻暗物质区域**：$m_\chi = 0.5 \text{ GeV}$ 的配置对应论文讨论的蒸发质量附近区域
- **散射靶为氢**：SD相互作用主要通过氢的核自旋 $J = 1/2$ 实现，代码中的 `DM_Scattering_Rate_Nucleus` 按同位素自旋加权

本文的关键创新是指出，即使在传统认为蒸发抹去暗物质信号的区域，太阳反射暗物质的速度谱仍然可以提供有意义的约束——因为反射不需要暗物质被捕获。

---

## 9. 总结

### 9.1 项目物理图景

DaMaSCUS-SUN 实现了暗物质与太阳相互作用的完整物理图景：

```
银河系暗物质晕                              地球探测器
     │                                         ▲
     │ 晕速度分布采样                           │ 反射速度谱
     ▼                                         │
  ┌──────┐  引力聚焦    ┌──────────┐          ┌─────────┐
  │ 初始  │ ─────────→ │ 太阳内部  │ ───────→ │ 1 AU处  │
  │ 条件  │  Kepler推进 │ 多次散射  │  逃逸后  │ 微分通量 │
  └──────┘             │          │  Kepler  └─────────┘
                       │  RK45积分 │
                       │  +MC散射  │
                       └────┬─────┘
                            │
                    ┌───────┴───────┐
                    │               │
               ┌────▼────┐    ┌────▼────┐
               │ 捕获粒子 │    │ 自由粒子 │
               │ (E < 0)  │    │ (N=0)   │
               └─────────┘    └─────────┘
```

### 9.2 核心技术特点总结

| 特性 | 实现 | 优势 |
|------|------|------|
| **蒙特卡洛轨迹模拟** | 逐粒子完整轨迹 | 精确处理多次散射，无解析近似 |
| **RK45自适应积分** | Fehlberg嵌入式方法 | 精度可控，步长自适应 |
| **AGSS09太阳模型** | 29种同位素密度+温度剖面 | 现代最佳太阳组成数据 |
| **暗光子模型** | 4种形因子 | 覆盖从接触到长程的完整参数空间 |
| **MPI并行化** | 各rank独立采样 + AllReduce/AllGather 汇总 | 线性可扩展（32+进程） |
| **2D散射率缓存** | $\Gamma(r,v)$ 内插表 | 数量级加速内循环 |
| **STA边界跟踪** | 排斥区域周长遍历 | 避免完整 $N^2$ 网格计算 |
| **在线径向bincount** | 每个RK45步累积 $\Delta t$ 与 $v^2\Delta t$ | 避免逐轨迹大文件输出，直接产出捕获/未捕获统计 |
| **理论-模拟交叉验证** | `evaporation_theory.py` ↔ 轨迹数据 | 确保数值结果的物理一致性 |

### 9.3 数据管线完整流程

```
     evaporation_theory.py                    config_Lingyu.cfg
     (理论蒸发率计算)                          (运行参数)
            │                                       │
            │ (历史产物 → 硬编码LUT)                 │
            ▼                                       ▼
  Get_Recording_Interval()                 DaMaSCUS-SUN (C++/MPI)
  (96条记录间隔查找表,                     (主模拟程序)
   已嵌入 Simulation_Trajectory.cpp)              │
            │                                       │
            └──────────┬────────────────────────────┘
                       ▼
              Simulation_Trajectory.cpp
              (轨迹模拟, 基于步数的记录
               + 运行时自校准)
                       │
          ┌────────────┼────────────┐
          ▼            ▼            ▼
   evaporation.py  bincount2.py  Reflection_Spectrum
   (捕获时间)      (捕获粒子      (反射速度谱)
                    径向分布)
          │            │            │
          └────────────┼────────────┘
                       ▼
              Parameter_Scan / 探测器约束
              (排斥极限曲线)
```

### 9.4 物理意义与创新点

本项目的核心贡献在于将**直接蒙特卡洛轨迹模拟**应用于轻暗物质的太阳物理：

1. **精确处理多次散射**：不同于单次散射近似或扩散方程方法，MCTS 追踪每次散射事件的完整运动学，自然包含了非平衡效应和高阶散射贡献。

2. **统一框架**：同一模拟同时产出捕获率、反射谱和蒸发信息，无需分别运行不同程序。

3. **蒸发边界以下的探索**：通过精确模拟轻暗物质（$m_\chi < m_\text{evap}$）的轨迹，验证理论蒸发率公式的准确性，并探索该区域的实验约束可能性。

4. **暗光子模型的系统实现**：完整的 $q$-依赖形因子处理使得代码可以研究从contact到long-range的所有介质子质量区域，覆盖论文 Emken (2022) 讨论的全部物理场景。

---

*本文档基于对项目全部源代码、头文件、Python脚本、配置文件、构建系统和参考文献的完整审阅编写而成。*
