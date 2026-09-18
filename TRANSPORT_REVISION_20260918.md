# Post-capture transport revision — 2026-09-18

本轮以用户提供的 ARTICLE_REVISION_PLAN.md 和 CODE_AND_FUTURE_WORK_PLAN.md 为修改依据。它们是研究建议；本文逐项区分已实现的接口、已执行的验证和仍需计算的研究工作。原始代码基线：`8d5dadb4a94ec1b9b5f3027fef1e8941e82667e9`。

## 当前计算定义

- 入射状态位于 2 R_sun；抽样仍是引力聚焦后会穿过太阳的入射集合，不能用普通局域速度抽样替代。
- 散射区为 r < R_sun，默认数值匹配面 1.1 R_sun。负能量外轨道解析返回；若远日点达到外边界，则在第一次向外到达该面时终止。默认外边界为 1100 R_sun，可设 3000 或 10000 R_sun。
- `physical_escape` 与 `outer_orbit_removed` 都是完整捕获历史，但分别计数。外边界只是轨道移除处方，未模拟行星扰动。
- `sample_size` 在 Parameter point 模式表示完整捕获历史数；在 Capture 模式表示固定入射尝试数。Capture 模式每条入射在首次捕获或未捕获离开时结束。
- 设置 `production_mode = true`，只要存在数值失败、wall-time/step/scattering 截断，或目标未满足，就停止领取新轨迹，等待已在执行的轨迹结束后以非零状态退出，并写 `production_accepted: false`。无效历史不贡献最终 captured residence。即使关闭退出门禁，后处理仍拒绝不合格 metadata。
- `thermal_validation_mode = true` 仅用于内部归一化形状；允许显式计算截断，绝对归一化永远不合格。不能同时启用 production_mode。
- 完整捕获 occupation 在匹配面验证逃逸时结束；其后自由飞出外边界的 dt/v2dt 单独存入 post_evap。未散射的太阳相交 transit 单独保存；它不是完整银河晕密度。
- 均不使用太阳年龄截断。稳态占据解释仍要求所考虑种群达到统计稳态，且湮灭耗尽可忽略。

## 输出及单位

`metadata.json` 是最后发布的完成标记，schema_version=6；含实际非零主随机种子、MPI 核数、源码 SHA256、编译器、构建选项、物理参数、原生径向边界和 acceptance。`input.cfg` 保留完整输入。源码哈希覆盖本项目 C++、头文件、明确列出的交付 Python 分析脚本、构建 CMake 模块和随库太阳表，包含相对文件名与内容哈希；忽略的本地脚本不参与指纹；不包含外部依赖及环境覆盖的太阳表，实际太阳参考表另行导出并在配对分析时检查。

| 文件 | 内容 |
| --- | --- |
| capture_summary.json | 固定入射标记、N_inj、N_capt、C_geom、64-block 计数 |
| radial_blocks.tsv | 64 个块的捕获、未散射 transit、post_evap 的 dt[s] 和 v2dt[km^2/s]；边界[km] |
| block_counts.tsv | 完整捕获历史的块计数 |
| trajectory_summary.tsv | 每条捕获历史的终止原因、捕获状态、驻留时间、远日点、散射/回捕次数等标量 |
| orbit_class_blocks.tsv | 按最大远日点分组的稀疏块驻留和；以完整捕获集合归一化 |
| termination_counts.tsv | 捕获和未捕获终止原因计数 |
| solar_reference.tsv | 实际太阳模型的半径[cm]、温度[K]、氢密度[cm^-3]、相对中心势[km^2/s^2] |

旧 `bincount.txt`、蒸发时间和 survival 输出保留兼容性；绝对论文结果应从 schema-6 产品及独立 capture run 重建，不能把旧输出中的内部捕获估计当作新的独立归一化。outer_removed 不进入 evaporation_times。

## 可复现运行

依赖和 CMake 构建方式见 README。Python 分析需要 numpy、scipy、matplotlib；可用 `python3 scripts/test_analyze_point.py` 检查几何。

```bash
python3 scripts/prepare_transport_runs.py runs/pilot --phase pilot
# manifest.json 中每一点有 capture 与 transport 两份独立种子配置。
# 在项目根目录将实际可执行文件依次用于所选配置：
./build-transport-revision/src/DaMaSCUS-SUN /absolute/path/to/capture.cfg
./build-transport-revision/src/DaMaSCUS-SUN /absolute/path/to/transport.cfg
python3 scripts/analyze_point.py /absolute/path/to/results_MASS_SIGMA --capture-dir /absolute/path/to/results_capture_MASS_SIGMA
python3 scripts/analyze_thermal.py /absolute/path/to/thermal_result
python3 scripts/summarize_transport_scan.py /absolute/path/to/point1/derived.json /absolute/path/to/point2/derived.json --output /absolute/path/to/scan_summary
```

生成器另有 `--phase thermal`、`production`、`cutoff`。生成配置不会启动任务。默认 pilot 为 500 条完整捕获历史；固定入射 capture 数量必须根据 pilot 捕获率调整。极低捕获率下不能仅凭固定 10^5 次入射认定统计充分。生产阶段代表点配置三个独立种子。cutoff 默认比较点需要根据实际 pilot 的 removal fraction 再选择；不能把模板点当作已完成敏感性测试。

`analyze_point.py` 检查两个 run 的物理配置、太阳参考表、源码哈希及互不重叠的 MPI 秩随机种子（含 uint32 回绕），逐块验证计数、太阳内/外驻留时间及轨道类别的闭合。输出占据、有效体积、湮灭事件率、完整球面伽马视线积分、角包含尺度、1 AU 外源比例、可见率、中微子直接 νν̄ 的无衰减全味参考源。独立 capture 与 transport 的 64-block delete-one 方差相加，逐次重算非线性量。轨道类别采用 mu_class × mu_total 分配二次贡献，保留类别间交叉项。该分析尚未对密度平方的有限样本偏差作修正；需要样本量和分箱收敛检查。

`analyze_thermal.py` 从同一 C++ 太阳模型导出的势和靶分布解质子能量平衡温度，并比较内部概率和二阶速度矩。`halo_focused_density.py` 需要显式 SHM 速度参数，计算透明太阳、无碰撞的球对称 halo 参考以及 halo^2、cross、captured^2 同域积分；没有银河晕散射耗减或观测背景拟合。中微子传播/探测器响应尚未接入。

## 已验证与待做

已加入确定性势连续、椭圆周期/速度矩、近抛物线、单程移除、双曲线矩、太阳遮挡和地球轨道外源几何测试；保留现有 RK、MPI 调度、失败输出传播和 replay 回归。初轮本地验证记录位于 `validation-results/transport_20260918/VALIDATION.md`；原始输出按仓库策略不提交。下面的推送前审查摘要随源码交付，供没有本地结果目录的读者查阅。

以下项目仍需研究计算或单独开发，不能由本轮小样本推断完成：

1. 2/3/4/10 GeV、10^-36 cm^2 的形状恢复以及散射预算/独立历史收敛。
2. 全质量/截面网格的 500-capture pilots、5×10^4–10^5 捕获生产样本、种子重复与稳定性检查。
3. 1100/3000/10000 R_sun 在最敏感点的统计比较；1100 只是默认处方。
4. 旧结果历史回归的物理差异归因、完整 transit 分布比较和插值精度评估。
5. 中微子传播、仪器响应/背景、线搜索似然及依赖其的排除线。
6. 精确 checkpoint/restart：当前周期 snapshot 仅供监控，不包含恢复所需的 RNG、MPI 调度和所有统计量状态。metadata 明确写 `restart_supported=false`；中断后必须使用新目录重新计算。

本轮没有把旧 finite-age 数据重新标成 complete-history production，也没有自动启动整个大规模网格。

## 端到端测试触发的额外修复

固定 10000 次入射测试暴露了旧 Capture 快速光学深度近似的问题：低入射能量的第 8535 条历史在首次碰撞前被插值为负能量。现在 Capture 与运输使用相同的碰撞定位精度，重放该状态正常达到首次物理捕获。Capture 模式同时保留并跨 MPI 收集 invalid replay ledger；重放工具读取配置中的 Capture 模式和外边界。回归测试对相同初态/种子检查首次捕获时间、半径和能量一致。两组最终样本均使用原种子重新执行，不以筛选随机种子绕过失败。

## 推送前审查与修复

本次检查包含当前工作区差异及待推送的基线提交。修复了以下具体问题：

- 几何入射率补入 `DM.fractional_density`，与组分密度一致；`DM_fraction=0.25` 的实际输出率严格缩放为默认值的四分之一。halo 参考也使用同一组分密度。
- `C/C_geom`、移除比例和内部二阶矩温度进入完整 delete-block 误差传播。汇总器拒绝缺失误差，避免把缺失值输出为零。
- 扫描汇总检查分析系数、太阳表、源代码及物理配置；同一点的独立重复必须同时具有独立 capture 和 transport 随机流。不同模型或复用捕获归一化的结果不能被静默当作独立重复平均。
- 近抛物线椭圆/双曲线矩使用稳定半角及小量级数；远日点使用同一套轨道不变量并显式处理端点。独立径向数值积分覆盖正、负比能偏离逃逸条件至 `1e-12` 的情况。
- 轨迹摘要的太阳内/外驻留时间从原生保守分箱累计，替代另一套直线穿越近似；兼容 `bincount.txt` 的最后一格边界也按移除半径裁剪。
- 从结果目录中的 `input.cfg` 重跑时，先完整读取再原子替换，保留输入；生产模式发生失败后停止领取新任务，避免无限补样。单进程和三秩 MPI 均覆盖这些路径。
- 外边界接受整数、64 位整数和浮点数配置；错误类型显式拒绝。修复了 `outer_removal_radius_rsun = 3000` 曾静默回退到 1100 的问题，生产/热形状开关也校验布尔类型。
- 传统 detector-limit `Parameter scan` 拒绝不受支持的生产/热形状模式和自定义移除半径，防止配置被静默忽略。运输扫描通过生成独立 Parameter point/Capture 配置执行。

最终验证包括 Release 构建、23/23 项 CTest（含 4 项 MPI）、7 项 transport contract 和 13 项 Python 几何/统计/数据完整性测试。实际单秩样本包含 10000 次入射（4480 次捕获）及独立的 500 条完整捕获历史（1067 次入射）；两者均无数值/计算失败，并通过逐块后处理闭合。该组输运样本源码指纹为 `4251662005d98a2d14a8e9d9da58eb27cead5b9160dae37542fcc35fa5470e15`。配置类型修复后的最终版本另完成三秩 MPI 配对：1000 次入射及 64 条完整捕获历史，使用整数形式的 3000 R_sun 边界；同样通过严格分析。最终源码指纹为 `25e2f17b4b948a7c035885589b3ddc9aa39d715b29caa1ec5bc76ea3c41d770f`。两组验证均为 0.01 GeV、10^-32 cm^2，capture/transport 主种子分别为 20260919/20260918。源代码清洁副本与构建指纹一致；本地另有四个被忽略的脚本，不影响指纹。完整数据和日志保存在本地 `validation-results/transport_review_20260918/`，不推送生成数据。

以上是实现及小样本验证，不等同于完成大网格、热形状恢复、外边界敏感性或尾部统计收敛；这些研究限制仍按前述待办保留。
