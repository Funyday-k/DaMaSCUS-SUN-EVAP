# Bash 使用手册

> 本手册基于当前工作环境（CHPC HPC 集群）定制，涵盖已安装的工具和常用技巧。

---

## 目录

1. [Prompt 说明](#1-prompt-说明)
2. [文件浏览 — eza](#2-文件浏览--eza)
3. [文件查看 — bat](#3-文件查看--bat)
4. [代码搜索 — ripgrep](#4-代码搜索--ripgrep)
5. [模糊搜索 — fzf](#5-模糊搜索--fzf)
6. [SLURM 任务管理](#6-slurm-任务管理)
7. [项目专用别名与函数](#7-项目专用别名与函数)
8. [历史记录技巧](#8-历史记录技巧)
9. [目录导航技巧](#9-目录导航技巧)
10. [Tab 补全增强](#10-tab-补全增强)
11. [常用 Bash 快捷键](#11-常用-bash-快捷键)
12. [文件与文本处理](#12-文件与文本处理)
13. [进程与系统监控](#13-进程与系统监控)
14. [Git 快捷操作](#14-git-快捷操作)

---

## 1. Prompt 说明

提示符格式：

```
┌─ ✔/✘  用户名@主机名  当前目录  (git分支)  [N jobs]
└─▶ 
```

| 标记 | 含义 |
|---|---|
| `✔`（绿色） | 上一条命令执行成功（退出码 0） |
| `✘ N`（红色） | 上一条命令失败，N 为退出码 |
| `(main)`（紫色） | 当前 Git 分支 |
| `[9 jobs]`（黄色） | 当前 SLURM 运行中任务数 |

---

## 2. 文件浏览 — eza

`eza` 是 `ls` 的现代替代品，支持彩色、Git 状态、文件图标。

### 别名速查

| 命令 | 说明 |
|---|---|
| `ll` | 详细列表（含权限、大小、日期、Git 状态） |
| `la` | 显示所有文件（含隐藏文件） |
| `lt` | 按修改时间倒序 |
| `lS` | 按文件大小倒序 |
| `tree` | 树形目录结构 |
| `tree -L 2` | 限制树形深度为 2 层 |

### 实用示例

```bash
# 查看当前目录树（只显示 2 层）
tree -L 2

# 只显示 .cpp 文件
ll src/*.cpp

# 按大小查看 Trajectoies 中的大文件
lS /project/kennyng/backup_DM/Trajectoies/results_0.000000_-33.000000/ | head -10

# 带 Git 状态的详细列表
eza -alF --git src/
```

---

## 3. 文件查看 — bat

`bat` 是 `cat` 的替代品，支持语法高亮、行号、Git diff。

### 别名速查

| 命令 | 说明 |
|---|---|
| `cat <file>` | 纯输出（无行号，适合管道） |
| `catn <file>` | 带行号 + 语法高亮 |
| `catf <file>` | 完整模式（行号 + 边框 + Git diff） |

### 实用示例

```bash
# 查看 C++ 源文件（语法高亮）
catn src/Simulation_Trajectory.cpp

# 只看第 90-100 行
bat -r 90:100 src/Simulation_Trajectory.cpp

# 查看配置文件
catn bin/config_Lingyu.cfg

# 查看 job 输出（语法高亮的日志）
catn bin/job-1842563

# 用 bat 查看 man 手册（已自动配置）
man squeue

# 对比两个文件的差异
bat --diff file1.cpp file2.cpp
```

---

## 4. 代码搜索 — ripgrep

`rg` 比 `grep` 快 10-100 倍，自动递归目录、自动忽略 `.git` 等。

### 别名速查

| 命令 | 说明 |
|---|---|
| `rg <pattern>` | 智能大小写搜索（有大写则区分大小写） |
| `rgi <pattern>` | 包含 `.gitignore` 中的文件 |

### 实用示例

```bash
# 在整个项目中搜索函数名
rg "max_scattering"

# 只搜索 .cpp 和 .hpp 文件
rg "Trajectory_Result" -t cpp

# 搜索并显示行号上下文（前后各 2 行）
rg "maximum_scatterings" -C 2 src/

# 搜索正则表达式（变量定义）
rg "long int.*scattering" --pcre2

# 统计每个文件中的匹配次数
rg "Scatter" --count src/

# 列出包含关键词的文件名（不显示具体行）
rg -l "DM_Scattering_Rate" src/

# 搜索并替换预览（不修改文件）
rg "old_function" --replace "new_function" src/ --passthru | head -20
```

---

## 5. 模糊搜索 — fzf

fzf 为 bash 提供交互式模糊搜索界面。

### 快捷键

| 快捷键 | 功能 |
|---|---|
| `Ctrl+R` | 模糊搜索历史命令（极其常用） |
| `Ctrl+T` | 模糊搜索当前目录下的文件（带 bat 预览） |
| `Alt+C` | 模糊搜索子目录并进入 |

### 实用示例

```bash
# 模糊搜索历史中含 "squeue" 的命令
# 按 Ctrl+R 然后输入 squeue

# 用 fzf 选择文件并用 bat 查看
bat $(fzf)

# 模糊搜索 job 文件
cat $(ls bin/job-* | fzf)

# 结合 rg 和 fzf 交互式搜索代码
rg "Scatter" src/ | fzf

# 选择 SLURM job 并取消
scancel $(squeue -u $USER -h | awk '{print $1}' | fzf)
```

---

## 6. SLURM 任务管理

### 别名速查

| 命令 | 说明 |
|---|---|
| `sq` | 格式化显示当前任务（含状态、运行时间） |
| `si` | 查看节点信息 |
| `sc <jobid>` | 取消任务 |
| `jout <jobid>` | 快速查看 job 输出最后 30 行 |

### 实用示例

```bash
# 查看所有正在运行的任务
sq

# 查看指定 job 的输出
jout 1842563

# 实时监控 job 输出（持续刷新）
watch -n 5 tail -20 bin/job-1842563

# 提交新任务
cd bin && sbatch job_Lingyu.sh

# 取消所有任务（危险！谨慎使用）
squeue -u $USER -h | awk '{print $1}' | xargs scancel

# 查看任务的详细信息
scontrol show job 1842563

# 查看节点使用情况
sinfo -o "%n %T %C" | head -20

# 排队中的任务预计等待时间
squeue -u $USER --start
```

---

## 7. 项目专用别名与函数

### 目录快捷跳转

| 命令 | 目标目录 |
|---|---|
| `cdp` | `/project/kennyng/backup_DM` |
| `cdb` | `/project/kennyng/backup_DM/bin` |
| `cds` | `/project/kennyng/backup_DM/src` |
| `cdr` | `/project/kennyng/backup_DM/Trajectoies` |

### `trajinfo` — 轨迹文件信息

```bash
# 快速查看轨迹文件概要（大小、数据点数、模拟时间）
trajinfo Trajectoies/results_0.000000_-33.000000/trajectory_1_task4.txt

# 输出示例：
# File:        ...trajectory_1_task4.txt
# Size:        27G
# Data points: 1 → 320929504
# Sim time:    6.42e+09 seconds
```

### `duf` — 目录大小统计

```bash
# 查看当前目录各子项大小（按大小倒序 Top 20）
duf

# 查看指定目录
duf Trajectoies/
```

### `extract` — 万能解压

```bash
extract archive.tar.gz
extract data.zip
extract file.tar.bz2
```

---

## 8. 历史记录技巧

历史记录已增强：5 万条、带时间戳、去重。

```bash
# 搜索历史（模糊，推荐）
Ctrl+R  → 输入关键词

# 查看最近 20 条历史（带时间戳）
history 20

# 执行上一条命令
!!

# 执行历史第 N 条命令
!N

# 执行上一条以 "squeue" 开头的命令
!squeue

# 上一条命令的最后一个参数
!$
# 例：vim src/main.cpp → catn !$  (等价于 catn src/main.cpp)

# 替换上一条命令中的字符串并执行
^old^new
# 例：上一条是 rg "max_scat"，改成
^max_scat^min_scat

# 清空当前输入行
Ctrl+C 或 Ctrl+U
```

---

## 9. 目录导航技巧

```bash
# 返回上一个目录（来回切换）
cd -

# 多级返回上级目录（已设别名）
..      # cd ..
...     # cd ../..
....    # cd ../../..

# 直接输入目录名进入（已启用 autocd）
src/    # 相当于 cd src/

# 目录栈操作
pushd /some/dir   # 进入目录并入栈
popd              # 返回栈顶目录
dirs -v           # 查看目录栈

# 快速查找大目录并进入
cd $(find . -type d -name "results*" | fzf)
```

---

## 10. Tab 补全增强

已安装 `bash-completion` 和 `complete-alias`。

```bash
# Git 命令补全
git ch<Tab>        → git checkout / cherry-pick / ...
git checkout <Tab> → 分支名自动补全

# gs (git status 别名) 也支持补全
gs <Tab>

# squeue 参数补全
squeue --<Tab>     → 显示所有参数选项

# 路径补全中自动纠错（已启用 cdspell）
cd src/Simualtion<Tab>  → 自动纠正为 Simulation_Trajectory.cpp
```

---

## 11. 常用 Bash 快捷键

### 光标移动

| 快捷键 | 功能 |
|---|---|
| `Ctrl+A` | 跳到行首 |
| `Ctrl+E` | 跳到行尾 |
| `Ctrl+←` / `Alt+B` | 向左跳一个词 |
| `Ctrl+→` / `Alt+F` | 向右跳一个词 |

### 编辑

| 快捷键 | 功能 |
|---|---|
| `Ctrl+U` | 删除光标到行首 |
| `Ctrl+K` | 删除光标到行尾 |
| `Ctrl+W` | 删除光标前一个词 |
| `Alt+D` | 删除光标后一个词 |
| `Ctrl+Y` | 粘贴（还原刚删除的内容） |
| `Ctrl+_` | 撤销 |

### 控制

| 快捷键 | 功能 |
|---|---|
| `Ctrl+C` | 终止当前命令 |
| `Ctrl+Z` | 挂起当前命令（后台） |
| `fg` | 恢复挂起的命令 |
| `Ctrl+L` | 清屏（等同 `clear`） |
| `Ctrl+D` | 退出当前 shell |

---

## 12. 文件与文本处理

```bash
# 查看大文件尾部（实时刷新）
tail -f bin/job-1842563
tail -100 Trajectoies/results_0.000000_-33.000000/trajectory_1_task4.txt

# 查看文件前几行
head -5 Trajectoies/results_0.000000_-33.000000/trajectory_1_task4.txt

# 统计文件行数（大文件慎用 wc -l，超慢）
wc -l src/Simulation_Trajectory.cpp

# 估算大文件行数（用文件大小除以平均行长）
# 27G 文件 ≈ 27*1024^3 / 85 bytes/行 ≈ 3.4 亿行

# 快速查看轨迹文件最后一个数据点
tail -1 Trajectoies/results_0.000000_-33.000000/trajectory_1_task4.txt | awk '{print "point="$1, "time="$2}'

# 统计目录下的文件数
ls Trajectoies/results_0.000000_-33.000000/ | wc -l

# 查看第 N 行
sed -n '100p' src/Simulation_Trajectory.cpp

# 统计各文件扩展名数量
find src/ -type f | sed 's/.*\.//' | sort | uniq -c | sort -rn

# AWK 提取列（例：提取轨迹文件radius列）
awk '{print $10}' trajectory.txt | head -20

# 对结果排序去重统计
cat file.txt | sort | uniq -c | sort -rn | head -10
```

---

## 13. 进程与系统监控

```bash
# 交互式进程查看（已安装）
htop

# 查看磁盘使用（友好格式）
df -h

# 查看当前目录各文件夹占用（友好格式）
duf

# 查看内存使用
free -h

# 查看当前 shell 进程
ps aux | rg $USER

# 后台运行命令并保留输出
nohup command > output.log 2>&1 &

# 查看后台任务
jobs -l

# 查看端口占用
ss -tulnp 2>/dev/null | head -20
```

---

## 14. Git 快捷操作

| 命令 | 说明 |
|---|---|
| `gs` | `git status` |
| `gl` | `git log --oneline --graph -20` |
| `gd` | `git diff`（查看未暂存的修改） |

```bash
# 查看修改了哪些文件
gs

# 查看具体修改内容
gd src/Simulation_Trajectory.cpp

# 提交修改
git add src/Simulation_Trajectory.cpp
git commit -m "fix: update max_scattering type to long int"

# 查看提交历史（图形化）
gl

# 查看某个文件的历史修改
git log --oneline src/Simulation_Trajectory.cpp

# 查看某次提交的内容
git show <commit-hash>

# 撤销未暂存的修改（危险！）
git checkout -- src/Simulation_Trajectory.cpp

# 暂存当前修改（临时保存）
git stash
git stash pop   # 恢复
```

---

## 附录：工具安装位置

| 工具 | 路径 | 版本 |
|---|---|---|
| `fzf` | `~/.fzf/bin/fzf` | 0.70.0 |
| `rg` | `~/.local/bin/rg` | 14.1.1 |
| `eza` | `~/.local/bin/eza` | 0.20.14 |
| `bat` | `~/.local/bin/bat` | 0.24.0 |
| `bash-completion` | `~/.bash-completion/` | latest |
| `complete-alias` | `~/.complete-alias/` | latest |
| bashrc 备份 | `~/.bashrc.bak.20260305` | — |

---

*最后更新：2026-03-05*
