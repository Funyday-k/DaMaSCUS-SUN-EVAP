#!/bin/bash
#SBATCH --job-name=build_damascus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

#SBATCH --time=01:00:00
#SBATCH --output=build_%j.log
#SBATCH --error=build_%j.err

set -e   # 出错立即停止

echo "===== 节点: $(hostname) ====="
echo "===== 开始时间: $(date) ====="

# 计算节点是干净环境，必须重新加载所有模块（这一步解决 libmpfr.so.4）
module purge
module load gcc/10.2.0
module load boost/1.77.0
module load openmpi/gcc/4.0.5
module load mpfr/3.1.4

echo "===== 已加载模块 ====="
module list

# 手装的 cmake 放进 PATH
export PATH=$HOME/cmake-3.12.4-Linux-x86_64/bin:$PATH

# 显式指定编译器（绝对路径，避免出岔子）
export CC=/opt/share/openmpi/4.0.5/bin/mpicc
export CXX=/opt/share/openmpi/4.0.5/bin/mpicxx

cd /lustre/project/kennyng/damascus-earth/DaMaSCUS-NewEarth-EVAP/build

# 彻底清缓存重新配置
rm -rf CMakeCache.txt CMakeFiles/

echo "===== 开始 CMake 配置 ====="
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=/opt/share/openmpi/4.0.5/bin/mpicc \
    -DCMAKE_CXX_COMPILER=/opt/share/openmpi/4.0.5/bin/mpicxx \
    -DLIBCONFIG_INCLUDE_DIRs=$HOME/libconfig-1.7.3/include \
    -DLIBCONFIGPP_LIBRARY=$HOME/libconfig-1.7.3/lib/libconfig++.so \
    -DBUILD_TESTING=OFF \
    -DCODE_COVERAGE=OFF

echo "===== 开始编译 ====="
make -j4

echo "===== 完成时间: $(date) ====="
