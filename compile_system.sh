#!/bin/bash
#SBATCH --job-name=Damascus_compile
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=compile_%j.out
#SBATCH --error=compile_%j.err

# 不加载任何模块，使用系统默认
module purge

export LD_LIBRARY_PATH=$HOME/libconfig-1.7.3/lib:$LD_LIBRARY_PATH

# 使用系统的默认gcc（通常是 /usr/bin/gcc）
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++

cd /lustre/project/kennyng/damascus-earth/DaMaSCUS-NewEarth-EVAP

# 清理旧的build
rm -rf build
mkdir build && cd build

# 配置CMake
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DCMAKE_C_COMPILER=/usr/bin/gcc \
         -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
         -DLIBCONFIGPP_LIBRARY=$HOME/libconfig-1.7.3/lib/libconfig++.so \
         -DBUILD_TESTING=OFF \
         -DCODE_COVERAGE=OFF

# 编译
make -j4

echo "编译完成，检查可执行文件："
ls -la src/DaMaSCUS-SUN 2>/dev/null || echo "可执行文件未找到"
