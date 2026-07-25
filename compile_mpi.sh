#!/bin/bash
#SBATCH --job-name=Damascus_compile
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=compile_%j.out
#SBATCH --error=compile_%j.err

# 加载 MPI
module purge
module load openmpi/gcc/3.0.1

# 环境变量
export LD_LIBRARY_PATH=$HOME/libconfig-1.7.3/lib:$LD_LIBRARY_PATH

cd /lustre/project/kennyng/damascus-earth/DaMaSCUS-NewEarth-EVAP

# 清理并创建build
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
