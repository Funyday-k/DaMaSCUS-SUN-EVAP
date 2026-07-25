#!/bin/bash
#SBATCH --job-name=Damascus_compile
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=compile_%j.out
#SBATCH --error=compile_%j.err

module load gcc/10.2.0
module load openmpi/gcc/4.0.5
module load boost/1.77.0

export LD_LIBRARY_PATH=$HOME/libconfig-1.7.3/lib:$LD_LIBRARY_PATH
export BOOST_ROOT=/opt/share/boost/1.77.0
export BOOST_INCLUDEDIR=/opt/share/boost/1.77.0/include

cd /lustre/project/kennyng/damascus-earth/DaMaSCUS-NewEarth-EVAP

# 清理旧的build
rm -rf build
mkdir build && cd build

# 配置CMake
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DLIBCONFIGPP_LIBRARY=$HOME/libconfig-1.7.3/lib/libconfig++.so \
         -DBUILD_TESTING=OFF \
         -DCODE_COVERAGE=OFF \
         -DBOOST_ROOT=/opt/share/boost/1.77.0 \
         -DBoost_INCLUDE_DIR=/opt/share/boost/1.77.0/include

# 编译
make -j4

echo "编译完成，检查可执行文件："
ls -la src/DaMaSCUS-SUN
