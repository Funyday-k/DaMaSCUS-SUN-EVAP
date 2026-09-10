# DaMaSCUS-SUN 在 CUHK CHPC 集群的构建笔记

## 环境模块（关键，缺一不可）
    module purge
    module load gcc/10.2.0
    module load boost/1.77.0
    module load openmpi/gcc/4.0.5
    module load mpfr/3.1.4          # 解决 cc1plus: libmpfr.so.4 not found

## 手动安装的依赖（装在 $HOME）
- libconfig-1.7.3 -> $HOME/libconfig-1.7.3 (include/ 和 lib/)
- cmake 3.12.4    -> $HOME/cmake-3.12.4-Linux-x86_64/bin

## 运行时库路径
    export LD_LIBRARY_PATH=$HOME/libconfig-1.7.3/lib:$LD_LIBRARY_PATH

## CMake 配置（注意 LIBCONFIG_INCLUDE_DIRs 结尾是小写 s）
    cd build
    rm -rf CMakeCache.txt CMakeFiles/
    cmake .. -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER=/opt/share/openmpi/4.0.5/bin/mpicc \
        -DCMAKE_CXX_COMPILER=/opt/share/openmpi/4.0.5/bin/mpicxx \
        -DLIBCONFIG_INCLUDE_DIRs=$HOME/libconfig-1.7.3/include \
        -DLIBCONFIGPP_LIBRARY=$HOME/libconfig-1.7.3/lib/libconfig++.so \
        -DBUILD_TESTING=OFF -DCODE_COVERAGE=OFF
    make -j4

## 产物
- 可执行文件: build/src/DaMaSCUS-SUN

## 运行（MPI 程序，必须用 mpirun）
    cd bin
    mpirun --mca btl ^openib -np <N> ../build/src/DaMaSCUS-SUN config.cfg

## 编译时踩过的坑与解法
1. /usr/bin/g++ 不存在 -> 用 -DCMAKE_CXX_COMPILER 指定绝对路径
2. 登录节点 OOM Killed -> 必须在计算节点编译（srun 或 sbatch）
3. libmpfr.so.4 not found -> module load mpfr/3.1.4
4. libconfig++ not found -> 变量名是 LIBCONFIG_INCLUDE_DIRs（结尾小写 s）
5. boost gauss.hpp not found -> module load boost/1.77.0
