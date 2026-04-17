#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH -o job-%j
#SBATCH --mail-type=end
#SBATCH --mail-user=lingyuxia@link.cuhk.edu.hk

#module load intel/2022.2
#module load cmake/3.21.2 boost/1.77.0

export LD_LIBRARY_PATH=/project/kennyng/backup_DM/lib:$LD_LIBRARY_PATH
mpirun -n 32 ./DaMaSCUS-SUN config_Lingyu.cfg
