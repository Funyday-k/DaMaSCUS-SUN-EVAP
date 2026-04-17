#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH -o job-%j
#SBATCH --mail-type=end --mail-user=lingyuxia@link.cuhk.edu.hk

#module load python/3.9

singularity exec /project/kennyng/eROSITA/esass-x64_latest_fix.sif python3 bincount2.py