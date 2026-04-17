#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -o job-%j


mpirun -n 4 ./DaMaSCUS-SUN config1gev_-31.cfg
