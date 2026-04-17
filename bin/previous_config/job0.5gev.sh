#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -o job-%j


mpirun -n 4 ./DaMaSCUS-SUN config0.5gev_-32.cfg
