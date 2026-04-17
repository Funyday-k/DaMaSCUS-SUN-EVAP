#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -o job-%j


python3 a.py
python3 a.py
python3 a.py
python3 a.py
python3 a.py
