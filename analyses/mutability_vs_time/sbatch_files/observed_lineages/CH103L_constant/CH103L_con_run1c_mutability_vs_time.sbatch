#!/bin/bash
#SBATCH --job-name=CH103L_con_run1c_mutability_vs_time
#SBATCH --output=CH103L_con_run1c_mutability_vs_time.out
#SBATCH --error=CH103L_con_run1c_mutability_vs_time.err
#SBATCH --time=100:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

cd ../../../
python mutability_vs_time.py CH103L_con_run1c
