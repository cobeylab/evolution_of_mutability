#!/bin/bash
#SBATCH --job-name=VRC26int_con_run1c_mutability_vs_time
#SBATCH --output=VRC26int_con_run1c_mutability_vs_time.out
#SBATCH --error=VRC26int_con_run1c_mutability_vs_time.err
#SBATCH --time=100:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

cd ../../../
python mutability_vs_time.py VRC26int_con_run1c
