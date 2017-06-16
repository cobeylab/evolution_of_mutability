#!/bin/bash
#SBATCH --job-name=VRC01_19_con_run1a_mutability_vs_time
#SBATCH --output=VRC01_19_con_run1a_mutability_vs_time.out
#SBATCH --error=VRC01_19_con_run1a_mutability_vs_time.err
#SBATCH --time=100:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

cd ../../../
python mutability_vs_time.py VRC01_19_con_run1a
