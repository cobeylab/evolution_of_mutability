#!/bin/bash
#SBATCH --job-name=CH103_exp_run1d_mutability_vs_time
#SBATCH --output=CH103_exp_run1d_mutability_vs_time.out
#SBATCH --error=CH103_exp_run1d_mutability_vs_time.err
#SBATCH --time=100:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000

cd ../../../
python mutability_vs_time.py CH103_exp_run1d
