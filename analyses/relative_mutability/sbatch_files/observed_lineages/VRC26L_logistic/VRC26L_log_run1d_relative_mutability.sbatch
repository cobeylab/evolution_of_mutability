#!/bin/bash
#SBATCH --job-name=VRC26L_log_run1d_relative_mutability
#SBATCH --output=VRC26L_log_run1d_relative_mutability.out
#SBATCH --error=VRC26L_log_run1d_relative_mutability.err
#SBATCH --time=600:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000

cd ../../../
python relative_mutability.py VRC26L_log_run1d
