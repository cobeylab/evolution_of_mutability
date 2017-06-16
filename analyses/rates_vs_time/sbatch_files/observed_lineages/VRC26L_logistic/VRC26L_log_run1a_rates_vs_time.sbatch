#!/bin/bash
#SBATCH --job-name=VRC26L_log_run1a_rates_vs_time
#SBATCH --output=VRC26L_log_run1a_rates_vs_time.out
#SBATCH --error=VRC26L_log_run1a_rates_vs_time.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

cd ../../../
python rates_vs_time.py VRC26L_log_run1a
