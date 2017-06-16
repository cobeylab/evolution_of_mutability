#!/bin/bash
#SBATCH --job-name=CH103_con_run1b_rates_vs_time
#SBATCH --output=CH103_con_run1b_rates_vs_time.out
#SBATCH --error=CH103_con_run1b_rates_vs_time.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

cd ../../../
python rates_vs_time.py CH103_con_run1b
