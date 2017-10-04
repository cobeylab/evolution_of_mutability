#!/bin/bash
#SBATCH --job-name=scenario4b_rep1_rates_vs_time
#SBATCH --output=scenario4b_rep1_rates_vs_time.out
#SBATCH --error=scenario4b_rep1_rates_vs_time.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

cd ../../../
python rates_vs_time.py scenario4b_rep1
