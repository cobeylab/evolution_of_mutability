#!/bin/bash
#SBATCH --job-name=scenario4e_rep4_rates_vs_time
#SBATCH --output=scenario4e_rep4_rates_vs_time.out
#SBATCH --error=scenario4e_rep4_rates_vs_time.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

cd ../../../
python rates_vs_time.py scenario4e_rep4
