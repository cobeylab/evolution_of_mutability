#!/bin/bash
#SBATCH --job-name=CH103L_log_run1b_contrasts
#SBATCH --output=CH103L_log_run1b_contrasts.out
#SBATCH --error=CH103L_log_run1b_contrasts.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000

cd ../../../
python contrasts.py CH103L_log_run1b
