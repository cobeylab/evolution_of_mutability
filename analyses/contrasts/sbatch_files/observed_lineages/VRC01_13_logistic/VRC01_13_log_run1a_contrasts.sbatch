#!/bin/bash
#SBATCH --job-name=VRC01_13_log_run1a_contrasts
#SBATCH --output=VRC01_13_log_run1a_contrasts.out
#SBATCH --error=VRC01_13_log_run1a_contrasts.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000

cd ../../../
python contrasts.py VRC01_13_log_run1a
