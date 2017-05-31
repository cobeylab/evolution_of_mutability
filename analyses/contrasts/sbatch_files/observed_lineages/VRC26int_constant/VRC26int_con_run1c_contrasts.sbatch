#!/bin/bash
#SBATCH --job-name=VRC26int_con_run1c_contrasts
#SBATCH --output=VRC26int_con_run1c_contrasts.out
#SBATCH --error=VRC26int_con_run1c_contrasts.err
#SBATCH --time=200:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000

cd ../../../
python contrasts.py VRC26int_con_run1c
