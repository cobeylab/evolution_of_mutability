#!/bin/bash
#SBATCH --job-name=VRC26_con_run1c
#SBATCH --output=VRC26_con_run1c.out
#SBATCH --error=VRC26_con_run1c.err
#SBATCH --time=2376:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8

cd ../../../../results/BEAST/observed_lineages/VRC26_constant/VRC26_con_run1c/

beast ../../../../../analyses/BEAST/observed_lineages/VRC26_constant/VRC26_con_run1.xml