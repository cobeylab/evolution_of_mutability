#!/bin/bash
#SBATCH --job-name=VRC26L_con_run1d
#SBATCH --output=VRC26L_con_run1d.out
#SBATCH --error=VRC26L_con_run1d.err
#SBATCH --time=2376:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8

cd ../../../../results/BEAST/observed_lineages/VRC26L_constant/VRC26L_con_run1d/

beast ../../../../../analyses/BEAST/observed_lineages/VRC26L_constant/VRC26L_con_run1.xml