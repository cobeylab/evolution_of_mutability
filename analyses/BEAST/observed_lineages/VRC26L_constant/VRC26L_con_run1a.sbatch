#!/bin/bash
#SBATCH --job-name=VRC26L_con_run1a
#SBATCH --output=VRC26L_con_run1a.out
#SBATCH --error=VRC26L_con_run1a.err
#SBATCH --time=2376:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8

cd ../../../../results/BEAST/observed_lineages/VRC26L_constant/VRC26L_con_run1a/

beast ../../../../../analyses/BEAST/observed_lineages/VRC26L_constant/VRC26L_con_run1.xml