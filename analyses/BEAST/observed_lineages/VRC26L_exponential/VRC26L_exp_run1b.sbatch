#!/bin/bash
#SBATCH --job-name=VRC26L_exp_run1b
#SBATCH --output=VRC26L_exp_run1b.out
#SBATCH --error=VRC26L_exp_run1b.err
#SBATCH --time=500:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8

cd ../../../../results/BEAST/observed_lineages/VRC26L_exponential/VRC26L_exp_run1b/

beast ../../../../../analyses/BEAST/observed_lineages/VRC26L_exponential/VRC26L_exp_run1.xml