#!/bin/bash
#SBATCH --job-name=CH103L_exp_run1a
#SBATCH --output=CH103L_exp_run1a.out
#SBATCH --error=CH103L_exp_run1a.err
#SBATCH --time=500:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8

cd ../../../../results/BEAST/observed_lineages/CH103L_exponential/CH103L_exp_run1a/

beast ../../../../../analyses/BEAST/observed_lineages/CH103L_exponential/CH103L_exp_run1.xml