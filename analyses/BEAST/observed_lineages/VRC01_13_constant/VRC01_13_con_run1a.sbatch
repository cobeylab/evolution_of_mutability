#!/bin/bash
#SBATCH --job-name=VRC01_13_con_run1a
#SBATCH --output=VRC01_13_con_run1a.out
#SBATCH --error=VRC01_13_con_run1a.err
#SBATCH --time=500:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8

cd ../../../../results/BEAST/observed_lineages/VRC01_13_constant/VRC01_13_con_run1a/

beast ../../../../../analyses/BEAST/observed_lineages/VRC01_13_constant/VRC01_13_con_run1.xml