#!/bin/bash
#SBATCH --job-name=VRC01_19_con_run1d
#SBATCH --output=VRC01_19_con_run1d.out
#SBATCH --error=VRC01_19_con_run1d.err
#SBATCH --time=500:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8

cd ../../../../results/BEAST/observed_lineages/VRC01_19_constant/VRC01_19_con_run1d/

beast ../../../../../analyses/BEAST/observed_lineages/VRC01_19_constant/VRC01_19_con_run1.xml