#!/bin/bash
#SBATCH --job-name=CH103_con_run1d
#SBATCH --output=CH103_con_run1d.out
#SBATCH --error=CH103_con_run1d.err
#SBATCH --time=2376:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8

cd ../../../../results/BEAST/observed_lineages/CH103_constant/CH103_con_run1d/

beast ../../../../../analyses/BEAST/observed_lineages/CH103_constant/CH103_con_run1.xml