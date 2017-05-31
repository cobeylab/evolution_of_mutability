#!/bin/bash
#SBATCH --job-name=CH103_con_run1b
#SBATCH --output=CH103_con_run1b.out
#SBATCH --error=CH103_con_run1b.err
#SBATCH --time=2376:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8

cd ../../../../results/BEAST/observed_lineages/CH103_constant/CH103_con_run1b/

beast ../../../../../analyses/BEAST/observed_lineages/CH103_constant/CH103_con_run1.xml