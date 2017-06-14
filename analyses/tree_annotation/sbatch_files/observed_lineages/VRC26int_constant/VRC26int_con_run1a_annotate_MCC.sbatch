#!/bin/bash
#SBATCH --job-name=VRC26int_con_run1a_annotate_MCC
#SBATCH --output=VRC26int_con_run1a_annotate_MCC.out
#SBATCH --error=VRC26int_con_run1a_annotate_MCC.err
#SBATCH --time=20:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000

cd ../../../
python annotate_MCC.py VRC26int_con_run1a
