#!/bin/bash
#SBATCH --job-name=CH103_exp_run1c_annotate_MCC
#SBATCH --output=CH103_exp_run1c_annotate_MCC.out
#SBATCH --error=CH103_exp_run1c_annotate_MCC.err
#SBATCH --time=20:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000

cd ../../../
python annotate_MCC.py CH103_exp_run1c
