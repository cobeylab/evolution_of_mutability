#!/bin/bash
#SBATCH --job-name=scenario6b_rep1_annotate_MCC
#SBATCH --output=scenario6b_rep1_annotate_MCC.out
#SBATCH --error=scenario6b_rep1_annotate_MCC.err
#SBATCH --time=30:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000

cd ../../../
python annotate_MCC.py scenario6b_rep1
