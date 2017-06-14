#!/bin/bash
#SBATCH --job-name=scenario3_rep4_annotate_MCC
#SBATCH --output=scenario3_rep4_annotate_MCC.out
#SBATCH --error=scenario3_rep4_annotate_MCC.err
#SBATCH --time=30:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000

cd ../../../
python annotate_MCC.py scenario3_rep4
