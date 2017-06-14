#!/bin/bash
#SBATCH --job-name=scenario4a_rep2_annotate_MCC
#SBATCH --output=scenario4a_rep2_annotate_MCC.out
#SBATCH --error=scenario4a_rep2_annotate_MCC.err
#SBATCH --time=30:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000

cd ../../../
python annotate_MCC.py scenario4a_rep2
